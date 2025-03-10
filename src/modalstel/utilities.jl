using StatsBase
using SoleModels: LeafModel, naturalconditions
using ModalDecisionTrees
using ModalDecisionTrees: DForest, apply

############################################################################################
################################# UTILITIES ################################################
############################################################################################

function from_to(n::Int; rng::AbstractRNG=default_rng(), one=false)
    from, to = begin
        if one
            idx = rand(rng,1:n)
            (idx,idx)
        else
            rand(rng,1:n,2)
        end
    end
    return from > to ? (to,from) : (from,to)
end

default_rng(k::Int=1) = MersenneTwister(k)

function combining(v1::Vector,v2::Vector)
    @assert length(v1) == length(v2) "Two vector must have same length"
    v = []
    for i in 1:length(v1)
        push!(v,v1[i])
        push!(v,v2[i])
    end

    return v
end

############################################################################################

nleaves(m::DecisionForest) = sum(nleaves.(trees(m)))
nleaves(m::DecisionTree) = nleaves(root(m))
nleaves(m::Branch) = nleaves(posconsequent(m)) + nleaves(negconsequent(m))
nleaves(m::LeafModel) = 1

nrules(m::DecisionList; kwargs...) = length(rulebase(m)) + 1

function kappa(
    m::AbstractModel,
    X,
    Y::AbstractVector{<:Label};
    kwargs...
)
    Y_pred = apply(m, X; kwargs...)
    return MLJ.Kappa()(Y_pred, Y)*100
end

global ApplyMemoStructure = Dict{Tuple{Union{DecisionForest,DecisionList},Any},Vector{Label}}()

function accuracy(
    m::AbstractModel,
    X,
    Y::AbstractVector{<:Label};
    kwargs...
)
    Y_pred = begin
        if haskey(ApplyMemoStructure, (m, X))
            ApplyMemoStructure[(m,X)]
        else
            ApplyMemoStructure[(m,X)] = apply(m, X; kwargs...)
        end
    end
    return (MLJ.Accuracy()(Y_pred, Y))*100
end

function _error(
    m::AbstractModel,
    X,
    Y::AbstractVector{<:Label};
    kwargs...
)
    return (100 - accuracy(m,X,Y; kwargs...))
end

function nsymbolscomplexity(m::Rule; kwargs...)
    f = antecedent(m)

    return f isa LeftmostLinearForm ? SoleLogics.nchildren(f) : 1
end

function symbolnumber(m::DecisionList; kwargs...)
    formulas = antecedent.(rulebase(m))

    return length(formulas) == 0 ? 0 : sum(ntokens.(formulas))
end

function meandelay(
    m::DecisionList,
    X;
    kwargs...
)
    newdl = apply!(m, X; compute_metrics=true, kwargs...)

    return meandelaydl(newdl)
end

function absnrulescomplexity(m::DecisionList, Y::AbstractVector{<:Label}; kwargs...)
    nclasses = length(unique(Y))
    nrls = nrules(m)
    #abs(nclasses - nrules(m))
    #
    # abs(log(nclasses,nrules(m)) - 1)*100
    return abs(log(nclasses,nrls) - 1)*exp((nclasses)^2/nrls)/2
end

function nrulescomplexity(m::DecisionList, Y::AbstractVector{<:Label}; kwargs...)
    nclasses = length(unique(Y))
    _nrules = nrules(m) - 1 # Note: Account for the default rule
    return (_nrules >= nclasses ? (1.0 - (nclasses/_nrules)) : 1.0)
end

function nsymbolscomplexity(m::DecisionList, Y::AbstractVector{<:Label}; kwargs...)
    nclasses = length(unique(Y))
    formulas = antecedent.(rulebase(m))
    nsymbols = length(formulas) == 0 ? 1 : sum(ntokens.(formulas))

    return (nsymbols >= nclasses ? (1.0 - (nclasses/nsymbols)) : 1.0)
end

############################################################################################
#################################### Final Metrics #########################################
############################################################################################

function finalmetrics(
    m::Vector{<:DecisionTree},
    X,
    Y::AbstractVector{<:Label};
    finalmetricsfuns::Vector{<:Function} = [
        kappa, _error, nrules, symbolnumber, meandelay, absnrulescomplexity,
    ],
    originalmodel::DecisionForest,
    kwargs...,
)
    dlstrees = extract_decision_lists(originalmodel, (X, Y))

    return [
        mean(nleaves.(m)),
        finalmetrics(dlstrees,X,Y)...,
    ]
end

function finalmetrics(
    m::DecisionList,
    X,
    Y::AbstractVector{<:Label};
    finalmetricsfuns::Vector{<:Function} = [
        kappa, _error, nrules, symbolnumber, meandelay, absnrulescomplexity
    ],
    kwargs...,
)
    #return map(f->f(m; kwargs...), finalmetricsfuns)

    Y_pred = begin
        if haskey(ApplyMemoStructure, (m, X))
            ApplyMemoStructure[(m,X)]
        else
            ApplyMemoStructure[(m,X)] = apply(m, X; kwargs...)
        end
    end

    return [
        (MLJ.Kappa()(Y_pred, Y)*100),
        (100 - (MLJ.Accuracy()(Y_pred, Y)*100)),
        nrules(m),
        symbolnumber(m, kwargs...),
        meandelay(m, X; kwargs...),
        absnrulescomplexity(m, Y; kwargs...),
    ]
end

function finalmetrics(
    m::Vector{<:DecisionList},
    X,
    Y::AbstractVector{<:Label};
    finalmetricsfuns::Vector{<:Function} = [
        kappa, _error, nrules, symbolnumber, meandelay, absnrulescomplexity,
    ],
    kwargs...
)
    Y_pred = map(i->begin
        if haskey(ApplyMemoStructure, (i, X))
            ApplyMemoStructure[(i,X)]
        else
            ApplyMemoStructure[(i,X)] = apply(i, X; kwargs...)
        end
        #apply(i, X; kwargs...)
    end, m)
    metrics = [
        map(j -> (MLJ.Kappa()(Y_pred[j], Y)*100), 1:length(m)),
        map(j -> 100 - (MLJ.Accuracy()(Y_pred[j], Y)*100), 1:length(m)),
        nrules.(m),
        map(j -> symbolnumber(j, kwargs...), m),
        map(j -> meandelay(j, X; kwargs...), m),
        map(j -> absnrulescomplexity(j, Y; kwargs...), m),
    ]

    return mean.(metrics)
end

function finalmetrics(
    m::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    finalmetricsfuns::Vector{<:Function} = [nleaves, kappa, _error],
    kwargs...
)
    Y_pred = begin
        if haskey(ApplyMemoStructure, (m, X))
            ApplyMemoStructure[(m,X)]
        else
            ApplyMemoStructure[(m,X)] = apply(m, X; kwargs...)
        end
    end

    return [
        nleaves(m),
        (MLJ.Kappa()(Y_pred, Y)*100),
        100-((MLJ.Accuracy()(Y_pred, Y))*100),
    ]
end

############################################################################################
################################# Modal Decision Trees #####################################
############################################################################################

function kappa(
    m::M,
    X,
    Y::AbstractVector{<:Label};
    suppress_parity_warning = false,
    kwargs...
) where {M<:Union{DTree,DForest}}
    Y_pred = begin
        if m isa DTree
            ModalDecisionTrees.apply(m, X)
        else
            ModalDecisionTrees.apply(m, X; suppress_parity_warning = suppress_parity_warning)
        end
    end

    return MLJ.Kappa()(Y_pred, Y)*100
end

function accuracy(
    m::M,
    X,
    Y::AbstractVector{<:Label};
    suppress_parity_warning = false,
    kwargs...
) where {M<:Union{DTree,DForest}}
    Y_pred = begin
        if m isa DTree
            ModalDecisionTrees.apply(m,X)
        else
            ModalDecisionTrees.apply(m,X; suppress_parity_warning=suppress_parity_warning)
        end
    end

    return (MLJ.Accuracy()(Y_pred, Y))*100
end

function _error(
    m::M,
    X,
    Y::AbstractVector{<:Label};
    suppress_parity_warning = false,
    kwargs...
) where {M<:Union{DTree,DForest}}
    return (100 - accuracy(
        m,X,Y; suppress_parity_warning=suppress_parity_warning, kwargs...,
    ))
end

function finalmetrics(
    m::M,
    X,
    Y::AbstractVector{<:Label};
    suppress_parity_warning = false,
    kwargs...
) where {M<:Union{DTree,DForest}}
    levels = unique(Y)

    Y_pred = begin
        if m isa DTree
            ModalDecisionTrees.apply(m,X)
        else
            ModalDecisionTrees.apply(m,X; suppress_parity_warning=suppress_parity_warning)
        end
    end
    Y_pred = categorical(Y_pred, ordered=true, levels = levels)

    Y_nonpred = categorical(Y, ordered=true, levels=levels)

    cm = ConfusionMatrix()(Y_pred, Y_nonpred)

    TP, FP, TN, FN = cm[1], cm[2], cm[4], cm[3]

    k = kappa(m,X,Y; suppress_parity_warning=suppress_parity_warning)
    acc = (TP+TN) / (TP+TN+FP+FN)
    sensitivity = TP / (TP+FN)
    specificity = TN / (TN+FP)

    return [k, acc*100.0, sensitivity*100.0, specificity*100.0]
end

############################################################################################
##################################### Model -> DecisionList ################################
############################################################################################

function extract_decision_lists(model::DecisionForest, args...; kwargs...)
    ts = trees(model)
    map(i -> extract_decision_lists(
            ts[i],
            args...;
            additional_info = (; tree_id = i),
            kwargs...), 1:length(ts)
        )

    #=map((i,t)->extract_decision_lists(
            t,
            args...;
            additional_info = (; tree_id = i),
            kwargs...), enumerate(trees(model)))=#
end

function extract_decision_lists(
    model::Union{MixedSymbolicModel,DecisionTree,DecisionList,Branch,Rule},
    args...;
    additional_info = (),
    kwargs...,
)
    rules = listrules(model)
    default_consequent = rand(consequent.(rules))
    extract_decision_lists(
        model,
        default_consequent,
        args...;
        additional_info = additional_info,
        kwargs...,
    )
end

function extract_decision_lists(
    model::Union{MixedSymbolicModel,DecisionTree,DecisionList,Branch,Rule},
    dataset::Union{AbstractVector,NTuple{2,Union{SoleBase.AbstractDataset,AbstractVector}}},
    args...;
    additional_info = (),
    kwargs...,
)
    Y = (dataset isa NTuple{2,Union{SoleBase.AbstractDataset,AbstractVector}} ?
                    last(dataset) : dataset)
    default_consequent = SoleModels.ConstantModel(bestguess(Y; suppress_parity_warning = true))
    extract_decision_lists(
        model,
        default_consequent,
        args...;
        additional_info = additional_info,
        kwargs...,
    )
end

function extract_decision_lists(
    model::Union{MixedSymbolicModel,DecisionTree,DecisionList,Branch,Rule},
    default_consequent::AbstractModel,
    args...;
    additional_info = (),
    kwargs...
)
    rules = listrules(model)
    i = merge(info(model), additional_info)
    #i = merge(i, (; originalrules = Ref(rules)))

    DecisionList(rules, default_consequent, i)
end

############################################################################################
############################################################################################
############################################################################################

function compute_dataset(dataset::String, isspatial::Bool, isstatic::Bool)
    println("Initializating dataset $(dataset) in ...")
    cube, Y = @time begin
        d = Serialization.deserialize(dataset)

        X, Y = begin
            if isspatial
                (X, Y) = d
                (X, Y) = round_dataset((X, Y), UInt16)
                X, Y = begin
                    if isstatic
                        #TODO check if correct
                        (X, Y), mixed_conditions, relations =
                            apply_multimodal_modes(
                                (X, Y),
                                [:static],
                                MixedCondition[StatsBase.mean],
                                AbstractRelation[],
                            )
                        (X, Y) = round_dataset((X, Y), UInt16)
                    end
                    (X, Y)
                end
                X = (X[1])

                (X, Y)
            elseif isstatic
                (X, Y), _ = d
                X = (X)
                (X, Y)
            else
                (X1, Y1), _ = d.train_n_test
                (X2, Y2), _ = d.only_training
                X1, Y1, X2, Y2 = begin
                    if isstatic
                        ((X1, Y1), (X2, Y2)), mixed_conditions, relations =
                            apply_multimodal_modes(
                                ((X1[1], Y1), (X2[1], Y2)),
                                [:static],
                                MixedCondition[SoleData.canonical_geq, SoleData.canonical_leq],
                                AbstractRelation[]
                            )
                    end
                    (X1, Y1, X2, Y2)
                end
                X, Y = concat_labeled_datasets(((X1, Y1)), ((X2, Y2)))
                @assert length(X) == 1
                X = (X[1])

                (X,Y)
            end
        end

        (X, Y)
    end

    features = MixedCondition[minimum, maximum]

    relations = begin
        if isstatic
            AbstractRelation[]
        elseif isspatial
            [globalrel, SoleLogics.RCC5Relations...]
        else
            [globalrel, SoleLogics.IARelations...]
        end
    end

    println("\nComputing scalar logiset in...")
    X = @time scalarlogiset(
        cube,
        features;
        use_full_memoization = true,
        use_onestep_memoization = false,
        #conditions = naturalconditions(cube, features),
        #relations = relations
    )
    println(SoleData.displaystructure(X))
    println()

    return (X, Y)
end
