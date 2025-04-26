using Evolutionary
using Plots
using Random

############################################################################################
############################################################################################
############################################################################################

"""
    samplingensemble(
        m::Union{DForest,DecisionForest},
        X::SoleBase.AbstractDataset,
        Y::AbstractVector{<:Label};
        policy = :increment,
        kwargs...,
    )

Function that applies a sampling method to input forest and returns a vector of subforests.
If the policy is not passed, then by default the incremental policy is used. It's possible
to pass a custom function to sampling the forest. If the number of subforests is one, then
this is wrapped in a vector because the function returns a vector.

A new sampling function must accept forest, dataset and label as arguments.
"""
function samplingensemble(
    m::DForest,
    policy::Function,
    X::SoleBase.AbstractDataset,
    Y::AbstractVector{<:Label};
    kwargs...,
)
    submodel = policy(m,X,Y; kwargs...)

    return submodel isa AbstractVector ? submodel : [submodel]
end

############################################################################################
####################################### POLICIES ###########################################
############################################################################################

function ranticsampling(
    forest::DForest,
    X::SoleBase.AbstractDataset,
    Y::AbstractVector{<:Label};
    nforests::Integer = 10,
    rng::AbstractRNG = default_rng(),
    fileplots::Union{Nothing,String} = nothing,
    kwargs...,
)
    ftrees = ModalDecisionTrees.trees(forest)
    interval = 5:round(Integer, 0.15*length(ftrees))
    paretoforests = []

    for i in 1:nforests
        ntrees = rand(rng, interval)
        rtrees = StatsBase.sample(rng, ftrees, ntrees; replace=false)

        push!(paretoforests, DForest(rtrees))
    end

    paretoaccuracy = accuracy.(paretoforests; X=X, Y=Y, suppress_parity_warning=true)

    if !isnothing(fileplots)
        paretontrees = map(f -> length(ModalDecisionTrees.trees(f)), paretoforests)
        idxsort = sortperm(paretontrees)

        p = plot(paretontrees[idxsort], paretoaccuracy[idxsort])
        savefig(p, fileplots)
    end

    return paretoforests[findmax(paretoaccuracy)[2]]
end

function incrementsampling(
    forest::DForest,
    X::SoleBase.AbstractDataset,
    Y::AbstractVector{<:Label};
    step::Integer = 1,
    fileplots::Union{Nothing,String} = nothing,
    kwargs...,
)
    ftrees = ModalDecisionTrees.trees(forest)
    fmetrics = ModalDecisionTrees.metrics(forest)
    facc = accuracy(forest; X=X, Y=Y, suppress_parity_warning=true)
    ntrees = length(ftrees)
    if step > ntrees
        return forest
    end

    paretoforests = []
    paretoaccuracy = []
    for i in collect(1:step:ntrees)
        for j in i:i+step
            #check indexes
            interval = begin
                j > ntrees ? [collect(i:ntrees)..., collect(1:(j-ntrees))...] : collect(i:j)
            end
            f = DForest(ftrees[interval])
            tmpacc = accuracy(f; X=X, Y=Y, suppress_parity_warning=true)

            if interval[end] == (i-1) || abs(facc-tmpacc) <= 5.0
                push!(paretoforests,f)
                push!(paretoaccuracy,tmpacc)
                break
            end
        end
    end

    if !isnothing(fileplots)
        paretontrees = map(f -> length(ModalDecisionTrees.trees(f)), paretoforests)
        idxsort = sortperm(paretontrees)

        p = plot(paretontrees[idxsort], paretoaccuracy[idxsort])
        savefig(p, fileplots)
    end

    return paretoforests[findmax(paretoaccuracy)[2]]
end

function geneticsampling(
    forest::Union{DForest,DecisionForest},
    X::SoleBase.AbstractDataset,
    Y::AbstractVector{<:Label};
    npopulation::Integer = 200,
    bounds::Vector{<:Vector{<:Integer}} = [[0,1],[100,10]],
    rng::AbstractRNG=MersenneTwister(1),
    kwargs...,
)
    function _evaluate(indiv::BitArray{1})
        return [
            _error(DForest(alltrees[indiv]); X=X, Y=Y, suppress_parity_warning = true),
            sum(indiv),
        ]
    end

    function _evaluate(indivs::Vector{<:BitArray{1}})
        return map(indiv->_evaluate(indiv), indivs)
    end

    alltrees = ModalDecisionTrees.trees(forest)
    population = [bitrand(rng,length(alltrees)) for _ in 1:npopulation]
    constrs = Evolutionary.BoxConstraints(bounds[1], bounds[2])
    alg = Evolutionary.NSGA2(; populationSize = npopulation)
    opts = Evolutionary.Options(rng=rng, iterations=500, parallelization=:thread)
    F = [0.0,1]
    res = @time Evolutionary.optimize(
        indiv->_evaluate(indiv), F, constrs, alg, population, opts,
    )

    bestforests = Evolutionary.minimizer(res)
    paretoforests = map(b->DForest(alltrees[b]), bestforests)
    paretoaccuracy = accuracy.(paretoforests)
    #idxsort = sortperm(accuracy.(forests; X=X,Y=Y,suppress_parity_warning=true); rev=true)
    #return forests[idxsort[1]]

    return paretoforests[findmax(paretoaccuracy)[2]]
end

############################################################################################
####################################### OLD POLICIES #######################################
############################################################################################

function besttrees(
    forest::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    ntrees::Integer = 10,
    memostruct,
    kwargs...,
)
    alltrees = trees(forest)
    macc = map(t->accuracy(t; X=X, Y=Y, memostruct=memostruct), alltrees)
    idxsbestacc = sortperm(macc; rev=true)

    return DecisionForest(alltrees[idxsbestacc][1:ntrees])
end

function randomtrees(
    forest::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    ntrees::Integer = 10,
    rng::AbstractRNG=MersenneTwister(1),
    kwargs...,
)
    alltrees = trees(forest)

    return DecisionForest(rand(rng, alltrees, ntrees))
end

function bestrandomsubforest(
    forest::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    ntrees::Integer = 10,
    nattempts::Integer = 100,
    rng::AbstractRNG=MersenneTwister(1),
    memostruct,
    kwargs...,
)
    best_forest = nothing
    best_accforest = -1
    alltrees = trees(forest)
    randcombs = [
        StatsBase.sample(rng, 1:length(alltrees), ntrees, replace=false)
        for _ in 1:nattempts
    ]

    for ts in randcombs
        f = DecisionForest(alltrees[ts])
        acc = accuracy(f; X=X, Y=Y, memostruct=memostruct, suppress_parity_warning = true)

        if acc > best_accforest
            best_forest = f
            best_accforest = acc
        end
    end

    return best_forest
end

function bestelementssampling(
    forest::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    metric::Function=accuracy,
    nbests::Integer=10,
    rev::Bool=true,
    kwargs...,
)
    v = trees(forest)
    nbests = nbests > length(v) ? length(v) : nbests

    bests = map(t->metric(t,X,Y; kwargs...), v)
    idxsbests = sortperm(bests; rev=rev)

    return DecisionForest(v[idxsbests][1:nbests])
end

function reservationsampling(
    forest::DecisionForest,
    args...;
    k::Integer=10,
    rng::AbstractRNG=MersenneTwister(1),
    kwargs...,
)
    v = trees(forest)
    k = k > length(v) ? length(v) : k

    reserve = v[1:k]
    for (i,elem) in enumerate(v)
        n = rand(rng, 1:i)
        if n <= k
            reserve[n] = elem
        end
    end

    return DecisionForest(reserve)
end

function bestsubforestsampling(
    forest::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    ntrees::Integer = 10,
    nattempts::Integer = 100,
    nforests::Integer = 1,
    memostruct,
    rng::AbstractRNG=MersenneTwister(1),
    kwargs...,
)
    alltrees = trees(forest)
    bestsforest = []
    minaccsforest = []
    randcombs = [
        StatsBase.sample(rng, 1:length(alltrees), ntrees, replace=false)
        for _ in 1:nattempts
    ]

    X1, Y1, X2, Y2 = begin
        (d1_ids,d2_ids) , _ = balanced_cv_dataset_slices(Y, 2, rng)
        X1,Y1 = @time slicedataset((X,Y), d1_ids; return_view = true)
        X2,Y2 = @time slicedataset((X,Y), d2_ids; return_view = true)

        (X1, Y1, X2, Y2)
    end

    for (i,ts) in enumerate(randcombs)
        f = DecisionForest(alltrees[ts])
        acc = mean([
            accuracy(f; X=X1, Y=Y1, memostruct=memostruct, suppress_parity_warning = true),
            accuracy(f; X=X2, Y=Y2, memostruct=memostruct, suppress_parity_warning = true),
        ])

        if length(bestsforest) < nforests
            push!(bestsforest, f)
            push!(minaccsforest, acc)
        elseif acc > min(minaccsforest...)
            minimum = min(minaccsforest...)
            idxmin = rand(rng, findall(minaccsforest .== minimum))
            bestsforest[idxmin] = f
            minaccsforest[idxmin] = acc
        end
    end

    return bestsforest
end
