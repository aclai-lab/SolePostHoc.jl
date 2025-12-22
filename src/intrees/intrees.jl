import DecisionTree as DT
using ComplexityMeasures
using SoleModels

# ---------------------------------------------------------------------------- #
#                                InTrees struct                                #
# ---------------------------------------------------------------------------- #
"""
    InTreesRuleExtractor(; kwargs...)

Create a rule extractor based on the InTrees method.

# Keyword Arguments
- `prune_rules::Bool=true`: access to prune or not
- `pruning_s::Union{Float64,Nothing}=nothing`: parameter that limits the denominator in the pruning metric calculation
- `pruning_decay_threshold::Union{Float64,Nothing}=nothing`: threshold used in pruning to remove or not a joint from the rule
- `rule_selection_method::Symbol=:CBC`: rule selection method. Currently only supports `:CBC`
- `rule_complexity_metric::Symbol=:natoms`: Metric to use for estimating a rule complexity measure
- `min_coverage::Union{Float64,Nothing}=nothing`: minimum rule coverage for stel
- `rng::AbstractRNG=Random.TaskLocalRNG()`: RNG used for any randomized steps (e.g., feature selection)

See also [`intrees`](@ref).
"""
struct InTreesRuleExtractor <: RuleExtractor
    prune_rules             :: Bool
    pruning_s               :: Float64
    pruning_decay_threshold :: Float64
    min_coverage            :: Float64
    max_rules               :: Int64
    rule_selection_method   :: Symbol
    rule_complexity_metric  :: Symbol
    rng                     :: AbstractRNG

    function InTreesRuleExtractor(;
        prune_rules             :: Bool=true,
        pruning_s               :: Float64=1.0e-6,
        pruning_decay_threshold :: Float64=0.05,
        min_coverage            :: Float64=0.01,
        max_rules               :: Int64=-1,
        rule_selection_method   :: Symbol=:CBC,
        rule_complexity_metric  :: Symbol=:natoms,
        # accuracy_rule_selection = nothing,
        rng                     :: AbstractRNG=Random.TaskLocalRNG()
    )
        new(
            prune_rules,
            pruning_s,
            pruning_decay_threshold,
            min_coverage,
            max_rules,
            rule_selection_method,
            rule_complexity_metric,
            rng
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
get_prune_rules(r::InTreesRuleExtractor)             = r.prune_rules
get_pruning_s(r::InTreesRuleExtractor)               = r.pruning_s
get_pruning_decay_threshold(r::InTreesRuleExtractor) = r.pruning_decay_threshold
get_min_coverage(r::InTreesRuleExtractor)            = r.min_coverage
get_max_rules(r::InTreesRuleExtractor)               = r.max_rules
get_rule_selection_method(r::InTreesRuleExtractor)   = r.rule_selection_method
get_rule_complexity_metric(r::InTreesRuleExtractor)  = r.rule_complexity_metric
get_rng(r::InTreesRuleExtractor)                     = r.rng

# ---------------------------------------------------------------------------- #
#                           Intrees pruning utility                            #
# ---------------------------------------------------------------------------- #
"""
    intrees_prunerule(rc::Rule)::Rule

Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
considering the error metric

See also
[`Rule`](@ref),
[`rulemetrics`](@ref).
"""
@inline intrees_prunerule(
    r::Rule,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label};
    kwargs...
) = _prune_rule(typeof(antecedent(r)), r, X, y; kwargs...)

@inline _prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where{O} =
    intrees_prunerule(Rule(LeftmostConjunctiveForm(
        [antecedent(r)]), consequent(r), info(r)), args...; kwargs...,)

function _prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r                 :: Rule{O},
    X                 :: AbstractInterpretationSet,
    y                 :: AbstractVector{<:SoleModels.Label};
    pruning_s         :: AbstractFloat,
    pruning_decay_thr :: AbstractFloat,
    kwargs...,
) where {O}
    nruleconjuncts = nconjuncts(r)
    e_zero         = rulemetrics(r, X, y)[:error]
    valid_idxs     = 1:nruleconjuncts
    antd, cons     = antecedent(r), consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) < 2) && break

        # Indices to be considered to evaluate the rule
        other_idxs = intersect!(vcat(1:(idx-1), (idx+1):nruleconjuncts), valid_idxs)
        rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

        # Return error of the rule without idx-th pair
        e_minus_i = rulemetrics(rule, X, y)[:error]
        decay_i   = (e_minus_i - e_zero) / max(e_zero, pruning_s)

        if decay_i ≤ pruning_decay_thr 
             # Remove the idx-th pair in the vector of decisions
            deleteat!([valid_idxs...], idx)
            e_zero = e_minus_i
        end
    end

    return Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]), cons)
end

function _prune_rule(::Type{<:MultiFormula}, r::Rule{O}, args...; kwargs...) where {O}
    @assert antecedent(r) isa MultiFormula "Cannot use this function on $(antecedent(r))"
    children = [
        MultiFormula(i_modality, modant) for (i_modality, modant) in modforms(antecedent(r))
    ]

    return  length(children) < 2 ? r :
            intrees_prunerule(
                Rule(LeftmostConjunctiveForm(children), consequent(r), info(r)),
                args...;
                kwargs...,
            )
end

# ---------------------------------------------------------------------------- #
#                                   InTrees                                    #
# ---------------------------------------------------------------------------- #
"""
    intrees(model::Union{AbstractModel,DecisionForest}, X, y::AbstractVector{<:Label}; kwargs...)::DecisionList

Return a decision list which approximates the behavior of the input `model` on the specified supervised dataset.
The set of relevant and
non-redundant rules in the decision list are obtained by means of rule selection, rule pruning,
and sequential covering (stel).

# References
- Deng, Houtao. "Interpreting tree ensembles with intrees." International Journal of Data Science and Analytics 7.4 (2019): 277-287.

# Keyword Arguments
- `prune_rules::Bool=true`: access to prune or not
- `pruning_s::Union{Float64,Nothing}=nothing`: parameter that limits the denominator in the pruning metric calculation
- `pruning_decay_threshold::Union{Float64,Nothing}=nothing`: threshold used in pruning to remove or not a joint from the rule
- `rule_selection_method::Symbol=:CBC`: rule selection method. Currently only supports `:CBC`
- `rule_complexity_metric::Symbol=:natoms`: Metric to use for estimating a rule complexity measure
- `max_rules::Int=-1`: maximum number of rules in the final decision list (excluding default rule). Use -1 for unlimited rules.
- `min_coverage::Union{Float64,Nothing}=nothing`: minimum rule coverage for stel
- See [`modalextractrules`](@ref) keyword arguments...

Although the method was originally presented for forests it is hereby extended to work with any symbolic models.

See also
[`AbstractModel`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
function intrees(
    extractor :: InTreesRuleExtractor;
    model     :: AbstractModel,
    X         :: AbstractInterpretationSet,
    y         :: AbstractVector{<:SoleModels.Label}
)
    """
        cfs()

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`Rule`](@ref),
    [`rulemetrics`](@ref).
    """
    function cfs(X::AbstractInterpretationSet, y::AbstractVector{<:Label})
        @show typeof(_x)
        @show typeof(_y)
        entropyd(_x) = ComplexityMeasures.entropy(probabilities(_x))
        midd(_x, _y) = -entropyd(collect(zip(_x, _y)))+entropyd(_x)+entropyd(_y)
        information_gain(f1, f2) = entropyd(f1) - (entropyd(f1) - midd(f1, f2))
        su(f1, f2) = (2.0 * information_gain(f1, f2) / (entropyd(f1) + entropyd(f2)))
        function merit_calculation(X, y::AbstractVector{<:Label})
            n_samples, n_features = size(X)
            rff = 0
            rcf = 0

            for i in collect(1:n_features)
                fi = X[:, i]
                rcf += su(fi, y)  # su is the symmetrical uncertainty of fi and y
                for j in collect(1:n_features)
                    if j > i
                        fj = X[:, j]
                        rff += su(fi, fj)
                    end
                end
            end

            rff *= 2
            merits = rcf / sqrt(n_features + rff)
        end

        n_samples, n_features = size(X)
        F = [] # vector of returned indices samples
        M = [] # vector which stores the merit values

        while true
            merit = -100000000000
            idx = -1

            for i in collect(1:n_features)
                if i ∉ F
                    append!(F, i)
                    idxs_column = F[findall(F .> 0)]
                    t = merit_calculation(X[:, idxs_column], y)

                    if t > merit
                        merit = t
                        idx = i
                    end

                    pop!(F)
                end
            end

            append!(F, idx)
            append!(M, merit)

            (length(M) > 5) &&
                (M[length(M)] <= M[length(M)-1]) &&
                (M[length(M)-1] <= M[length(M)-2]) &&
                (M[length(M)-2] <= M[length(M)-3]) &&
                (M[length(M)-3] <= M[length(M)-4]) &&
                break
        end

        valid_idxs = findall(F .> 0)
        # @show F
        # @show valid_idxs
        return F[valid_idxs]
    end

    @inline starterruleset(model::AbstractModel; kwargs...) =
        unique!(
            reduce(vcat, [listrules(subm; kwargs...) for subm in SoleModels.models(model)]),
        )

    # Extract rules from each tree, obtain full ruleset
    listrules_kwargs = (use_shortforms=true, normalize=true)
    ruleset =
        isensemble(model) ?
            starterruleset(model; listrules_kwargs...) :
            listrules(model; listrules_kwargs...)
    fnames = featurenames(model)

    # Prune rules with respect to a dataset
    if get_prune_rules(extractor)
        ruleset = begin
            afterpruningruleset = Vector{Rule}(undef, length(ruleset))

            Threads.@threads for i in axes(ruleset, 1)
                if ruleset[i].antecedent isa SoleLogics.BooleanTruth
                    # this case happens with XgBoost: the rule is a simply BooleanTruth
                    afterpruningruleset[i] = ruleset[i]
                else
                    afterpruningruleset[i] =
                        intrees_prunerule(
                            ruleset[i], X, y;
                            pruning_s=get_pruning_s(extractor),
                            pruning_decay_thr=get_pruning_decay_threshold(extractor)
                        )
                end
            end
            afterpruningruleset
        end
    end

    # Rule selection to obtain the best rules
    ruleset = begin
        if get_rule_selection_method(extractor) == :CBC
            matrixrulemetrics     = Matrix{Float64}(undef, length(ruleset), 3)
            afterselectionruleset = Vector{BitVector}(undef, length(ruleset))

            Threads.@threads for i in axes(ruleset, 1)
                eval_result = rulemetrics(ruleset[i], X, y)
                afterselectionruleset[i] = eval_result[:checkmask,]
                matrixrulemetrics[i, 1]  = eval_result[:coverage]
                matrixrulemetrics[i, 2]  = eval_result[:error]
                matrixrulemetrics[i, 3]  = eval_result[get_rule_complexity_metric(extractor)]
            end

            rf = DT.build_forest(y, hcat(afterselectionruleset...), 2, 50, 0.7, -1; rng=get_rng(extractor))
            importances = begin
                importance = DT.impurity_importance(rf)
                importance/max(importance...)
            end

            best_idxs = begin
                selected_features = findall(importances .> 0.01)

                ruleset_pruned_rrf = hcat(
                    matrixrulemetrics[selected_features, :],
                    importances[selected_features],
                    selected_features,
                )
                finalmatrix = sortslices(
                    ruleset_pruned_rrf,
                    dims = 1,
                    by = x->(x[4], x[2], x[3]),
                    rev = true,
                )

                # Get all selected rules indices or limit if max_rules is specified
                get_max_rules(extractor) > 0 ?
                    Int64.(finalmatrix[1:min(max_rules, size(finalmatrix, 1)), 5]) :
                    Int64.(finalmatrix[:, 5])
            end

            ruleset[best_idxs]
        else
            error("Unexpected rule selection method specified: $(rule_selection_method)")
        end
    end

    # Construct a rule-based model from the set of best rules
    stel(
        ruleset, X, y;
        max_rules=get_max_rules(extractor),
        min_coverage=get_min_coverage(extractor), 
        rule_complexity_metric=get_rule_complexity_metric(extractor),
        rng=get_rng(extractor)
    )
end





function stel(
    ruleset,
    X,
    y;
    max_rules::Int = -1,
    min_coverage,
    rule_complexity_metric = :natoms,
    rng::AbstractRNG = MersenneTwister(1)
)
    D = deepcopy(X) # Copy of the original dataset
    L = deepcopy(y)
    R = Rule[]      # Ordered rule list
    # S is the vector of rules left
    S = [deepcopy(ruleset)..., Rule(bestguess(y; suppress_parity_warning = true))]
    # Rules with a frequency less than min_coverage
    S = begin
        rules_coverage = Vector{Float64}(undef, length(S))

        Threads.@threads for (i, s) in collect(enumerate(S))
            if isa(antecedent(s), BooleanTruth)
                # Se l'antecedente è BooleanTruth, verifica il valore
                if antecedent(s) == BooleanTruth(true) || string(antecedent(s)) == "⊤"
                    rules_coverage[i] = 1.0
                else
                    rules_coverage[i] = 0.0
                end
            else
                rules_coverage[i] = rulemetrics(s, X, y)[:coverage]
            end
        end

        idxs_undeleted = findall(rules_coverage .>= min_coverage)
        S[idxs_undeleted]
    end

    # Metrics update based on remaining instances
    rules_coverage = Vector{Float64}(undef, length(S))
    rules_error = Vector{Float64}(undef, length(S))
    rules_length = Vector{Int}(undef, length(S))

    while true
        # Check if we've reached the maximum number of rules (excluding the default rule)
        if max_rules > 0 && length(R) >= max_rules - 1
            return DecisionList(R, bestguess(L; suppress_parity_warning = true))
        end

        resize!(rules_coverage, length(S))
        resize!(rules_error, length(S))
        resize!(rules_length, length(S))

        Threads.@threads for (i, s) in collect(enumerate(S))
            if isa(antecedent(s), BooleanTruth)                                             #TODO ASK MICHI {}
                # Gestione speciale per BooleanTruth
                if antecedent(s) == BooleanTruth(true) || string(antecedent(s)) == "⊤"
                    rules_coverage[i] = 1.0
                    # Calcola l'errore basato sulla distribuzione delle classi
                    pred_model = consequent(s)
                    # Estrai il valore effettivo dal modello
                    if isa(pred_model, ConstantModel)
                        pred_class = pred_model.outcome  # Assumo che il valore sia accessibile con .value
                    else
                        pred_class = string(pred_model)
                    end
                    rules_error[i] = sum(L .!= pred_class) / length(L)
                else
                    rules_coverage[i] = 0.0
                    rules_error[i] = 1.0  # Errore massimo per false
                end
                rules_length[i] = 1  # Lunghezza minima per BooleanTruth
            else # TODO ASK MICHI }
                metrics = rulemetrics(s, D, L)
                rules_coverage[i] = metrics[:coverage]
                rules_error[i] = metrics[:error]
                rules_length[i] = metrics[rule_complexity_metric]
            end # TODO ASK MICHI ~
        end

        # Best rule index
        idx_best = begin
            idx_best = nothing
            # First: find the rule with minimum error
            m = minimum(filter(!isnan, rules_error))
            best_idxs = findall(rules_error .== m)
            (length(best_idxs) == 1) && (idx_best = best_idxs[1])

            # If not one, find the rule with maximum coverage
            if isnothing(idx_best)
                # @show rules_coverage[best_idxs]
                m = maximum(filter(!isnan, rules_coverage[best_idxs]))
                idx_coverage = findall(rules_coverage[best_idxs] .== m)
                best_idxs = best_idxs[idx_coverage]
                (length(best_idxs) == 1) && (idx_best = best_idxs[1])
            end

            # If not one, find the rule with minimum length
            if isnothing(idx_best)
                m = minimum(filter(!isnan, rules_length[best_idxs]))
                idx_length = findall(rules_length[best_idxs] .== m)
                best_idxs = best_idxs[idx_length]
                (length(best_idxs) == 1) && (idx_best = best_idxs[1])
            end

            # Final case: more than one rule with minimum length
            # Randomly choose a rule
            if isnothing(idx_best)
                idx_best = rand(rng, best_idxs)
            end

            idx_best
        end

        # Add at the end the best rule
        push!(R, S[idx_best])

        # Indices of the remaining instances
        idx_remaining = begin
            if isa(antecedent(S[idx_best]), BooleanTruth)                                                       # TODO ASK MICHI {
                if antecedent(S[idx_best]) == BooleanTruth(true) ||
                   string(antecedent(S[idx_best])) == "⊤"
                    # Se l'antecedente è true, nessuna istanza rimane
                    Int[]
                else
                    # Se l'antecedente è false, tutte le istanze rimangono
                    collect(1:length(L))
                end
            else                                                                                                # TODO ASK MICHI  }
                sat_unsat = evaluaterule(S[idx_best], D, L)[:checkmask,]
                # Remain in D the rule that not satisfying the best rule's pruning_s condition
                findall(sat_unsat .== false)
            end                                                                                                 # TODO ASK MICHI ~
        end

        # Exit condition
        if idx_best == length(S)
            # @show S
            # @show R
            # @show DecisionList(R[1:end-1], consequent(R[end]))
            return DecisionList(R[1:(end-1)], consequent(R[end]))
        elseif length(idx_remaining) == 0
            return DecisionList(R, bestguess(y; suppress_parity_warning = true))
        end

        D = slicedataset(D, idx_remaining; return_view = true)
        L = L[idx_remaining]
        # Delete the best rule from S
        deleteat!(S, idx_best)
        # Update of the default rule
        S[end] = Rule(bestguess(L; suppress_parity_warning = true))
    end
    error("Unexpected error.")
end
