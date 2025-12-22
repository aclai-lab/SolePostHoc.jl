import DecisionTree as DT
using  ComplexityMeasures
using  SoleModels
const  SM = SoleModels

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
function _prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where{O}
    r = Rule(LeftmostConjunctiveForm([antecedent(r)]), consequent(r), info(r))
    _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
end

function _prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r                 :: Rule{O},
    X                 :: AbstractInterpretationSet,
    y                 :: AbstractVector{<:SM.Label};
    pruning_s         :: AbstractFloat,
    pruning_decay_thr :: AbstractFloat,
    kwargs...,
) where {O}
    nruleconjuncts = SM.nconjuncts(r)
    e_zero         = SM.rulemetrics(r, X, y)[:error]
    valid_idxs     = 1:nruleconjuncts
    antd, cons     = SM.antecedent(r), SM.consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) < 2) && break

        # indices to be considered to evaluate the rule
        other_idxs = setdiff(valid_idxs, idx)
        rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

        # return error of the rule without idx-th pair
        e_minus_i = SM.rulemetrics(rule, X, y)[:error]
        decay_i   = (e_minus_i - e_zero) / max(e_zero, pruning_s)

        if decay_i ≤ pruning_decay_thr 
            # remove the idx-th pair in the vector of decisions
            valid_idxs = setdiff(valid_idxs, idx)
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

    return  length(children) < 2 ? r : begin
            r = Rule(LeftmostConjunctiveForm(children), consequent(r), info(r))
            _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
    end
end

function _prune_ruleset(
    ruleset :: Vector{<:Rule},
    X :: AbstractInterpretationSet,
    y :: AbstractVector{<:SM.Label},
    extractor :: InTreesRuleExtractor
)
    pruned = similar(ruleset)
    
    @inbounds Threads.@threads for i in eachindex(ruleset)
        pruned[i] = if ruleset[i].antecedent isa SoleLogics.BooleanTruth
            ruleset[i]  # keep BooleanTruth rules as-is (e.g., from XGBoost)
        else
            _prune_rule(
                typeof(antecedent(ruleset[i])), ruleset[i], X, y;
                pruning_s=get_pruning_s(extractor),
                pruning_decay_thr=get_pruning_decay_threshold(extractor)
            )
        end
    end
    
    return pruned
end

# ---------------------------------------------------------------------------- #
#                                    Stel                                      #
# ---------------------------------------------------------------------------- #
@inline _is_true_antecedent(ant) =
    isa(ant, BooleanTruth) && (ant == BooleanTruth(true) || string(ant) == "⊤")

@inline _get_pred_class(pred_model) =
    isa(pred_model, ConstantModel) ? pred_model.outcome : string(pred_model)

function _compute_rule_metrics(s, X, y, rule_complexity_metric)
    return if _is_true_antecedent(antecedent(s))
        pred_class = _get_pred_class(consequent(s))
        (coverage = 1.0, error = sum(y .!= pred_class) / length(y), length = 1)
    else
        metrics = SM.rulemetrics(s, X, y)
        (coverage = metrics[:coverage], error = metrics[:error], length = metrics[rule_complexity_metric])
    end
end

function _select_best_rule(rules_error, rules_coverage, rules_length, rng)
    # filter out NaN values and find candidates
    valid_mask = .!isnan.(rules_error)
    valid_idxs = findall(valid_mask)
    
    isempty(valid_idxs) && return first(valid_idxs)
    
    # minimum error
    min_error = minimum(rules_error[valid_idxs])
    candidates = findall(rules_error .== min_error)
    
    length(candidates) == 1 && return candidates[1]
    
    # maximum coverage among candidates
    max_coverage = maximum(rules_coverage[candidates])
    candidates = candidates[rules_coverage[candidates] .== max_coverage]
    
    length(candidates) == 1 && return candidates[1]
    
    # minimum length among candidates
    min_length = minimum(rules_length[candidates])
    candidates = candidates[rules_length[candidates] .== min_length]
    
    length(candidates) == 1 && return candidates[1]
    
    # random selection
    return rand(rng, candidates)
end

function stel(
    r                      :: Vector{<:Rule},
    X                      :: AbstractInterpretationSet,
    y                      :: AbstractVector{<:SM.Label};
    max_rules              :: Int64=-1,
    min_coverage           :: Float64=0.01,
    rule_complexity_metric :: Symbol=:natoms,
    rng                    :: AbstractRNG=Random.TaskLocalRNG()
)
    rules = Rule[]
    ruleset = [r..., Rule(SM.bestguess(y; suppress_parity_warning = true))]

    # filter rules by minimum coverage
    ruleset = filter(ruleset) do s
        return _is_true_antecedent(antecedent(s)) ? true : SM.rulemetrics(s, X, y)[:coverage] ≥ min_coverage
    end

    nrules = length(ruleset)
    rules_coverage = Vector{Float64}(undef, nrules)
    rules_error    = Vector{Float64}(undef, nrules)
    rules_length   = Vector{Int64}(undef, nrules)

    while true
        # check max rules limit
        if max_rules > 0 && length(rules) ≥ max_rules - 1
            return DecisionList(rules, SM.bestguess(y; suppress_parity_warning=true))
        end

        nrules = length(ruleset)
        resize!(rules_coverage, nrules)
        resize!(rules_error, nrules)
        resize!(rules_length, nrules)

        Threads.@threads for i in eachindex(ruleset)
            metrics = _compute_rule_metrics(ruleset[i], X, y, rule_complexity_metric)
            rules_coverage[i] = metrics.coverage
            rules_error[i] = metrics.error
            rules_length[i] = metrics.length
        end

        # select best rule
        idx_best = _select_best_rule(rules_error, rules_coverage, rules_length, rng)
        push!(rules, ruleset[idx_best])

        # compute remaining instances
        idx_remaining = _is_true_antecedent(antecedent(ruleset[idx_best])) ?
            Int64[] :
            findall(.!evaluaterule(ruleset[idx_best], X, y)[:checkmask,])

        # exit condition
        idx_best == nrules && return DecisionList(rules[1:(end-1)], consequent(rules[end]))
        length(idx_remaining) == 0 && DecisionList(rules, SM.bestguess(y; suppress_parity_warning = true))

        # update for next iteration
        @views begin
            X = X[idx_remaining, :]
            y = y[idx_remaining]
        end
        
        idx_best != length(ruleset) && (ruleset[idx_best] = ruleset[end])
        resize!(ruleset, length(ruleset) - 1)
        ruleset[end] = Rule(SM.bestguess(y; suppress_parity_warning = true))
    end
    error("Unexpected error.")
end

# ---------------------------------------------------------------------------- #
#                                   InTrees                                    #
# ---------------------------------------------------------------------------- #
@inline starterruleset(model::AbstractModel; kwargs...) = unique!(
    reduce(vcat, [listrules(subm; kwargs...) for subm in SM.models(model)]),
)

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
    y         :: AbstractVector{<:SM.Label}
)
    # Extract rules from model
    listrules_kwargs = (use_shortforms=true, normalize=true)
    ruleset = isensemble(model) ?
        starterruleset(model; listrules_kwargs...) :
        listrules(model; listrules_kwargs...)

    # prune rules if enabled
    get_prune_rules(extractor) && (ruleset = _prune_ruleset(ruleset, X, y, extractor))

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
