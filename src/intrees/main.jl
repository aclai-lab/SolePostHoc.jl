module InTrees

import DecisionTree as DT
using ComplexityMeasures
using SoleModels
using DataFrames
using Random

include("config.jl")

export intrees, InTreesConfig

# ---------------------------------------------------------------------------- #
#                           Intrees pruning utility                            #
# ---------------------------------------------------------------------------- #
"""
    _prune_rule(::Type{<:Atom}, r::Rule, args...; kwargs...) -> Rule

Wrap a single-atom antecedent in a `LeftmostConjunctiveForm` and delegate to
the conjunctive-form pruning method.
"""
function _prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where {O}
    r = Rule(LeftmostConjunctiveForm([antecedent(r)]), consequent(r), info(r))
    _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
end

"""
    _prune_rule(
        ::Type{<:LeftmostConjunctiveForm},
        r::Rule,
        X::AbstractInterpretationSet,
        y::AbstractVector{<:SoleModels.Label};
        pruning_s::AbstractFloat,
        pruning_decay_thr::AbstractFloat,
        kwargs...
    ) -> Rule

Greedily remove conjuncts from `r`'s antecedent whose removal does not
increase the rule's error by more than `pruning_decay_thr` (normalized by
`pruning_s`), as in the original InTrees pruning procedure.
"""
function _prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r::Rule{O},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label};
    pruning_s::AbstractFloat,
    pruning_decay_thr::AbstractFloat,
    kwargs...,
) where {O}
    nruleconjuncts = SoleModels.nconjuncts(r)
    e_zero = SoleModels.rulemetrics(r, X, y)[:error]
    valid_idxs = 1:nruleconjuncts
    antd, cons = SoleModels.antecedent(r), SoleModels.consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) < 2) && break

        # indices to be considered to evaluate the rule
        other_idxs = setdiff(valid_idxs, idx)
        rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

        # return error of the rule without idx-th pair
        e_minus_i = SoleModels.rulemetrics(rule, X, y)[:error]
        decay_i = (e_minus_i - e_zero) / max(e_zero, pruning_s)

        if decay_i ≤ pruning_decay_thr
            # remove the idx-th pair in the vector of decisions
            valid_idxs = setdiff(valid_idxs, idx)
            e_zero = e_minus_i
        end
    end

    return Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]), cons)
end

"""
    _prune_rule(::Type{<:MultiFormula}, r::Rule, args...; kwargs...) -> Rule

Flatten a `MultiFormula` antecedent into its per-modality children and
delegate to the conjunctive-form pruning method.
"""
function _prune_rule(::Type{<:MultiFormula}, r::Rule{O}, args...; kwargs...) where {O}
    @assert antecedent(r) isa MultiFormula "Cannot use this function on $(antecedent(r))"
    children = [
        MultiFormula(i_modality, modant) for (i_modality, modant) in modforms(antecedent(r))
    ]

    return length(children) < 2 ? r : begin
        r = Rule(LeftmostConjunctiveForm(children), consequent(r), info(r))
        _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
    end
end

"""
    _prune_ruleset(
        ruleset::Vector{<:Rule},
        X::AbstractInterpretationSet,
        y::AbstractVector{<:SoleModels.Label},
        config::InTreesConfig
    ) -> Vector{<:Rule}

Prune every rule in `ruleset` in parallel via [`_prune_rule`](@ref), using the
pruning parameters stored in `config`. Rules with a `BooleanTruth` antecedent
(e.g. produced by XGBoost) are kept unchanged.
"""
function _prune_ruleset(
    ruleset::Vector{<:Rule},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label},
    config::InTreesConfig
)
    pruned = similar(ruleset)

    @inbounds Threads.@threads for i in eachindex(ruleset)
        pruned[i] = if ruleset[i].antecedent isa SoleLogics.BooleanTruth
            ruleset[i]  # keep BooleanTruth rules as-is (e.g., from XGBoost)
        else
            _prune_rule(
                typeof(antecedent(ruleset[i])), ruleset[i], X, y;
                pruning_s=get_pruning_s(config),
                pruning_decay_thr=get_pruning_decay_threshold(config)
            )
        end
    end

    return pruned
end

# ---------------------------------------------------------------------------- #
#                                    Cbc                                       #
# ---------------------------------------------------------------------------- #
"""
    _select_rules_cbc(ruleset, X, y, config::InTreesConfig) -> Vector{<:Rule}

Select a compact subset of `ruleset` via the Complexity-guided Boolean
Combination (CBC) procedure: rule-coverage checkmasks are used as features to
fit a random forest against `y`, and rules whose normalized impurity
importance exceeds `get_cbc_threshold(config)` are kept, sorted by
importance (desc), error (asc), and complexity (asc).
"""
function _select_rules_cbc(ruleset, X, y, config::InTreesConfig)
    n_rules = length(ruleset)
    metrics = Matrix{Float64}(undef, n_rules, 3)
    checkmasks = Vector{BitVector}(undef, n_rules)

    Threads.@threads for i in eachindex(ruleset)
        eval_result = rulemetrics(ruleset[i], X, y)
        checkmasks[i] = eval_result[:checkmask,]
        metrics[i, 1] = eval_result[:coverage]
        metrics[i, 2] = eval_result[:error]
        metrics[i, 3] = eval_result[get_rule_complexity_metric(config)]
    end

    # build random forest for feature importance
    rf = DT.build_forest(
        y, hcat(checkmasks...),
        get_n_subfeatures(config),
        get_n_trees(config),
        get_partial_sampling(config),
        get_max_depth(config);
        rng=get_rng(config))
    importance = DT.impurity_importance(rf)
    importances = importance ./ maximum(importance)

    # select features with sufficient importance
    selected_idxs = findall(importances .> get_cbc_threshold(config))
    isempty(selected_idxs) && return ruleset

    # combine metrics with importance and original indices
    combined = hcat(
        metrics[selected_idxs, :],
        importances[selected_idxs],
        selected_idxs
    )

    # sort by importance (desc), error (asc), complexity (asc)
    sorted = sortslices(combined, dims=1, by=x->(-x[4], x[2], x[3]))

    # extract final indices, limiting if max_rules is set
    max_rules = get_max_rules(config)
    n_selected = max_rules > 0 ? min(max_rules, size(sorted, 1)) : size(sorted, 1)
    final_idxs = Int64.(sorted[1:n_selected, 5])

    return ruleset[final_idxs]
end

# ---------------------------------------------------------------------------- #
#                                    Stel                                      #
# ---------------------------------------------------------------------------- #
@inline _is_true_antecedent(ant) =
    isa(ant, BooleanTruth) && (ant == BooleanTruth(true) || string(ant) == "⊤")

@inline _get_pred_class(pred_model) =
    isa(pred_model, ConstantModel) ? pred_model.outcome : string(pred_model)

"""
    _compute_rule_metrics(s, X, y, rule_complexity_metric)
        -> NamedTuple{(:coverage, :error, :length)}

Compute coverage, error and complexity for rule `s`. The tautological
(`⊤`) default rule is handled as a special case (full coverage, length 1).
"""
function _compute_rule_metrics(s, X, y, rule_complexity_metric)
    return if _is_true_antecedent(antecedent(s))
        pred_class = _get_pred_class(consequent(s))
        (coverage=1.0, error=sum(y .!= pred_class) / length(y), length=1)
    else
        metrics = SoleModels.rulemetrics(s, X, y)
        (coverage=metrics[:coverage], error=metrics[:error], length=metrics[rule_complexity_metric])
    end
end

"""
    _select_best_rule(rules_error, rules_coverage, rules_length, rng) -> Int

Return the index of the best rule among candidates, breaking ties in order by
minimum error, then maximum coverage, then minimum length, then at random
using `rng`.
"""
function _select_best_rule(
    rules_error::Vector{T},
    rules_coverage::Vector{T},
    rules_length::Vector{Int64},
    rng::AbstractRNG
) where {T<:Float64}
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

"""
    _stel(
        r::AbstractVector{<:Rule},
        X::AbstractInterpretationSet,
        y::AbstractVector{<:SoleModels.Label};
        max_rules::Int64=-1,
        min_coverage::Float64=0.01,
        rule_complexity_metric::Symbol=:natoms,
        rng::AbstractRNG=Random.TaskLocalRNG()
    ) -> DecisionList

Sequential covering (STEL): repeatedly pick the best remaining rule (by error,
coverage, then complexity), append it to the decision list, and restrict `X`/`y`
to the instances it does not cover, until no instances remain, the tautological
default rule is picked, or `max_rules` is reached.
"""
function _stel(
    r::AbstractVector{<:Rule},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label};
    max_rules::Int64=-1,
    min_coverage::Float64=0.01,
    rule_complexity_metric::Symbol=:natoms,
    rng::AbstractRNG=Random.TaskLocalRNG()
)
    rules = Rule[]
    ruleset = [r..., Rule(SoleModels.bestguess(y; suppress_parity_warning=true))]

    # filter rules by minimum coverage
    ruleset = filter(ruleset) do s
        return _is_true_antecedent(antecedent(s)) ?
               true :
               SoleModels.rulemetrics(s, X, y)[:coverage] ≥ min_coverage
    end

    nrules = length(ruleset)
    rules_coverage = Vector{Float64}(undef, nrules)
    rules_error = Vector{Float64}(undef, nrules)
    rules_length = Vector{Int64}(undef, nrules)

    while true
        # check max rules limit
        if max_rules > 0 && length(rules) ≥ max_rules - 1
            return DecisionList(rules, SoleModels.bestguess(y; suppress_parity_warning=true))
        end

        nrules = length(ruleset)
        resize!(rules_coverage, nrules)
        resize!(rules_error, nrules)
        resize!(rules_length, nrules)

        @inbounds Threads.@threads for i in eachindex(ruleset)
            r = ruleset[i]
            if _is_true_antecedent(antecedent(r))
                rules_coverage[i] = 1.0
                pred_class = _get_pred_class(consequent(r))
                rules_error[i] = sum(@. y != pred_class) / length(y)
                rules_length[i] = 1
            else
                metrics = rulemetrics(r, X, y)
                rules_coverage[i] = metrics.coverage
                rules_error[i] = metrics.error
                rules_length[i] = metrics[rule_complexity_metric]
            end
        end

        # select best rule
        idx_best = _select_best_rule(rules_error, rules_coverage, rules_length, rng)
        push!(rules, ruleset[idx_best])

        # compute remaining instances
        idx_remaining = _is_true_antecedent(antecedent(ruleset[idx_best])) ?
                        Int64[] :
                        findall(.!evaluaterule(ruleset[idx_best], X, y)[:checkmask,])

        # exit condition
        if idx_best == length(ruleset)
            return DecisionList(rules[1:(end-1)], consequent(rules[end]))
        elseif length(idx_remaining) == 0
            return DecisionList(rules, bestguess(y; suppress_parity_warning=true))
        end

        # update for next iteration
        @views begin
            X = X[idx_remaining, :]
            y = y[idx_remaining]
        end

        deleteat!(ruleset, idx_best)
        ruleset[end] = Rule(bestguess(y; suppress_parity_warning=true))
    end
    error("Unexpected error.")
end

# ---------------------------------------------------------------------------- #
#                                   InTrees                                    #
# ---------------------------------------------------------------------------- #
@inline _starterruleset(model::AbstractModel; kwargs...) = unique!(
    reduce(vcat, [listrules(subm; kwargs...) for subm in SoleModels.models(model)]),
)

"""
    intrees(config::InTreesConfig, model::AbstractModel, X, y::AbstractVector{<:Label})
        -> DecisionList

Return a decision list which approximates the behavior of the input `model` on
the specified supervised dataset. The set of relevant and non-redundant rules
in the decision list is obtained by means of rule extraction, rule pruning,
CBC rule selection, and sequential covering (STEL), using the parameters
encoded in `config`.

# References
- Deng, Houtao. "Interpreting tree ensembles with intrees." International
  Journal of Data Science and Analytics 7.4 (2019): 277-287.

# Pipeline
1. Extract the starting ruleset from `model` (per-tree rules for ensembles,
   `listrules` otherwise).
2. If `get_prune_rules(config)`, prune every rule's antecedent
   ([`_prune_ruleset`](@ref)).
3. Select a compact subset of rules via CBC ([`_select_rules_cbc`](@ref)).
4. Build the final decision list via sequential covering ([`_stel`](@ref)).

Although the method was originally presented for forests it is hereby extended
to work with any symbolic model.

# Arguments
- `config::InTreesConfig`: Algorithm configuration (pruning, CBC, STEL
  parameters).
- `model::AbstractModel`: A single (possibly ensemble) symbolic model.
- `X::AbstractInterpretationSet`: The dataset used to evaluate rules.
- `y::AbstractVector{<:SoleModels.Label}`: Ground-truth labels for `X`.

# Returns
- `DecisionList`: The extracted decision list.

---

    intrees(model::AbstractModel, X, y::AbstractVector{<:Label}; kwargs...)
        -> DecisionList

Convenience method that builds an [`InTreesConfig`](@ref) internally from
`kwargs` and forwards the call to `intrees(config, model, X, y)`. Any keyword
argument accepted by `InTreesConfig` (e.g. `prune_rules`, `pruning_s`,
`pruning_decay_threshold`, `rule_selection_method`, `rule_complexity_metric`,
`min_coverage`, `max_rules`, `n_subfeatures`, `n_trees`, `partial_sampling`,
`max_depth`, `rng`, `dns`, `cbc_threshold`) can be passed here.

`X` can be an `AbstractInterpretationSet` or an `AbstractDataFrame` (in the
latter case it is converted via `SoleData.scalarlogiset`).

# Examples
```julia
# Default configuration
dl = intrees(model, X, y)

# Explicit config object
config = InTreesConfig(prune_rules=true, max_rules=20)
dl = intrees(config, model, X, y)

# Custom CBC / pruning parameters via keyword arguments
dl = intrees(model, X, y; n_trees=100, pruning_decay_threshold=0.1)
```

See also
[`InTreesConfig`](@ref),
[`AbstractModel`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
function intrees(
    config::InTreesConfig,
    model::AbstractModel,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label}
)
    # Extract rules from model
    listrules_kwargs = (use_shortforms=true, normalize=true)
    set = isensemble(model) ?
          _starterruleset(model; listrules_kwargs...) :
          listrules(model; listrules_kwargs...)

    # prune rules if enabled
    get_prune_rules(config) && (set = _prune_ruleset(set, X, y, config))

    # rule selection
    ruleset = if get_rule_selection_method(config) == :CBC
        _select_rules_cbc(set, X, y, config)
    else
        error("Unexpected rule selection method: $(get_rule_selection_method(config))")
    end

    # construct final decision list via sequential covering
    _stel(
        ruleset, X, y;
        max_rules=get_max_rules(config),
        min_coverage=get_min_coverage(config),
        rule_complexity_metric=get_rule_complexity_metric(config),
        rng=get_rng(config)
    )
end

intrees(
    config::InTreesConfig,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label},
    model::AbstractModel
) = intrees(config, model, X, y)

intrees(config::InTreesConfig, m, X::AbstractDataFrame, y) =
    intrees(config, m, SoleData.scalarlogiset(X; allow_propositional=true), y)

function intrees(
    model::AbstractModel,
    X,
    y::AbstractVector{<:SoleModels.Label};
    kwargs...
)
    config = InTreesConfig(; kwargs...)
    return intrees(config, model, X, y)
end

end