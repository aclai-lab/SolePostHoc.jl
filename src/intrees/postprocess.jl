# ---------------------------------------------------------------------------- #
#                                   utils                                      #
# ---------------------------------------------------------------------------- #
@inline _is_true_antecedent(ant) =
    isa(ant, BooleanTruth) && (ant == BooleanTruth(true) || string(ant) == "⊤")

@inline _get_pred_class(pred_model) =
    isa(pred_model, ConstantModel) ? pred_model.outcome : string(pred_model)

# """
#     _compute_rule_metrics(s, X, y, rule_complexity_metric)
#         -> NamedTuple{(:coverage, :error, :length)}

# Compute coverage, error and complexity for rule `s`. The tautological
# (`⊤`) default rule is handled as a special case (full coverage, length 1).
# """
# function _compute_rule_metrics(s, X, y, rule_complexity_metric)
#     return if _is_true_antecedent(antecedent(s))
#         pred_class = _get_pred_class(consequent(s))
#         (coverage=1.0, error=sum(y .!= pred_class) / length(y), length=1)
#     else
#         metrics = SoleModels.rulemetrics(s, X, y)
#         (coverage=metrics[:coverage], error=metrics[:error], length=metrics[rule_complexity_metric])
#     end
# end

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

# ---------------------------------------------------------------------------- #
#                          stel sequential covering                            #
# ---------------------------------------------------------------------------- #
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
function postprocess(
    config::InTreesConfig{T,S},
    set::Vector{<:Rule{O}},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label}
) where {T,S<:STEL,O}
    max_rules = get_max_rules(config)
    rules = Rule[]
    ruleset = [set..., Rule(SoleModels.bestguess(y; suppress_parity_warning=true))]

    # filter rules by minimum coverage
    ruleset = filter(ruleset) do s
        return _is_true_antecedent(antecedent(s)) ?
               true :
               SoleModels.rulemetrics(s, X, y)[:coverage] ≥ get_min_coverage(config)
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
                rules_length[i] = metrics[get_complexity_metric(config)]
            end
        end

        # select best rule
        idx_best = _select_best_rule(rules_error, rules_coverage, rules_length, get_rng(config))
        push!(rules, ruleset[idx_best])

        # compute remaining instances
        idx_remaining = _is_true_antecedent(antecedent(ruleset[idx_best])) ?
            Int[] :
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