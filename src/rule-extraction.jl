module RuleExtraction

export intrees

using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet
using SoleFeatures: findcorrelation
using SoleModels
using SoleModels: AbstractModel
using SoleModels: Rule, antecedent, consequent, rule_metrics
using SoleModels: FinalModel, Branch, DecisionForest, DecisionList
using SoleModels: RuleCascade, antecedents, convert, unroll_rules_cascade
using SoleModels: best_guess, Label, evaluate_rule
using SoleData: slice_dataset
using Statistics: cor
# using ModalDecisionTrees: MultiFrameModalDataset

############################################################################################
# Rule extraction from random forest
############################################################################################

"""
    intrees(args...; kwargs...)::DecisionList

Extracts rules from a model, reduces the length of each rule (number of variable value
pairs), and applies a sequential coverage approach to obtain a set of relevant and
non-redundant rules

# Arguments
- `model::Union{AbstractModel,DecisionForest}`: input model
- `X::Any`: dataset
- `Y::AbstractVector{<:Label}`: label vector

# Keywords
- `prune_rules::Bool=true`: access to prune or not
- `pruning_s::Float=nothing`: parameter that limits the denominator in the pruning metric calculation
- `pruning_decay_threshold::Float=nothing`: threshold used in pruning to remove or not a joint from the rule
- `method_rule_selection::Symbol=:CBC`: method of rule selection
- `accuracy_rule_selection::Float=nothing`: percentage of rules that rule selection must follow
- `min_frequency::Float=nothing`: minimum frequency that the rules must have at the beginning of the definition of the learner

# Returns
- `DecisionList`: decision list that represent a new learner

TODO cite paper

See also
[`AbstractModel`](@ref),
[`DecisionForest`](@ref),
[`DecisionList`](@ref),
[`unroll_rules_cascade`](@ref),
[`rule_metrics`](@ref).
"""
# Extract rules from a forest, with respect to a dataset
function intrees(
    model::Union{AbstractModel,DecisionForest},
    # X::Union{AbstractInterpretationSet,MultiFrameModalDataset}, # TODO
    X::Any,
    Y::AbstractVector{<:Label};
    #
    prune_rules = true,
    pruning_s = nothing,
    pruning_decay_threshold = nothing,
    #
    method_rule_selection = :CBC,
    accuracy_rule_selection = nothing,
    min_frequency = nothing,
)

    isnothing(pruning_s) && !isnothing(pruning_decay_threshold) && (prune_rules = false)
    isnothing(pruning_decay_threshold) && !isnothing(pruning_s) && (prune_rules = false)
    isnothing(pruning_s) && (pruning_s = 1.0e-6)
    isnothing(pruning_decay_threshold) && (pruning_decay_threshold = 0.05)
    isnothing(accuracy_rule_selection) && (accuracy_rule_selection = 0.0)
    isnothing(min_frequency) && (min_frequency = 0.01)

    @assert model isa DecisionForest || SoleModels.issymbolic(model) "Cannot extract rules for model of type $(typeof(model))."

    """
        prune_rc(rc::RuleCascade)::RuleCascade

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`RuleCascade`](@ref),
    [`rule_metrics`](@ref).
    """
    function prune_rc(rc::RuleCascade)
        E_zero = rule_metrics(SoleModels.convert(Rule,rc),X,Y)[:error]
        valid_idxs = collect(1:length(rc))

        for idx in reverse(valid_idxs)
            (length(valid_idxs) < 2) && break

            # Indices to be considered to evaluate the rule
            other_idxs = intersect!(vcat(1:(idx-1),(idx+1):length(rc)),valid_idxs)
            rule = SoleModels.convert(Rule,rc[other_idxs])

            # Return error of the rule without idx-th pair
            E_minus_i = rule_metrics(rule,X,Y)[:error]

            decay_i = (E_minus_i - E_zero) / max(E_zero, pruning_s)

            if decay_i < pruning_decay_threshold
                # Remove the idx-th pair in the vector of decisions
                deleteat!(valid_idxs, idx)
                E_zero = E_minus_i
            end
        end

        return rc[valid_idxs]
    end

    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    rcset = begin
        if model isa DecisionForest
            unique([unroll_rules_cascade(tree) for tree in trees(model)])
            # TODO maybe also sort?
        else
            unroll_rules_cascade(model)
        end
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    ########################################################################################
    if prune_rules
        rcset = prune_rc.(rcset)
    end
    ########################################################################################
    # Convert set of rules cascade in a set of rules
    ########################################################################################
    ruleset = convert.(Rule,rcset)

    ########################################################################################
    # Rule selection to obtain the best rules
    ########################################################################################
    best_rules = begin
        if method_rule_selection == :CBC
            M = hcat([evaluate_rule(rule, X, Y)[:ant_sat] for rule in ruleset]...)
            best_idxs = findcorrelation(cor(M), threshold = accuracy_rule_selection)
            ruleset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    ########################################################################################
    # Construct a rule-based model from the set of best rules
    ########################################################################################

    D = deepcopy(X) # Copy of the original dataset
    R = Rule[]      # Ordered rule list
    S = [deepcopy(best_rules)..., Rule(best_guess(Y))] # Vector of rules left

    # Rules with a frequency less than min_frequency
    S = begin
        rules_support = [rule_metrics(s,X,Y)[:support] for s in S]
        idxs_undeleted = findall(rules_support .>= min_frequency)
        S[idxs_undeleted]
    end

    while true
        # Metrics update based on remaining instances
        metrics = [rule_metrics(s,D,Y) for s in S]
        println(metrics[1])
        rules_support = [metrics[i][:support] for i in eachindex(metrics)]
        rules_error = [metrics[i][:error] for i in eachindex(metrics)]
        rules_length = [metrics[i][:length] for i in eachindex(metrics)]

        # Best rule index
        idx_best = begin
            idx_best = nothing
            # First: find the rule with minimum error
            idx = findall(rules_error .== min(rules_error...))
            (length(idx) == 1) && (idx_best = idx)

            # If not one, find the rule with maximum frequency
            idx_support = findall(rules_support .== max(rules_support[idx]...))
            (length(intersect!(idx, idx_support)) == 1) && (idx_best = idx)

            # If not one, find the rule with minimum length
            idx_length = findall(rules_length .== min(rules_length[idx]...))
            (length(intersect!(idx, idx_length)) == 1) && (idx_best = idx)

            # Final case: more than one rule with minimum length
            # Randomly choose a rule
            isnothing(idx_best) && (idx_best = rand(idx))

            length(idx_best) > 1 && error("More than one best indexes")
            idx_best[1]
        end

        # Add at the end the best rule
        push!(R, S[idx_best])

        # Indices of the remaining instances
        idx_remaining = begin
            eval_result = evaluate_rule(S[idx_best], D, Y)
            sat_unsat = eval_result[:ant_sat]
            # Remain in D the rule that not satisfying the best rule'pruning_s condition
            findall(sat_unsat .== false)
        end
        D = length(idx_remaining) > 0 ? slice_dataset(D,idx_remaining) : nothing

        if idx_best == length(S)
            return DecisionList(R[1:end-1],consequent(R[end]))
        elseif isnothing(D) || nsamples(D) == 0
            return DecisionList(R,best_guess(Y))
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        S[end] = Rule(best_guess(Y[idx_remaining]))
    end

    return error("Unexpected error in intrees!")
end
