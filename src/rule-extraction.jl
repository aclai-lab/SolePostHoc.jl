module RuleExtraction

export intrees

using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet, nconjuncts, LeftmostConjunctiveForm
using SoleFeatures: findcorrelation
using SoleModels
using SoleModels: AbstractModel
using SoleModels: Rule, antecedent, consequent, rulemetrics
using SoleModels: LeafModel, Branch, DecisionForest, DecisionList
using SoleModels: listrules, LogicalTruthCondition
using SoleModels: bestguess, Label, evaluaterule
using SoleData: slicedataset
using Statistics: cor

############################################################################################
# Rule extraction from random forest
############################################################################################

"""
    intrees(args...; kwargs...)::DecisionList

Extract rules from a model, reduces the length of each rule (number of variable value
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

TODO cite paper and specify that the method is for forests, but was extended to work with any other model

See also
[`AbstractModel`](@ref),
[`DecisionForest`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
# Extract rules from a forest, with respect to a dataset
function intrees(
    model::Union{AbstractModel,DecisionForest},
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
        prune_rule(rc::Rule)::Rule

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`Rule`](@ref),
    [`rulemetrics`](@ref).
    """
    function prune_rule(
        r::Rule{O,<:LogicalTruthCondition{<:LeftmostConjunctiveForm}} # TODO add TrueCondition and nchildren
    ) where {O}
        E_zero = rulemetrics(r,X,Y)[:error]
        valid_idxs = collect(1:nconjuncts(r))

        for idx in reverse(valid_idxs)
            (length(valid_idxs) < 2) && break

            # Indices to be considered to evaluate the rule
            other_idxs = intersect!(vcat(1:(idx-1),(idx+1):nconjuncts(r)),valid_idxs)
            rule = r[other_idxs]

            # Return error of the rule without idx-th pair
            E_minus_i = rulemetrics(rule,X,Y)[:error]

            decay_i = (E_minus_i - E_zero) / max(E_zero, pruning_s)

            if decay_i < pruning_decay_threshold
                # Remove the idx-th pair in the vector of decisions
                deleteat!(valid_idxs, idx)
                E_zero = E_minus_i
            end
        end

        return r[valid_idxs]
    end

    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    ruleset = begin
        if model isa DecisionForest
            unique([listrules(tree) for tree in trees(model)])
            # TODO maybe also sort?
        else
            listrules(model)
        end
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    ########################################################################################
    if prune_rules
        ruleset = prune_rule.(ruleset)
    end

    ########################################################################################
    # Rule selection to obtain the best rules
    ########################################################################################
    best_rules = begin
        if method_rule_selection == :CBC
            M = hcat([evaluaterule(rule, X, Y)[:antsat] for rule in ruleset]...)
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
    S = [deepcopy(best_rules)..., Rule(bestguess(Y))] # Vector of rules left

    # Rules with a frequency less than min_frequency
    S = begin
        rules_support = [rulemetrics(s,X,Y)[:support] for s in S]
        idxs_undeleted = findall(rules_support .>= min_frequency)
        S[idxs_undeleted]
    end

    while true
        # Metrics update based on remaining instances
        metrics = [rulemetrics(s,D,Y) for s in S]
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
            sat_unsat = evaluaterule(S[idx_best], D, Y)[:antsat]
            # Remain in D the rule that not satisfying the best rule'pruning_s condition
            findall(sat_unsat .== false)
        end
        D = length(idx_remaining) > 0 ? slicedataset(D,idx_remaining) : nothing

        if idx_best == length(S)
            return DecisionList(R[1:end-1],consequent(R[end]))
        elseif isnothing(D) || ninstances(D) == 0
            return DecisionList(R,bestguess(Y))
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        S[end] = Rule(bestguess(Y[idx_remaining]))
    end

    return error("Unexpected error in intrees!")
end

end
