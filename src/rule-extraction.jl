using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet
using SoleFeatures: findcorrelation
using SoleModels
using SoleModels: Rule, antecedent, consequent, rule_metrics
using SoleModels: FinalModel, Branch, DecisionForest, DecisionList
using SoleModels: RuleCascade, antecedents, convert, unroll_rules_cascade
using SoleModels: majority_vote, Label
using SoleModels.ModalLogic: _slice_dataset
using Statistics: cor
# using ModalDecisionTrees: MultiFrameModalDataset

############################################################################################
# Rule extraction from random forest
############################################################################################

# Extract rules from a forest, with respect to a dataset
function extract_rules(
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
        prune_rc(rc::RuleCascade)
            -> RuleCascade

        Prune the paths in rc with error metric

    # Arguments
    - `rc::RuleCascade`: rule to prune

    # Returns
    - `RuleCascade`: rule after the prune
    """
    function prune_rc(rc::RuleCascade)
        E_zero = rule_metrics(SoleModels.convert(Rule,rc),X,Y)[:error]
        valid_idxs = collect(1:length(rc))

        for idx in reverse(valid_idxs)
            (length(valid_idxs) < 2) && break #TODO: check

            # Indices to be considered to evaluate the rule
            other_idxs = intersect!(vcat(1:(idx-1),(idx+1):length(rc)),valid_idxs)

            # Return error of the rule without idx-th pair
            E_minus_i =
                rule_metrics(SoleModels.convert(Rule,rc[other_idxs]),X,Y)[:error]

            decay_i = (E_minus_i - E_zero) / max(E_zero, pruning_s)

            if decay_i < pruning_decay_threshold
                # Remove the idx-th pair in the vector of decisions
                deleteat!(valid_idxs, idx)
                E_zero = E_minus_i
            end
        end

        rc[valid_idxs]
    end

    ########################################################################################
    # Extract rules from each tree
    ########################################################################################
    # Obtain full ruleset
    rcset = begin
        if model isa DecisionForest
            rcset = []
            for tree in trees(model)
                tree_paths = unroll_rules_cascade(tree) #list_paths(tree)
                append!(rcset, tree_paths)
            end
            unique(rcset) # TODO maybe also sort (which requires a definition of isless(formula1, formula2))
        else model
            unroll_rules_cascade(model)
        end
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    if prune_rules
        rcset = prune_rc.(rcset)
    end
    ########################################################################################

    #Vector{Union{Branch,FinalOutcome}} -> Vector{Union{Condition,FinalOutcome}} -> RuleNest -> Rule
    ruleset = convert.(Rule,rcset)

    println("Numero di regole: $(length(ruleset))")
    ########################################################################################
    # Obtain the best rules
    best_rules = begin
        if method_rule_selection == :CBC

            M = hcat([evaluate_antecedent(rule, X) for rule in ruleset]...)

            # correlation() -> function in SoleFeatures
            #best_idxs = 1:5
            best_idxs = findcorrelation(cor(M), threshold = accuracy_rule_selection)
            #M = M[:, best_idxs]
            ruleset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    ########################################################################################
    # Construct a rule-based model from the set of best rules

    D = deepcopy(X) # Copy of the original dataset
    # Ordered rule list
    R = Rule[]
    # Vector of rules left
    #S = deepcopy(best_rules)
    #push!(S,Rule(majority_vote(Y)))
    S = [deepcopy(best_rules)..., Rule(majority_vote(Y))]

    # Rules with a frequency less than min_frequency
    S = begin
        metrics = [rule_metrics(s,X,Y) for s in S]
        rules_support = [metrics[i][:support] for i in eachindex(metrics)]
        idxs_undeleted = findall(rules_support .>= min_frequency) # Undeleted rule indexes
        S[idxs_undeleted]
    end

    while true
        # Metrics update based on remaining instances
        metrics = [rule_metrics(s,D,Y) for s in S]
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
        D = _slice_dataset(D,idx_remaining)

        if idx_best == length(S)
            return DecisionList(R[1:end-1],consequent(R[end]))
        elseif nsamples(D) == 0
            return DecisionList(R,majority_vote(Y))
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        S[end] = Rule(majority_vote(Y[idx_remaining]))
    end

    return error("Unexpected error in extract_rules!")
end
