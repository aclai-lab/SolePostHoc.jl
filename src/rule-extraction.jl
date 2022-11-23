############################################################################################
# Rule extraction from random forest
############################################################################################

using SoleLearning: Consequent,
            Rule, antecedent, consequent, rule_metrics
            DecisionList, list_rules

# Patch single-frame _-> multi-frame
extract_rules(model::Any, X::ModalDataset, args...; kwargs...) =
    extract_rules(model, MultiFrameModalDataset(X), args...; kwargs...)

# Extract rules from a forest, with respect to a dataset
# TODO avoid Formula{L}()
function extract_rules(
        forest::DecisionForest,
        X::MultiFrameModalDataset,
        Y::AbstractVector{<:Consequent};
        prune_rules = false,
        s = nothing,
        decay_threshold = nothing,
        #
        method = :CBC,
        min_frequency = nothing,
)

    isnothing(s) && (s = 1.0e-6)
    isnothing(decay_threshold) && (decay_threshold = 0.05)
    isnothing(min_frequency) && (min_frequency = 0.01)

    """
        prune_ruleset(ruleset::AbstractVector{<:Rule}) -> DecisionList

        Prune the rules in a ruleset according to the error metric

    If `s` and `decay_threshold` is unspecified, their values are set to nothing and the
    first two rows of the function set s and decay_threshold with their default values

    # Arguments
    - `ruleset::AbstractVector{<:Rule}`: rules to prune

    # Returns
    - `DecisionList`: rules after the prune
    """
    function prune_ruleset(
        ruleset::AbstractVector{<:Rule}
    )
        [begin
            # Extract decisions from rule
            decs = extract_decisions(antecedent(rule))
            cons = consequent(rule)

            E_zero = rule_metrics(decs, cons, X, Y)[:error]

            for idx in length(decs):1
                # Indices to be considered to evaluate the rule
                other_idxs = vcat(1:(idx-1), (idx+1):length(decs))
                # Return error of the rule without idx-th pair
                E_minus_i = rule_metrics(decs[other_idxs], cons, X, Y)[:error]
                decay_i = (E_minus_i - E_zero) / max(E_zero, s)
                if decay_i < decay_threshold
                    # Remove the idx-th pair in the vector of decisions
                    deleteat!(decs, idx)
                    E_zero = rule_metrics(decs, cons, X, Y)[:error]
                end
            end
            # Assemble formula from vector of decisions (decs)
            #TODO: formula_update(antecedent(rule),nodes_deleted)
            antecedent = TODO
            #TODO check if this works:
            Rule(antecedent, cons)
            # Rule{typeof(antecedent),typeof(cons)}(antecedent, cons)
        end for rule in ruleset]
    end

    ########################################################################################
    # Extract rules from each tree
    ########################################################################################
    # Obtain full ruleset
    ruleset = begin
        ruleset = []
        for every tree in the forest
            tree_rules = list_rules(tree) # TODO implement
            append!(ruleset, tree_rules)
        end
        unique(ruleset) # TODO maybe also sort (which requires a definition of isless(formula1, formula2))
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    if prune_rules
        ruleset = prune_ruleset(ruleset)
    end
    ########################################################################################

    ########################################################################################
    # Obtain the best rules
    best_rules = begin
        if method == :CBC
            # Extract antecedents
            antset = antecedent.(ruleset)
            # Build the binary satisfuction matrix (m × j, with m instances and j antecedents)
            M = hcat([evaluate_antecedent(antecedent(rule), X) for rule in antset]...)
            # correlation() -> function in SoleFeatures
            best_idxs = findcorrelation(M)
            M = M[:, best_idxs]
            ruleset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    ########################################################################################
    # Construct a rule-based model from the set of best rules

    D = copy(X) # Copy of the original dataset
    # TODO @Michele: R and S should be Vector{Rule}; only at the end, you return DecisionList(R)
    R = DecisionList()  # Vector of ordered rules
    S = DecisionList()  # Vector of rules left
    append!(S, best_rules)
    #TODO: remove Formula{L}()
    push!(S, Rule{L}(Formula{L}(), majority_vote(Y)))

    # Rules with a frequency less than min_frequency
    S = begin
        metrics = rule_metrics.(S, X, Y)
        rules_support = [metrics[i][:support] for i in eachindex(metrics)]
        idxs_undeleted = findall(rules_support .>= min_frequency) # Undeleted rule indexes
        S[idxs_undeleted]
    end

    while true
        # Metrics update based on remaining instances
        metrics = rule_metrics.(S, D, Y)
        rules_support = [metrics[i][:support] for i in eachindex(metrics)]
        rules_error = [metrics[i][:error] for i in eachindex(metrics)]
        rules_length = [metrics[i][:length] for i in eachindex(metrics)]

        # Best rule index
        idx_best = begin
            # First: find the rule with minimum error
            idx = findall(rules_error .== min(rules_error...))
            (length(idx) == 1) && (return idx)

            # If not one, find the rule with maximum frequency
            idx_support = findall(rules_support .== max(rules_support[idx]...))
            (length(intersect!(idx, idx_support)) == 1) && (return idx)

            # If not one, find the rule with minimum length
            idx_length = findall(rules_length .== min(rules_length[idx]...))
            (length(intersect!(idx, idx_length)) == 1) && (return idx)

            # Final case: more than one rule with minimum length
            # Randomly choose a rule
            rand(idx)
        end

        # Add at the end the best rule
        push!(R, S[idx_best])

        # Indices of the remaining instances
        idx_remaining = begin
            eval_result = evaluate_rule(S[idx_best], D, Y)
            sat_unsat = eval_result[:ant_sat]
            # Remain in D the rule that not satisfying the best rule's condition
            findall(sat_unsat .== false)
        end
        D = D[idx_remaining,:]

        if idx_best == length(S)
            return R
        elseif size(D, 1) == 0
            # TODO: remove Formula{L}()
            push!(R, Rule{L}(Formula{L}(),majority_vote(Y)))
            return R
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        # TODO: remove Formula{L}()
        # TODO fix S[length(S)] into into S[end] maybe?
        S[length(S)] = Rule{L}(Formula{L}(), majority_vote(Y[idx_remaining]))
    end
    ########################################################################################
end
