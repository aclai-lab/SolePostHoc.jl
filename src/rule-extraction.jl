using SoleLogics
using SoleModels: Consequent,
    Rule, antecedent, consequent, rule_metrics,
    Branch,
    AbstractDecisionTree, DecisionTreeNode, convert_path_rule,
    DecisionList, list_paths

abstract type AbstractDataset end

############################################################################################
# Convert path to rule
############################################################################################
# TODO move to SoleModels
# Convert path to rule
# TODO: Fix the function
# TODO: Rename into convert(::Rule, Type{RuleNest}) ...
"""
    Convert a rule cascade in a rule
"""
convert(::Rule,rule_cascade::RuleCascade) =
    convert(::Rule, antecedent::AbstractVector{<:Formula}, consequent::AbstractModel)

function convert(::Rule,antecedent::AbstractVector{<:Formula}, consequent::AbstractModel)
    # Building antecedent
    function _build_formula(conjuncts::AbstractVector{<:Formula{L}}) where {L<:Logic}
        if length(conjuncts) > 2
            SoleLogics.CONJUNCTION(conjuncts[1],_build_formula(conjuncts[2:end]))
        else
            SoleLogics.CONJUNCTION(conjuncts[1],conjuncts[2])
        end
    end

    # Antecedent of the rule
    ant = begin

        root = begin
            # Number of internal nodes in the antecedent
            n_internal = length(antecedent)

            if n_internal == 0
                FNode(SoleLogics.TOP)
            elseif n_internal == 1
                tree(antecedent[1])
            else
                _build_formula(antecedent)
            end
        end

        Formula(root)
    end

    Rule(ant, consequent)
end

############################################################################################
# List paths that represent rules of DecisionTree
############################################################################################
# TODO: Move to SoleLoearning/SoleModels
# TODO: fix the prediction, each model have a specific field for the prediction

"""
    List all paths of a decision tree by performing a tree traversal
"""

function list_paths(tree::DecisionTree)
    # tree(f) [where f is a Formula object] is used to
    # retrieve the root FNode of the formula(syntax) tree
    pathset = list_paths(root(tree))

    (length(pathset) == 1) && (return [RuleCascade(SoleLogics.TOP,pathset[1])])

    return [RuleCascade(path[1:end-1],path[end]) for path in pathset]
end

function list_paths(node::Branch)
    # NOTE: antecedent(node) or tree(antecedent(node)) to obtain a FNode?
    positive_path  = [antecedent(node)]
    negative_path = [NEGATION(antecedent(node))]
    return [
        list_paths(positive_consequent(node),  positive_path)...,
        list_paths(negative_consequent(node), negative_path)...,
    ]
end

function list_paths(node::AbstractModel)
    return [node]
end

function list_paths(node::Branch, this_path::AbstractVector)
    # NOTE: antecedent(node) or tree(antecedent(node)) to obtain a FNode?
    positive_path  = [this_path..., antecedent(node)]
    negative_path = [this_path..., NEGATION(antecedent(node))]
    return [
        list_paths(positive_consequent(node),  positive_path)...,
        list_paths(negative_consequent(node), negative_path)...,
    ]
end

function list_paths(node::AbstractModel,this_path::AbstractVector)
    return [[this_path..., node], ]
end

############################################################################################
# Rule evaluation
############################################################################################

# Evaluation for an antecedent
function evaluate_antecedent(rule::Rule, X::AbstractDataset)
    evaluate_antecedent(antecedent(rule), X)
end

function evaluate_antecedent(antecedent::Formula{L}, X::AbstractDataset)
    check(antecedent, X)
end

function evaluate_rule(
    rule::Rule,
    X::AbstractDataset,
    Y::AbstractVector{<:Consequent}
)
    # Antecedent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable.
    ant_sat = evaluate_antecedent(antecedent(rule),X)

    # Indices of satisfiable instances
    idxs_sat = findall(ant_sat .== true)

    # Consequent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable,
    #  - `nothing` when antecedent does not hold.
    cons_sat = begin
        cons_sat = Vector{Union{Bool, Nothing}}(fill(nothing, length(Y)))
        idxs_true = begin
            idx_cons = findall(consequent(rule) .== Y)
            intersect(idxs_sat,idx_cons)
        end
        idxs_false = begin
            idx_cons = findall(consequent(rule) .!= Y)
            intersect(idxs_sat,idx_cons)
        end
        cons_sat[idxs_true]  .= true
        cons_sat[idxs_false] .= false
        cons_sat
    end

    y_pred = begin
        y_pred = Vector{Union{Consequent, Nothing}}(fill(nothing, length(Y)))
        y_pred[idxs_sat] .= consequent(rule)
        y_pred
    end

    return (;
        ant_sat   = ant_sat,
        idxs_sat  = idxs_sat,
        cons_sat  = cons_sat,
        y_pred    = y_pred,
    )
end

############################################################################################
# Rule length
############################################################################################

rule_length(rule::Rule) = rule_length(antecedent(rule))
rule_length(formula::Formula) = rule_length(root(formula))

function rule_length(node::FNode)
    if isleaf(node)
        return 1
    else
        return ((isdefined(node,:leftchild) ? rule_length(leftchild(node)) : 0)
                    + (isdefined(node,:rightchild) ? rule_length(rightchild(node)) : 0))
    end
end

############################################################################################
# Rule metrics
############################################################################################

function rule_metrics(
    rule::Rule,
    X::MultiFrameModalDataset,
    Y::AbstractVector{<:Consequent}
)

    eval_result = evaluate_rule(rule, X, Y)
    n_instances = size(X, 1)
    n_satisfy = sum(eval_result[:ant_sat])

    # Support of the rule
    rule_support =  n_satisfy / n_instances

    # Error of the rule
    rule_error = begin
        if outcome_type(consequent(rule)) <: CLabel
            # Number of incorrectly classified instances divided by number of instances
            # satisfying the rule condition.
            misclassified_instances = length(findall(eval_result[:y_pred] .== Y))
            misclassified_instances / n_satisfy
        elseif outcome_type(consequent(rule)) <: RLabel
            # Mean Squared Error (mse)
            idxs_sat = eval_result[:idxs_sat]
            mse(eval_result[:y_pred][idxs_sat], Y[idxs_sat])
        end
    end

    return (;
        support   = rule_support,
        error     = rule_error,
        length    = rule_length(rule),
    )
end

############################################################################################
# Rule extraction from random forest
############################################################################################

# Patch single-frame _-> multi-frame
extract_rules(model::Any, X::ModalDataset, args...; kwargs...) =
    extract_rules(model, AbstractDataset(X), args...; kwargs...)

# Extract rules from a forest, with respect to a dataset
# TODO: SoleLogics.True
function extract_rules(
        forest::DecisionForest,
        X::AbstractDataset,
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
        prune_pathset(pathset::AbstractVector{<:AbstractVector{<:DecisionTreeNode}})
            -> AbstractVector{<:AbstractVector{<:DecisionTreeNode}}

        Prune the paths in pathset with error metric

    # Arguments
    - `pathset::AbstractVector{<:AbstractVector{<:DecisionTreeNode}}`: paths to prune

    # Returns
    - `AbstractVector{<:AbstractVector{<:DecisionTreeNode}}`: paths after the prune
    """
    function prune_pathset(pathset::AbstractVector{<:RuleCascade})
        [begin
            ant = antecedent(path)
            cons = consequents(path)

            E_zero = rule_metrics(convert(Rule,ant,cons), X, Y)[:error]
            valid_idxs = collect(length(ant):-1:1)

            for idx in reverse(valid_idxs)
                # Indices to be considered to evaluate the rule
                other_idxs = intersect!(vcat(1:(idx-1),(idx+1):length(ant)),valid_idxs)

                # Return error of the rule without idx-th pair
                E_minus_i = rule_metrics(convert(Rule,ant[other_idxs],cons), X, Y)[:error]

                decay_i = (E_minus_i - E_zero) / max(E_zero, s)

                if decay_i < decay_threshold
                    # Remove the idx-th pair in the vector of decisions
                    deleteat!(valid_idxs, idx)
                    E_zero = E_minus_i
                end
            end

            RuleCascade(ant[valid_idxs],cons)
        end for path in pathset]
    end

    ########################################################################################
    # Extract rules from each tree
    ########################################################################################
    # Obtain full ruleset
    pathset = begin
        pathset = []
        for tree in forest
            tree_paths = list_paths(tree)
            append!(pathset, tree_paths)
        end
        unique(pathset) # TODO maybe also sort (which requires a definition of isless(formula1, formula2))
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    if prune_rules
        pathset = prune_pathset(pathset)
    end
    ########################################################################################

    #Vector{Union{Branch,FinalOutcome}} -> Vector{Union{Condition,FinalOutcome}} -> RuleNest -> Rule
    ruleset = convert.(Rule,pathset)

    ########################################################################################
    # Obtain the best rules
    best_rules = begin
        if method == :CBC

            M = hcat([evaluate_antecedent(rule, X) for rule in ruleset]...)

            # correlation() -> function in SoleFeatures
            best_idxs = findcorrelation(M)
            #M = M[:, best_idxs]
            ruleset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    ########################################################################################
    # Construct a rule-based model from the set of best rules

    #TODO: fix majority_vote

    D = copy(X) # Copy of the original dataset
    # Ordered rule list
    R = Rule[]
    # Vector of rules left
    S = copy(best_rules)
    #TODO: SoleLogics.TOP
    #TODO: Fix Default Rule
    push!(S,Rule(Formula(FNode(SoleLogics.TOP)),majority_vote(Y)))

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
            #TODO: fix default field; majority_vote(Y)?
            return DecisionList(R[end-1],consequent(R[end]))
        elseif size(D, 1) == 0
            return DecisionList(R,majority_vote(Y))
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        S[end] = Rule(Formula(FNode(SoleLogics.TOP)),majority_vote(Y[idx_remaining]))
    end

    return error("Unexpected error in extract_rules!")
end
