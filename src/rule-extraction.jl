############################################################################################

const Consequent = Any

struct Rule{L<:Logic,C<:Consequent}
    antecedent :: Formula{L}
    consequent :: C
end

antecedent(rule::Rule) = rule.antecedent
consequent(rule::Rule) = rule.consequent

const ClassificationRule = Rule{L,CLabel} where {L}
const RegressionRule     = Rule{L,RLabel} where {L}

struct RuleBasedModel{L<:Logic,C<:Consequent}
    rules :: Vector{<:Rule{L,C}}
end

rules(model::RuleBasedModel) = model.rules

const RuleBasedClassifier = RuleBasedModel{L,CLabel} where {L}
const RuleBasedRegressor  = RuleBasedModel{L,RLabel} where {L}

"""
TODO document
"""
function list_rules(tree::DTree)
    list_rules(tree.root) # TODO what about world_types and init_conditions?
end

function list_rules(node::DTInternal)
    pos_lambda = get_lambda(node.decision)
    neg_lambda = get_inverse_lambda(node.decision)
    return [
        [advance_rule(rule, pos_lambda) for rule in list_rules(node.left)]...,
        [advance_rule(rule, neg_lambda) for rule in list_rules(node.right)]...,
    ]
end

function list_rules(leaf::DTLeaf{L})
    Rule{TODO, L}(TODO: SoleLogics.True, prediction(leaf))
end

function list_rules(leaf::NSDTLeaf{L})
    Rule{TODO, PredictingFunction{L}}(TODO: SoleLogics.True, predicting_function(leaf))
end


function get_lambda(decision::Decision)
    Formula( # TODO formula of any logic, or formula of a specific logic?
        if is_propositional_decision(decision)
            return PropositionalDimensionalLetter( # TODO
                decision.feature,
                decision.test_operator,
                decision.threshold
            )
        else
            ModalExistentialOperator{decision.relation}(
                PropositionalDimensionalLetter( # TODO
                    decision.feature,
                    decision.test_operator,
                    decision.threshold
                )
            )
        end
    )
end
function get_inverse_lambda(decision::Decision)
    Formula( # TODO formula of any logic, or formula of a specific logic?
        if is_propositional_decision(decision)
            return Operator{negation...}(
                PropositionalDimensionalLetter( # TODO
                    decision.feature,
                    decision.test_operator,
                    decision.threshold
                )
            )
        else
            ModalUniversalOperator{decision.relation}(
                Operator{negation...}(
                    PropositionalDimensionalLetter( # TODO
                        decision.feature,
                        decision.test_operator,
                        decision.threshold
                    )
                )
            )
        end
    )
end

function advance_rule(rule::Rule{L,C}, lambda::Formula{L}) where {L,C}
    Rule{L,C}(advance_rule(antecedent(rule), lambda), consequent(rule))
end

function advance_rule(rule_antecedent::Formula{L}, lambda::Formula{L}) where {L}
    conjuncts = begin
        Formula{L}[TODO derive conjuncts from rule_antecedent...]
        una sorta di `conjuncts = split(formula, :(\wedge))`
    end
    Formula{L}(tipo join(advance_conjunct.(conjuncts, lambda), :(\wedge)))
end

function advance_conjunct(conjunct::Formula{L}, lambda::Formula{L}) where {L}
    # TODO
    is_positive(conjunct)
    is_left(lambda)
end

function is_positive(conjunct::Formula{L}) where {L}
    # TODO
end

function is_left(lambda::Formula{L}) where {L}
    # TODO
end

# Evaluation for single decision
# TODO
function evaluate_decision(dec::Decision, X::MultiFrameModalDataset) end

# TODO remove
logic = extract_logic

############################################################################################
############################################################################################
############################################################################################

############################################################################################
# path --> Rule
############################################################################################

# Convert path to rule
function path2rule(path::AbstractVector{<:AbstractVector{<:DTNode}})
    function _union_family(p_node::AbstractOperator,l_node::DTNode,r_node::Union{Fnode,DTNode})
        p_node = FNode(p_node)
        l_node = FNode(decision(l_node))
        r_node = begin
            if l_node <: DTNode
                return Fnode(decision(r_node))
            end
            l_node
        end

        leftchild!(p_node,l_node)
        rightchild!(p_node,r_node)
        parent!(l_node,p_node)
        parent!(r_node,p_node)

        return p_node
    end

    # Building antecedent
    function _build_formula(conjuncts::AbstractVector{<:DTNode})
        if length(conjuncts) > 2
            _union_family(SoleLogics.CONJUNCTION,conjuncts[1],_build_formula(conjuncts[2:end]))
        else
            _union_family(SoleLogics.CONJUNCTION,conjuncts[1],conjuncts[2])
        end
    end

    # Antecedent of the rule
    ant = begin
        root = begin
            # Number of internal nodes in the path
            n_internal = length(path) - 1

            (n_internal == 0) && (return FNode(SoleLogics.True))
            (n_internal == 1) && (return FNode(decision(path[1])))
            _build_formula(path[1:(end-1)])
        end

        Formula(root)
    end

    # Consequent of the rule
    cons = prediction(path[end])

    Rule(ant,cons)
end

############################################################################################
# Rule evaluation
############################################################################################

# Evaluation for an antecedent
function evaluate_antecedent(decs::AbstractVector{<:Decision}, X::MultiFrameModalDataset)
    D = hcat([evaluate_decision(d, X) for d in decs]...)
    # If all values in a row is true, then true (and logical)
    return map(all, eachrow(D))
end

# Evaluation for a rule
function evaluate_rule(
    path::AbstractVector{<:DTNode},
    X::MultiFrameModalDataset,
    Y::AbstractVector{<:Consequent}
)
    # Antecedent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable.
    ant_sat = evaluate_antecedent(decision.(path[1:(end-1)]),X)

    # Indices of satisfiable instances
    idxs_sat = findall(ant_sat .== true)

    # Consequent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable,
    #  - `nothing` when antecedent does not hold.
    cons_sat = begin
        cons_sat = Vector{Union{Bool, Nothing}}(fill(nothing, length(Y)))
        idxs_true = begin
            idx_cons = findall(prediction(path[end]) .== Y)
            intersect(idxs_sat,idx_cons)
        end
        idxs_false = begin
            idx_cons = findall(prediction(path[end]) .!= Y)
            intersect(idxs_sat,idx_cons)
        end
        cons_sat[idxs_true]  .= true
        cons_sat[idxs_false] .= false
        cons_sat
    end

    y_pred = begin
        y_pred = Vector{Union{Consequent, Nothing}}(fill(nothing, length(Y)))
        y_pred[idxs_sat] .= prediction(path[end])
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
# Rule extraction from random forest
############################################################################################

# Patch single-frame _-> multi-frame
extract_rules(model::Any, X::ModalDataset, args...; kwargs...) =
    extract_rules(model, MultiFrameModalDataset(X), args...; kwargs...)

# Extract rules from a forest, with respect to a dataset
# TODO: SoleLogics.True
function extract_rules(
        forest::DForest,
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
        rule_metrics(args...) -> AbstractVector

        Compute frequency, error and length of the rule

    # Arguments
    - `path::AbstractVector{<:DTNode}`: path of the rule
    - `X::MultiFrameModalDataset`: dataset
    - `Y::AbstractVector{<:Consequent}`: target values of X

    # Returns
    - `AbstractVector`: metrics values vector of the rule
    """
    function rule_metrics(
        path::AbstractVector{<:DTNode},
        X::MultiFrameModalDataset,
        Y::AbstractVector{<:Consequent}
    )
        eval_result = evaluate_rule(path, X, Y)
        n_instances = size(X, 1)
        n_satisfy = sum(eval_result[:ant_sat])

        # Support of the rule
        rule_support =  n_satisfy / n_instances

        # Error of the rule
        rule_error = begin
            cons = prediction(path[end])
            if typeof(cons) <: CLabel
                # Number of incorrectly classified instances divided by number of instances
                # satisfying the rule condition.
                misclassified_instances = length(findall(eval_result[:y_pred] .!= Y))
                misclassified_instances / n_satisfy
            elseif typeof(cons) <: RLabel
                # Mean Squared Error (mse)
                idxs_sat = eval_result[:idxs_sat]
                mse(eval_result[:y_pred][idxs_sat], Y[idxs_sat])
            end
        end

        # Length of the rule
        rule_length = length(path)

        return (;
            support   = rule_support,
            error     = rule_error,
            length    = rule_length,
        )
    end

    """
        prune_pathset(pathset::AbstractVector{<:AbstractVector{<:DTNode}})
            -> AbstractVector{<:AbstractVector{<:DTNode}}

        Prune the paths in pathset with error metric

    # Arguments
    - `pathset::AbstractVector{<:AbstractVector{<:DTNode}}`: paths to prune

    # Returns
    - `AbstractVector{<:AbstractVector{<:DTNode}}`: paths after the prune
    """
    function prune_pathset(
        pathset::AbstractVector{<:AbstractVector{<:DTNode}}
    )
        [begin
            E_zero = rule_metrics(path, X, Y)[:error]

            for idx in (length(path)-1):1
                # Indices to be considered to evaluate the rule
                other_idxs = vcat(1:(idx-1), (idx+1):length(path))
                # Return error of the rule without idx-th pair
                E_minus_i = rule_metrics(path[other_idxs], X, Y)[:error]
                decay_i = (E_minus_i - E_zero) / max(E_zero, s)
                if decay_i < decay_threshold
                    # Remove the idx-th pair in the vector of decisions
                    deleteat!(path, idx)
                    E_zero = E_minus_i #rule_metrics(path, X, Y)[:error]
                end
            end
        end for path in pathset]
    end

    ########################################################################################
    # Extract rules from each tree
    ########################################################################################
    # Obtain full ruleset
    pathset = begin
        pathset = []
        for tree in forest
            tree_paths = list_paths(tree) # TODO implement
            append!(pathset, tree_paths)
        end
        unique(pathset) # TODO maybe also sort (which requires a definition of isless(formula1, formula2))
    end
    ########################################################################################

    ########################################################################################
    # Prune rules according to the error metric (with respect to a dataset)
    #  (and similar metrics: support, error, and length)
    if prune_rules
        pathset = prune_pathset(pathset)
    end
    ########################################################################################

    ########################################################################################
    # Obtain the best rules
    best_rules = begin
        if method == :CBC
            # Extract antecedents
            antset = [decision.(path) for path in pathset]
            # Build the binary satisfuction matrix (m × j, with m instances and j antecedents)
            M = hcat([evaluate_antecedent(ant, X) for ant in antset]...)
            # correlation() -> function in SoleFeatures
            best_idxs = findcorrelation(M)
            #M = M[:, best_idxs]
            pathset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    ########################################################################################
    # Construct a rule-based model from the set of best rules

    D = copy(X) # Copy of the original dataset
    # Ordered rule list
    R = RuleBasedModel([])
    # Vector of rules left
    S = copy(best_rules)
    #TODO: SoleLogics.True
    push!(S, [DTNode(SoleLogics.True), DTLeaf(majority_vote(Y))])

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
        push!(rules(R), path2rule(S[idx_best]))

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
            #TODO: SoleLogics.True
            push!(rules(R), [DTNode(SoleLogics.True), DTLeaf(majority_vote(Y))])
            return R
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        # TODO: SoleLogics.True
        S[length(S)] = [DTNode(SoleLogics.True), DTLeaf(majority_vote(Y[idx_remaining]))]
    end
    ########################################################################################
end
