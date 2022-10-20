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


############################################################################################
# MDTv1 translation
############################################################################################
function translate_mdtv1(tree:DTree)
    formula = Formula{L}(SoleLogics.TOP) # Formula wrapped in the root
    # Maybe it's more realistic to translate each node as follows
    # pure_node = deepcopy(non_pure_node)
    # formula!(pure_root) = formula
    # In this sketch of implementation, I only create nodes containing a formula.
    # This is why the variable `tree` here is never used.
    pure_root = DTInternal(formula)
    pure_tree = DTree(pure_root)

    translate_mdtv1(pure_root, leftchild(tree))
    translate_mdtv1(pure_root, rightchild(tree))

    return pure_tree
end

function translate_mdtv1(prev_node::DTInternal, node::Union{DTInternal, DTLeaf, NSFTLeaf})
    # `node` comes from a non-pure DTree, while `prev_node` is part of a pure DTree

    # formula in `node` is modified, depending on its parent formula
    # accordingly to ModalDecisionTree paper (pg. 17)
    formula  = compose_pureformula(node, prev_node)
    new_node = DTInternal(formula)  # new pure-tree node

    # Linking
    # if node was a left child
    if is_lchild(node)
        # then the new pure-tree node has to be a left child too
        leftchild!(prev_node, new_node)
    else
        # like above
        rightchild!(prev_node, new_node)
    end
    # the pure-tree has now, officialy, the `new_node`
    parent!(new_node, prev_node)

    # if `node` is a DTLeaf or NSFTLeaf, stop the recursion
    # else continue linking some nodes as `new_node` descendents
    if !(typeof(node) <: Union{DTLeaf, NSFTleaf})
        translate_mdtv1(leftchild(new_node), new_node)
        translate_mdtv1(rightchild(new_node), new_node)
    end
end

function compose_pureformula(node::DTInternal, lambda::DTInternal)
    # ModalDecisionTree paper, pg.17
    # base case (h=1) and "λ ∈ S⁻" case ("p ∧ φπ" is returned)
    if height(lambda) == 1 || is_rchild(lambda)
        conj = FNode{L}(CONJUNCTION)
        p    = FNode{L}(formula(lambda))
        φ    = FNode{L}(formula(node))
        link_nodes(conj, l, r)
        return Formula{L}(conj)
    else
        # other cases dispatch
        return compose_pureformula(
            node::DTInternal,
            lambda::DTInternal,
            Val(agreement(node, parent(node)))
        )
    end
end

# first and second case
function compose_pureformula(node::DTInternal, lambda::DTInternal, ::Val{true})
    conj = FNode{L}(CONJUNCTION)
    p    = FNode{L}(formula(lambda))
    φ    = FNode{L}(formula(node))
    link_nodes(conj, p, φ)

    if is_propositional_decision(lambda)
        # first case (p ∧ φπ)
        return Formula{L}(conj)
    else
        # second case (⟨X⟩(p ∧ φπ))
        diamond = FNode{L}(EXMODOP(relation(lambda)))
        link_nodes(diamond, conj)
        return Formula{L}(diamond)
    end
end

# third and fourth case
function compose_pureformula(node::DTInternal, lambda::DTInternal, ::Val{false})
    p = FNode{L}(formula(lambda))
    φ = FNode{L}(formula(node))
    impl = FNode{L}(IMPLICATION)
    conj = FNode{L}(CONJUNCTION)

    # NOTE: p2 is identical to p, but we need double the memory
    # to correctly represent every link in the tree.
    # See if-else statement to better understand p2 role.
    p2 = FNode{L}(formula(lambda))

    # In both cases the fragment p => φπ is common
    link_nodes(impl, p, φ)

    if is_propositional_decision(lambda)
        # third case p2 ∧ (p => φπ)
        link_nodes(conj, p2, impl)
    else
        # fourth case ⟨X⟩p2 ∧ [X](p => φπ)
        diamond = FNode{L}(EXMODOP(relation(lambda)))
        box = FNode{L}(UNIVMODOP(relation(lambda)))

        link_nodes(box, impl)
        link_nodes(diamond, p2)
        link_nodes(conj, diamond, box)
    end

    # `conj` is the formula(syntax)-tree root
    return Formula{L}(conj)
end

# utility to link one parent to exactly one children
function link_nodes(p::FNode{L}, r::FNode{L}) where {L <: AbstractLogic}
    rightchild!(p, r)
    parent!(r, p)
end

# utility to link one parent to its left and right childrens
function link_nodes(p::FNode{L}, l::FNode{L}, r::FNode{L}) where {L <: AbstractLogic}
    leftchild!(p, l)
    rightchild!(p, r)
    parent!(l, p)
    parent!(r, p)
end


############################################################################################
# List rules
############################################################################################

# utility to list all the rules composed by traversing a pure-tree
# from root to its leaves.
function list_rules(tree::DTree)
    # tree(f) [where f is a Formula object] is used to
    # retrieve the root FNode of the formula(syntax) tree
    return [
        list_rules(leftchild(root(tree)), tree(formula(root(tree)))),
        list_rules(rightchild(root(tree)), tree(formula(root(tree))))
    ]
end

function list_rules(node::DTInternal, conjuncts::FNode{L}) where {L<:AbstractLogic}
    # (conjuncts retrieved until now) ∧ (formula in current node)
    conj = FNode{L}(CONJUNCTION)
    link_nodes(conj, conjuncts, tree(formula(node)))
    return [list_rules(leftchild(node), conj), list_rules(rightchild(node), conj)]
end

function list_rules(
    node::Union{DTLeaf, NSFTLeaf},
    conjuncts::FNode{L}
    ) where {L<:AbstractLogic}
    conj = Fnode{L}(CONJUNCTION)
    link_nodes(conj, conjuncts, tree(formula(node)))
    formula = Formula{L}(conj)
    return Rule{L,C}(formula, prediction(node))
end

"""
Previous idea to extract rules (code commented by Mauro, gio 20 oct)

TODO document

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
"""

# Evaluation for single decision
# TODO
function evaluate_decision(dec::Decision, X::MultiFrameModalDataset) end

# TODO remove
logic = extract_logic

############################################################################################
############################################################################################
############################################################################################

# Extract decisions from rule
function extract_decisions(formula::Formula{L}) where {L<:Logic}
    # TODO remove in favor of operators_set = operators(L)
    operators_set = operators(logic(antecedent(rule)))
    function _extract_decisions(node::FNode, decs::AbstractVector{<:Decision})
        # Leaf or internal node
        if !isdefined(node, :leftchild) && !isdefined(node, :rightchild)
            if token(node) in operators_set
                return decs
            else
                return push!(decs, token(node))
            end
        else
            isdefined(node, :leftchild)  && _extract_decisions(leftchild(node),  decs)
            isdefined(node, :rightchild) && _extract_decisions(rightchild(node), decs)

            if !(token(node) in operators_set)
                return push!(decs, token(node))
            end
            decs
        end
    end
    _extract_decisions(tree(formula), [])
end

############################################################################################
# Formula Update
############################################################################################

function formula_update(formula::Formula{L},nodes_deleted::AbstractVector)
    root = tree(formula)

    function _transplant(u::FNode,v::FNode)
        #u è radice
        u == root ? root = v : nothing

        #u è figlio sx
        u == leftchild(parent(u)) ? leftchild!(parent(u),v) : nothing

        #u è figlio dx
        u == rightchild(parent(u)) ? rightchild!(parent(u),v) : nothing

        #v definito
        isdefined(v,:token) ? parent!(v,parent(u)) : nothing

        return nothing
    end

    function _formula_update(node::FNode,node_deleted::FNode)

        #è il nodo da eliminare
        if node == node_deleted
            if leftchild(parent(node)) == node
                return _transplant(parent(node),rightchild(parent(node)))
            else
                return _transplant(parent(node),leftchild(parent(node)))
            end
        end

        #non è il nodo da eliminare

        #se non sono in una foglia, passo ai rami
        isdefined(node, :leftchild)  && _formula_update(leftchild(node), node_deleted)
        isdefined(node, :rightchild) && _formula_update(rightchild(node), node_deleted)

        return nothing
    end

    for node in nodes_deleted
        _formula_update(root,node)
    end

    return Formula{L}(root)
end

############################################################################################
# Rule evaluation
############################################################################################

# Evaluation for an antecedent

evaluate_antecedent(antecedent::Formula{L}, X::MultiFrameModalDataset) where {L<:Logic} =
    evaluate_antecedent(extract_decisions(antecedent), X)

function evaluate_antecedent(decs::AbstractVector{<:Decision}, X::MultiFrameModalDataset)
    D = hcat([evaluate_decision(d, X) for d in decs]...)
    # If all values in a row is true, then true (and logical)
    return map(all, eachrow(D))
end

# Evaluation for a rule

# From rule to antecedent and consequent
evaluate_rule(rule::Rule, X::MultiFrameModalDataset, Y::AbstractVector{<:Consequent}) =
    evaluate_rule(antecedent(rule), consequent(rule), X, Y)

# From antecedent to decision
evaluate_rule(
    ant::Formula{L},
    cons::Consequent,
    X::MultiFrameModalDataset,
    Y::AbstractVector{<:Consequent}
) where {L<:Logic} = evaluate_rule(extract_decisions(ant),cons,X,Y)

# Use decision and consequent
function evaluate_rule(
    decs::AbstractVector{<:Decision},
    cons::Consequent,
    X::MultiFrameModalDataset,
    Y::AbstractVector{<:Consequent}
)
    # Antecedent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable.
    ant_sat = evaluate_antecedent(decs,X)

    # Indices of satisfiable instances
    idxs_sat = findall(ant_sat .== true)

    # Consequent satisfaction. For each instances in X:
    #  - `false` when not satisfiable,
    #  - `true` when satisfiable,
    #  - `nothing` when antecedent does not hold.
    cons_sat = begin
        cons_sat = Vector{Union{Bool, Nothing}}(fill(nothing, length(Y)))
        idxs_true = begin
            idx_cons = findall(cons .== Y)
            intersect(idxs_sat,idx_cons)
        end
        idxs_false = begin
            idx_cons = findall(cons .!= Y)
            intersect(idxs_sat,idx_cons)
        end
        cons_sat[idxs_true]  .= true
        cons_sat[idxs_false] .= false
    end

    y_pred = begin
        y_pred = Vector{Union{Consequent, Nothing}}(fill(nothing, length(Y)))
        y_pred[idxs_sat] .= C
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
# TODO avoid Formula{L}()
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

    # """
    #     length_rule(node::FNode, operators::Operators) -> Int

    #     Computer the number of pairs in a rule (length of the rule)

    # # Arguments
    # - `node::FNode`: node on which you refer
    # - `operators::Operators`: set of operators of the considered logic

    # # Returns
    # - `Int`: number of pairs
    # """
    # function length_rule(node::FNode, operators::Operators)
    #     left_size = 0
    #     right_size = 0

    #     if !isdefined(node, :leftchild) && !isdefined(node, :rightchild)
    #         # Leaf
    #         if token(node) in operators
    #             return 0
    #         else
    #             return 1
    #         end
    #     end

    #     isdefined(node, :leftchild) && (left_size = length_rule(leftchild(node), operators))
    #     isdefined(node, :rightchild) && (right_size = length_rule(rightchild(node), operators))

    #     if token(node) in operators
    #         return left_size + right_size
    #     else
    #         return 1 + left_size + right_size
    #     end
    # end

    rule_metrics(rule::Rule{L,C}, X::MultiFrameModalDataset, Y::AbstractVector{<:Consequent}) =
        rule_metrics(extract_decisions(antecedent(rule)),cons,X,Y)

    """
        rule_metrics(args...) -> AbstractVector

        Compute frequency, error and length of the rule

    # Arguments
    - `decs::AbstractVector{<:Decision}`: vector of decisions
    - `cons::Consequent`: rule's consequent
    - `X::MultiFrameModalDataset`: dataset
    - `Y::AbstractVector{<:Consequent}`: target values of X

    # Returns
    - `AbstractVector`: metrics values vector of the rule
    """
    function rule_metrics(
        decs::AbstractVector{<:Decision},
        cons::Consequent,
        X::MultiFrameModalDataset,
        Y::AbstractVector{<:Consequent}
    )
        eval_result = evaluate_rule(decs, cons, X, Y)
        n_instances = size(X, 1)
        n_satisfy = sum(eval_result[:ant_sat])

        # Support of the rule
        rule_support =  n_satisfy / n_instances

        # Error of the rule
        rule_error = begin
            if typeof(cons) <: CLabel
                # Number of incorrectly classified instances divided by number of instances
                # satisfying the rule condition.
                misclassified_instances = length(findall(eval_result[:y_pred] .== Y))
                misclassified_instances / n_satisfy
            elseif typeof(cons) <: RLabel
                # Mean Squared Error (mse)
                idxs_sat = eval_result[:idxs_sat]
                mse(eval_result[:y_pred][idxs_sat], Y[idxs_sat])
            end
        end

        # Length of the rule
        rule_length = length(decs)

        return (;
            support   = rule_support,
            error     = rule_error,
            length    = rule_length,
        )
    end

    """
        prune_ruleset(ruleset::AbstractVector{<:Rule}) -> RuleBasedModel

        Prune the rules in ruleset with error metric

    If `s` and `decay_threshold` is unspecified, their values are set to nothing and the
    first two rows of the function set s and decay_threshold with their default values

    # Arguments
    - `ruleset::AbstractVector{<:Rule}`: rules to prune

    # Returns
    - `RuleBasedModel`: rules after the prune
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
    # Prune rules according to the confidence metric (with respect to a dataset)
    #  (and similar metrics: support, confidence, and length)
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
    R = RuleBasedModel()  # Vector of ordered list
    S = RuleBasedModel()  # Vector of rules left
    append!(rules(S), best_rules)
    #TODO: Formula{L}()
    push!(rules(S), Rule{L}(Formula{L}(), majority_vote(Y)))

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
        push!(rules(R), S[idx_best])

        # Indices of the remaining instances
        idx_remaining = begin
            eval_result = evaluate_rule(S[idx_best], D, Y)
            sat_unsat = eval_result[:ant_sat]
            # Remain in D the rule that not satisfying the best rule's condition
            findall(sat_unsat .== false)
        end
        D = D[idx_remaining,:]

        if idx_best == length(rules(S))
            return R
        elseif size(D, 1) == 0
            push!(R, Rule{L}(Formula{L}(),majority_vote(Y)))
            return R
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        # TODO: Formula{L}()
        rules(S)[length(rules(S))] = Rule{L}(Formula{L}(), majority_vote(Y[idx_remaining]))
    end
    ########################################################################################
end
