using SoleModels: DecisionList, rulebase, defaultconsequent, info, Rule
using SoleModels: metacond, feature
using SoleLogics: SyntaxTree, height, children, atoms, noperators, value

"""
    GENERATE_MUTATION_FUNCTION(;
        nmaxmutations::Integer = 1,
        decisionlistlevel::Bool = false,
        rulelevel::Bool = false,
        conjunctlevel::Bool = false,
        leaflevel::Bool = false,
    )::Function

It generates a mutation function following the indicated kwargs:
 - `nmaxmutations` is a parameter that indicates whether the individual should be mutated
once or a random number of times
 - `decisionlistlevel, rulelevel, conjunctlevel, leaflevel` indicate the depth level of the
mutation to be achieved; if all are set to false, then one of the possible mutations is
chosen
"""
function GENERATE_MUTATION_FUNCTION(;
    nmaxmutations::Integer = 1,
    decisionlistlevel::Bool = false,
    rulelevel::Bool = false,
    conjunctlevel::Bool = false,
    leaflevel::Bool = false,
)::Function
    function mutation(indiv::DLIndividual; rng::AbstractRNG = default_rng())::DLIndividual

        dl_level = [swaprules, removerules, updatedefault, phoenix]
        r_level = [removeconjunct] #updateconsequent,
        #c_level = [transplant] #negateconjunct, negationsubtree, logicalchange
        l_level = [randomatom] #inverseatom

        user_choices = [decisionlistlevel, rulelevel, leaflevel] #conjunctlevel,
        all_level = [dl_level, r_level, l_level] #c_level,
        possible_mutations = begin
            if user_choices == [false, false, false] #false,
                reduce(vcat, all_level)
            else
                println("user_choices: $(user_choices)")
                reduce(vcat, all_level[user_choices])
            end
        end

        nmutations = rand(rng, 1:nmaxmutations)
        mutation_methods = rand(rng, possible_mutations, nmutations)

        declist, alp, cls, adls = dl(indiv), alphabet(indiv), classes(indiv), alldls(indiv)
        new_declist = begin
            new_declist = declist
            for m in mutation_methods
                new_declist = m(declist, alp, cls, adls; rng = rng)
            end
            new_declist
        end

        return DLIndividual(new_declist, alp, cls, adls)
    end

    return mutation
end

############################################################################################
############################################################################################
############################################################################################

############################################################################################
# Dispatch mutation
############################################################################################

function apply_mutation_to_rule_conjunct(
    r::Rule,
    method::Function,
    args...;
    rng::AbstractRNG = default_rng(),
)
    form = antecedent(r)
    ants = form isa LeftmostLinearForm ? children(form) : [form]
    u = rand(rng, 1:length(ants))
    new_ant = method(ants[u], args...; rng = rng)

    ants = [ants[1:(u-1)]..., new_ant, ants[(u+1):end]...]
    ant = begin
        if length(ants) == 1
            first(ants)
        else
            LeftmostLinearForm(∧, ants)
        end
    end

    return Rule(ant, consequent(r), info(r))
end

function apply_mutation_to_dl_rule(
    dl::DecisionList,
    method::Function,
    args...;
    rng::AbstractRNG = default_rng(),
)
    if method ∈ [swaprules, removerules, updatedefault]
        return method(dl, args...; rng = rng)
    end

    d_rules = rulebase(dl)
    return length(d_rules) == 0 ? dl :
    begin
        u = rand(rng, 1:length(d_rules))
        new_rule = method(d_rules[u], args...; rng = rng)

        rules = [d_rules[1:(u-1)]..., new_rule, d_rules[(u+1):end]...]

        DecisionList(rules, defaultconsequent(dl), info(dl))
    end
end

############################################################################################
# DecisionList level
############################################################################################
# Move a block of rules within a single decision list
function swaprules(dl::DecisionList, args...; rng::AbstractRNG = default_rng())
    d_rules = rulebase(dl)
    return length(d_rules) == 0 ? dl :
    begin
        from, to = from_to(length(d_rules); rng = rng)
        swapped = d_rules[from:to]
        d_rules = deleteat!(d_rules, from:to)

        u = length(d_rules) == 0 ? 1 : rand(rng, 1:length(d_rules))
        d_rules = [d_rules[1:(u-1)]..., swapped..., d_rules[u:end]...]

        DecisionList(d_rules, defaultconsequent(dl), info(dl))
    end
end

# Remove a block of rules
function removerules(dl::DecisionList, args...; rng::AbstractRNG = default_rng())
    d_rules = rulebase(dl)
    return length(d_rules) == 0 ? dl :
    begin
        from, to = from_to(length(d_rules); rng = rng)
        d_rules = deleteat!(d_rules, from:to)

        DecisionList(d_rules, defaultconsequent(dl), info(dl))
    end
end

function updatedefault(
    dl::DecisionList,
    alp::Base.RefValue{<:AbstractConditionalAlphabet},
    classes::Base.RefValue{<:AbstractVector{<:Label}},
    args...;
    rng::AbstractRNG = default_rng(),
)
    return DecisionList(rulebase(dl), rand(rng, unique(classes[])), info(dl))
end

function phoenix(
    dl::DecisionList,
    alp::Base.RefValue{<:AbstractConditionalAlphabet},
    classes::Base.RefValue{<:AbstractVector{<:Label}},
    alldls::Base.RefValue{<:AbstractVector{<:DecisionList}},
    args...;
    rng::AbstractRNG = default_rng(),
)
    d_rules = rulebase(dl)
    return length(d_rules) > 0 ? dl : rand(rng, alldls[])
end

############################################################################################
# Rule level
############################################################################################
function updateconsequent(
    r::Rule,
    alp::Base.RefValue{<:AbstractConditionalAlphabet},
    classes::Base.RefValue{<:AbstractVector{<:Label}},
    args...;
    rng::AbstractRNG = default_rng(),
)
    return Rule(antecedent(r), rand(rng, unique(classes[])), info(r))
end
function updateconsequent(dl::DecisionList, args...; rng::AbstractRNG = default_rng())
    return apply_mutation_to_dl_rule(dl, updateconsequent, args...; rng = rng)
end

function removeconjunct(r::Rule, args...; rng::AbstractRNG = default_rng())
    form_rule = antecedent(r)
    return length(form_rule) == 1 ? r :
    begin
        ants = children(form_rule)
        u = rand(rng, 1:length(ants))
        deleteat!(ants, u)
        ant = begin
            if length(ants) == 1
                first(ants)
            else
                LeftmostLinearForm(connective(form_rule), ants)
            end
        end

        Rule(ant, consequent(r), info(r))
    end
end
function removeconjunct(dl::DecisionList, args...; rng::AbstractRNG = default_rng())
    return apply_mutation_to_dl_rule(dl, removeconjunct, args...; rng = rng)
end

############################################################################################
############################################################################################
############################################################################################

############################################################################################
# Conjunct level
############################################################################################
negateconjunct(st::SyntaxTree, args...; kwargs...) = ¬st
negateconjunct(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, negateconjunct, args...; rng = rng)
negateconjunct(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, negateconjunct, args...; rng = rng)

_transplant(st::SyntaxTree, rng::AbstractRNG = default_rng(), args...; kwargs...) =
    treewalk(st; rng = rng, returnnode = true)
transplant(st::SyntaxTree, args...; rng::AbstractRNG = default_rng()) =
    treewalk(st, rng; rng = rng, transformnode = _transplant)
transplant(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, transplant, args...; rng = rng)
transplant(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, transplant, args...; rng = rng)

# Using for testing:
# SyntaxTree: ⟨G⟩(max[V28] <= 7.245112655929639 ∧ (⟨L̅⟩(min[V20] >= 4222.6591159789605 ∧ (⟨L̅⟩(max[V11] <= 0.0038141608453366675 ∧ (⟨A̅⟩(max[V29] <= 178.31522392540964)))))))
# st = antecedent(posconsequent(posconsequent(posconsequent(root(tree2))))).formula
# SyntaxTree: ([G](max[V28] > 7.245112655929639)) ∧ ((⟨G⟩(min[V8] >= 494.33421895459713)) ∧ ([G](min[V8] >= 494.33421895459713 → (([L̅](min[V27] < 87446.39318797569)) ∧ (⟨=⟩(max[V2] <= 42.36525041432014))))))
# st_row = antecedent(negconsequent(posconsequent(negconsequent(root(tree2))))).formula

_negationsubtree(st::SyntaxTree, args...; kwargs...) = ¬st
negationsubtree(st::SyntaxTree, args...; rng::AbstractRNG = default_rng()) =
    treewalk(st; rng = rng, transformnode = _negationsubtree)
negationsubtree(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, negationsubtree, args...; rng = rng)
negationsubtree(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, negationsubtree, args...; rng = rng)

function _logicalchange(
    st::SyntaxTree,
    opers::Vector{<:Operator},
    rng::AbstractRNG = default_rng(),
    args...;
    kwargs...,
)
    tok = token(st)
    chs = children(st)
    opdiff = setdiff(opers, [tok])

    return isempty(opdiff) ? st : begin
        op = rand(rng, opdiff)
        connective(chs[1], chs[2])
    end
end
function logicalchange(st::SyntaxTree, args...; rng::AbstractRNG = default_rng())
    opers = filter(op->arity(op)==2, unique(operators(st)))

    return treewalk(
        st,
        opers,
        rng;
        rng = rng,
        criterion = tok->(tok isa Operator && arity(tok)==2),
        transformnode = _logicalchange,
    )
end
logicalchange(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, logicalchange, args...; rng = rng)
logicalchange(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, logicalchange, args...; rng = rng)

############################################################################################
# Leaf level
############################################################################################
function _randomatom(
    st::SyntaxTree,
    alp::Base.RefValue{<:AbstractConditionalAlphabet},
    rng::AbstractRNG = default_rng(),
    args...;
    kwargs...,
)
    a = alp[]
    featcond = value(token(st))
    orderedoperators = [<, ≤, ≥, >]
    plusminusone = rand(rng, [-1, 1])
    grouped_featconditions = a.grouped_featconditions
    policy = rand(rng, [:feature, :operator, :threshold, :stdthreshold])

    newfeatcond = begin
        if policy == :feature
            rand(rng, a; metaconditions = metacond(featcond))
        elseif policy == :operator
            idx = findall(orderedoperators .== test_operator(featcond))[]
            if (idx == 1 && plusminusone == -1) ||
               (idx == length(orderedoperators) && plusminusone == 1)
                featcond
            else
                ScalarCondition(
                    feature(featcond),
                    orderedoperators[idx+plusminusone],
                    threshold(featcond),
                )
            end
        else
            metacondition = metacond(featcond)
            filtered_featcond = filter(
                mc_thresholds->first(mc_thresholds) in [metacondition, dual(metacondition)],
                grouped_featconditions,
            )
            @assert length(filtered_featcond) == 1 "Metacondition $(metacondition) not in " *
                                                   "the alphabet\n filtered metaconditions: $(filtered_featcond)\n" *
                                                   "grouped_featconditions: $(map(first,grouped_featconditions))"

            metacondsupport = last(first(filtered_featcond))

            if policy == :threshold
                choices =
                    plusminusone == -1 ?
                    findall(metacondsupport .< threshold(featcond)) :
                    findall(metacondsupport .> threshold(featcond))
                if length(choices) != 0
                    plusminusone == -1 ?
                    ScalarCondition(metacondition, metacondsupport[last(choices)]) :
                    ScalarCondition(metacondition, metacondsupport[first(choices)])
                else
                    policy = :stdthreshold
                end
            end

            if policy == :stdthreshold
                k = rand(rng, [-0.1, 0.1])
                stdsupport = std(metacondsupport)
                ScalarCondition(metacondition, threshold(featcond) + (k*stdsupport))
            end
        end
    end

    return SyntaxTree(Atom(newfeatcond))
end
function randomatom(
    st::SyntaxTree,
    alp::Base.RefValue{<:AbstractConditionalAlphabet},
    classes::Base.RefValue{<:AbstractVector{<:Label}},
    args...;
    rng::AbstractRNG = default_rng(),
)
    return treewalk(
        st,
        alp,
        rng;
        rng = rng,
        criterion = tok->(tok isa Atom && arity(tok)==0),
        transformnode = _randomatom,
    )
end
randomatom(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, randomatom, args...; rng = rng)
randomatom(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, randomatom, args...; rng = rng)

#TODO: don't stop in an atom
_inverseatom(st::SyntaxTree, args...; kwargs...) = SyntaxTree(dual(token(st)))
inverseatom(st::SyntaxTree, args...; rng::AbstractRNG = default_rng()) = treewalk(
    st;
    rng = rng,
    criterion = tok->(tok isa Atom && arity(tok)==0),
    transformnode = _inverseatom,
)
inverseatom(r::Rule, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_rule_conjunct(r, inverseatom, args...; rng = rng)
inverseatom(dl::DecisionList, args...; rng::AbstractRNG = default_rng()) =
    apply_mutation_to_dl_rule(dl, inverseatom, args...; rng = rng)
