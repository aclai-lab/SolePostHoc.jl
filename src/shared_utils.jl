# ---------------------------------------------------------------------------- #
#                              element to string                               #
# ---------------------------------------------------------------------------- #
@inline get_op_str(op) = op == (<)  ? "<" : op == (<=) ? "≤" : op == (>)  ? ">" : op == (>=) ? "≥" : string(op)

function antecedent_to_string(antecedent)
    atoms = antecedent.grandchildren
    parts = String[]
    for atom in atoms
        cond   = SoleLogics.value(atom)
        i_name = SoleData.featurename(SoleData.feature(cond))
        op_str = get_op_str(SoleData.test_operator(cond))
        thr    = SoleData.threshold(cond)

        push!(parts, "([$i_name] $op_str $thr)")
    end
    return join(parts, " ∧ ")
end

# recursive function that converts an element (Atom, SyntaxBranch, or similar) to a string
function _element_to_string(x::Atom)
    cond   = SoleLogics.value(x)
    i_name = SoleData.featurename(SoleData.feature(cond))
    digits = SoleData.get_threshold_display_method(nothing, 2)

    return if cond isa SoleData.RangeScalarCondition
        lower_op = SoleData.minincluded(cond) ? "≥" : ">"
        upper_op = SoleData.maxincluded(cond) ? "≤" : "<"
        
        "(" *
            join(
                [
                    "[$(i_name)] $(lower_op) $(digits(SoleData.minval(cond)))",
                    "[$(i_name)] $(upper_op) $(digits(SoleData.maxval(cond)))",
                ],
                " ∧ ",
            ) *
        ")"
    elseif hasproperty(cond, :metacond)
        # standard scalar condition
        op_str = get_op_str(SoleData.test_operator(cond))
        
        "[$(i_name)] $(op_str) $(digits(SoleData.threshold(cond)))"
    else
        syntaxstring(x; threshold_digits=2)
    end
end

function _element_to_string(x::SyntaxBranch)
    t = SoleLogics.token(x)
    children_strs = map(_element_to_string, SoleLogics.children(x))
    return if t == ¬
        # for denial we assume only one child
        "¬ " * children_strs[1]
    else
        # for other operators (e.g. "∧") join the children, enclosing the expression in parentheses
        "(" * join(children_strs, " " * string(t) * " ") * ")"
    end
end

function _element_to_string(x::Any)
    if hasproperty(x, :grandchildren)
        # if the element has the grandchildren field (e.g. LeftmostConjunctiveForm), process it with the specific function
        _element_to_string(grandchildren(x))
    else
        syntaxstring(x; threshold_digits=2)
    end
end

# Function to convert a LeftmostConjunctiveForm to a readable string,
# using _element_to_string on each element of cf.grandchildren.
@inline lf_to_string(cf::SoleLogics.LeftmostConjunctiveForm) =
    "(" * join(map(_element_to_string, SoleLogics.grandchildren(cf)), " ∧ ") * ")"

# ---------------------------------------------------------------------------- #
#                                  atom parser                                 #
# ---------------------------------------------------------------------------- #
# Custom Atom parser for SoleLogics.parseformula.
# Since conditions are now written as scalar conditions, we use the parser for ScalarCondition here.
atom_parser = function (a::String)
    #println("Parsing atom: ", a)
    return Atom(
        parsecondition(
            SoleData.ScalarCondition,
            a;
            featuretype = SoleData.VariableValue,
            featvaltype = Real,
        ),
    )
end

# ---------------------------------------------------------------------------- #
#                        convert classification rules                          #
# ---------------------------------------------------------------------------- #
# function that, given a vector of ClassificationRule (ll),
# groups the antecedents by outcome and creates a new Rule for each outcome.
# the antecedent of the new Rule is obtained by concatenating (with " ∨ ")
# the strings corresponding to each rule (obtained with _lf_to_string).
function convert_classif_rules(
    ::SoleModels.DecisionList,
    ll::AbstractVector{<:SoleModels.ClassificationRule},
)::Vector{Rule}
    # group antecedent strings by outcome
    grouped = Dict{String,Vector{String}}()

    for rule in ll
        outcome = consequent(rule).outcome
        antecedent_str = hasproperty(antecedent(rule), :grandchildren) ?
            lf_to_string(antecedent(rule)) :
            string(antecedent(rule))
        push!(get!(grouped, outcome, String[]), antecedent_str)
    end

    rules = Rule[]

    for (outcome, antecedent_list) in grouped
        dnf_string = join(antecedent_list, " ∨ ")

        # parsing and apply dnf transformation
        φ = SoleLogics.parseformula(dnf_string; atom_parser)
        φ = SoleLogics.dnf(φ)

        # convert dnf form back to SyntaxBranch representation
        φ = dnf_to_syntaxbranch(φ)

        # create new rule with the SyntaxBranch representation
        push!(rules, Rule(φ, outcome))
    end

    return rules

    # φ = SoleModels.joinrules(ll)
    # Rule.(SoleLogics.dnf.(SoleModels.antecedent.(φ)), SoleModels.consequent.(φ))
end

# ---------------------------------------------------------------------------- #
#                           syntaxbranch conversions                           #
# ---------------------------------------------------------------------------- #
function _to_syntaxbranch(
    formula    :: Union{LeftmostLinearForm, LeftmostConjunctiveForm},
    connective :: NamedConnective
)
    # this is a disjunction of conjunctions (typical DNF structure)
    disjuncts = formula.grandchildren

    return if length(disjuncts) == 1
        # if there's only one disjunct, we just need the conjunction
        conjunction_to_syntaxbranch(disjuncts[1])
    elseif length(disjuncts) == 2
        # binary disjunction of two conjunctions
        child1 = conjunction_to_syntaxbranch(disjuncts[1])
        child2 = conjunction_to_syntaxbranch(disjuncts[2])
        SyntaxBranch(connective, (child1, child2))
    else
        # for more than two disjuncts, create nested binary disjunctions
        current_branch = conjunction_to_syntaxbranch(disjuncts[1])

        # combine with each remaining disjunct using binary disjunction
        foldl(2:length(disjuncts); init=conjunction_to_syntaxbranch(disjuncts[1])) do branch, i
            SyntaxBranch(connective, (branch, conjunction_to_syntaxbranch(disjuncts[i])))
        end
        current_branch
    end
end

@inline literal_to_syntaxbranch(literal::Literal) =
    literal.ispos ? literal.atom : SyntaxBranch(NamedConnective{:¬}(), (literal.atom,))

@inline dnf_to_syntaxbranch(dnf_formula::LeftmostLinearForm) = _to_syntaxbranch(dnf_formula, NamedConnective{:∨}())
@inline dnf_to_syntaxbranch(dnf_formula::LeftmostConjunctiveForm) = conjunction_to_syntaxbranch(dnf_formula)
@inline dnf_to_syntaxbranch(dnf_formula::Literal) = literal_to_syntaxbranch(dnf_formula)
@inline dnf_to_syntaxbranch(dnf_formula::Atom) = dnf_formula

@inline conjunction_to_syntaxbranch(conjunction::LeftmostConjunctiveForm) =
    _to_syntaxbranch(conjunction, NamedConnective{:∧}())
@inline conjunction_to_syntaxbranch(conjunction::Literal) = literal_to_syntaxbranch(conjunction)

# ---------------------------------------------------------------------------- #
#                                decision set                                  #
# ---------------------------------------------------------------------------- #
function build_dnf_rules(rules)
    # group antecedent strings (conjunctions) by class
    class_to_antecedents = Dict{String,Vector{String}}()
    for r in rules
        c = r.consequent.outcome
        ant_str = antecedent_to_string(r.antecedent)
        push!(get!(class_to_antecedents, c, String[]), ant_str)
    end

    # create a vector of rules
    minimized_rules = Rule[]

    # sort classes to have a repeatable order
    sorted_classes = sort(collect(keys(class_to_antecedents)))
    for c in sorted_classes
        # all conjunctions for class c
        all_conjunctions = class_to_antecedents[c]
        # join with " ∨ "
        big_dnf_str = join(all_conjunctions, " ∨ ")

        # parse the string into a SoleLogics formula
        φ = SoleLogics.parseformula(
            big_dnf_str;
            atom_parser = a -> Atom(
                parsecondition(
                    ScalarCondition,
                    a;
                    featuretype = VariableValue,
                    featvaltype = Real,
                ),
            ),
        )
        # create the rule (formula + outcome as string)
        push!(minimized_rules, Rule(φ, c))
    end

    return minimized_rules
end

"""
    make_decisionset(dl::DecisionList)
    make_decisionset(de::DecisionEnsemble)

Converts a DecisionList or a DecisionEnsemble to a DecisionSet
by grouping rules by outcome.
Returns a new DecisionSet with one rule per class in DNF form.
"""
function make_decisionset(dl::DecisionList)::DecisionSet
    rules_list = listrules(dl, use_shortforms = true)

    # group antecedents by outcome
    class_to_antecedents = Dict{String,Vector{String}}()

    for r in rules_list
        c = r.consequent.outcome
        ant_str = el_to_string(r.antecedent)
        push!(get!(class_to_antecedents, c, String[]), ant_str)
    end

    # Create minimized rules
    minimized_rules = Rule[]

    # Sort classes for repeatable order
    sorted_classes = sort(collect(keys(class_to_antecedents)))
    for c in sorted_classes
        all_conjunctions = class_to_antecedents[c]
        big_dnf_str = join(all_conjunctions, " ∨ ")

        φ = SoleLogics.parseformula(
            big_dnf_str;
            atom_parser = a -> Atom(
                parsecondition(
                    ScalarCondition,
                    a;
                    featuretype = VariableValue,
                    featvaltype = Real,
                ),
            ),
        )
        println("pre : ", φ)
        φ = dnf(φ)
        push!(minimized_rules, Rule(φ, c))
        println("post :", φ)
    end

    return DecisionSet(minimized_rules)
end

function make_decisionset(de::DecisionEnsemble)::DecisionSet
    ll = listrules(de, use_shortforms = true)
    minimized_rules = build_dnf_rules(ll)
    ds = DecisionSet(minimized_rules)
end