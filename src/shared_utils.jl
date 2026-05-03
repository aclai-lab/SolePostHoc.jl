#===========================================================================================
                            used by InTrees, Trepan and Refne methods
===========================================================================================#
# struct MyRule
#     formula::Formula   # Here the parsed formula is saved
#     outcome::String    # Now storing the class label directly as string
# end

# ---------------------------------------------------------------------------- #
#                              element to string                               #
# ---------------------------------------------------------------------------- #
@inline get_op_str(op) = op == (<) ?
                         "<" : op == (<=) ?
                         "≤" : op == (>) ?
                         ">" : op == (>=) ?
                         "≥" : string(op)

function antecedent_to_string(antecedent)
    atoms = antecedent.grandchildren
    parts = String[]
    for atom in atoms
        cond = SoleLogics.value(atom)
        i_name = SoleData.featurename(SoleData.feature(cond))
        op_str = get_op_str(SoleData.test_operator(cond))
        thr = SoleData.threshold(cond)

        push!(parts, "([$i_name] $op_str $thr)")
    end
    return join(parts, " ∧ ")
end

# recursive function that converts an element (Atom, SyntaxBranch, or similar) to a string
function _element_to_string(x::Atom)
    cond = SoleLogics.value(x)
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
        # for other operators (e.g. "∧") join the children,
        # enclosing the expression in parentheses
        "(" * join(children_strs, " " * string(t) * " ") * ")"
    end
end

function _element_to_string(x::Any)
    if hasproperty(x, :grandchildren)
        # if the element has the grandchildren field (e.g. LeftmostConjunctiveForm),
        # process it with the specific function
        _element_to_string(grandchildren(x))
    else
        syntaxstring(x; threshold_digits=2)
    end
end


# Function to convert a LeftmostConjunctiveForm to a readable string,
# using _element_to_string on each element of cf.grandchildren.
@inline lf_to_string(cf::SoleLogics.LeftmostConjunctiveForm) =
    "(" * join(map(_element_to_string, SoleLogics.grandchildren(cf)), " ∧ ") * ")"

function element_to_string(x)
    if x isa Atom
        # Returns a string representing the range as two conditions
        cond = x.value
        if cond isa SoleData.RangeScalarCondition
            i_var = cond.feature.i_variable
            lower_op = cond.minincluded ? "≥" : ">"
            upper_op = cond.maxincluded ? "≤" : "<"
            # Returns a string representing the range as two conditions
            return "(" *
                   join(
                       [
                           "V$(i_var) $(lower_op) $(cond.minval)",
                           "V$(i_var) $(upper_op) $(cond.maxval)",
                       ],
                       " ∧ ",
                   ) *
                   ")"
        elseif hasproperty(cond, :metacond)
            # Standard scalar condition
            i_var = cond.metacond.feature.i_variable
            op_fun = cond.metacond.test_operator
            thr = cond.threshold
            op_str =
                op_fun === (<) ? "<" :
                op_fun === (<=) ? "≤" :
                op_fun === (>) ? ">" : op_fun === (>=) ? "≥" : string(op_fun)
            return "V$(i_var) $(op_str) $(thr)"
        else
            return string(cond)
        end
    elseif x isa SyntaxBranch
        # Process a SyntaxBranch: we convert the token to a string and recursively process
        # the children
        t = string(x.token)
        children_strs = map(element_to_string, x.children)
        if t == "¬"
            # For denial we assume only one child
            return "¬ " * children_strs[1]
        else
            # For other operators (e.g. "∧") join the children, enclosing the expression
            # in parentheses
            return "(" * join(children_strs, " " * t * " ") * ")"
        end
    elseif hasproperty(x, :grandchildren)
        # If the element has the grandchildren field (e.g. LeftmostConjunctiveForm),
        # process it with the specific function
        return leftmost_conjunctive_form_to_string(x)
    else
        return string(x)
    end
end

function leftmost_conjunctive_form_to_string(cf)
    return "(" * join(map(element_to_string, cf.grandchildren), " ∧ ") * ")"
end

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
)
    # Group antecedent strings by outcome
    grouped = Dict{String,Vector{String}}()
    for rule in ll
        outcome = rule.consequent.outcome
        antecedent_str =
            hasproperty(rule.antecedent, :grandchildren) ?
            leftmost_conjunctive_form_to_string(rule.antecedent) : string(rule.antecedent)
        push!(get!(grouped, outcome, String[]), antecedent_str)
    end

    # For each outcome, join the strings (each one is a conjunction) with " ∨ "
    # and create a new Rule parsing the DNF string
    rules = Rule[]
    for (outcome, antecedent_list) in grouped
        dnf_string = join(antecedent_list, " ∨ ")
        #println("DNF per outcome ", outcome, ": ", dnf_string)

        # Parsing DNF string into a formula
        φ = SoleLogics.parseformula(dnf_string; atom_parser = atom_parser)

        #Apply DNF transformation
        dnf_result = dnf(φ)

        # Convert DNF form back to SyntaxBranch representation
        syntax_branch = dnf_to_syntaxbranch(dnf_result)

        # Create new rule with the SyntaxBranch representation
        new_rule = Rule(syntax_branch, outcome)
        push!(rules, new_rule)
    end
    return rules
end

# ---------------------------------------------------------------------------- #
#                           syntaxbranch conversions                           #
# ---------------------------------------------------------------------------- #
# function to convert a formula back to SyntaxBranch representation
function _to_syntaxbranch(
    lf::Union{LeftmostLinearForm,LeftmostConjunctiveForm},
    connective::NamedConnective,
    func::Base.Callable
)
    # this is a disjunction of conjunctions (typical DNF structure)
    formula = lf.grandchildren

    if length(formula) == 1
        # if there's only one formula, we just need the conjunction
        return func(formula[1])
    elseif length(formula) == 2
        # Binary disjunction of two conjunctions
        child1 = func(formula[1])
        child2 = func(formula[2])
        return SyntaxBranch(connective, (child1, child2))
    else
        # for more than two formula, create nested binary disjunctions
        current_branch = func(formula[1])
        current_branch = foldl(formula[2:end]; init=func(formula[1])) do branch, f
            SyntaxBranch(connective, (branch, func(f)))
        end
    end

    return current_branch
end

@inline dnf_to_syntaxbranch(dnf_formula::LeftmostLinearForm) = _to_syntaxbranch(dnf_formula, NamedConnective{:∨}(), conjunction_to_syntaxbranch)
@inline dnf_to_syntaxbranch(dnf_formula::LeftmostConjunctiveForm) = conjunction_to_syntaxbranch(dnf_formula)
@inline dnf_to_syntaxbranch(dnf_formula::Literal) = literal_to_syntaxbranch(dnf_formula)
@inline dnf_to_syntaxbranch(dnf_formula::Atom) = dnf_formula

# convert a conjunction (LeftmostConjunctiveForm) to SyntaxBranch
@inline conjunction_to_syntaxbranch(conjunction::LeftmostConjunctiveForm) = _to_syntaxbranch(conjunction, NamedConnective{:∧}(), literal_to_syntaxbranch)

# convert a Literal to SyntaxBranch
@inline literal_to_syntaxbranch(literal) = return literal.ispos ? literal.atom : SyntaxBranch(NamedConnective{:¬}(), (literal.atom,))

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

# # ---------------------------------------------------------------------------- #
# # Function to convert a LeftmostConjunctiveForm to a readable string,
# # using _element_to_string on each element of cf.grandchildren.
# @inline lf_to_string(cf::SoleLogics.LeftmostConjunctiveForm) =
#     "(" * join(map(_element_to_string, SoleLogics.grandchildren(cf)), " ∧ ") * ")"

# # Function to convert a LeftmostConjunctiveForm to a readable string,
# # using element_to_string on each element of cf.grandchildren.
# function leftmost_conjunctive_form_to_string(cf)
#     return "(" * join(map(element_to_string, cf.grandchildren), " ∧ ") * ")"
# end

# function build_dnf_rules(rules)
#     # 1) Group antecedent strings (conjunctions) by class
#     class_to_antecedents = Dict{String,Vector{String}}()
#     for r in rules
#         c = r.consequent.outcome   # now is string =to "Iris-setosa"
#         ant_str = antecedent_to_string(r.antecedent)
#         push!(get!(class_to_antecedents, c, String[]), ant_str)
#     end

#     # 2) Create a vector of rules (Rule)
#     minimized_rules = Rule[]  # an empty "Vector{Rule}"

#     # Sort the classes to have a repeatable order (optional)
#     sorted_classes = sort(collect(keys(class_to_antecedents)))
#     for c in sorted_classes
#         # All conjunctions for class c
#         all_conjunctions = class_to_antecedents[c]
#         # Join with " ∨ "
#         big_dnf_str = join(all_conjunctions, " ∨ ")

#         # Parse the string into a SoleLogics formula
#         φ = SoleLogics.parseformula(
#             big_dnf_str;
#             atom_parser = a -> Atom(
#                 parsecondition(
#                     ScalarCondition,
#                     a;
#                     featuretype = VariableValue,
#                     featvaltype = Real,
#                 ),
#             ),
#         )
#         # Create the rule (formula + outcome as string)
#         push!(minimized_rules, Rule(φ, c))
#     end

#     return minimized_rules
# end

function convertApi(f)
    ll = listrules(f, use_shortforms = true)
    minimized_rules = build_dnf_rules(ll)
    ds = DecisionSet(minimized_rules)
end

# #==========================================================================================#
# # Function that, given a vector of ClassificationRule (ll),
# # groups the antecedents by outcome and creates a new Rule for each outcome.
# # The antecedent of the new Rule is obtained by concatenating (with " ∨ ")
# # the strings corresponding to each rule
# # (obtained with leftmost_conjunctive_form_to_string).
# function convert_classification_rules(
#     ::SoleModels.DecisionList,
#     ll::AbstractVector{<:SoleModels.ClassificationRule},
# )
#     # Group antecedent strings by outcome
#     grouped = Dict{String,Vector{String}}()
#     for rule in ll
#         outcome = rule.consequent.outcome
#         antecedent_str =
#             hasproperty(rule.antecedent, :grandchildren) ?
#             leftmost_conjunctive_form_to_string(rule.antecedent) : string(rule.antecedent)
#         push!(get!(grouped, outcome, String[]), antecedent_str)
#     end

#     # For each outcome, join the strings (each one is a conjunction) with " ∨ "
#     # and create a new Rule parsing the DNF string
#     rules = Rule[]
#     for (outcome, antecedent_list) in grouped
#         dnf_string = join(antecedent_list, " ∨ ")
#         #println("DNF per outcome ", outcome, ": ", dnf_string)

#         # Parsing DNF string into a formula
#         φ = SoleLogics.parseformula(dnf_string; atom_parser = atom_parser)

#         #println("pre : ", φ)
#         #dump(φ)

#         #Apply DNF transformation
#         dnf_result = dnf(φ)
#         #println("post :", dnf_result)
#         #dump(dnf_result)

#         # Convert DNF form back to SyntaxBranch representation
#         syntax_branch = dnf_to_syntaxbranch(dnf_result)

#         # Create new rule with the SyntaxBranch representation
#         new_rule = Rule(syntax_branch, outcome)
#         push!(rules, new_rule)
#     end
#     return rules
# end

# function convert_classification_rules(
#     ::SoleModels.DecisionEnsemble,
#     ll::AbstractVector{<:SoleModels.ClassificationRule},
# )
#     # Group antecedent strings by outcome
#     grouped = Dict{String,Vector{String}}()
#     for rule in ll
#         outcome = rule.consequent.outcome
#         antecedent_str =
#             hasproperty(rule.antecedent, :grandchildren) ?
#             leftmost_conjunctive_form_to_string(rule.antecedent) : string(rule.antecedent)
#         push!(get!(grouped, outcome, String[]), antecedent_str)
#     end

#     # For each outcome, join the strings (each one is a conjunction) with " ∨ "
#     # and create a new Rule parsing the DNF string
#     rules = Rule[]
#     for (outcome, antecedent_list) in grouped
#         dnf_string = join(antecedent_list, " ∨ ")
#         println("DNF per outcome ", outcome, ": ", dnf_string)
#         # Parsing DNF string into a formula
#         φ = SoleLogics.parseformula(dnf_string; atom_parser = atom_parser)
#         println("pre : ", φ)
#         φ = dnf(φ, reduce_negations = false, allow_atom_flipping = true)
#         new_rule = Rule(φ, outcome)
#         println("post :", φ)
#         push!(rules, new_rule)
#     end
#     return rules
# end
