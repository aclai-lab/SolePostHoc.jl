# ---------------------------------------------------------------------------- #
#                  used by InTrees, Trepan and Refne methods                   #
# ---------------------------------------------------------------------------- #
struct MyRule
    formula::Formula   # Here the parsed formula is saved
    outcome::String    # Now storing the class label directly as string
end

function antecedent_to_string(antecedent)
    atoms = antecedent.grandchildren
    parts = String[]
    for atom in atoms
        cond = atom.value
        feat = cond.metacond.feature.i_variable
        op = cond.metacond.test_operator
        thr = cond.threshold

        op_str =
            op === (<) ? "<" :
            op === (<=) ? "≤" : op === (>) ? ">" : op === (>=) ? "≥" : string(op)

        push!(parts, "(V$feat $op_str $thr)")
    end
    return join(parts, " ∧ ")
end

# ---------------------------------------------------------------------------- #
#                      used by InTrees and Refne methods                       #
# ---------------------------------------------------------------------------- #
# Recursive function that converts an element (Atom, SyntaxBranch, or similar) to a string
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
        # Process a SyntaxBranch: we convert the token to a string and recursively process the children
        t = string(x.token)
        children_strs = map(element_to_string, x.children)
        if t == "¬"
            # For denial we assume only one child
            return "¬ " * children_strs[1]
        else
            # For other operators (e.g. "∧") join the children, enclosing the expression in parentheses
            return "(" * join(children_strs, " " * t * " ") * ")"
        end
    elseif hasproperty(x, :grandchildren)
        # If the element has the grandchildren field (e.g. LeftmostConjunctiveForm), process it with the specific function
        return leftmost_conjunctive_form_to_string(x)
    else
        return string(x)
    end
end

# Function to convert a LeftmostConjunctiveForm to a readable string,
# using element_to_string on each element of cf.grandchildren.
function leftmost_conjunctive_form_to_string(cf)
    return "(" * join(map(element_to_string, cf.grandchildren), " ∧ ") * ")"
end

function build_dnf_rules(rules)
    # 1) Group antecedent strings (conjunctions) by class
    class_to_antecedents = Dict{String,Vector{String}}()
    for r in rules
        c = r.consequent.outcome   # now is string =to "Iris-setosa"
        ant_str = antecedent_to_string(r.antecedent)
        push!(get!(class_to_antecedents, c, String[]), ant_str)
    end

    # 2) Create a vector of rules (Rule)
    minimized_rules = Rule[]  # an empty "Vector{Rule}"

    # Sort the classes to have a repeatable order (optional)
    sorted_classes = sort(collect(keys(class_to_antecedents)))
    for c in sorted_classes
        # All conjunctions for class c
        all_conjunctions = class_to_antecedents[c]
        # Join with " ∨ "
        big_dnf_str = join(all_conjunctions, " ∨ ")

        # Parse the string into a SoleLogics formula
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
        # Create the rule (formula + outcome as string)
        push!(minimized_rules, Rule(φ, c))
    end

    return minimized_rules
end

function convertApi(f)
    ll = listrules(f, use_shortforms = true)
    minimized_rules = build_dnf_rules(ll)
    ds = DecisionSet(minimized_rules)
end

# ---------------------------------------------------------------------------- #
#                            Convert Symbolic Rules                            #
# ---------------------------------------------------------------------------- #
# function that, given a vector of ClassificationRule (ll),
# groups the antecedents by outcome and creates a new Rule for each outcome.
# The antecedent of the new Rule is obtained by concatenating (with " ∨ ")
# the strings corresponding to each rule (obtained with leftmost_conjunctive_form_to_string).
function convert_symbolic_rules(
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

        #println("pre : ", φ)
        #dump(φ)

        #Apply DNF transformation
        dnf_result = dnf(φ)
        #println("post :", dnf_result)
        #dump(dnf_result)

        # Convert DNF form back to SyntaxBranch representation
        syntax_branch = dnf_to_syntaxbranch(dnf_result)

        # Create new rule with the SyntaxBranch representation
        new_rule = Rule(syntax_branch, outcome)
        push!(rules, new_rule)
    end
    return rules
end

function convert_symbolic_rules(
    ::SoleModels.DecisionEnsemble,
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
        println("DNF per outcome ", outcome, ": ", dnf_string)
        # Parsing DNF string into a formula
        φ = SoleLogics.parseformula(dnf_string; atom_parser = atom_parser)
        println("pre : ", φ)
        φ = dnf(φ, reduce_negations = false, allow_atom_flipping = true)
        new_rule = Rule(φ, outcome)
        println("post :", φ)
        push!(rules, new_rule)
    end
    return rules
end
