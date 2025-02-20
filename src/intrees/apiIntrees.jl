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
            return "(" * join(["V$(i_var) $(lower_op) $(cond.minval)", "V$(i_var) $(upper_op) $(cond.maxval)"], " ∧ ") * ")"
        elseif hasproperty(cond, :metacond)
            # Standard scalar condition
            i_var = cond.metacond.feature.i_variable
            op_fun = cond.metacond.test_operator
            thr = cond.threshold
            op_str = op_fun === (<)  ? "<"  :
                     op_fun === (<=) ? "≤"  :
                     op_fun === (>)  ? ">"  :
                     op_fun === (>=) ? "≥"  : string(op_fun)
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

# Custom Atom parser for SoleLogics.parseformula.
# Since conditions are now written as scalar conditions, we use the parser for ScalarCondition here.
atom_parser = function(a::String)
    println("Parsing atom: ", a)
    return Atom(parsecondition(SoleData.ScalarCondition, a;
             featuretype = SoleData.VariableValue,
             featvaltype = Real))
end

# Function that, given a vector of ClassificationRule (ll),
# groups the antecedents by outcome and creates a new Rule for each outcome.
# The antecedent of the new Rule is obtained by concatenating (with " ∨ ")
# the strings corresponding to each rule (obtained with leftmost_conjunctive_form_to_string).
function convert_classification_rules(ll::Vector)
    # Group antecedent strings by outcome
    grouped = Dict{String, Vector{String}}()
    for rule in ll
        outcome = rule.consequent.outcome
        antecedent_str = hasproperty(rule.antecedent, :grandchildren) ? 
            leftmost_conjunctive_form_to_string(rule.antecedent) :
            string(rule.antecedent)
        push!(get!(grouped, outcome, String[]), antecedent_str)
    end

    # For each outcome, join the strings (each one is a conjunction) with " ∨ "
    # and create a new Rule parsing the DNF string
    rules = Rule[]
    for (outcome, antecedent_list) in grouped
         dnf_string = join(antecedent_list, " ∨ ")
         println("DNF per outcome ", outcome, ": ", dnf_string)
         # Parsing DNF string into a formula
         φ = SoleLogics.parseformula(
             dnf_string;
             atom_parser = atom_parser
         )
         new_rule = Rule(φ, outcome)
         push!(rules, new_rule)
    end
    return rules
end