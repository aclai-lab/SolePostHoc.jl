# Recursive function that converts an element (Atom, SyntaxBranch, or similar) to a string

#### already defined in apiREFNESole.jl

# function element_to_string(x)
#     if x isa Atom
#         # Returns a string representing the range as two conditions
#         cond = x.value
#         if cond isa SoleData.RangeScalarCondition
#             i_var = cond.feature.i_variable
#             lower_op = cond.minincluded ? "≥" : ">"
#             upper_op = cond.maxincluded ? "≤" : "<"
#             # Returns a string representing the range as two conditions
#             return "(" * join(["V$(i_var) $(lower_op) $(cond.minval)", "V$(i_var) $(upper_op) $(cond.maxval)"], " ∧ ") * ")"
#         elseif hasproperty(cond, :metacond)
#             # Standard scalar condition
#             i_var = cond.metacond.feature.i_variable
#             op_fun = cond.metacond.test_operator
#             thr = cond.threshold
#             op_str = op_fun === (<)  ? "<"  :
#                      op_fun === (<=) ? "≤"  :
#                      op_fun === (>)  ? ">"  :
#                      op_fun === (>=) ? "≥"  : string(op_fun)
#             return "V$(i_var) $(op_str) $(thr)"
#         else
#             return string(cond)
#         end
#     elseif x isa SyntaxBranch
#         # Process a SyntaxBranch: we convert the token to a string and recursively process the children
#         t = string(x.token)
#         children_strs = map(element_to_string, x.children)
#         if t == "¬"
#             # For denial we assume only one child
#             return "¬ " * children_strs[1]
#         else
#             # For other operators (e.g. "∧") join the children, enclosing the expression in parentheses
#             return "(" * join(children_strs, " " * t * " ") * ")"
#         end
#     elseif hasproperty(x, :grandchildren)
#         # If the element has the grandchildren field (e.g. LeftmostConjunctiveForm), process it with the specific function
#         return leftmost_conjunctive_form_to_string(x)
#     else
#         return string(x)
#     end
# end

# Function to convert a LeftmostConjunctiveForm to a readable string,
# using element_to_string on each element of cf.grandchildren.

#### already defined in apiREFNESole.jl

# function leftmost_conjunctive_form_to_string(cf)
#     return "(" * join(map(element_to_string, cf.grandchildren), " ∧ ") * ")"
# end

# Custom Atom parser for SoleLogics.parseformula.
# Since conditions are now written as scalar conditions, we use the parser for ScalarCondition here.
atom_parser = function(a::String)
    #println("Parsing atom: ", a)
    return Atom(parsecondition(SoleData.ScalarCondition, a;
             featuretype = SoleData.VariableValue,
             featvaltype = Real))
end

# Function to convert a DNF formula back to SyntaxBranch representation
function dnf_to_syntaxbranch(dnf_formula)
    if dnf_formula isa LeftmostLinearForm
        # This is a disjunction of conjunctions (typical DNF structure)
        disjuncts = dnf_formula.grandchildren
        
        if length(disjuncts) == 1
            # If there's only one disjunct, we just need the conjunction
            return conjunction_to_syntaxbranch(disjuncts[1])
        elseif length(disjuncts) == 2
            # Binary disjunction of two conjunctions
            child1 = conjunction_to_syntaxbranch(disjuncts[1])
            child2 = conjunction_to_syntaxbranch(disjuncts[2])
            return SyntaxBranch(NamedConnective{:∨}(), (child1, child2))
        else
            # For more than two disjuncts, create nested binary disjunctions
            current_branch = conjunction_to_syntaxbranch(disjuncts[1])
            
            # Combine with each remaining disjunct using binary disjunction
            for i in 2:length(disjuncts)
                next_disjunct = conjunction_to_syntaxbranch(disjuncts[i])
                current_branch = SyntaxBranch(NamedConnective{:∨}(), (current_branch, next_disjunct))
            end
            
            return current_branch
        end
    elseif dnf_formula isa LeftmostConjunctiveForm
        # This is just a conjunction
        return conjunction_to_syntaxbranch(dnf_formula)
    elseif dnf_formula isa Literal
        # This is a single literal
        return literal_to_syntaxbranch(dnf_formula)
    elseif dnf_formula isa Atom
        # This is a single atom
        return dnf_formula
    else
        error("Unexpected formula type: $(typeof(dnf_formula))")
    end
end

# Convert a conjunction (LeftmostConjunctiveForm) to SyntaxBranch
function conjunction_to_syntaxbranch(conjunction)
    literals = conjunction.grandchildren
    
    if length(literals) == 1
        # If there's only one literal, we don't need a conjunction
        return literal_to_syntaxbranch(literals[1])
    elseif length(literals) == 2
        # Create a conjunction of two literals (binary operation)
        child1 = literal_to_syntaxbranch(literals[1])
        child2 = literal_to_syntaxbranch(literals[2])
        return SyntaxBranch(NamedConnective{:∧}(), (child1, child2))
    else
        # For more than two literals, create nested binary conjunctions
        # Convert first literal to SyntaxBranch
        current_branch = literal_to_syntaxbranch(literals[1])
        
        # Combine with each remaining literal using binary conjunction
        for i in 2:length(literals)
            next_literal = literal_to_syntaxbranch(literals[i])
            current_branch = SyntaxBranch(NamedConnective{:∧}(), (current_branch, next_literal))
        end
        
        return current_branch
    end
end

# Convert a Literal to SyntaxBranch

#### already defined in apiREFNESole.jl

# function literal_to_syntaxbranch(literal)
#     if literal.ispos
#         # Positive literal, just return the atom
#         return literal.atom
#     else
#         # Negative literal, wrap the atom in a negation
#         return SyntaxBranch(NamedConnective{:¬}(), (literal.atom,))
#     end
# end

# Function that, given a vector of ClassificationRule (ll),
# groups the antecedents by outcome and creates a new Rule for each outcome.
# The antecedent of the new Rule is obtained by concatenating (with " ∨ ")
# the strings corresponding to each rule (obtained with leftmost_conjunctive_form_to_string).

#### already defined in apiREFNESole.jl

# function convert_classification_rules(ll::Vector)
#     # Group antecedent strings by outcome
#     grouped = Dict{String, Vector{String}}()
#     for rule in ll
#         outcome = rule.consequent.outcome
#         antecedent_str = hasproperty(rule.antecedent, :grandchildren) ? 
#             leftmost_conjunctive_form_to_string(rule.antecedent) :
#             string(rule.antecedent)
#         push!(get!(grouped, outcome, String[]), antecedent_str)
#     end

#     # For each outcome, join the strings (each one is a conjunction) with " ∨ "
#     # and create a new Rule parsing the DNF string
#     rules = Rule[]
#     for (outcome, antecedent_list) in grouped
#         dnf_string = join(antecedent_list, " ∨ ")
#         #println("DNF per outcome ", outcome, ": ", dnf_string)
        
#         # Parsing DNF string into a formula
#         φ = SoleLogics.parseformula(
#             dnf_string;
#             atom_parser = atom_parser
#         )

#         #println("pre : ", φ)
#         #dump(φ)
        
#         # Apply DNF transformation
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


#==============================================================================================#
#### already defined in apiREFNESole.jl

# struct MyRule
#     formula::Formula   # Here the parsed formula is saved
#     outcome::String    # Now storing the class label directly as string
# end

function antecedent_to_string(antecedent)
    atoms = antecedent.grandchildren
    parts = String[]
    for atom in atoms
        cond = atom.value
        feat = cond.metacond.feature.i_variable
        op   = cond.metacond.test_operator
        thr  = cond.threshold

        op_str = op === (<)  ? "<" :
                 op === (<=) ? "≤" :
                 op === (>)  ? ">" :
                 op === (>=) ? "≥" : 
                 string(op)

        push!(parts, "(V$feat $op_str $thr)")
    end
    return join(parts, " ∧ ")
end

function build_dnf_rules(rules)
    # 1) Group antecedent strings (conjunctions) by class
    class_to_antecedents = Dict{String, Vector{String}}()
    for r in rules
        c = r.consequent.outcome   # ora è già una stringa come "Iris-setosa"
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
                    featvaltype = Real
                )
            )
        )

        # Create the rule (formula + outcome as string)
        push!(minimized_rules, Rule(φ, c))
    end
    
    return minimized_rules
end

#### already defined in apiREFNESole.jl

# function convertApi(f)
#     ll = listrules(f)
#     minimized_rules = build_dnf_rules(ll)
#     ds = DecisionSet(minimized_rules)
# end