using DecisionTree
using SoleModels
using SoleLogics


# struct MyRule
#     formula::Formula   # Here the parsed formula is saved
#     outcome::String    # Now storing the class label directly as string
# end

# function antecedent_to_string(antecedent)
#     atoms = antecedent.grandchildren
#     parts = String[]
#     for atom in atoms
#         cond = atom.value
#         feat = cond.metacond.feature.i_variable
#         op   = cond.metacond.test_operator
#         thr  = cond.threshold

#         op_str = op === (<)  ? "<" :
#                  op === (<=) ? "≤" :
#                  op === (>)  ? ">" :
#                  op === (>=) ? "≥" : 
#                  string(op)

#         push!(parts, "(V$feat $op_str $thr)")
#     end
#     return join(parts, " ∧ ")
# end

# function build_dnf_rules(rules)
#     # 1) Group antecedent strings (conjunctions) by class
#     class_to_antecedents = Dict{String, Vector{String}}()
#     for r in rules
#         c = r.consequent.outcome   # ora è già una stringa come "Iris-setosa"
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
#                     featvaltype = Real
#                 )
#             )
#         )

#         # Create the rule (formula + outcome as string)
#         push!(minimized_rules, Rule(φ, c))
#     end
    
#     return minimized_rules
# end

# function convertApi(f)
#     ll = listrules(f)
#     minimized_rules = build_dnf_rules(ll)
#     ds = DecisionSet(minimized_rules)
# end