using SoleLogics
using SoleData

"""
    condition_to_string(cond)

Converts a scalar condition to a string representation.
Returns a string in the format "(Vx op val)" where:
- x is the variable number
- op is the operator (>, ≥, <, ≤, ==)
- val is the threshold value
"""
function condition_to_string(cond)
    if cond isa SoleData.RangeScalarCondition
        i_var = cond.feature.i_variable
        # Remove any existing 'V' prefix to avoid duplication
        i_var = replace(string(i_var), r"^V+" => "")
        lower_op = cond.minincluded ? "≥" : ">"
        upper_op = cond.maxincluded ? "≤" : "<"
        return "(" * join(["V$(i_var) $(lower_op) $(cond.minval)",
                       "V$(i_var) $(upper_op) $(cond.maxval)"], " ∧ ") * ")"
    elseif hasproperty(cond, :metacond)
        i_var = cond.metacond.feature.i_variable
        # Remove any existing 'V' prefix to avoid duplication
        i_var = replace(string(i_var), r"^V+" => "")
        op_fun = cond.metacond.test_operator
        thr = cond.threshold
        op_str = op_fun === (<) ? "<" :
                 op_fun === (<=) ? "≤" :
                 op_fun === (>) ? ">" :
                 op_fun === (>=) ? "≥" : string(op_fun)
        return "(V$(i_var) $(op_str) $(thr))"
    else
        return string(cond)
    end
end

# [Il resto del codice rimane invariato...]

"""
    element_to_string(x)

Recursively converts a decision tree element to a string representation.
Handles Atoms, SyntaxBranches, and structures with grandchildren.
"""

#### already defined in apiREFNESole.jl

function element_to_string(x)
    if x isa Atom
        return condition_to_string(x.value)
    elseif x isa SyntaxBranch
        t = string(x.token)
        children_strs = map(element_to_string, x.children)
        if t == "¬"
            return "¬ " * children_strs[1]
        else
            return "(" * join(children_strs, " " * t * " ") * ")"
        end
    elseif hasproperty(x, :grandchildren)
        return "(" * join(map(element_to_string, x.grandchildren), " ∧ ") * ")"
    else
        return string(x)
    end
end

"""
    rule_to_string(rule)

Converts a classification rule to a string representation.
Returns a string in the format "antecedent ↣ outcome"
"""
function rule_to_string(rule)
    antecedent_str = element_to_string(rule.antecedent)
    return "$(antecedent_str) ↣ $(rule.consequent.outcome)"
end

"""
    decision_list_to_string(dl)

Converts an entire decision list to a formatted string representation.
Includes rule numbering and default consequent.
"""
function decision_list_to_string(dl)
    result = "▣\n"
    for (i, rule) in enumerate(dl.rulebase)
        result *= "├[$i/$(length(dl.rulebase))]┐ $(rule_to_string(rule)) : NamedTuple()\n"
    end
    result *= "└✘ $(dl.defaultconsequent.outcome) : NamedTuple()"
    return result
end

"""
    parse_rule_string(rule_str)

Parses a rule string back into a formula.
Accepts strings in the format created by rule_to_string.
"""
function parse_rule_string(rule_str)
    parts = split(rule_str, "↣")
    if length(parts) != 2
        error("Invalid rule format: $rule_str")
    end

    antecedent_str = strip(parts[1])
    outcome = strip(replace(parts[2], r": .*$" => ""))

    φ = SoleLogics.parseformula(
        antecedent_str;
        atom_parser=a -> Atom(
            parsecondition(
                ScalarCondition,
                a;
                featuretype=VariableValue,
                featvaltype=Real
            )
        )
    )

    return φ, outcome
end

"""
    convertApi(ds::DecisionList)

Converts a DecisionList to a DecisionSet by grouping rules by outcome.
Returns a new DecisionSet with one rule per class in DNF form.
"""
function convertApi(ds::DecisionList)
    # Get list of rules
    rules_list = listrules(ds, use_shortforms=true)

    # Group antecedents by outcome
    class_to_antecedents = Dict{String,Vector{String}}()

    for r in rules_list
        c = r.consequent.outcome
        ant_str = element_to_string(r.antecedent)
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
            atom_parser=a -> Atom(
                parsecondition(
                    ScalarCondition,
                    a;
                    featuretype=VariableValue,
                    featvaltype=Real
                )
            )
        )
        println("pre : ", φ)
        φ = dnf(φ)
        push!(minimized_rules, Rule(φ, c))
        println("post :", φ)
    end

    return DecisionSet(minimized_rules)
end
