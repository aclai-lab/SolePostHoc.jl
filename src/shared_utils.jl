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
        thr    = SoleData.(cond)

        push!(parts, "($i_name $op_str $thr)")
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
#                                  dnf rules                                   #
# ---------------------------------------------------------------------------- #
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

        φ = SoleLogics.parseformula(dnf_string; atom_parser)

        # parsing and apply dnf transformation
        φ = SoleLogics.dnf(φ)

        # convert dnf form back to SyntaxBranch representation
        φ = dnf_to_syntaxbranch(φ)

        # create new rule with the SyntaxBranch representation
        new_rule = Rule(φ, outcome)
        push!(rules, new_rule)
    end

    return rules
end
