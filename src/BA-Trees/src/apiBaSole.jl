module TradModule
export BAinDS
export BAinSoleTree

using DecisionTree
using SoleModels
using SoleLogics

# Structure to represent an input tree node
mutable struct TreeNode
    id::Int
    node_type::String
    left_child::Int
    right_child::Int
    feature::Int
    threshold::Float64
    depth::Int
    majority_class::Int
end

function parse_tree_file(filename::String)
    nodes = Dict{Int,TreeNode}()
    current_tree = -1

    open(filename) do file
        for line in eachline(file)
            if startswith(line, "DATASET_NAME:") ||
               startswith(line, "ENSEMBLE:") ||
               startswith(line, "NB_TREES:") ||
               startswith(line, "NB_FEATURES:") ||
               startswith(line, "NB_CLASSES:") ||
               startswith(line, "MAX_TREE_DEPTH:") ||
               startswith(line, "Format:") ||
               isempty(strip(line))
                continue
            end

            if startswith(line, "[TREE")
                current_tree = parse(Int, split(line)[2][1:(end-1)])
                continue
            end

            if startswith(line, "NB_NODES:")
                continue
            end

            # Parse node line
            parts = split(line)
            id = parse(Int, parts[1])
            node_type = parts[2]
            left_child = parse(Int, parts[3])
            right_child = parse(Int, parts[4])
            feature = parse(Int, parts[5])
            threshold = parse(Float64, parts[6])
            depth = parse(Int, parts[7])
            majority_class = parse(Int, parts[8])

            nodes[id] = TreeNode(
                id,
                node_type,
                left_child,
                right_child,
                feature,
                threshold,
                depth,
                majority_class,
            )
        end
    end
    return nodes
end

# Recursive function to convert the node into DecisionTree.jl format
function convert_to_decision_tree_node(nodes::Dict{Int,TreeNode}, node_id::Int)
    node = nodes[node_id]

    if node.node_type == "LN"  # Leaf node
        # Create a vector with a single element for the majority class
        return Leaf(node.majority_class, [node.majority_class])
    else  # Internal node
        left = convert_to_decision_tree_node(nodes, node.left_child)
        right = convert_to_decision_tree_node(nodes, node.right_child)

        # In DecisionTree.jl, feature indices are 1-based
        feature = node.feature + 1

        return Node(feature, node.threshold, left, right)
    end
end

# Main function that reads the file and converts the tree
function convert_tree(filename::String)
    nodes = parse_tree_file(filename)
    # Convert starting from the root node (id = 0)
    root = convert_to_decision_tree_node(nodes, 0)
    return root
end

# Function to test the converted tree
function test_converted_tree(tree::Node, X::Matrix{Float64})
    n_samples = size(X, 1)
    predictions = zeros(Int, n_samples)

    for i = 1:n_samples
        predictions[i] = apply_tree(tree, X[i, :])
    end

    return predictions
end

using SoleLogics
using SoleData

# Define a custom "container" for formula + outcome
struct MyRule
    formula::Formula   # Here the parsed formula is saved
    outcome::Int       # Here the class/label is saved
end

"""
Converts the antecedent into a string like
(V4 < 0.75) ∧ (V3 < 4.85) ∧ (V4 < 0.7)
"""
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


function build_dnf_rules(rules, class_map)
    # Map: from integer (0,1,2) to string with the name of the iris

    # 1) Group antecedent strings (conjunctions) by class
    class_to_antecedents = Dict{Int,Vector{String}}()
    for r in rules
        c = r.consequent.outcome   # es. 0,1,2
        ant_str = antecedent_to_string(r.antecedent)
        push!(get!(class_to_antecedents, c, String[]), ant_str)
    end

    # 2) Create a vector of rules (Rule)
    minimized_rules = Rule[]

    # Sort the classes to have a repeatable order
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

        # NOTE: Convert the numeric class c to the corresponding string name
        class_label = class_map[c]

        # Create the rule (formula + outcome as string)
        push!(minimized_rules, Rule(φ, class_label))
    end

    return minimized_rules
end



function BAinSoleTree()
    println(
        "Converting the tree in : ",
        joinpath(@__DIR__, "temp_ba_trees", "result.txt.tree"),
    )
    tree = convert_tree(joinpath(@__DIR__, "temp_ba_trees", "result.txt.tree"))
    t = solemodel(tree)
    return t
end


function BAinDS(class_map; silent = true)
    silent || println(
        "Converting the tree in : ",
        joinpath(@__DIR__, "temp_ba_trees", "result.txt.tree"),
    )
    tree = convert_tree(joinpath(@__DIR__, "temp_ba_trees", "result.txt.tree"))
    t = solemodel(tree)

    ll = listrules(t)
    silent || println("Rules: ", ll)
    silent || println("Class map: ", class_map)
    inverted_map = Dict(value => key for (key, value) in class_map)
    silent || println("Inverted map: ", inverted_map)
    minimized_rules = build_dnf_rules(ll, inverted_map)
    ds = DecisionSet(minimized_rules)
    return ds
end

end
