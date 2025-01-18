using DecisionTree
using SoleModels

# Struttura per rappresentare un nodo dell'albero di input
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
    nodes = Dict{Int, TreeNode}()
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
                current_tree = parse(Int, split(line)[2][1:end-1])
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
            
            nodes[id] = TreeNode(id, node_type, left_child, right_child, 
                               feature, threshold, depth, majority_class)
        end
    end
    return nodes
end

# Funzione ricorsiva per convertire il nodo nel formato DecisionTree.jl
function convert_to_decision_tree_node(nodes::Dict{Int, TreeNode}, node_id::Int)
    node = nodes[node_id]
    
    if node.node_type == "LN"  # Leaf node
        # Creiamo un vettore con un singolo elemento per la classe maggioritaria
        return Leaf(node.majority_class, [node.majority_class])
    else  # Internal node
        left = convert_to_decision_tree_node(nodes, node.left_child)
        right = convert_to_decision_tree_node(nodes, node.right_child)
        
        # In DecisionTree.jl, feature indices are 1-based
        feature = node.feature + 1
        
        return Node(feature, node.threshold, left, right)
    end
end

# Funzione principale che legge il file e converte l'albero
function convert_tree(filename::String)
    nodes = parse_tree_file(filename)
    # Converti partendo dal nodo radice (id = 0)
    root = convert_to_decision_tree_node(nodes, 0)
    return root
end

# Funzione per testare l'albero convertito
function test_converted_tree(tree::Node, X::Matrix{Float64})
    n_samples = size(X, 1)
    predictions = zeros(Int, n_samples)
    
    for i in 1:n_samples
        predictions[i] = apply_tree(tree, X[i, :])
    end
    
    return predictions
end





# Salva il tuo file come "tree.txt"
tree = convert_tree(joinpath(@__DIR__,"..", "temp_ba_trees", "result.txt.tree"))
t = solemodel(tree)
ll = listrules(t)
