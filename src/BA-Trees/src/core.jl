module RunModule

using DelimitedFiles
using Random
using CSV
using DataFrames
using Statistics

using DecisionTree
using SoleModels

# <-- PATH FIX: usa @__DIR__ per includere "other.jl"
include(joinpath(@__DIR__, "utils/minor.jl"))

export WRAP_batrees

"""
Calculate the mode (most frequent value) of an array.

# Arguments
- `x`: An array of elements for which the mode is to be calculated.

# Returns
- The mode of the array, which is the first most frequent value in case of ties.

# Throws
- `ArgumentError`: If the input array is empty.
"""
function mode(x)
    if isempty(x)
        throw(ArgumentError("array must be non-empty"))
    end

    # Count occurrences of each value
    counts = Dict{eltype(x),Int}()
    for val in x
        counts[val] = get(counts, val, 0) + 1
    end

    # Find value with maximum count
    max_count = maximum(values(counts))
    modes = [k for (k, v) in counts if v == max_count]

    # Return first modal value (in case of multi-modality)
    return first(modes)
end

"""
Convert an input feature to a numeric value.

# Arguments
- `x`: The feature to be converted, which can be a number, string, or symbol.

# Returns
- A `Float64` representation of the input feature.
"""
function convert_feature(x)
    if x isa Number
        return Float64(x)
    elseif x isa String || x isa Symbol
        # For strings/symbols, convert to a numeric index
        return 1.0  # Simplified for this example
    else
        return Float64(x)
    end
end

"""
Get the majority class from an array of target values.

# Arguments
- `targets::AbstractArray`: An array of target values.

# Returns
- The majority class as an integer.
"""
function get_majority_class(targets)
    if isempty(targets)
        return 0
    end
    targets_int = convert(Vector{Int}, targets)
    return mode(targets_int)
end

"""
Prepare dataset for analysis by converting features and targets to numeric values.

# Arguments
- `dataset::AbstractArray`: The dataset to be processed.

# Returns
- Tuple of numeric data matrix and target vector.
"""
function prepare_dataset(dataset)
    n_features = size(dataset, 2) - 1
    n_samples = size(dataset, 1)

    # Convert all features to numbers
    numeric_data = zeros(Float64, n_samples, n_features)
    for i = 1:n_features
        numeric_data[:, i] = map(x -> convert_feature(x), dataset[:, i])
    end

    # Convert target classes to 0-based indices
    unique_classes = unique(dataset[:, end])
    class_map = Dict(class => i - 1 for (i, class) in enumerate(unique_classes))
    targets = [class_map[x] for x in dataset[:, end]]

    return numeric_data, convert(Vector{Int}, targets)
end

"""
Create a decision tree node with random splits.

# Arguments
- `data`: Dataset matrix
- `targets`: Target vector
- `depth`: Current depth
- `max_depth`: Maximum allowed depth
- `min_samples`: Minimum samples for split
- `node_id`: Current node ID

# Returns
- Array of strings representing tree nodes
"""
function create_tree_node(
    data,
    targets,
    depth = 0,
    max_depth = 3,
    min_samples = 5,
    node_id = 0,
)
    nodes = String[]
    n_samples, n_features = size(data)

    # Stop conditions
    if depth >= max_depth || n_samples < min_samples
        majority_class = get_majority_class(targets)
        push!(nodes, "$node_id LN -1 -1 -1 -1 $depth $majority_class")
        return nodes
    end

    # Select feature and threshold randomly
    feature = rand(0:(n_features-1))
    feature_values = data[:, feature+1]
    if length(unique(feature_values)) > 1
        threshold =
            rand() * (maximum(feature_values) - minimum(feature_values)) +
            minimum(feature_values)
    else
        threshold = feature_values[1]
    end

    # Split data
    left_mask = feature_values .<= threshold
    if !any(left_mask) || all(left_mask)
        majority_class = get_majority_class(targets)
        push!(nodes, "$node_id LN -1 -1 -1 -1 $depth $majority_class")
        return nodes
    end

    # Create child nodes
    left_child = node_id + 1
    left_subtree = create_tree_node(
        data[left_mask, :],
        targets[left_mask],
        depth + 1,
        max_depth,
        min_samples,
        left_child,
    )
    right_child = left_child + length(left_subtree)
    right_subtree = create_tree_node(
        data[.!left_mask, :],
        targets[.!left_mask],
        depth + 1,
        max_depth,
        min_samples,
        right_child,
    )

    # Add current node
    push!(
        nodes,
        "$node_id IN $left_child $right_child $feature $(round(threshold, digits=3)) $depth -1",
    )

    # Combine all nodes
    append!(nodes, left_subtree)
    append!(nodes, right_subtree)

    return nodes
end

"""
Create random forest input content.

# Arguments
- `dataset`: Input dataset
- `num_trees`: Number of trees to generate
- `max_depth`: Maximum tree depth

# Returns
- Formatted string for random forest input
"""
function create_random_forest_input(dataset, num_trees = 10, max_depth = 3)
    data, targets = prepare_dataset(dataset)
    n_features = size(data, 2)
    n_classes = length(unique(targets))

    lines = String[]
    push!(lines, "DATASET_NAME: dataset.train.csv")
    push!(lines, "ENSEMBLE: RF")
    push!(lines, "NB_TREES: $num_trees")
    push!(lines, "NB_FEATURES: $n_features")
    push!(lines, "NB_CLASSES: $n_classes")
    push!(lines, "MAX_TREE_DEPTH: $max_depth")
    push!(
        lines,
        "Format: node / node type (LN - leave node, IN - internal node) left child / right child / feature / threshold / node_depth / majority class (starts with index 0)",
    )
    push!(lines, "")

    for tree_idx = 1:(num_trees-1)
        sample_indices = rand(1:size(data, 1), size(data, 1))
        bootstrap_data = data[sample_indices, :]
        bootstrap_targets = targets[sample_indices]

        push!(lines, "[TREE $tree_idx]")
        tree_nodes = create_tree_node(bootstrap_data, bootstrap_targets, 0, max_depth)
        push!(lines, "NB_NODES: $(length(tree_nodes))")
        append!(lines, tree_nodes)
        push!(lines, "")
    end

    return join(lines, "\n")
end



"""
Utility function that calculates the maximum depth
of a SoleModels tree (Branch/ConstantModel).
"""
function compute_tree_depth(node, depth::Int=0)::Int
    if isa(node, ConstantModel{String})
        return depth
    elseif isa(node, Branch{String})
        left_d  = compute_tree_depth(node.posconsequent, depth+1)
        right_d = compute_tree_depth(node.negconsequent, depth+1)
        return max(left_d, right_d)
    else
        throw(ArgumentError("Unrecognized node type: $(typeof(node))"))
    end
end

"""
Given a node (Branch or ConstantModel), recursively visit in pre-order DFS,
creating BA-Trees lines in a Vector{String}.

If the node is a leaf (ConstantModel), print "LN", if it is a Branch, print "IN".
"""
function create_tree_node_from_branch(
    node::Union{Branch{String}, ConstantModel{String}},
    node_id::Ref{Int},       # counter for the next node ID
    depth::Int,
    class_map::Dict{String,Int},
    feature_offset::Int=1
)::Vector{String}

    lines = String[]

    # ID assigned to this node
    current_node = node_id[]
    node_id[] += 1  # increment for the next node

    if isa(node, ConstantModel{String})
        # Nodo foglia
        leaf = node::ConstantModel{String}
        leaf_class_idx = class_map[leaf.outcome]
        push!(lines, "$current_node LN -1 -1 -1 -1 $depth $leaf_class_idx")
        return lines
    else
        # Internal node
        branch = node::Branch{String}
        feature_raw = branch.antecedent.value.metacond.feature.i_variable
        threshold   = branch.antecedent.value.threshold

        # Shift from 1-based to 0-based (if needed)
        feature_0 = feature_raw - feature_offset
        threshold_str = round(threshold, digits=3)

        # Left subtree
        left_child_id = node_id[]
        left_lines = create_tree_node_from_branch(
            branch.posconsequent,
            node_id,
            depth+1,
            class_map,
            feature_offset
        )

        # Right subtree
        right_child_id = node_id[]
        right_lines = create_tree_node_from_branch(
            branch.negconsequent,
            node_id,
            depth+1,
            class_map,
            feature_offset
        )

        # Line for the internal node
        push!(lines,
              "$current_node IN $left_child_id $right_child_id $feature_0 $threshold_str $depth -1"
        )

        # Combine the lines of the two children
        append!(lines, left_lines)
        append!(lines, right_lines)

        return lines
    end
end

"""
Create a BA-Trees file similar to create_random_forest_input, but using
your model f::DecisionEnsemble{String} (SoleModels).

- Does not write the class map to the file (to avoid interfering with BA-Trees),
    prints it to the terminal with "println" so you can see it.
- Leaf classes are 0-based numbers (0 -> the first encountered class, 1 -> the second, etc.).
"""
function create_random_forest_input_from_model(
    f::DecisionEnsemble{String};
    output_file::String="forest.txt",
    feature_offset::Int=1
)::String

    # 1) Number of trees
    num_trees = length(f.models)

    # 2) Calculate NB_FEATURES
    max_feat_used = 0
    for tree in f.models
        if isa(tree, Branch{String})
            feat = tree.antecedent.value.metacond.feature.i_variable
            max_feat_used = max(max_feat_used, feat)
        end
    end
    nb_features = max_feat_used - (feature_offset - 1)

    # 3) NB_CLASSES
    all_labels = f.info.supporting_labels
    unique_labels = unique(all_labels)
    nb_classes = length(unique_labels)

    # 4) MAX_TREE_DEPTH
    max_depth = num_trees == 0 ? 0 :
        maximum(tree -> compute_tree_depth(tree), f.models)

    # 5) class_map
    class_map = Dict{String,Int}()
    for (i, cl) in enumerate(unique_labels)
        class_map[cl] = i - 1
    end

    println("\nClass mapping (0-based) per il modello:")
    for (lbl, idx) in sort(collect(class_map), by=x->x[2])
        println("   $idx -> $lbl")
    end
    println()

    # 6) File rows
    lines = String[]
    push!(lines, "DATASET_NAME: dataset.train.csv")
    push!(lines, "ENSEMBLE: RF")
    push!(lines, "NB_TREES: $num_trees")
    push!(lines, "NB_FEATURES: $nb_features")
    push!(lines, "NB_CLASSES: $nb_classes")
    push!(lines, "MAX_TREE_DEPTH: $max_depth")
    push!(lines, "Format: node / node type (LN - leaf node, IN - internal node) left child / right child / feature / threshold / node_depth / majority class (starts with index 0)")
    push!(lines, "")

    # 7) Convert any tree in the model
    for (i, tree) in enumerate(f.models)
        push!(lines, "[TREE $(i)]")

        node_id = Ref(0)
        tree_lines = create_tree_node_from_branch(
            tree,
            node_id,
            0,
            class_map,
            feature_offset
        )

        push!(lines, "NB_NODES: $(length(tree_lines))")
        append!(lines, tree_lines)
        push!(lines, "")
    end

    # 8) Save
    output_str = join(lines, "\n")
    open(output_file, "w") do io
        write(io, output_str)
    end

    return output_str
end

function convert_tree_structure(node, depth=0, node_counter=Ref(0), node_map=Dict())::Vector{String}
    lines = String[]
    current_node = node_counter[] += 1
    
    if isa(node, ConstantModel)
        # Leaf node
        push!(lines, "$current_node / LN / -1 / -1 / -1 / -1.0 / $depth / $(node.outcome)")
        node_map[node] = current_node
        return lines
    end
    
    # Internal node
    feature_idx = node.antecedent.value.metacond.feature.i_variable
    threshold = node.antecedent.value.threshold
    
    # Process children recursively
    left_subtree = convert_tree_structure(node.posconsequent, depth + 1, node_counter, node_map)
    right_subtree = convert_tree_structure(node.negconsequent, depth + 1, node_counter, node_map)
    
    # Get children node numbers
    left_child = node_map[node.posconsequent]
    right_child = node_map[node.negconsequent]
    
    # Add current node
    majority_class = get_majority_class(node.info.supporting_labels)
    push!(lines, "$current_node / IN / $left_child / $right_child / $feature_idx / $threshold / $depth / $majority_class")
    node_map[node] = current_node
    
    # Combine all nodes
    append!(lines, left_subtree)
    append!(lines, right_subtree)
    
    return lines
end

function get_majority_class(labels::Vector{String})::String
    counts = Dict{String, Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    return maximum(p -> p[2], collect(pairs(counts)))[1]
end

"""
Prepare and run BA-Trees algorithm.

# Arguments
- `dataset_name`: Name of dataset file
- `num_trees`: Number of trees to generate
- `max_depth`: Maximum tree depth
"""
function prepare_and_run_ba_trees_hardcoded(; dataset_name = "iris", num_trees = 3, max_depth = 3)
    # <-- PATH FIX: usa @__DIR__ invece di pwd()
    base_dir = @__DIR__
    # Se il CSV è nella stessa cartella di questo file:
    dataset_path = joinpath(base_dir, dataset_name * ".csv")
    @show dataset_path

    if !isfile(dataset_path)
        println("ERROR: Dataset file not found: $dataset_path")
        println("Current directory: ", pwd())
        println("Content of src directory:")
        if isdir("src")
            for file in readdir("src")
                println(" - $file")
            end
        else
            println("The 'src' directory doesn't exist!")
        end
        return
    end

    dataset = DataFrame(CSV.File(dataset_path))

    # <-- PATH FIX: cartella temporanea creata accanto a main.jl
    temp_dir = joinpath(base_dir, "temp_ba_trees")
    isdir(temp_dir) || mkdir(temp_dir)

    input_file = joinpath(temp_dir, "forest.txt")
    output_base = joinpath(temp_dir, "result.txt")
    output_stats = output_base * ".out"
    output_tree = output_base * ".tree"

    # <-- PATH FIX: se bornAgain è nella stessa cartella di main.jl
    executable_path = joinpath(base_dir, "bornAgain")

    try
        println("Generating random forest...")
        f, model, start_time = learn_and_convert(num_trees, "iris", 3)
        forest_content = create_random_forest_input_from_model(f)   
        write(input_file, forest_content)

        cmd = `$executable_path $input_file $output_base -trees $num_trees -obj 4`
        println("Executing command: ", cmd)
        run(cmd)

        if isfile(output_stats)
            println("\nBorn-Again Tree Analysis Results:")
            println(read(output_stats, String))
            
            if isfile(output_tree)
                println("\nBorn-Again Tree Structure:")
                println(read(output_tree, String))
            else
                println("Tree structure file not found: $output_tree")
            end
        else
            println("Statistics file not found: $output_stats")
        end

    catch e
        println("\nError during execution:")
        println(e)
        println("\nDebug information:")
        println("Current directory: ", pwd())
        println("Executable exists? ", isfile(executable_path))
        println("Temporary directory contents:")
        for (root, dirs, files) in walkdir(temp_dir)
            println("Files: ", join(files, ", "))
        end
        if isfile(input_file)
            println("\nFirst 10 lines of input file:")
            println(join(collect(Iterators.take(eachline(input_file), 10)), "\n"))
        end
    finally
        if !@isdefined(e)
            #rm(temp_dir, recursive = true, force = true) 
        end
    end
end

# Main execution
function WRAP_batrees_hardcoded()
    println("Born-Again Tree Analysis")
    println("="^30)

    # Configurable parameters
    dataset_name = "iris"
    num_trees = 3
    max_depth = 3

    println("Dataset: $dataset_name")
    println("Number of trees: $num_trees")
    println("Maximum depth: $max_depth")
    println("="^30)

    # Run analysis
    prepare_and_run_ba_trees_hardcoded(
        dataset_name = dataset_name,
        num_trees = num_trees,
        max_depth = max_depth,
    )
end

"""
Display BA-Trees analysis results.

# Arguments
- `output_file`: Path to output file
"""
function display_results(output_file)
    output_file = endswith(output_file, ".out") ? output_file : output_file * ".out"
    
    if isfile(output_file)
        println("\nBorn-Again Tree Analysis:")
        println("-"^40)
        for line in eachline(output_file)
            if startswith(line, "TIME") ||
               startswith(line, "DEPTH") ||
               startswith(line, "LEAVES") ||
               startswith(line, "ACCURACY")
                parts = split(line, ":")
                if length(parts) == 2
                    key = rpad(strip(parts[1]), 20)
                    value = strip(parts[2])
                    println("$key: $value")
                end
            end
        end
        println("-"^40)
    end
end

using DecisionTree
using DataFrames

"""
Node structure for BA-Trees format
"""
struct BANode
    id::Int
    node_type::String  # "IN" or "LN"
    left_child::Int
    right_child::Int
    feature::Int
    threshold::Float64
    depth::Int
    class::Int
end

"""
Parse a BA-Trees node line into a BANode struct

# Arguments
- `line::String`: A line from the BA-Trees output file

# Returns
- `BANode`: Structured representation of the node
"""
function parse_ba_node(line::String)
    parts = split(line)
    return BANode(
        parse(Int, parts[1]),     # id
        parts[2],                 # node_type
        parse(Int, parts[3]),     # left_child
        parse(Int, parts[4]),     # right_child
        parse(Int, parts[5]),     # feature
        parse(Float64, parts[6]), # threshold
        parse(Int, parts[7]),     # depth
        parse(Int, parts[8])      # class
    )
end

"""
Convert BA-Trees node to DecisionTree.jl node recursively

# Arguments
- `nodes::Dict{Int,BANode}`: Dictionary of all BA nodes indexed by ID
- `current_id::Int`: ID of the current node to convert
- `features::Vector{String}`: Feature names for the tree

# Returns
- `Node`: DecisionTree.jl node structure
"""
function convert_to_dt_node(nodes::Dict{Int,BANode}, current_id::Int, features::Vector{String})
    node = nodes[current_id]
    if node.node_type == "LN"
        # Leaf node
        return Leaf(node.class)
    else
        # Internal node
        feature_name = features[node.feature + 1]  # BA-Trees uses 0-based indexing
        left = convert_to_dt_node(nodes, node.left_child, features)
        right = convert_to_dt_node(nodes, node.right_child, features)
        return Node(feature_name, node.threshold, left, right)
    end
end

"""
Convert BA-Trees format to DecisionTree.jl tree

# Arguments
- `ba_tree_file::String`: Path to the BA-Trees output file
- `features::Vector{String}`: Names of features in the dataset
- `classes::Vector`: Possible class values

# Returns
- `DecisionTreeClassifier`: DecisionTree.jl compatible tree
"""
function convert_ba_to_dt(ba_tree_file::String, features::Vector{String}, classes::Vector)
    # Read and parse BA-Trees nodes
    nodes = Dict{Int,BANode}()
    open(ba_tree_file) do f
        for line in eachline(f)
            # Skip non-node lines
            if !occursin(r"^\d+\s+(IN|LN)", line)
                continue
            end
            node = parse_ba_node(line)
            nodes[node.id] = node
        end
    end

    # Find root node (ID 0)
    root = convert_to_dt_node(nodes, 0, features)
    
    # Create DecisionTreeClassifier
    n_classes = length(classes)
    dt = DecisionTreeClassifier(
        root,
        features,
        classes,
        n_classes
    )
    
    return dt
end

"""
Example usage function to demonstrate the conversion

# Arguments
- `ba_tree_file::String`: Path to the BA-Trees output file
- `dataset::DataFrame`: Original dataset used to create the tree
"""
function demonstrate_conversion(ba_tree_file::String, dataset::DataFrame)
    # Get feature names and classes from dataset
    features = names(dataset)[1:end-1] |> Vector{String}
    classes = unique(dataset[!, end]) |> Vector
    
    # Convert the tree
    dt = convert_ba_to_dt(ba_tree_file, features, classes)
    
    # Print some information about the converted tree
    println("Converted Decision Tree:")
    println("Features: ", dt.features)
    println("Classes: ", dt.classes)
    
    # You can now use the tree for predictions
    # Example: predict on first row of dataset
    x = Array(dataset[1, 1:end-1])
    prediction = DecisionTree.predict(dt, x)
    println("\nPrediction for first row: ", prediction)
    println("Actual class: ", dataset[1, end])
    
    return dt
end

function prepare_and_run_ba_trees(; 
    dataset_name::String = "iris", 
    num_trees::Int = 10, 
    max_depth::Int = 3, 
    forest = nothing
 )
    # <-- PATH FIX: usa @__DIR__ invece di pwd()
    base_dir = @__DIR__
    dataset_path = joinpath(base_dir, dataset_name * ".csv")
    @show dataset_path
 
    # Check if born_again_db exists
    db_path = joinpath(base_dir, "born_again_db")
    if !isdir(db_path)
        println("""
        ERROR: The 'born_again_db' directory is missing!
        Please follow these steps:
        1. Clone the BA-Trees repository: git clone https://github.com/vidalt/BA-Trees.git
        2. Follow the installation instructions in the README
        3. Make sure the 'born_again_db' directory is present in the same folder as this script
        
        For more details, visit: https://github.com/vidalt/BA-Trees
        """)
        return
    end
 
    # Check if bornAgain executable exists
    executable_path = joinpath(base_dir, "bornAgain")
    if !isfile(executable_path)
        println("""
        ERROR: The 'bornAgain' executable is missing!
        Please follow these steps:
        1. If you haven't already, clone the BA-Trees repository: git clone https://github.com/vidalt/BA-Trees.git
        2. Follow the compilation instructions in the README to build the executable
        3. Make sure the 'bornAgain' executable is present in the same folder as this script
        
        For more details, visit: https://github.com/vidalt/BA-Trees
        """)
        return
    end
 
    # Check if the dataset exists
    if !isfile(dataset_path)
        println("ERROR: Dataset file not found: $dataset_path")
        println("Current directory: ", pwd())
        println("Content of src directory:")
        if isdir("src")
            for file in readdir("src")
                println(" - $file")
            end
        else
            println("The 'src' directory doesn't exist!")
        end
        return
    end
 
    # Carica il dataset
    dataset = DataFrame(CSV.File(dataset_path))
 
    # <-- PATH FIX: cartella temporanea affiancata a main.jl
    temp_dir = joinpath(base_dir, "temp_ba_trees")
    isdir(temp_dir) || mkdir(temp_dir)
 
    input_file   = joinpath(temp_dir, "forest.txt")
    output_base  = joinpath(temp_dir, "result.txt")
    output_stats = output_base * ".out"
    output_tree  = output_base * ".tree"
 
    # <-- PATH FIX: if "bornAgain" is in the same folder
    executable_path = joinpath(base_dir, "bornAgain")
 
    try
        println("Preparing random forest data...")
 
        forest_content = ""
        if forest === nothing
            forest, model, start_time = learn_and_convert(num_trees, "iris", max_depth)
            forest_content = create_random_forest_input_from_model(forest)
        else
            forest_content = create_random_forest_input_from_model(forest)
        end
 
        # Write the content to file
        write(input_file, forest_content)
 
        # number of tree you are interesad
        num_trees = length(forest.models)
        # Comando per eseguire l'analisi Born-Again
        #=
        Usage:
            ./bornAgain input_ensemble_path output_BAtree_path [list of options]
              Available options:
            -obj X	       Objective used in the algorithm: 0 = Depth ; 1 = NbLeaves ; 2 = Depth then NbLeaves ; 4 = Heuristic BA-Tree (defaults to 4)
            -trees X      Limits the number of trees read by the algorithm from the input file (X in 3 to 10, defaults to 10)
            -seed X       Defines the random seed (defaults to 1)
 
            TODO more parameterization
        =#
        cmd = `$executable_path $input_file $output_base -trees $num_trees -obj 0`
        println("Executing command: ", cmd)
        run(cmd)
 
        # Read results
        if isfile(output_stats)
            println("\nBorn-Again Tree Analysis Results:")
            println(read(output_stats, String))
            
            if isfile(output_tree)
                println("\nBorn-Again Tree Structure:")
                println(read(output_tree, String))
            else
                println("Tree structure file not found: $output_tree")
            end
        else
            println("Statistics file not found: $output_stats")
        end
 
    catch e
        println("\nError during execution:")
        println(e)
        println("\nDebug information:")
        println("Current directory: ", pwd())
        println("Executable exists? ", isfile(executable_path))
        println("Temporary directory contents:")
        for (root, dirs, files) in walkdir(temp_dir)
            println("Files: ", join(files, ", "))
        end
        if isfile(input_file)
            println("\nFirst 10 lines of input file:")
            println(join(collect(Iterators.take(eachline(input_file), 10)), "\n"))
        end
    finally
        if !@isdefined(e)
            # Example: final cleanup (if desired)
            # rm(temp_dir, recursive = true, force = true) its better if not do that for now 
        end
    end
 end

function WRAP_batrees(f, max_depth=10; dataset_name="iris", num_trees=10)
    println("Born-Again Tree Analysis")
    println("="^30)

    println("Dataset: $dataset_name")
    println("Number of trees: $num_trees")
    println("Maximum depth: $max_depth")
    println("="^30)

    # is indifferent if i have to create a new f ora f is passed
    if (isnothing(f))
        prepare_and_run_ba_trees(
            dataset_name = dataset_name,
            num_trees    = num_trees,
            max_depth    = max_depth,
            forest = nothing,
        )
    else
        prepare_and_run_ba_trees(
            dataset_name = dataset_name, # is indifferent
            num_trees    = length(f.models),
            max_depth    = max_depth,    # is indifferent ? 
            forest       = f,
        )
    end
end

end


