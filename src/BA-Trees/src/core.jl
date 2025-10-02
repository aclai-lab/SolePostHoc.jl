module RunModule

using DelimitedFiles
using Random
using CSV
using DataFrames
using Statistics
using CategoricalArrays

const TreeType = Union{String, CategoricalArrays.CategoricalValue{String, UInt32}}

using DecisionTree
using SoleBase: Label
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
    if isa(node, ConstantModel{<:TreeType})
        return depth
    elseif isa(node, Branch{<:TreeType})
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
    node::Union{Branch, ConstantModel{<:TreeType}},
    node_id::Ref{Int},       # counter for the next node ID
    depth::Int,
    class_map::Dict{String,Int},
    feature_offset::Int=1
)::Vector{String}

    lines = String[]

    # ID assigned to this node
    current_node = node_id[]
    node_id[] += 1  # increment for the next node

    if isa(node, ConstantModel{<:TreeType})
        # Nodo foglia
        leaf = node::ConstantModel{<:TreeType}
        leaf_class_idx = class_map[leaf.outcome]
        push!(lines, "$current_node LN -1 -1 -1 -1 $depth $leaf_class_idx")
        return lines
    else
        # Internal node
        branch = node::Branch{<:TreeType}
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


function get_number_of_features(models::Vector{<:Branch{<:Label}})
    function get_feature_from_model(model::ConstantModel; maxvalue=0)
        return maxvalue
    end
    function get_feature_from_model(model::Branch{<:Label}; maxvalue=0)
        return maximum([
            get_feature_from_model(model.posconsequent, maxvalue=maxvalue), 
            get_feature_from_model(model.negconsequent, maxvalue=maxvalue), 
            get_feature_from_model(model.antecedent, maxvalue=maxvalue)
        ]);
    end
    function get_feature_from_model(model::Atom; maxvalue=0)
        return max(maxvalue, model.value.metacond.feature.i_variable)
    end

    return maximum([get_feature_from_model(model) for model in models])
end

"""
Create a BA-Trees file similar to create_random_forest_input, but using
your model f::DecisionEnsemble{String} (SoleModels).

- Does not write the class map to the file (to avoid interfering with BA-Trees),
    prints it to the terminal with "println" so you can see it.
- Leaf classes are 0-based numbers (0 -> the first encountered class, 1 -> the second, etc.).
"""
function create_random_forest_input_from_model(
    f::Union{DecisionEnsemble, SoleModels.DecisionXGBoost};
    output_file::String="forest.txt",
    feature_offset::Int=1
)
    class_map = nothing
    # 1) Number of trees
    num_trees = length(f.models)

    # 2) Calculate NB_FEATURES
    max_feat_used = get_number_of_features(f.models)
    nb_features = max_feat_used - feature_offset + 1
    
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
    println()
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

    return output_str,class_map
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
    base_dir = @__DIR__

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

"""
Creates a simplified Makefile that works cross-platform without CPLEX dependencies.
"""
function create_simple_makefile(born_again_dp_path::String)

    makefile_content = """
CXX = g++
CXXFLAGS = -O2 -std=c++11

SRCS = main.cpp BornAgainDecisionTree.cpp FSpace.cpp RandomForest.cpp
OBJS = \$(SRCS:.cpp=.o)
TARGET = bornAgain

all: \$(TARGET)

\$(TARGET): \$(OBJS)
\t\$(CXX) \$(CXXFLAGS) -o \$@ \$(OBJS)

%.o: %.cpp
\t\$(CXX) \$(CXXFLAGS) -c \$< -o \$@

clean:
\trm -f \$(OBJS) \$(TARGET)

.PHONY: all clean
"""
    
    makefile_path = joinpath(born_again_dp_path, "Makefile.simple")
    write(makefile_path, makefile_content)
    println("✓ Created simplified Makefile at: $makefile_path")
    return makefile_path
end


"""
Automatically clones, compiles and sets up BA-Trees if not present.
Returns true if setup is successful, false otherwise.
"""
function setup_ba_trees(base_dir::String)
    
    born_again_dp_path = joinpath(base_dir, "born_again_dp")
    executable_path = joinpath(base_dir, "bornAgain")
    
    # Check if born_again_dp directory exists, if not clone it
    if !isdir(born_again_dp_path)
        println("BA-Trees repository not found. Cloning...")
        try
            # Clone the repository
            repo_url = "https://github.com/vidalt/BA-Trees.git"
            temp_clone_dir = joinpath(base_dir, "BA-Trees_temp")
            
            run(`git clone $repo_url $temp_clone_dir`)
            
            # Move born_again_dp directory to the correct location
            cloned_dp_path = joinpath(temp_clone_dir, "born_again_dp")
            if isdir(cloned_dp_path)
                mv(cloned_dp_path, born_again_dp_path)
                println("✓ BA-Trees repository cloned successfully")
            else
                println("ERROR: born_again_dp directory not found in cloned repository")
                rm(temp_clone_dir, recursive=true, force=true)
                return false
            end
            
            # Clean up temporary clone directory
            rm(temp_clone_dir, recursive=true, force=true)
            
        catch e
            println("ERROR: Failed to clone BA-Trees repository")
            println("Error details: ", e)
            println("Please ensure git is installed and you have internet connection")
            return false
        end
    end
    
    # Check if executable exists, if not compile it
    if !isfile(executable_path)
        println("bornAgain executable not found. Compiling...")
        try
            # Detect OS
            os_type = Sys.iswindows() ? "windows" : (Sys.isapple() ? "macos" : "linux")
            println("Detected OS: $os_type")
            
            # Change to born_again_dp directory for compilation
            cd(born_again_dp_path) do
                # Create simplified Makefile for reliable compilation
                simple_makefile = create_simple_makefile(born_again_dp_path)
                
                # Check if original makefile exists (backup for reference)
                original_makefile = nothing
                if isfile("makefile")
                    original_makefile = "makefile"
                elseif isfile("Makefile")
                    original_makefile = "Makefile"
                end
                
                if original_makefile !== nothing
                    println("ℹ Original Makefile found: $original_makefile (keeping as backup)")
                    println("ℹ Using simplified Makefile for compilation")
                end
                
                # Try compilation with simplified Makefile first
                compilation_success = false
                
                println("Running make with simplified Makefile...")
                try
                    if Sys.iswindows()
                        # For Windows, try different make commands
                        try
                            run(`make -f Makefile.simple`)
                            compilation_success = true
                        catch
                            println("Trying mingw32-make...")
                            run(`mingw32-make -f Makefile.simple`)
                            compilation_success = true
                        end
                    else
                        # For Unix-like systems (Linux/Mac)
                        run(`make -f Makefile.simple`)
                        compilation_success = true
                    end
                catch e
                    println("⚠ Simplified Makefile compilation failed: ", e)
                    
                    # Fallback: try original Makefile without CPLEX
                    if original_makefile !== nothing
                        println("Attempting compilation with original Makefile (without CPLEX)...")
                        try
                            if Sys.iswindows()
                                try
                                    run(`make -f $original_makefile`)
                                    compilation_success = true
                                catch
                                    run(`mingw32-make -f $original_makefile`)
                                    compilation_success = true
                                end
                            else
                                run(`make -f $original_makefile`)
                                compilation_success = true
                            end
                        catch e2
                            println("⚠ Original Makefile compilation also failed: ", e2)
                        end
                    end
                end
                
                if !compilation_success
                    println("ERROR: All compilation attempts failed")
                    println("\nPlease check:")
                    println("1. C++ compiler is installed (g++ or clang)")
                    println("2. Make is installed")
                    println("3. All source files are present")
                    return false
                end
                
                # Check if executable was created
                possible_names = ["bornAgain", "bornAgain.exe", "main", "main.exe"]
                compiled_exec = nothing
                
                for name in possible_names
                    if isfile(name)
                        compiled_exec = name
                        println("✓ Found compiled executable: $name")
                        break
                    end
                end
                
                if compiled_exec === nothing
                    println("ERROR: Compilation completed but executable not found")
                    println("Looking for: ", join(possible_names, ", "))
                    println("Directory contents:")
                    for file in readdir()
                        println("  ", file)
                    end
                    return false
                end
                
                # Move executable to base directory
                target_name = Sys.iswindows() ? "bornAgain.exe" : "bornAgain"
                target_path = joinpath(base_dir, target_name)
                cp(compiled_exec, target_path, force=true)
                
                # Make executable on Unix systems
                if !Sys.iswindows()
                    run(`chmod +x $target_path`)
                end
                
                println("✓ bornAgain compiled successfully and moved to: $target_path")
            end
            
        catch e
            println("ERROR: Failed to compile bornAgain")
            println("Error details: ", e)
            println("\nPlease check:")
            println("1. C++ compiler is installed (g++ or clang)")
            println("2. Make is installed")
            println("3. All dependencies are available")
            return false
        end
    else
        println("✓ bornAgain executable already exists")
    end
    
    println("✓ BA-Trees setup complete")
    return true
end


function prepare_and_run_ba_trees(; 
    dataset_name::String = "", 
    num_trees::Int = nothing, 
    max_depth::Int = -1, 
    forest = nothing,
    mode_obj = 0
)
    class_map = nothing
    
    # Check if mode_obj is valid
    if !(mode_obj in [0, 1, 2, 4])
        println("""
        ERROR: Invalid mode_obj value: $mode_obj
        Valid options for -obj parameter are:
        0 = Depth
        1 = NbLeaves
        2 = Depth then NbLeaves
        4 = Heuristic BA-Tree
        """)
        return
    end

    base_dir = @__DIR__
    
    # Automatic setup: clone and compile if needed
    println("Checking BA-Trees setup...")
    if !setup_ba_trees(base_dir)
        println("ERROR: Failed to setup BA-Trees. Please check the error messages above.")
        return
    end
    
    # Verify paths after setup
    db_path = joinpath(base_dir, "born_again_dp")
    executable_name = Sys.iswindows() ? "bornAgain.exe" : "bornAgain"
    executable_path = joinpath(base_dir, executable_name)
    
    if !isdir(db_path) || !isfile(executable_path)
        println("ERROR: Setup verification failed")
        println("born_again_dp exists: ", isdir(db_path))
        println("bornAgain exists: ", isfile(executable_path))
        return
    end
 
    # Create temporary directory
    temp_dir = joinpath(base_dir, "temp_ba_trees")
    isdir(temp_dir) || mkdir(temp_dir)
 
    input_file   = joinpath(temp_dir, "forest.txt")
    output_base  = joinpath(temp_dir, "result.txt")
    output_stats = output_base * ".out"
    output_tree  = output_base * ".tree"
 
    try
        println("Preparing random forest data...")
 
        forest_content = ""
        if forest === nothing
            if num_trees <= 0 || dataset_name == ""
                println("ERROR: Invalid parameters for random forest generation.")
                return
            end
            forest, model, start_time = learn_and_convert(num_trees, dataset_name, max_depth)
            forest_content, class_map = create_random_forest_input_from_model(forest)
        else
            forest_content, class_map = create_random_forest_input_from_model(forest)
        end
 
        # Write the content to file
        write(input_file, forest_content)
 
        # Number of trees
        num_trees = length(forest.models)
        
        # Run Born-Again
        cmd = `$executable_path $input_file $output_base -trees $num_trees -obj $mode_obj`
        println("Executing command: ", cmd)
        run(cmd)
 
        # Read results
        if isfile(output_stats)
            println("\n" * "="^60)
            println("Born-Again Tree Analysis Results:")
            println("="^60)
            println(read(output_stats, String))
            
            if isfile(output_tree)
                println("\n" * "="^60)
                println("Born-Again Tree Structure:")
                println("="^60)
                println(read(output_tree, String))
            else
                println("⚠ Tree structure file not found: $output_tree")
            end
        else
            println("⚠ Statistics file not found: $output_stats")
        end
 
    catch e
        println("\n" * "="^60)
        println("ERROR during execution:")
        println("="^60)
        println(e)
        println("\nDebug information:")
        println("Current directory: ", pwd())
        println("Executable exists? ", isfile(executable_path))
        println("Executable path: ", executable_path)
        println("\nTemporary directory contents:")
        if isdir(temp_dir)
            for (root, dirs, files) in walkdir(temp_dir)
                println("Directory: ", root)
                println("Files: ", join(files, ", "))
            end
        end
        if isfile(input_file)
            println("\nFirst 10 lines of input file:")
            println(join(collect(Iterators.take(eachline(input_file), 10)), "\n"))
        end
        rethrow(e)
    finally
        if !@isdefined(e)
            return class_map
            # Optional cleanup (commented out for debugging)
            # rm(temp_dir, recursive=true, force=true)
        end
    end
end

function WRAP_batrees(f, max_depth=10; dataset_name="", num_trees=nothing, mod = 0)
    println("Born-Again Tree Analysis")
    println("="^30)

    println("Dataset: $dataset_name")
    println("Number of trees: $num_trees")
    println("Maximum depth: $max_depth")
    println("="^30)

    # is indifferent if i have to create a new f ora f is passed
    if (isnothing(f))
        class_map = prepare_and_run_ba_trees(
            dataset_name = dataset_name,
            num_trees    = num_trees,
            max_depth    = max_depth,
            forest = nothing,
            mode_obj = mod
        )
    else
        class_map = prepare_and_run_ba_trees(
            num_trees    = length(f.models),
            forest       = f,
            mode_obj = mod
        )
    end
    return class_map
end

end