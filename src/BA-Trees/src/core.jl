"""
RunModule - Born-Again Trees Integration Module

This module provides functionality to convert machine learning models (particularly
decision ensembles from SoleModels) into BA-Trees format and execute the Born-Again
algorithm for model compression and analysis.

Key Features:
- Automatic BA-Trees repository setup (clone and compile)
- Model-to-BA-Trees format conversion
- Cross-platform support (Linux, macOS, Windows)
- Flexible tree generation and analysis
"""
module RunModule

# ============================================================================
# Dependencies
# ============================================================================

using DelimitedFiles
using Random
using CSV
using DataFrames
using Statistics
using CategoricalArrays
using DecisionTree
using SoleBase: Label
using SoleModels

# ============================================================================
# Type Definitions
# ============================================================================

const TreeType = Union{String, CategoricalArrays.CategoricalValue{String, UInt32}}

"""
    BANode

Structure representing a node in BA-Trees format.

# Fields
- `id::Int`: Unique node identifier
- `node_type::String`: Node type ("IN" for internal, "LN" for leaf)
- `left_child::Int`: ID of left child (-1 if leaf)
- `right_child::Int`: ID of right child (-1 if leaf)
- `feature::Int`: Feature index for splitting (-1 if leaf)
- `threshold::Float64`: Split threshold (-1.0 if leaf)
- `depth::Int`: Depth of node in tree
- `class::Int`: Majority class for the node
"""
struct BANode
    id::Int
    node_type::String
    left_child::Int
    right_child::Int
    feature::Int
    threshold::Float64
    depth::Int
    class::Int
end

# ============================================================================
# Utility Functions - Path Management
# ============================================================================

include(joinpath(@__DIR__, "utils/minor.jl"))

# ============================================================================
# Utility Functions - Statistics
# ============================================================================

"""
    mode(x::AbstractArray) -> eltype(x)

Calculate the mode (most frequent value) of an array.

Returns the first modal value in case of ties.

# Throws
- `ArgumentError`: If the input array is empty.

# Examples
```julia
mode([1, 2, 2, 3, 3, 3]) # Returns 3
mode(["a", "b", "a", "c"]) # Returns "a"
```
"""
function mode(x::AbstractArray)
    isempty(x) && throw(ArgumentError("array must be non-empty"))
    
    counts = Dict{eltype(x), Int}()
    for val in x
        counts[val] = get(counts, val, 0) + 1
    end
    
    max_count = maximum(values(counts))
    modes = [k for (k, v) in counts if v == max_count]
    
    return first(modes)
end

"""
    get_majority_class(targets::AbstractArray) -> Int

Determine the majority class from target values.

Returns the most frequent class as an integer, or 0 if the array is empty.

# Examples
```julia
get_majority_class([1, 2, 2, 3, 2]) # Returns 2
```
"""
function get_majority_class(targets::AbstractArray)
    isempty(targets) && return 0
    return mode(convert(Vector{Int}, targets))
end

"""
    get_majority_class(labels::Vector{String}) -> String

Determine the majority class from string labels.

# Examples
```julia
get_majority_class(["cat", "dog", "cat", "cat"]) # Returns "cat"
```
"""
function get_majority_class(labels::Vector{String})
    counts = Dict{String, Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    return maximum(p -> p[2], collect(pairs(counts)))[1]
end

# ============================================================================
# Tree Structure Analysis
# ============================================================================

"""
    compute_tree_depth(node, depth::Int=0) -> Int

Recursively compute the maximum depth of a SoleModels tree.

# Arguments
- `node`: Tree node (Branch or ConstantModel)
- `depth::Int=0`: Current depth in the tree

# Returns
- Maximum depth of the tree
"""
function compute_tree_depth(node, depth::Int=0)::Int
    if isa(node, ConstantModel{<:TreeType})
        return depth
    elseif isa(node, Branch{<:TreeType})
        left_depth = compute_tree_depth(node.posconsequent, depth + 1)
        right_depth = compute_tree_depth(node.negconsequent, depth + 1)
        return max(left_depth, right_depth)
    else
        throw(ArgumentError("Unrecognized node type: $(typeof(node))"))
    end
end

"""
    get_number_of_features(models::Vector{<:Branch{<:Label}}) -> Int

Determine the maximum feature index used across all models in an ensemble.

This function recursively traverses the tree structure to find the highest
feature index referenced, accounting for the feature offset.

# Arguments
- `models::Vector{<:Branch{<:Label}}`: Collection of tree models

# Returns
- Maximum feature index found in the ensemble
"""
function get_number_of_features(models::Vector{<:Branch{<:Label}})
    function get_feature_from_model(model::ConstantModel; maxvalue=0)
        return maxvalue
    end
    
    function get_feature_from_model(model::Branch{<:Label}; maxvalue=0)
        return maximum([
            get_feature_from_model(model.posconsequent, maxvalue=maxvalue), 
            get_feature_from_model(model.negconsequent, maxvalue=maxvalue), 
            get_feature_from_model(model.antecedent, maxvalue=maxvalue)
        ])
    end
    
    function get_feature_from_model(model::Atom; maxvalue=0)
        return max(maxvalue, model.value.metacond.feature.i_variable)
    end

    return maximum([get_feature_from_model(model) for model in models])
end

# ============================================================================
# BA-Trees Format Conversion
# ============================================================================

"""
    create_tree_node_from_branch(
        node::Union{Branch, ConstantModel{<:TreeType}},
        node_id::Ref{Int},
        depth::Int,
        class_map::Dict{String,Int},
        feature_offset::Int=1
    ) -> Vector{String}

Recursively convert a SoleModels tree node to BA-Trees format.

Performs pre-order depth-first traversal, generating BA-Trees formatted strings
for each node. Leaf nodes (ConstantModel) are marked as "LN", while internal
nodes (Branch) are marked as "IN".

# Arguments
- `node`: Current tree node to convert
- `node_id::Ref{Int}`: Mutable counter for node IDs
- `depth::Int`: Current depth in the tree
- `class_map::Dict{String,Int}`: Mapping from class labels to 0-based indices
- `feature_offset::Int=1`: Offset to convert 1-based to 0-based indexing

# Returns
- Vector of strings representing nodes in BA-Trees format
"""
function create_tree_node_from_branch(
    node::Union{Branch, ConstantModel{<:TreeType}},
    node_id::Ref{Int},
    depth::Int,
    class_map::Dict{String,Int},
    feature_offset::Int=1
)::Vector{String}
    lines = String[]
    current_node = node_id[]
    node_id[] += 1

    if isa(node, ConstantModel{<:TreeType})
        # Leaf node: format = "ID LN -1 -1 -1 -1 depth class"
        leaf = node::ConstantModel{<:TreeType}
        leaf_class_idx = class_map[leaf.outcome]
        push!(lines, "$current_node LN -1 -1 -1 -1 $depth $leaf_class_idx")
        return lines
    else
        # Internal node: format = "ID IN left_child right_child feature threshold depth -1"
        branch = node::Branch{<:TreeType}
        feature_raw = branch.antecedent.value.metacond.feature.i_variable
        threshold = branch.antecedent.value.threshold

        # Convert from 1-based to 0-based indexing
        feature_0 = feature_raw - feature_offset
        threshold_str = round(threshold, digits=3)

        # Process left subtree (positive consequence)
        left_child_id = node_id[]
        left_lines = create_tree_node_from_branch(
            branch.posconsequent,
            node_id,
            depth + 1,
            class_map,
            feature_offset
        )

        # Process right subtree (negative consequence)
        right_child_id = node_id[]
        right_lines = create_tree_node_from_branch(
            branch.negconsequent,
            node_id,
            depth + 1,
            class_map,
            feature_offset
        )

        # Add current internal node
        push!(
            lines,
            "$current_node IN $left_child_id $right_child_id $feature_0 $threshold_str $depth -1"
        )

        # Append child nodes
        append!(lines, left_lines)
        append!(lines, right_lines)

        return lines
    end
end

"""
    create_random_forest_input_from_model(
        f::Union{DecisionEnsemble, SoleModels.DecisionXGBoost};
        output_file::String="forest.txt",
        feature_offset::Int=1
    ) -> Tuple{String, Dict{String,Int}}

Convert a SoleModels ensemble to BA-Trees input format.

Generates a formatted text file compatible with the BA-Trees algorithm,
including metadata about the ensemble and detailed tree structures.

# Arguments
- `f`: Decision ensemble or XGBoost model from SoleModels
- `output_file::String="forest.txt"`: Path for the output file
- `feature_offset::Int=1`: Offset for feature indexing conversion

# Returns
- Tuple of (formatted_string, class_mapping)

# Notes
- The class mapping is also printed to stdout for reference
- Leaf classes use 0-based indexing
"""
function create_random_forest_input_from_model(
    f::Union{DecisionEnsemble, SoleModels.DecisionXGBoost};
    output_file::String="forest.txt",
    feature_offset::Int=1
)
    # Extract ensemble metadata
    num_trees = length(f.models)
    max_feat_used = get_number_of_features(f.models)
    nb_features = max_feat_used - feature_offset + 1
    
    # Determine unique classes
    all_labels = f.info.supporting_labels
    unique_labels = unique(all_labels)
    nb_classes = length(unique_labels)

    # Calculate maximum tree depth
    max_depth = num_trees == 0 ? 0 :
        maximum(tree -> compute_tree_depth(tree), f.models)

    # Create 0-based class mapping
    class_map = Dict{String,Int}()
    for (i, cl) in enumerate(unique_labels)
        class_map[cl] = i - 1
    end

    # Display class mapping for user reference
    println("\n" * "="^60)
    println("Class Mapping (0-based indexing):")
    println("="^60)
    for (lbl, idx) in sort(collect(class_map), by=x->x[2])
        println("  $idx => $lbl")
    end
    println("="^60 * "\n")

    # Build BA-Trees format output
    lines = String[]
    push!(lines, "DATASET_NAME: dataset.train.csv")
    push!(lines, "ENSEMBLE: RF")
    push!(lines, "NB_TREES: $num_trees")
    push!(lines, "NB_FEATURES: $nb_features")
    push!(lines, "NB_CLASSES: $nb_classes")
    push!(lines, "MAX_TREE_DEPTH: $max_depth")
    push!(lines, "Format: node / node type (LN - leaf node, IN - internal node) left child / right child / feature / threshold / node_depth / majority class (starts with index 0)")
    push!(lines, "")

    # Convert each tree in the ensemble
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

    # Write to file
    output_str = join(lines, "\n")
    open(output_file, "w") do io
        write(io, output_str)
    end

    return output_str, class_map
end

# ============================================================================
# BA-Trees Setup and Compilation
# ============================================================================

"""
    create_simple_makefile(born_again_dp_path::String) -> String

Generate a simplified, cross-platform Makefile without CPLEX dependencies.

This Makefile is designed to work across different platforms and doesn't
require the CPLEX optimization library, making it more portable.

# Arguments
- `born_again_dp_path::String`: Path to the born_again_dp directory

# Returns
- Path to the created Makefile
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
    println("Created simplified Makefile at: $makefile_path")
    return makefile_path
end

"""
    setup_ba_trees(base_dir::String) -> Bool

Automatically clone, compile, and set up BA-Trees if not present.

This function handles the complete setup process:
1. Checks if BA-Trees repository exists
2. Clones the repository if needed
3. Compiles the executable if not present
4. Handles cross-platform differences

# Arguments
- `base_dir::String`: Base directory for BA-Trees installation

# Returns
- `true` if setup is successful, `false` otherwise

# Notes
- Requires git, g++ (or clang), and make to be installed
- Supports Linux, macOS, and Windows
- Creates a simplified Makefile for reliable compilation
"""
function setup_ba_trees(base_dir::String)
    born_again_dp_path = joinpath(base_dir, "born_again_dp")
    executable_path = joinpath(base_dir, "bornAgain")
    
    # ========================================================================
    # Step 1: Clone repository if needed
    # ========================================================================
    if !isdir(born_again_dp_path)
        println("BA-Trees repository not found. Cloning...")
        try
            repo_url = "https://github.com/vidalt/BA-Trees.git"
            temp_clone_dir = joinpath(base_dir, "BA-Trees_temp")
            
            run(`git clone $repo_url $temp_clone_dir`)
            
            # Extract born_again_dp from cloned repository
            cloned_dp_path = joinpath(temp_clone_dir, "src", "born_again_dp")
            
            if isdir(cloned_dp_path)
                mv(cloned_dp_path, born_again_dp_path)
                println("BA-Trees repository cloned successfully")
                println("Moved born_again_dp to: $born_again_dp_path")
            else
                println("ERROR: born_again_dp directory not found in cloned repository")
                println("Expected location: $cloned_dp_path")
                _print_directory_structure(temp_clone_dir)
                rm(temp_clone_dir, recursive=true, force=true)
                return false
            end
            
            # Clean up temporary directory
            rm(temp_clone_dir, recursive=true, force=true)
            
        catch e
            println("ERROR: Failed to clone BA-Trees repository")
            println("Error details: ", e)
            println("Please ensure git is installed and you have internet connection")
            return false
        end
    else
        println(" born_again_dp directory already exists")
    end
    
    # ========================================================================
    # Step 2: Compile executable if needed
    # ========================================================================
    if !isfile(executable_path)
        println("bornAgain executable not found. Compiling...")
        try
            os_type = Sys.iswindows() ? "windows" : (Sys.isapple() ? "macos" : "linux")
            println("Detected OS: $os_type")
            
            cd(born_again_dp_path) do
                # Create simplified Makefile
                simple_makefile = create_simple_makefile(born_again_dp_path)
                
                # Check for original makefile (keep as backup)
                original_makefile = nothing
                if isfile("makefile")
                    original_makefile = "makefile"
                elseif isfile("Makefile")
                    original_makefile = "Makefile"
                end
                
                if original_makefile !== nothing
                    println("Original Makefile found: $original_makefile (keeping as backup)")
                    println("Using simplified Makefile for compilation")
                end
                
                # Attempt compilation
                compilation_success = _compile_ba_trees(simple_makefile, original_makefile)
                
                if !compilation_success
                    println("ERROR: All compilation attempts failed")
                    _print_compilation_help()
                    return false
                end
                
                # Locate and move executable
                compiled_exec = _find_compiled_executable()
                if compiled_exec === nothing
                    _print_compilation_error()
                    return false
                end
                
                # Move to base directory and set permissions
                target_name = Sys.iswindows() ? "bornAgain.exe" : "bornAgain"
                target_path = joinpath(base_dir, target_name)
                cp(compiled_exec, target_path, force=true)
                
                if !Sys.iswindows()
                    run(`chmod +x $target_path`)
                end
                
                println("bornAgain compiled successfully and moved to: $target_path")
            end
            
        catch e
            println("ERROR: Failed to compile bornAgain")
            println("Error details: ", e)
            _print_compilation_help()
            return false
        end
    else
        println(" bornAgain executable already exists")
    end
    
    println(" BA-Trees setup complete")
    return true
end

# Helper function to attempt compilation with different strategies
function _compile_ba_trees(simple_makefile::String, original_makefile::Union{String,Nothing})
    println("Running make with simplified Makefile...")
    
    try
        if Sys.iswindows()
            # Try different make commands on Windows
            try
                run(`make -f Makefile.simple`)
                return true
            catch
                println("Trying mingw32-make...")
                run(`mingw32-make -f Makefile.simple`)
                return true
            end
        else
            # Unix-like systems
            run(`make -f Makefile.simple`)
            return true
        end
    catch e
        println(" Simplified Makefile compilation failed: ", e)
        
        # Fallback to original Makefile
        if original_makefile !== nothing
            println("Attempting compilation with original Makefile (without CPLEX)...")
            try
                if Sys.iswindows()
                    try
                        run(`make -f $original_makefile`)
                        return true
                    catch
                        run(`mingw32-make -f $original_makefile`)
                        return true
                    end
                else
                    run(`make -f $original_makefile`)
                    return true
                end
            catch e2
                println("Original Makefile compilation also failed: ", e2)
            end
        end
    end
    
    return false
end

# Helper function to find compiled executable
function _find_compiled_executable()
    possible_names = ["bornAgain", "bornAgain.exe", "main", "main.exe"]
    
    for name in possible_names
        if isfile(name)
            println("Found compiled executable: $name")
            return name
        end
    end
    
    return nothing
end

# Helper functions for error reporting
function _print_directory_structure(path::String)
    println("\nRepository structure:")
    for (root, dirs, files) in walkdir(path)
        level = count(c -> c == '/', replace(root, path => ""))
        indent = "  " ^ level
        println("$(indent)$(basename(root))/")
        if level < 2
            for dir in dirs
                println("$(indent)  $(dir)/")
            end
        end
    end
end

function _print_compilation_help()
    println("""
    
    Please check:
    1. C++ compiler is installed (g++ or clang)
    2. Make is installed
    3. All dependencies are available
    """)
end

function _print_compilation_error()
    println("ERROR: Compilation completed but executable not found")
    println("Looking for: ", join(["bornAgain", "bornAgain.exe", "main", "main.exe"], ", "))
    println("\nDirectory contents:")
    for file in readdir()
        println("  ", file)
    end
end

# ============================================================================
# BA-Trees Execution
# ============================================================================

"""
    prepare_and_run_ba_trees(;
        dataset_name::String = "",
        num_trees::Int = nothing,
        max_depth::Int = -1,
        forest = nothing,
        mode_obj = 0
    ) -> Union{Dict{String,Int}, Nothing}

Prepare data and execute the BA-Trees algorithm.

This is the main entry point for running BA-Trees analysis. It handles:
- Automatic setup (clone and compile if needed)
- Forest generation or conversion
- BA-Trees execution
- Results display

# Keyword Arguments
- `dataset_name::String`: Name of the dataset (without extension)
- `num_trees::Int`: Number of trees to generate (if forest is nothing)
- `max_depth::Int`: Maximum depth for generated trees
- `forest`: Existing DecisionEnsemble to analyze (optional)
- `mode_obj::Int`: Optimization objective (0=Depth, 1=Leaves, 2=Depth+Leaves, 4=Heuristic)

# Returns
- Dictionary mapping class labels to indices, or nothing if execution fails

# Examples
```julia
# Generate and analyze a new forest
class_map = prepare_and_run_ba_trees(
    dataset_name="iris",
    num_trees=10,
    max_depth=5,
    mode_obj=4
)

# Analyze an existing forest
class_map = prepare_and_run_ba_trees(
    forest=my_forest,
    mode_obj=4
)
```
"""
function prepare_and_run_ba_trees(; 
    dataset_name::String = "", 
    num_trees::Int = nothing, 
    max_depth::Int = -1, 
    forest = nothing,
    mode_obj = 0
)
    # Validate mode_obj parameter
    valid_modes = [0, 1, 2, 4]
    if !(mode_obj in valid_modes)
        println("""
         ERROR: Invalid mode_obj value: $mode_obj
        Valid options for -obj parameter are:
          0 = Minimize Depth
          1 = Minimize Number of Leaves
          2 = Minimize Depth then Number of Leaves
          4 = Heuristic BA-Tree
        """)
        return nothing
    end

    base_dir = @__DIR__
    
    # ========================================================================
    # Step 1: Setup BA-Trees (automatic clone and compile)
    # ========================================================================
    println("\n" * "="^60)
    println("Checking BA-Trees setup...")
    println("="^60)
    
    if !setup_ba_trees(base_dir)
        println("ERROR: Failed to setup BA-Trees. Please check the error messages above.")
        return nothing
    end
    
    # Verify paths after setup
    db_path = joinpath(base_dir, "born_again_dp")
    executable_name = Sys.iswindows() ? "bornAgain.exe" : "bornAgain"
    executable_path = joinpath(base_dir, executable_name)
    
    if !isdir(db_path) || !isfile(executable_path)
        println("ERROR: Setup verification failed")
        println("born_again_dp exists: ", isdir(db_path))
        println("bornAgain exists: ", isfile(executable_path))
        return nothing
    end
 
    # ========================================================================
    # Step 2: Prepare temporary directory and file paths
    # ========================================================================
    temp_dir = joinpath(base_dir, "temp_ba_trees")
    isdir(temp_dir) || mkdir(temp_dir)
 
    input_file = joinpath(temp_dir, "forest.txt")
    output_base = joinpath(temp_dir, "result.txt")
    output_stats = output_base * ".out"
    output_tree = output_base * ".tree"
 
    try
        # ====================================================================
        # Step 3: Generate or convert forest to BA-Trees format
        # ====================================================================
        println("\n" * "="^60)
        println("Preparing random forest data...")
        println("="^60)
 
        forest_content = ""
        class_map = nothing
        
        if forest === nothing
            # Validate parameters for new forest generation
            if num_trees <= 0 || dataset_name == ""
                println("ERROR: Invalid parameters for random forest generation.")
                println("Required: num_trees > 0 and non-empty dataset_name")
                return nothing
            end
            
            # Generate new forest
            forest, model, start_time = learn_and_convert(num_trees, dataset_name, max_depth)
            forest_content, class_map = create_random_forest_input_from_model(forest)
        else
            # Use existing forest
            forest_content, class_map = create_random_forest_input_from_model(forest)
        end
 
        # Write forest data to file
        write(input_file, forest_content)
        println(" Forest data written to: $input_file")
 
        # ====================================================================
        # Step 4: Execute BA-Trees algorithm
        # ====================================================================
        num_trees_actual = length(forest.models)
        
        println("\n" * "="^60)
        println("Executing BA-Trees algorithm...")
        println("="^60)
        println("Trees: $num_trees_actual")
        println("Objective mode: $mode_obj")
        
        cmd = `$executable_path $input_file $output_base -trees $num_trees_actual -obj $mode_obj`
        println("Command: $cmd")
        run(cmd)
 
        # ====================================================================
        # Step 5: Display results
        # ====================================================================
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
                println("Tree structure file not found: $output_tree")
            end
        else
            println("Statistics file not found: $output_stats")
        end
        
        return class_map
 
    catch e
        println("\n" * "="^60)
        println("ERROR during execution:")
        println("="^60)
        println(e)
        println("\nDebug information:")
        println("Current directory: ", pwd())
        println("Executable exists: ", isfile(executable_path))
        println("Executable path: ", executable_path)
        
        if isdir(temp_dir)
            println("\nTemporary directory contents:")
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
    end
end

# ============================================================================
# High-Level API
# ============================================================================

"""
    WRAP_batrees(
        f,
        max_depth=10;
        dataset_name="",
        num_trees=nothing,
        mod=0
    ) -> Union{Dict{String,Int}, Nothing}

High-level wrapper function for BA-Trees analysis.

This is the main user-facing function that provides a simple interface
for running BA-Trees analysis on either new or existing forests.

# Arguments
- `f`: Existing DecisionEnsemble to analyze (or nothing to generate new)
- `max_depth::Int=10`: Maximum depth for tree generation

# Keyword Arguments
- `dataset_name::String=""`: Dataset name for new forest generation
- `num_trees::Int=nothing`: Number of trees for new forest
- `mod::Int=0`: Optimization mode (0=Depth, 1=Leaves, 2=Both, 4=Heuristic)

# Returns
- Dictionary mapping class labels to 0-based indices

# Examples
```julia
# Analyze existing forest
class_map = WRAP_batrees(my_forest, mod=4)

# Generate and analyze new forest
class_map = WRAP_batrees(
    nothing,
    10,
    dataset_name="iris",
    num_trees=20,
    mod=4
)
```
"""
function WRAP_batrees(
    f,
    max_depth=10;
    dataset_name="",
    num_trees=nothing,
    mod=0
)
    println("\n" * "="^60)
    println("Born-Again Tree Analysis")
    println("="^60)
    println("Dataset: $dataset_name")
    println("Number of trees: $(isnothing(f) ? num_trees : length(f.models))")
    println("Maximum depth: $max_depth")
    println("="^60 * "\n")

    if isnothing(f)
        # Generate new forest and analyze
        return prepare_and_run_ba_trees(
            dataset_name = dataset_name,
            num_trees = num_trees,
            max_depth = max_depth,
            forest = nothing,
            mode_obj = mod
        )
    else
        # Analyze existing forest
        return prepare_and_run_ba_trees(
            num_trees = length(f.models),
            forest = f,
            mode_obj = mod
        )
    end
end

# ============================================================================
# Legacy/Deprecated Functions
# ============================================================================

"""
    convert_feature(x) -> Float64

Convert an input feature to a numeric value.

# Arguments
- `x`: The feature to be converted (number, string, or symbol)

# Returns
- Float64 representation of the input

# Note
This is a simplified conversion function. For production use, consider
implementing more sophisticated categorical encoding strategies.
"""
function convert_feature(x)
    if x isa Number
        return Float64(x)
    elseif x isa String || x isa Symbol
        # Simplified conversion for strings/symbols
        # In production, use proper categorical encoding
        return 1.0
    else
        return Float64(x)
    end
end

"""
    prepare_dataset(dataset::AbstractArray) -> Tuple{Matrix{Float64}, Vector{Int}}

Prepare dataset for analysis by converting features and targets to numeric values.

This function converts all features to Float64 and maps target classes to
0-based integer indices.

# Arguments
- `dataset::AbstractArray`: Input dataset with features and target column

# Returns
- Tuple of (numeric_data_matrix, target_vector)

# Example
```julia
data, targets = prepare_dataset(my_dataset)
```
"""
function prepare_dataset(dataset::AbstractArray)
    n_features = size(dataset, 2) - 1
    n_samples = size(dataset, 1)

    # Convert features to numeric matrix
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
    create_tree_node(
        data,
        targets,
        depth=0,
        max_depth=3,
        min_samples=5,
        node_id=0
    ) -> Vector{String}

Create a decision tree node with random splits.

This function generates random decision trees for testing purposes.
For production use, prefer proper tree learning algorithms.

# Arguments
- `data`: Dataset matrix
- `targets`: Target vector
- `depth::Int=0`: Current depth
- `max_depth::Int=3`: Maximum allowed depth
- `min_samples::Int=5`: Minimum samples for split
- `node_id::Int=0`: Current node ID

# Returns
- Vector of strings representing tree nodes in BA-Trees format
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

    # Stopping conditions
    if depth >= max_depth || n_samples < min_samples
        majority_class = get_majority_class(targets)
        push!(nodes, "$node_id LN -1 -1 -1 -1 $depth $majority_class")
        return nodes
    end

    # Random feature and threshold selection
    feature = rand(0:(n_features-1))
    feature_values = data[:, feature+1]
    
    if length(unique(feature_values)) > 1
        threshold = rand() * (maximum(feature_values) - minimum(feature_values)) +
                   minimum(feature_values)
    else
        threshold = feature_values[1]
    end

    # Split data based on threshold
    left_mask = feature_values .<= threshold
    
    if !any(left_mask) || all(left_mask)
        majority_class = get_majority_class(targets)
        push!(nodes, "$node_id LN -1 -1 -1 -1 $depth $majority_class")
        return nodes
    end

    # Create child nodes recursively
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

    # Add current internal node
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
    create_random_forest_input(
        dataset,
        num_trees=10,
        max_depth=3
    ) -> String

Create random forest input content for testing.

Generates a random forest with bootstrap sampling for testing purposes.
For production, use proper ensemble learning methods.

# Arguments
- `dataset`: Input dataset
- `num_trees::Int=10`: Number of trees to generate
- `max_depth::Int=3`: Maximum tree depth

# Returns
- Formatted string for BA-Trees input

# Note
This is a legacy function primarily for testing. Use
`create_random_forest_input_from_model` for production use.
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
        "Format: node / node type (LN - leaf node, IN - internal node) left child / right child / feature / threshold / node_depth / majority class (starts with index 0)",
    )
    push!(lines, "")

    # Generate trees with bootstrap sampling
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

# ============================================================================
# BA-Trees Node Parsing (for result interpretation)
# ============================================================================

"""
    parse_ba_node(line::String) -> BANode

Parse a BA-Trees node line into a BANode struct.

# Arguments
- `line::String`: A line from the BA-Trees output file

# Returns
- BANode structure with parsed node information

# Example
```julia
node = parse_ba_node("0 IN 1 2 0 5.1 0 -1")
```
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
    convert_tree_structure(
        node,
        depth=0,
        node_counter=Ref(0),
        node_map=Dict()
    ) -> Vector{String}

Convert tree structure to BA-Trees format with slashes.

This is an alternative formatting function that uses '/' as separators.

# Arguments
- `node`: Tree node to convert
- `depth::Int=0`: Current depth
- `node_counter::Ref{Int}`: Counter for node IDs
- `node_map::Dict`: Mapping from nodes to IDs

# Returns
- Vector of formatted strings representing the tree structure
"""
function convert_tree_structure(
    node,
    depth=0,
    node_counter=Ref(0),
    node_map=Dict()
)::Vector{String}
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
    left_subtree = convert_tree_structure(
        node.posconsequent,
        depth + 1,
        node_counter,
        node_map
    )
    right_subtree = convert_tree_structure(
        node.negconsequent,
        depth + 1,
        node_counter,
        node_map
    )
    
    # Get children node IDs
    left_child = node_map[node.posconsequent]
    right_child = node_map[node.negconsequent]
    
    # Add current node
    majority_class = get_majority_class(node.info.supporting_labels)
    push!(
        lines,
        "$current_node / IN / $left_child / $right_child / $feature_idx / $threshold / $depth / $majority_class"
    )
    node_map[node] = current_node
    
    # Combine all nodes
    append!(lines, left_subtree)
    append!(lines, right_subtree)
    
    return lines
end

# ============================================================================
# DecisionTree.jl Conversion (for interoperability)
# ============================================================================

"""
    convert_to_dt_node(
        nodes::Dict{Int,BANode},
        current_id::Int,
        features::Vector{String}
    ) -> Node

Convert BA-Trees node to DecisionTree.jl node recursively.

This enables using BA-Trees results with the DecisionTree.jl package.

# Arguments
- `nodes::Dict{Int,BANode}`: Dictionary of all BA nodes indexed by ID
- `current_id::Int`: ID of the current node to convert
- `features::Vector{String}`: Feature names for the tree

# Returns
- DecisionTree.jl Node or Leaf structure
"""
function convert_to_dt_node(
    nodes::Dict{Int,BANode},
    current_id::Int,
    features::Vector{String}
)
    node = nodes[current_id]
    
    if node.node_type == "LN"
        # Leaf node
        return Leaf(node.class)
    else
        # Internal node
        feature_name = features[node.feature + 1]  # Convert from 0-based
        left = convert_to_dt_node(nodes, node.left_child, features)
        right = convert_to_dt_node(nodes, node.right_child, features)
        return Node(feature_name, node.threshold, left, right)
    end
end

"""
    convert_ba_to_dt(
        ba_tree_file::String,
        features::Vector{String},
        classes::Vector
    ) -> DecisionTreeClassifier

Convert BA-Trees format to DecisionTree.jl tree.

This function reads a BA-Trees output file and converts it to a
DecisionTree.jl compatible classifier that can be used for predictions.

# Arguments
- `ba_tree_file::String`: Path to the BA-Trees output file
- `features::Vector{String}`: Names of features in the dataset
- `classes::Vector`: Possible class values

# Returns
- DecisionTreeClassifier compatible with DecisionTree.jl

# Example
```julia
features = ["sepal_length", "sepal_width", "petal_length", "petal_width"]
classes = ["setosa", "versicolor", "virginica"]
dt = convert_ba_to_dt("result.txt.tree", features, classes)
predictions = DecisionTree.predict(dt, test_data)
```
"""
function convert_ba_to_dt(
    ba_tree_file::String,
    features::Vector{String},
    classes::Vector
)
    # Parse BA-Trees file
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

    # Build tree starting from root (ID 0)
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
    demonstrate_conversion(
        ba_tree_file::String,
        dataset::DataFrame
    ) -> DecisionTreeClassifier

Demonstrate BA-Trees to DecisionTree.jl conversion.

This example function shows how to convert BA-Trees output and use it
for predictions.

# Arguments
- `ba_tree_file::String`: Path to the BA-Trees output file
- `dataset::DataFrame`: Original dataset used to create the tree

# Returns
- Converted DecisionTreeClassifier

# Example
```julia
dt = demonstrate_conversion("result.txt.tree", iris_dataset)
```
"""
function demonstrate_conversion(ba_tree_file::String, dataset::DataFrame)
    # Extract metadata from dataset
    features = names(dataset)[1:end-1] |> Vector{String}
    classes = unique(dataset[!, end]) |> Vector
    
    # Convert the tree
    dt = convert_ba_to_dt(ba_tree_file, features, classes)
    
    # Display information
    println("\n" * "="^60)
    println("Converted Decision Tree Information:")
    println("="^60)
    println("Features: ", dt.features)
    println("Classes: ", dt.classes)
    println("="^60)
    
    # Example prediction on first row
    x = Array(dataset[1, 1:end-1])
    prediction = DecisionTree.predict(dt, x)
    println("\nExample Prediction:")
    println("Input: ", x)
    println("Predicted class: ", prediction)
    println("Actual class: ", dataset[1, end])
    println("="^60 * "\n")
    
    return dt
end

"""
    display_results(output_file::String)

Display BA-Trees analysis results in a formatted manner.

This function reads and pretty-prints the key metrics from a BA-Trees
output file.

# Arguments
- `output_file::String`: Path to output file (with or without .out extension)

# Example
```julia
display_results("result.txt")
```
"""
function display_results(output_file::String)
    output_file = endswith(output_file, ".out") ? output_file : output_file * ".out"
    
    if !isfile(output_file)
        println("Output file not found: $output_file")
        return
    end
    
    println("\n" * "="^60)
    println("Born-Again Tree Analysis Results:")
    println("="^60)
    
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
    
    println("="^60 * "\n")
end

# ============================================================================
# Hardcoded Test Functions (for development/testing)
# ============================================================================

"""
    prepare_and_run_ba_trees_hardcoded(;
        dataset_name="iris",
        num_trees=3,
        max_depth=3
    )

Run BA-Trees with hardcoded parameters for testing.

This is a convenience function for quick testing and development.
For production use, prefer the more flexible `prepare_and_run_ba_trees`.

# Keyword Arguments
- `dataset_name::String="iris"`: Name of dataset file (without extension)
- `num_trees::Int=3`: Number of trees to generate
- `max_depth::Int=3`: Maximum tree depth
"""
function prepare_and_run_ba_trees_hardcoded(;
    dataset_name = "iris",
    num_trees = 3,
    max_depth = 3
)
    base_dir = @__DIR__
    dataset_path = joinpath(base_dir, dataset_name * ".csv")
    
    if !isfile(dataset_path)
        println("ERROR: Dataset file not found: $dataset_path")
        println("Current directory: ", pwd())
        println("\nContent of src directory:")
        if isdir("src")
            for file in readdir("src")
                println("  - $file")
            end
        else
            println("The 'src' directory doesn't exist!")
        end
        return
    end

    dataset = DataFrame(CSV.File(dataset_path))

    temp_dir = joinpath(base_dir, "temp_ba_trees")
    isdir(temp_dir) || mkdir(temp_dir)

    input_file = joinpath(temp_dir, "forest.txt")
    output_base = joinpath(temp_dir, "result.txt")
    output_stats = output_base * ".out"
    output_tree = output_base * ".tree"

    executable_path = joinpath(base_dir, "bornAgain")

    try
        println("Generating random forest...")
        f, model, start_time = learn_and_convert(num_trees, dataset_name, max_depth)
        forest_content, class_map = create_random_forest_input_from_model(f)
        write(input_file, forest_content)

        cmd = `$executable_path $input_file $output_base -trees $num_trees -obj 4`
        println("Executing command: ", cmd)
        run(cmd)

        display_results(output_stats)
        
        if isfile(output_tree)
            println("\n" * "="^60)
            println("Born-Again Tree Structure:")
            println("="^60)
            println(read(output_tree, String))
        else
            println("Tree structure file not found: $output_tree")
        end

    catch e
        println("\nERROR during execution:")
        println(e)
        println("\nDebug information:")
        println("Current directory: ", pwd())
        println("Executable exists: ", isfile(executable_path))
        
        if isdir(temp_dir)
            println("\nTemporary directory contents:")
            for (root, dirs, files) in walkdir(temp_dir)
                println("Files: ", join(files, ", "))
            end
        end
        
        if isfile(input_file)
            println("\nFirst 10 lines of input file:")
            println(join(collect(Iterators.take(eachline(input_file), 10)), "\n"))
        end
    end
end

"""
    WRAP_batrees_hardcoded()

Execute BA-Trees analysis with hardcoded test parameters.

This is a convenience function for quick testing during development.
"""
function WRAP_batrees_hardcoded()
    println("="^60)
    println("Born-Again Tree Analysis (Hardcoded Test)")
    println("="^60)

    # Test parameters
    dataset_name = "iris"
    num_trees = 3
    max_depth = 3

    println("Dataset: $dataset_name")
    println("Number of trees: $num_trees")
    println("Maximum depth: $max_depth")
    println("="^60 * "\n")

    prepare_and_run_ba_trees_hardcoded(
        dataset_name = dataset_name,
        num_trees = num_trees,
        max_depth = max_depth,
    )
end

# ============================================================================
# Module Exports
# ============================================================================

export WRAP_batrees,
       WRAP_batrees_hardcoded,
       prepare_and_run_ba_trees,
       create_random_forest_input_from_model,
       convert_ba_to_dt,
       display_results,
       BANode,
       parse_ba_node

end # module RunModule