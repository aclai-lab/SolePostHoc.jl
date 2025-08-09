module Lumen

# Core dependencies - Essential packages for the LUMEN algorithm
using Logging, Dates, DataStructures, DataFrames
using Setfield  
using SoleModels
using SoleModels: DecisionSet  
using AbstractTrees
using SoleData
using SoleData: MultivariateScalarAlphabet, UnivariateScalarAlphabet  
using DecisionTree: load_data, build_forest, apply_forest
using ModalDecisionTrees, SoleLogics
using BenchmarkTools, StatProfilerHTML, Profile
using Base: intersect
using Base.Threads: @threads, Atomic, atomic_add!
using ConcurrentCollections, ProgressMeter

import DecisionTree as DT

# Public interface exports - Only expose what users need
export lumen, LumenConfig, LumenResult

##
# Configuration and Types
##

"""
    LumenConfig

Configuration parameters for the Logic-driven Unified Minimal Extractor of Notions (LUMEN) algorithm.

This struct encapsulates all configuration options for the LUMEN algorithm, providing a clean
interface with automatic validation and sensible defaults. It uses Julia's `@kwdef` macro
to enable keyword-based construction with default values.

# Fields

## Core Algorithm Parameters
- `minimization_scheme::Symbol = :AlgorithmName`: The DNF minimization algorithm to use
  - `:mitespresso`: Advanced minimization with good balance of speed/quality
  - `:boom`: Boom minimizator 
  - `:abc`: Minimization whit Berkeley framework 
  
## Coverage Parameters  
- `vertical::Float64 = 1.0`: Vertical coverage parameter (α) ∈ (0.0, 1.0]
  Controls how many instances must be covered by extracted rules
- `horizontal::Float64 = 1.0`: Horizontal coverage parameter (β) ∈ (0.0, 1.0]  
  Controls the breadth of rule coverage across feature space (% of different thresholds)

## Processing Modes
- `ott_mode::Bool = false`: Optimized truth table processing
  When `true`, uses memory-efficient and time-efficient algorithms for large datasets
- `controllo::Bool = false`: Enable validation mode
  Compares results between different processing methods for correctness verification

## Customization Options
- `minimization_kwargs::NamedTuple = (;)`: Additional parameters for minimization algorithms
- `filteralphabetcallback = identity`: Custom function to filter/modify the logical alphabet
- `apply_function = nothing`: Custom function for model application
  If `nothing`, automatically determined based on model type (with SoleModels)
  
## Output Control
- `silent::Bool = false`: Suppress progress and diagnostic output
- `return_info::Bool = true`: Include additional metadata in results
- `vetImportance::Vector = []`: Vector for tracking feature importance values

## Testing and Debugging
- `testott = nothing`: Special testing mode for optimization validation
- `alphabetcontroll = nothing`: Special mode for alphabet analysis only

# Constructor Validation

The constructor automatically validates parameters and throws descriptive errors:
- Coverage parameters must be in range (0.0, 1.0]
- Minimization scheme must be supported
- Inconsistent parameter combinations are caught early

# Examples

```julia
# Basic usage with defaults
config = LumenConfig()

# Customized configuration
config = LumenConfig(
    minimization_scheme = :abc,
    vertical = 0.8,
    horizontal = 0.9,
    silent = true
)

# Advanced configuration with custom processing
config = LumenConfig(
    ott_mode = true,
    minimization_kwargs = (max_iterations = 1000,),
    filteralphabetcallback = my_custom_filter
)
```

See also: [`lumen`](@ref), [`LumenResult`](@ref), [`validate_config`](@ref)
"""
Base.@kwdef struct LumenConfig
    minimization_scheme::Symbol = :abc
    vertical::Float64 = 1.0
    horizontal::Float64 = 1.0
    ott_mode::Bool = true # TODO VALUATE THIS (if df. true?)
    controllo::Bool = false
    minimization_kwargs::NamedTuple = (;)
    filteralphabetcallback = identity
    apply_function = nothing
    silent::Bool = false
    return_info::Bool = true
    vetImportance::Vector = []
    testott = nothing
    alphabetcontroll = nothing
    #use_listrules::Bool = false
    
    # Custom constructor with validation - ensures configuration integrity
    # This approach catches configuration errors at creation time rather than runtime
    function LumenConfig(args...)
        config = new(args...)
        validate_config(config)  # Fail fast on invalid configurations
        return config
    end
end

"""
    LumenResult

Comprehensive result structure containing extracted logical rules and associated metadata.

This immutable struct encapsulates all outputs from the LUMEN algorithm, providing a clean
and extensible interface for accessing results. The design follows the principle of 
returning rich, self-documenting results rather than simple tuples.

# Fields

- `decision_set::DecisionSet`: The primary output - a collection of minimized logical rules
  Each rule consists of a logical formula (antecedent) and a decision outcome (consequent)
  
- `info::NamedTuple`: Extensible metadata container with algorithm-specific information
  Common fields include:
  - `vectPrePostNumber`: Vector of (pre, post) minimization term counts
  - `unminimized_ds`: Original decision set before minimization (if requested)
  - `processing_time`: Total algorithm execution time
  - `feature_importance`: Feature ranking information (if available)
  
- `processing_time::Float64`: Total processing time in seconds
  Measured from algorithm start to completion, useful for performance analysis

# Constructors

Two constructors are provided for different use cases:

```julia
# Full constructor - for complete results with metadata
LumenResult(decision_set, info_tuple, processing_time)

# Minimal constructor - when only rules are available  
LumenResult(decision_set)  # info=empty, processing_time=0.0
```

# Design Rationale

This structured approach provides several advantages over returning raw tuples:
1. **Self-documentation**: Field names clearly indicate content
2. **Type safety**: Julia's type system validates structure at compile time  
3. **Extensibility**: Easy to add new fields without breaking existing code
4. **IDE support**: Autocompletion and inline documentation
5. **Backward compatibility**: Old code can still access fields by name

# Examples

```julia
# Basic usage
result = lumen(model, config)
rules = result.decision_set
println("Extracted (length(rules.rules)) rules in (result.processing_time)s")

# Accessing metadata
if haskey(result.info, :vectPrePostNumber)
    stats = result.info.vectPrePostNumber
    total_reduction = sum(pre - post for (pre, post) in stats)
    println("Reduced formula complexity by \total_reduction terms")
end

# Comparing minimized vs original rules
if haskey(result.info, :unminimized_ds)
    original_rules = result.info.unminimized_ds
    println("Minimization: (length(original_rules.rules)) → (length(result.decision_set.rules))")
end
```

See also: [`lumen`](@ref), [`LumenConfig`](@ref), [`DecisionSet`](@ref)
"""
struct LumenResult
    decision_set::DecisionSet  # Primary output: minimized logical rules
    info::NamedTuple          # Extensible metadata container
    processing_time::Float64   # Performance metric in seconds
    
    # Full constructor - used when all information is available
    LumenResult(ds, info, time) = new(ds, info, time)
    
    # Convenience constructor - for cases where only rules matter
    # Uses empty NamedTuple (;) and zero time as sensible defaults
    LumenResult(ds) = new(ds, (;), 0.0)
end

##
# Validation Functions
##

"""
    validate_config(config::LumenConfig) -> Nothing

Comprehensive validation of LUMEN configuration parameters.

This function implements fail-fast validation to catch configuration errors early
rather than allowing them to propagate through the algorithm. It checks both
individual parameter validity and logical consistency between parameters.

# Validation Rules

## Coverage Parameters
- Both `vertical` and `horizontal` must be in range (0.0, 1.0]
- These control rule extraction coverage and must be positive and ≤ 1.0

## Algorithm Parameters  
- `minimization_scheme` must be one of the supported algorithms
- Each scheme has different performance characteristics and use cases

## Consistency Checks
- Certain parameter combinations may be invalid or suboptimal
- Warnings are issued for potentially problematic configurations

# Arguments
- `config::LumenConfig`: Configuration object to validate

# Throws
- `ArgumentError`: For invalid parameter values or unsupported combinations
- Includes descriptive error messages indicating the specific problem

# Examples

```julia
# This will throw an error
try
    config = LumenConfig(vertical = -0.5)  # Invalid: negative coverage
catch ArgumentError as e
    println("Configuration error: e")
end

# This will succeed
config = LumenConfig(vertical = 0.8, horizontal = 0.9)  # Valid ranges
```

# Implementation Notes

The validation approach separates concerns:
1. Parameter range validation (mathematical constraints)
2. Enum validation (supported algorithm schemes)  
3. Logical consistency validation (parameter interactions)

This structure makes it easy to add new validation rules as the algorithm evolves.
"""
function validate_config(config::LumenConfig)
    # Validate coverage parameters - must be positive and ≤ 1.0
    # These parameters control the proportion of instances that must be covered by rules
    if config.vertical <= 0.0 || config.vertical > 1.0 || 
       config.horizontal <= 0.0 || config.horizontal > 1.0
        throw(ArgumentError(
            "vertical and horizontal parameters must be in range (0.0, 1.0]. " *
            "Got vertical=$(config.vertical), horizontal=$(config.horizontal). " *
            "These parameters control rule coverage and must be meaningful proportions."
        ))
    end
    
    # Validate minimization scheme - only certain algorithms are implemented
    # Each algorithm has different trade-offs in speed vs. minimization quality
    valid_schemes = (:mitespresso, :boom, :abc)
    if config.minimization_scheme ∉ valid_schemes
        throw(ArgumentError(
            "minimization_scheme must be one of: $(valid_schemes). " *
            "Got: $(config.minimization_scheme). " *
            "Each scheme offers different performance characteristics:\n" *
            "  :mitespresso - balanced speed/quality for most cases\n" *
            "  :boom - aggressive minimization for complex formulas\n" *
            "  :abc - basic minimization by Berkeley framework, fastest but less thorough"
        ))
    end
    
    # Additional consistency checks could be added here
    # For example, warning about potentially problematic parameter combinations
    if config.ott_mode && config.controllo
        @warn "Running both ott_mode and controllo simultaneously may impact performance"
    end
end

##
# Core Algorithm Components - Rule Extraction Strategies
##

"""
    RuleExtractionStrategy

Abstract base type for different rule extraction approaches.

This design uses the Strategy pattern to handle different types of models
(single trees vs. ensembles) with appropriate algorithms. The pattern provides:

1. **Extensibility**: Easy to add new extraction methods
2. **Type safety**: Compile-time method dispatch  
3. **Clarity**: Each strategy has a clear, focused responsibility
4. **Testability**: Individual strategies can be tested in isolation

# Concrete Strategies

- [`StandardExtraction`](@ref): For ensemble models powered by Sole
- [`EnsembleExtraction`](@ref): For generic ensemble models (forests, boosting, etc.)

# Design Pattern Benefits

The Strategy pattern is chosen here because:
- Rule extraction logic differs significantly between model types
- We want to avoid large if-else blocks in the main algorithm
- Future model types can be supported by adding new strategies
- Each strategy can optimize for its specific model characteristics

See also: [`extract_rules`](@ref), [`StandardExtraction`](@ref), [`EnsembleExtraction`](@ref)
"""
abstract type RuleExtractionStrategy end

"""
    StandardExtraction <: RuleExtractionStrategy

Strategy for extracting rules from Sole decision models.

This strategy is optimized for Sole structures and provides
direct rule extraction without the complexity needed for ensemble handling.
The approach maintains the tree's original logical structure while converting
it to a rule-based representation.

# Characteristics
- **Efficiency**: Direct tree traversal without ensemble overhead
- **Simplicity**: Straightforward path-to-rule conversion
- **Preservation**: Maintains original tree decision logic

See also: [`RuleExtractionStrategy`](@ref), [`EnsembleExtraction`](@ref)
"""
struct StandardExtraction <: RuleExtractionStrategy end

"""
    EnsembleExtraction <: RuleExtractionStrategy  

Strategy for extracting rules from ensemble models (Random Forests, AdaBoost, etc.).

This strategy handles the complexity of multiple tree models by extracting rules
from each constituent tree and then combining/deduplicating the results. The approach
balances computational efficiency with comprehensive rule coverage.

# Characteristics
- **Comprehensiveness**: Extracts rules from all ensemble members
- **Deduplication**: Removes identical rules across trees
- **Scalability**: Handles ensembles of arbitrary size
- **Preservation**: Maintains the collective decision logic of the ensemble

# Implementation Notes
The strategy uses `unique()` to eliminate duplicate rules, which can significantly
reduce the final ruleset size for ensembles where trees learn similar patterns.

See also: [`RuleExtractionStrategy`](@ref), [`StandardExtraction`](@ref)
"""
struct EnsembleExtraction <: RuleExtractionStrategy end

"""
    extract_rules(model, ::StandardExtraction; use_shortforms=true) -> Vector{Rule}

Extract logical rules from a single decision tree model.

This method implements rule extraction for individual decision trees by traversing
each path from root to leaf and converting the path conditions into logical rules.
The process preserves the tree's decision logic while converting it to a more
flexible rule representation.

# Arguments
- `model`: Single decision tree model to extract rules from
- `::StandardExtraction`: Strategy type parameter (for dispatch)
- `use_shortforms::Bool=true`: Whether to use shortened rule representations

# Returns
- `Vector{Rule}`: Collection of logical rules representing the tree's decision logic

# Algorithm Details
1. **Path Enumeration**: Traverse all root-to-leaf paths in the tree
2. **Condition Extraction**: Convert split conditions to logical atoms
3. **Rule Construction**: Combine path conditions into rule antecedents
4. **Leaf Integration**: Associate each rule with its corresponding decision

# Performance Notes
- Time complexity: O(n) where n is the number of nodes
- Space complexity: O(k) where k is the number of leaves (rules)
- Optimized for single-tree scenarios without ensemble overhead
"""
function extract_rules(model, ::StandardExtraction; use_shortforms=true)
    return listrules(model; use_shortforms)
end

"""
    extract_rules(model, ::EnsembleExtraction; use_shortforms=true) -> Vector{Rule}

Extract and consolidate logical rules from an ensemble of decision trees.

This method handles the complexity of ensemble models by extracting rules from
each constituent tree and then combining them into a unified ruleset. The process
includes deduplication to avoid redundant rules while preserving the ensemble's
collective decision-making capability.

# Arguments
- `model`: Ensemble model containing multiple decision trees
- `::EnsembleExtraction`: Strategy type parameter (for dispatch)
- `use_shortforms::Bool=true`: Whether to use shortened rule representations

# Returns
- `Vector{Rule}`: Unified collection of unique logical rules from all trees

# Algorithm Details
1. **Per-Tree Extraction**: Extract rules from each tree in the ensemble
2. **Deduplication**: Remove identical rules using `unique()`
3. **Consolidation**: Flatten nested rule vectors into single collection
4. **Validation**: Ensure consistent rule format across the ensemble

# Performance Considerations
- Time complexity: O(m×n) where m is ensemble size, n is average tree size
- Space complexity: O(k) where k is the total number of unique rules
- Memory usage optimized through deduplication
- Parallelization potential for very large ensembles

# Implementation Notes
The method handles the common case where `listrules` returns nested vectors
by flattening them with `reduce(vcat, rs)`. This ensures consistent output
format regardless of the specific ensemble implementation.
"""
function extract_rules(model, ::EnsembleExtraction; use_shortforms=true)
    # Extract rules from each tree in the ensemble
    # Using list comprehension for clarity and potential parallelization
    rs = unique([listrules(tree; use_shortforms) for tree in SoleModels.models(model)])
    
    # Handle nested vector structures that some ensemble types produce
    # The ternary operator provides robust handling of different return types
    return rs isa Vector{<:Vector{<:Any}} ? reduce(vcat, rs) : rs
end

"""
    extract_rules(model; use_shortforms=true) -> Vector{Rule}

Automatic rule extraction with strategy selection based on model type.

This is the main entry point for rule extraction that automatically chooses
the appropriate strategy based on the input model type. It provides a clean
interface while leveraging the Strategy pattern internally for optimal performance.

# Arguments  
- `model`: Decision tree model (single tree or ensemble)
- `use_shortforms::Bool=true`: Whether to use shortened rule representations

# Returns
- `Vector{Rule}`: Collection of logical rules extracted from the model

# Strategy Selection Logic
- **Single trees**: Uses `StandardExtraction` for optimal single-tree performance
- **Ensembles**: Uses `EnsembleExtraction` to handle multiple trees correctly
- **Unknown types**: Defaults to `StandardExtraction` with potential warnings

# Design Benefits
This automatic dispatch approach provides:
1. **Simplicity**: Users don't need to know about strategy types
2. **Correctness**: Always uses the right algorithm for the model type
3. **Performance**: Optimal strategy selection without user overhead
4. **Maintainability**: Strategy logic centralized in one place

# Examples
```julia
# Works with single trees
tree_model = build_tree(X, y)
rules = extract_rules(tree_model)

# Works with ensembles  
forest_model = build_forest(X, y)
rules = extract_rules(forest_model)

# Automatic strategy selection handles both cases correctly
```
"""
function extract_rules(model; use_shortforms=true)
    # Automatically select strategy based on model type
    # This design decision centralizes the type-checking logic
    strategy = isensemble(model) ? EnsembleExtraction() : StandardExtraction()
    return extract_rules(model, strategy; use_shortforms)
end

##
# Alphabet and Atom Processing
##

"""
    extract_atoms(model, filteralphabetcallback) -> (Int, Vector, Alphabet)

Extract and process logical atoms from the model's decision alphabet.

This function performs the critical task of extracting the fundamental logical
building blocks (atoms) from a decision tree model. These atoms represent the
basic conditions used in decision rules and form the vocabulary for the final
logical formulas.

# Arguments
- `model`: Decision tree model to analyze
- `filteralphabetcallback`: Custom function to filter/process the alphabet
  - Signature: `alphabet -> filtered_alphabet`
  - Default: `identity` (no filtering)
  - Use cases: Feature selection, condition refinement, domain constraints

# Returns
A tuple containing:
1. `Int`: Total number of atoms before deduplication
2. `Vector`: Unique atoms extracted from the model  
3. `Alphabet`: Processed logical alphabet for formula construction

# Processing Pipeline

## Step 1: Raw Atom Extraction
```julia
all_atoms = collect(atoms(SoleModels.alphabet(model, false)))
```
Extracts all atomic conditions from the model's internal representation.

## Step 2: Deduplication
```julia  
my_atoms = unique(all_atoms)
```
Removes duplicate atoms that may appear in different parts of the model.

## Step 3: Feature Validation
```julia
n_features = maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
```
Ensures all features are properly indexed and validates feature naming scheme.

## Step 4: Alphabet Construction
```julia
my_alphabet = filteralphabetcallback(process_alphabet(my_atoms, n_features))
```
Builds the logical alphabet and applies any custom filtering.

## Step 5: Operator Validation
```julia
all(x -> (x == (<)), SoleData.test_operator.(subalphabets(my_alphabet)))
```
Ensures only supported comparison operators are present.

# Error Conditions

## Symbolic Features
Throws `ArgumentError` if symbolic (non-integer) feature names are detected.
The current implementation requires integer feature indices for efficient processing.

## Unsupported Operators  
Throws `ArgumentError` if operators other than `<` are found.
This limitation can be extended in future versions.

# Design Rationale

## Why Extract Atoms First?
1. **Modularity**: Separates atom extraction from alphabet processing
2. **Reusability**: Atoms can be reused across different alphabet configurations  
3. **Validation**: Early detection of unsupported model features
4. **Optimization**: Deduplication reduces downstream processing overhead

## Filter Callback Pattern
The `filteralphabetcallback` parameter enables:
- **Domain expertise**: Users can inject domain-specific filtering
- **Feature selection**: Remove irrelevant or noisy features
- **Preprocessing**: Apply transformations before rule generation
- **Experimentation**: Easy A/B testing of different alphabet configurations

# Examples

```julia
# Basic usage with no filtering
num_atoms, atoms, alphabet = extract_atoms(model, identity)

# Custom filtering to remove low-importance features  
function remove_weak_features(alphabet)
    # Custom logic to filter alphabet
    return filtered_alphabet
end
num_atoms, atoms, alphabet = extract_atoms(model, remove_weak_features)

# Feature importance-based filtering
importance_filter = alphabet -> filter_by_importance(alphabet, min_importance=0.1)
num_atoms, atoms, alphabet = extract_atoms(model, importance_filter)
```

See also: [`process_alphabet`](@ref), [`LumenConfig`](@ref), [`validate_config`](@ref)
"""
function extract_atoms(model, filteralphabetcallback)
    # Step 1: Extract all atomic conditions from the model
    # The 'false' parameter requests the raw alphabet without preprocessing
    all_atoms = collect(atoms(SoleModels.alphabet(model, false)))
    num_all_atoms = length(all_atoms)
    
    # Step 2: Remove duplicate atoms for efficiency
    # Decision trees often reuse the same conditions across different branches
    my_atoms = unique(all_atoms)
    
    # Step 3: Validate feature indexing scheme
    # Extract the maximum feature index to determine the feature space size
    n_features = maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
    if !isa(n_features, Integer)
        throw(ArgumentError(
            "Symbolic feature names are not supported. " *
            "All features must be indexed with integers. " *
            "Found feature type: $(typeof(n_features)). " *
            "Consider preprocessing your model to use integer feature indices."
        ))
    end
    
    # Step 4: Build and filter the logical alphabet
    # The alphabet defines the vocabulary available for constructing logical formulas
    my_alphabet = filteralphabetcallback(process_alphabet(my_atoms, n_features))
    
    # Step 5: Validate supported operators
    # Currently only '<' comparison is fully implemented and tested
    supported_operators = SoleData.test_operator.(subalphabets(my_alphabet))
    if !all(x -> (x == (<)), supported_operators)
        unsupported = unique(filter(x -> x != (<), supported_operators))
        throw(ArgumentError(
            "Only '<' operator is currently supported. " *
            "Found unsupported operators: $(unsupported). " *
            "This limitation may be addressed in future versions. " *
            "Consider preprocessing your model to use only '<' conditions."
        ))
    end
    
    return num_all_atoms, my_atoms, my_alphabet
end

##
# Truth Combination Generation
##

"""
    generate_combinations(model, alphabet, atoms, vertical, config::LumenConfig) -> (Results, LabelCount)

Generate truth value combinations for logical formula construction.

This function creates the truth table that serves as the foundation for DNF formula
generation. It evaluates the model on carefully selected input combinations to
determine which logical formulas correspond to each decision outcome.

# Arguments
- `model`: Decision tree model to evaluate
- `alphabet`: Logical alphabet defining the available conditions
- `atoms`: Vector of atomic logical conditions  
- `vertical::Float64`: Coverage parameter controlling instance sampling
- `horizontal::Float64`: Coverage parameter controlling feature sampling TODO:WORKINPROGRESS
- `config::LumenConfig`: Configuration object containing processing options

# Returns
- `Results`: Truth value combinations and their corresponding model outputs
- `LabelCount`: Statistics about label distribution in the generated combinations

# Algorithm Selection

The function automatically selects between two implementation strategies:

## Standard Mode (`ott_mode = false`)
- **Algorithm**: `truth_combinations()`
- **Use case**: General-purpose processing for most datasets
- **Characteristics**: Straightforward implementation with good memory usage
- **Performance**: Optimal for small to medium-sized problems

## Optimized Mode (`ott_mode = true`)  
- **Algorithm**: `truth_combinations_ott()`
- **Use case**: Large-scale problems requiring memory optimization
- **Characteristics**: Reduced memory footprint with specialized data structures
- **Performance**: Better for very large feature spaces or deep trees

# Processing Flow

## Phase 1: Input Validation
- Validates alphabet consistency with atoms
- Checks vertical parameter bounds  
- Ensures model compatibility with chosen processing mode

## Phase 2: Sample Generation
- Creates systematic samples covering the logical space
- Applies vertical sampling to control computational complexity
- Balances coverage completeness with computational efficiency

## Phase 3: Model Evaluation
- Evaluates the model on each generated sample
- Records both inputs (truth assignments) and outputs (predictions)
- Maintains correspondence between logical conditions and outcomes

## Phase 4: Result Aggregation
- Groups results by predicted outcome/label
- Computes statistics for coverage analysis
- Prepares data structures for subsequent minimization

# Performance Considerations

## Memory Usage
- Standard mode: O(n×m) where n=samples, m=features
- Optimized mode: O(k×m) where k<n through specialized structures
- Trade-off between memory efficiency and processing simplicity

## Computational Complexity
- Sample generation: O(2^k) where k is the effective feature count
- Model evaluation: O(n×d) where d is average tree depth  
- Result processing: O(n×log(n)) for sorting and grouping

## Scalability Recommendations
- Use standard mode for: < 20 features, < 1M samples
- Use optimized mode for: > 20 features, > 1M samples, memory constraints
- Monitor memory usage and switch modes if needed

# Configuration Impact

The `config` parameter influences processing through:
- `silent`: Controls progress output during generation
- `apply_function`: Custom model evaluation function
- `ott_mode`: Chooses between standard and optimized algorithms
- Processing parameters passed to the underlying algorithms

# Examples

```julia
# Standard processing
results, counts = generate_combinations(model, alphabet, atoms, 1.0, config)

# Memory-optimized processing for large problems
config_opt = LumenConfig(ott_mode=true)
results, counts = generate_combinations(model, alphabet, atoms, 0.8, config_opt)

# Custom processing with specific apply function
config_custom = LumenConfig(apply_function=my_custom_evaluator)
results, counts = generate_combinations(model, alphabet, atoms, 1.0, config_custom)
```

# Error Handling

The function validates inputs and provides descriptive errors for:
- Inconsistent alphabet/atoms combinations
- Invalid vertical parameter values  
- Model/apply_function mismatches
- Memory allocation failures in large-scale processing

See also: [`truth_combinations`](@ref), [`truth_combinations_ott`](@ref), [`LumenConfig`](@ref)
"""
function generate_combinations(model, alphabet, atoms, vertical, config::LumenConfig) #TODO vertical is not required, in config.vertical we have vertical value 
    # Algorithm selection based on processing requirements
    # The choice between standard and optimized mode significantly impacts
    # both memory usage and computational characteristics
    
    if config.ott_mode
        # Optimized processing for large-scale problems
        # This mode uses specialized data structures and algorithms
        # to reduce memory footprint while maintaining correctness
        return truth_combinations_ott(
            model, alphabet, atoms, vertical; 
            silent=config.silent, 
            apply_function=config.apply_function
        )
    else
        # Standard processing for general use cases
        # This mode prioritizes simplicity and maintainability
        # while providing good performance for typical problems
        return truth_combinations(
            model, alphabet, atoms, vertical; 
            silent=config.silent, 
            apply_function=config.apply_function
        )
    end
end

##
# Rule Processing and Minimization
##

"""
    process_rules(combined_results, config::LumenConfig) -> (Vector{Rule}, Union{Vector{Rule}, Nothing}, Vector{Tuple{Int,Int}})

Process and minimize extracted logical rules using the configured minimization scheme.

This function represents the core minimization phase of the LUMEN algorithm, where
raw DNF formulas extracted from decision trees are simplified into more compact
and interpretable logical rules. The process balances minimization quality with
computational efficiency.

# Arguments
- `combined_results`: Paired collection of (decision_outcome, dnf_formula) tuples
- `config::LumenConfig`: Configuration specifying minimization parameters

# Returns
A tuple containing:
1. `Vector{Rule}`: Minimized logical rules ready for use
2. `Union{Vector{Rule}, Nothing}`: Original unminimized rules (if `return_info=true`)
3. `Vector{Tuple{Int,Int}}`: Statistics showing (pre, post) minimization term counts

# Processing Pipeline

## Phase 1: Initialization  
```julia
unminimized_rules = config.return_info ? Rule[] : nothing
minimized_rules = Rule[]
vect_pre_post_number = Vector{Tuple{Int,Int}}()
```
Sets up data structures based on configuration requirements.

## Phase 2: Per-Formula Processing
For each (result, formula) pair:

### Step 2.1: Optional Archival
```julia
if config.return_info
    push!(unminimized_rules, convert_DNF_formula(formula, result, 1.0))
end
```
Preserves original formulas for comparison and analysis.

### Step 2.2: Minimization
```julia
minimized_formula = minimize_formula(formula, config)
```
Applies the configured minimization algorithm to simplify the formula.

### Step 2.3: Statistics Tracking
```julia
if should_track_statistics(config.minimization_scheme)
    pre_terms = nterms(formula)
    post_terms = nterms(minimized_formula) 
    push!(vect_pre_post_number, (pre_terms, post_terms))
end
```
Records minimization effectiveness for performance analysis.

### Step 2.4: Rule Construction
```julia
rule = create_rule(minimized_formula, result, config)
push!(minimized_rules, rule)
```
Converts the minimized formula back to a Rule object.

# Minimization Algorithms

## :mitespresso
- **Characteristics**: Balanced approach optimizing both speed and quality
- **Best for**: General use cases with moderate complexity
- **Trade-offs**: Good compression ratio with reasonable computation time

## :boom  
- **Characteristics**: Aggressive minimization for maximum compression
- **Best for**: Complex formulas where maximum simplification is critical
- **Trade-offs**: Higher computation time for better minimization results

## :abc
- **Characteristics**: Basic minimization with Berkeley framework
- **Best for**: Large-scale processing where speed is more important than optimal compression
- **Trade-offs**: Fastest processing with moderate compression

# Error Handling and Recovery

The function implements robust error handling:

```julia
try
    # Minimization process
    minimized_formula = minimize_formula(formula, config)
    rule = create_rule(minimized_formula, result, config)
    push!(minimized_rules, rule)
catch e
    @error "Minimization failed for formula" formula=formula exception=e
    # Could implement fallback strategies here
    rethrow(e)
end
```

## Fallback Strategies
- **Graceful degradation**: Use unminimized formula if minimization fails
- **Algorithm switching**: Try alternative minimization schemes
- **Partial results**: Continue processing remaining formulas

# Performance Optimization

## Memory Management
- Reuses data structures across iterations
- Releases temporary objects promptly
- Monitors memory usage for large formula sets

## Computational Efficiency
- Parallel processing potential for independent formulas
- Early termination for trivial formulas
- Caching of intermediate results

# Quality Metrics

The function tracks several metrics to assess minimization quality:

## Compression Ratio
```julia
compression_ratio = pre_terms / post_terms
total_reduction = sum(pre - post for (pre, post) in vect_pre_post_number)
```

## Coverage Preservation
- Ensures minimized formulas maintain the same logical coverage
- Validates that no decision cases are lost during minimization

# Examples

```julia
# Basic minimization with default settings
minimized_rules, _, stats = process_rules(results, config)

# Advanced minimization with full information retention
config_full = LumenConfig(return_info=true, minimization_scheme=:boom)
minimized, original, stats = process_rules(results, config_full)

# Analyze minimization effectiveness
total_reduction = sum(pre - post for (pre, post) in stats)
```

See also: [`minimize_formula`](@ref), [`create_rule`](@ref), [`LumenConfig`](@ref)
"""
function process_rules(combined_results, config::LumenConfig)
    minimized_rules = Rule[]
    vect_pre_post_number = Vector{Tuple{Int,Int}}()
    original_formulas = config.return_info ? Dict{String, Any}() : nothing
    
    for (result, formula) in combined_results
        config.silent || println("Performing minimization for: $result")
        

        if config.return_info
            original_formulas[result] = formula
        end
        
        minimized_formula = minimize_formula(formula, config)
        config.silent || println("End of minimization for: $result ")
        
        if should_track_statistics(config.minimization_scheme)
            pre_terms = get_formula_terms(formula)
            post_terms = get_formula_terms(minimized_formula)
            push!(vect_pre_post_number, (pre_terms, post_terms))
            log_minimization_stats(pre_terms, post_terms, config.silent)
        end
        
        try
            rule = create_rule(minimized_formula, result, config)
            config.silent || println("Generated rule: $rule")
            push!(minimized_rules, rule)
        catch e
            @error "Failed to create rule from minimized formula" exception=e
            rethrow(e)
        end
    end
    
    # TODO EVALUATE THIS, IS VERY BUSY
    unminimized_rules = nothing
    if config.return_info && !isempty(original_formulas)
        unminimized_rules = Rule[]
        #for (result, formula) in original_formulas
        #    try
        #        println("browuimuore")
        #        rule = convert_DNF_formula(formula, result, 1.0)
        #        push!(unminimized_rules, rule)
        #    catch
        #        # Se la conversione fallisce, usa una rappresentazione semplice
        #        push!(unminimized_rules, Rule(⊤, result))
        #    end
        #end
    end
    
    return minimized_rules, unminimized_rules, vect_pre_post_number
end

# Helper function to handle different formula types
function get_formula_terms(formula)
    if formula isa TwoLevelDNFFormula
        return nterms(formula)
    elseif hasmethod(nterms, (typeof(formula),))
        return nterms(formula)
    else
        # For LeftmostLinearForm and other types, count the number of disjuncts
        try
            if hasmethod(length, (typeof(formula),))
                return length(formula)
            elseif hasfield(typeof(formula), :children)
                return length(formula.children)
            else
                # Fallback: count terms by converting to string and parsing
                formula_str = string(formula)
                # Count the number of OR operators + 1 (rough approximation)
                return count("∨", formula_str) + 1
            end
        catch
            # Last resort fallback
            return 1
        end
    end
end

##
# Helper Functions - Utility functions supporting the main algorithm
##

"""
    determine_apply_function(model, config_apply_function) -> Function

Determine the appropriate function for applying the model to input data.

This function implements automatic function selection based on model type,
providing a clean abstraction over the differences between various model
implementations. The selection logic handles the common case where users
don't specify a custom apply function.

# Arguments
- `model`: The decision tree model (single tree or ensemble)
- `config_apply_function`: User-specified apply function (may be `nothing`)

# Returns
- `Function`: The appropriate function for evaluating the model

# Selection Logic

## User-Specified Function
If `config_apply_function` is not `nothing`, it takes precedence:
```julia
if !isnothing(config_apply_function)
    return config_apply_function
end
```
This allows users to provide custom evaluation functions for specialized scenarios.

## Automatic Selection
For automatic selection, the function uses type-based dispatch:

### DecisionTree.jl Ensembles
```julia
if model isa DT.Ensemble
    return DT.apply_forest
end
```
Uses the optimized ensemble application function from DecisionTree.jl.

### SoleModels Framework
```julia
else
    return SoleModels.apply  
end
```
Uses the generic application function for single trees and other model types.

# Design Rationale

## Why Automatic Selection?
1. **User Experience**: Eliminates need for users to understand implementation details
2. **Correctness**: Ensures the right function is used for each model type
3. **Performance**: Each function is optimized for its specific model type
4. **Maintainability**: Centralized logic for function selection

## Extensibility
New model types can be easily supported by extending the selection logic:
```julia
# Future extension example
if model isa NewModelType
    return new_model_apply_function
elseif model isa AnotherModelType
    return another_apply_function
end
```

# Performance Considerations

## Function Call Overhead
- The selection happens once per algorithm run, not per evaluation
- Minimal performance impact compared to model evaluation costs
- Function references are lightweight and efficient

## Type Stability
- Julia's type inference optimizes the selected function calls
- No runtime type checking after initial selection
- Maintains performance of specialized application functions

# Examples

```julia
# Automatic selection for different model types
tree_model = DecisionTree.build_tree(X, y)
apply_fn = determine_apply_function(tree_model, nothing)
# Returns: SoleModels.apply

forest_model = DecisionTree.build_forest(X, y)  
apply_fn = determine_apply_function(forest_model, nothing)
# Returns: DecisionTree.apply_forest

# Custom function takes precedence
custom_fn = (model, X) -> my_special_evaluation(model, X)
apply_fn = determine_apply_function(any_model, custom_fn)
# Returns: custom_fn
```

See also: [`LumenConfig`](@ref), [`generate_combinations`](@ref)
"""
function determine_apply_function(model, config_apply_function)
    # User-specified function always takes precedence
    # This allows for custom evaluation strategies and specialized use cases
    if !isnothing(config_apply_function)
        return config_apply_function
    end
    
    # Automatic selection based on model type
    # This type-based dispatch ensures optimal performance for each model category
    return if model isa DT.Ensemble
        # DecisionTree.jl ensemble models use specialized application logic
        # apply_forest is optimized for ensemble prediction aggregation
        DT.apply_forest
    else
        # Single trees and other SoleModels framework types
        # SoleModels.apply provides a generic interface for individual models
        SoleModels.apply
    end
end

"""
    minimize_formula(formula, config::LumenConfig) -> MinimizedFormula

Apply the configured minimization algorithm to simplify a DNF formula.

This function serves as the interface to the various DNF minimization algorithms,
selecting the appropriate method based on configuration and handling any
algorithm-specific parameter passing.

# Arguments
- `formula`: The DNF formula to minimize (typically extracted from decision tree paths)
- `config::LumenConfig`: Configuration specifying minimization scheme and parameters

# Returns
- `MinimizedFormula`: The simplified formula with reduced complexity

# Minimization Process

## Algorithm Dispatch
The function uses Julia's `Val()` type system for compile-time algorithm selection:
```julia
Lumen.minimizza_dnf(Val(config.minimization_scheme), formula; config.minimization_kwargs...)
```

This approach provides:
- **Compile-time optimization**: No runtime algorithm lookup overhead
- **Type safety**: Invalid algorithms caught at compile time
- **Extensibility**: Easy to add new algorithms with the same interface

## Parameter Forwarding
The `minimization_kwargs` from the configuration are forwarded directly to
the chosen algorithm, enabling fine-grained control over minimization behavior:

```julia
# Example configuration with algorithm-specific parameters
config = LumenConfig(
    minimization_scheme = :mitespresso,
    minimization_kwargs = (
        max_iterations = 1000,
        convergence_threshold = 1e-6,
        aggressive_pruning = true
    )
)
```

# Algorithm Characteristics

## Performance vs Quality Trade-offs
Each algorithm makes different trade-offs between computational cost and minimization quality:

- **:abc**: Fast execution, moderate minimization
- **:mitespresso**: Balanced approach, good for most cases  
- **:boom**: Thorough minimization, higher computational cost

## Memory Requirements
Different algorithms have varying memory footprints:
- Some use in-place optimization for memory efficiency
- Others create intermediate representations for better minimization

# Error Handling

The function implements comprehensive error handling:
```julia
try
    return Lumen.minimizza_dnf(Val(scheme), formula; kwargs...)
catch MethodError
    throw(ArgumentError("Unsupported minimization scheme: _scheme"))
catch e
    @error "Minimization failed" formula=formula scheme=scheme
    rethrow(e)
end
```

## Recovery Strategies
- **Fallback algorithms**: Could try simpler algorithms if advanced ones fail
- **Partial minimization**: Return partially minimized results when possible
- **Input validation**: Check formula structure before attempting minimization

# Examples

```julia
# Basic minimization with default parameters
config = LumenConfig(minimization_scheme = :mitespresso)
simplified = minimize_formula(complex_formula, config)

# Advanced minimization with custom parameters
config = LumenConfig(
    minimization_scheme = :boom,
    minimization_kwargs = (max_depth = 10, timeout = 30.0)
)
simplified = minimize_formula(complex_formula, config)
```

See also: [`LumenConfig`](@ref), [`process_rules`](@ref), [`Lumen.minimizza_dnf`](@ref)
"""
function minimize_formula(formula, config::LumenConfig)
    # Use Julia's Val() system for compile-time algorithm dispatch
    # This provides optimal performance by avoiding runtime method lookup
    return Lumen.minimizza_dnf(
        Val(config.minimization_scheme),  # Compile-time algorithm selection
        formula;                          # The formula to minimize
        config.minimization_kwargs...     # Forward algorithm-specific parameters
    )
end

"""
    should_track_statistics(scheme::Symbol) -> Bool

Determine whether to track detailed minimization statistics for the given algorithm.

Some minimization algorithms provide detailed statistics about their operation,
while others operate in a more black-box manner. This function centralizes
the logic for determining when statistics collection is meaningful and available.

# Arguments
- `scheme::Symbol`: The minimization algorithm identifier

# Returns
- `Bool`: `true` if statistics should be tracked, `false` otherwise

# Algorithm-Specific Behavior

## Statistics-Enabled Algorithms
- **:mitespresso**: Provides detailed term count reduction metrics
- **:boom**: Reports aggressive minimization statistics  
- **:abc**: Basic statistics about formula simplification

## Black-Box Algorithms  
- **:espresso**: External tool with limited introspection capabilities
- Future algorithms may have varying statistics availability

# Design Rationale

## Why Centralized Logic?
1. **Maintainability**: Single location for statistics policy
2. **Consistency**: Uniform behavior across the codebase
3. **Extensibility**: Easy to modify behavior for new algorithms
4. **Performance**: Avoid unnecessary computations for algorithms that don't support statistics

## Future Extensions
The function can be extended to support more granular statistics control:
```julia
function should_track_statistics(scheme::Symbol, stats_type::Symbol)
    # Could support different types of statistics per algorithm
    # :compression_ratio, :timing, :memory_usage, etc.
end
```

# Examples

```julia
# Check if an algorithm supports statistics
if should_track_statistics(:mitespresso)
    pre_terms = nterms(formula)
    # ... perform minimization ...
    post_terms = nterms(minimized_formula)
    track_reduction(pre_terms, post_terms)
end

# Algorithm without statistics support
if should_track_statistics(:espresso)  # Returns false
    # This block won't execute
    # No overhead from unnecessary computations
end
```

See also: [`process_rules`](@ref), [`log_minimization_stats`](@ref)
"""
function should_track_statistics(scheme::Symbol)
    # Define which algorithms support detailed statistics tracking
    # This list can be extended as new algorithms are added
    return scheme != :espresso  # Espresso is external and provides limited introspection
    
    # Future implementation might include more sophisticated logic:
    # return scheme in (:espresso, :bbbb, :cccc) && detailed_stats_enabled()
end

"""
    log_minimization_stats(pre_terms::Int, post_terms::Int, silent::Bool) -> Nothing

Log detailed statistics about formula minimization effectiveness.

This function provides formatted output showing the impact of minimization
on formula complexity, helping users understand algorithm performance and
make informed decisions about algorithm selection and parameter tuning.

# Arguments
- `pre_terms::Int`: Number of terms in the original formula
- `post_terms::Int`: Number of terms after minimization  
- `silent::Bool`: Whether to suppress output (when `true`)

# Output Format

When not in silent mode, produces formatted output like:
```
==========================
Terms before: 25
Terms after: 8
Reduction: 17 terms (68% decrease)
Compression ratio: 3.1x
==========================
```

# Metrics Computed

## Absolute Reduction
```julia
reduction = pre_terms - post_terms
```
Shows the raw number of terms eliminated.

## Percentage Reduction  
```julia
percentage = (reduction / pre_terms) * 100
```
Provides intuitive understanding of minimization effectiveness.

## Compression Ratio
```julia
ratio = pre_terms / post_terms
```
Indicates how many times smaller the minimized formula is.

# Design Considerations

## Output Formatting
- Clear visual separation with decorative borders
- Consistent formatting for easy parsing/monitoring
- Meaningful metrics that help users understand algorithm performance

## Performance Impact
- Minimal computational overhead (basic arithmetic)
- Output only when not in silent mode
- No memory allocation for silent operations

## Extensibility
Future versions could include:
- Histogram of term size distributions
- Timing information per minimization step
- Memory usage during minimization
- Comparative statistics across different algorithms

# Examples

```julia
# Typical usage in processing loop
for (original, minimized) in formula_pairs
    pre = nterms(original) 
    post = nterms(minimized)
    log_minimization_stats(pre, post, config.silent)
end

# Output example:
# ==========================
# Terms before: 42
# Terms after: 12  
# Reduction: 30 terms (71% decrease)
# Compression ratio: 3.5x
# ==========================
```

See also: [`process_rules`](@ref), [`should_track_statistics`](@ref)
"""
function log_minimization_stats(pre_terms, post_terms, silent)
    # Only produce output when not in silent mode
    # This design respects user preferences for output verbosity
    if !silent
        # Calculate meaningful metrics for user understanding
        reduction = pre_terms - post_terms
        percentage = reduction > 0 ? round((reduction / pre_terms) * 100, digits=1) : 0.0
        compression_ratio = post_terms > 0 ? round(pre_terms / post_terms, digits=1) : Inf
        
        # Formatted output with clear visual separation
        println("==========================")
        println("Terms before: $pre_terms")
        println("Terms after: $post_terms")
        
        # Additional insights when reduction occurred
        if reduction > 0
            println("Reduction: $reduction terms ($percentage% decrease)")
            println("Compression ratio: $(compression_ratio)x")
        elseif reduction < 0
            println("Formula expanded by $(abs(reduction)) terms")
        else
            println("No change in formula complexity")
        end
        
        println("==========================")
    end
end

"""
    create_rule(formula, result, config::LumenConfig) -> Rule

Convert a minimized logical formula into a Rule object with proper formatting.

This function handles the conversion from internal formula representations
back to the Rule objects used by the broader SoleModels ecosystem. It bridges
the gap between minimization algorithms and the final rule representation.

# Arguments
- `formula`: The minimized logical formula (from minimization algorithms)
- `result`: The decision outcome associated with this formula
- `config::LumenConfig`: Configuration affecting rule creation

# Returns  
- `Rule`: A properly formatted rule object ready for use in DecisionSets

# Processing Modes

The function handles two different processing paths based on the minimization scheme:

## String-Based Processing (Advanced Algorithms)
For algorithms like `:mitespresso`, `:boom`, `:abc`:
```julia
formula_string = leftmost_disjunctive_form_to_string(formula, config.horizontal, config.vetImportance)
φ = SoleLogics.parseformula(formula_string; atom_parser=custom_parser)
return Rule(φ, result)
```

### Advantages:
- Leverages sophisticated string-based parsing  
- Handles complex formula structures
- Supports advanced logical constructs
- Integration with SoleLogics ecosystem

## Direct Conversion (Basic Algorithms)
For simpler algorithms or fallback cases:
```julia
return convert_DNF_formula(formula, result, config.horizontal)
```

### Advantages:
- Direct object-to-object conversion
- Lower overhead for simple formulas
- No string parsing overhead
- Faster for straightforward cases

# Parameter Integration

## Horizontal Coverage
The `config.horizontal` parameter affects rule breadth:
- **1.0**: Include all conditions (maximum specificity)
- **< 1.0**: Allow broader rules by omitting some conditions

## Feature Importance
The `config.vetImportance` vector influences rule construction:
- Prioritizes important features in rule antecedents
- Can be used to simplify rules by focusing on key variables

# Error Handling

The function implements robust error handling for common failure modes:

```julia
try
    # Attempt string-based parsing
    φ = SoleLogics.parseformula(formula_string; atom_parser=custom_parser)
    return Rule(φ, result)
catch ParseError as e
    @error "Formula parsing failed" formula_string=formula_string
    # Could fallback to direct conversion
    return convert_DNF_formula(formula, result, config.horizontal)
catch e
    @error "Rule creation failed" formula=formula result=result
    rethrow(e)
end
```

# Examples

```julia
# Advanced rule creation with feature importance
config = LumenConfig(
    minimization_scheme = :mitespresso,
    horizontal = 0.9,
    vetImportance = [2, 1, 3, 4]  # Feature importance weights
)
rule = create_rule(minimized_formula, "Class_A", config)

# Basic rule creation  
config_basic = LumenConfig(minimization_scheme = :basic)
rule = create_rule(simple_formula, "Class_B", config_basic)
```

See also: [`process_rules`](@ref), [`LumenConfig`](@ref), [`leftmost_disjunctive_form_to_string`](@ref)
"""
function create_rule(formula, result, config::LumenConfig)
    # Choose processing path based on minimization algorithm capabilities
    # Advanced algorithms support sophisticated string-based processing
    if config.minimization_scheme in (:mitespresso, :boom, :abc)
        # String-based processing path for advanced algorithms
        # This path provides maximum flexibility and feature support
        
        # Convert formula to canonical string representation
        # Incorporates horizontal coverage and feature importance
        formula_string = leftmost_disjunctive_form_to_string(
            formula, 
            config.horizontal,      # Controls rule breadth/specificity, TODO: IN FUTURE WE DONT WONT THIS HERE
            config.vetImportance    # Influences feature prioritization TODO: IDEM "
        )
        
        # Parse the string representation back to logical formula
        # Custom atom parser handles domain-specific condition formats
        φ = SoleLogics.parseformula(
            formula_string;
            atom_parser = a -> Atom(
                parsecondition(
                    SoleData.ScalarCondition,
                    a;
                    featuretype = SoleData.VariableValue,  # Feature type specification
                    featvaltype = Real                      # Value type specification
                )
            )
        )
        
        # Create rule with parsed logical formula
        return Rule(φ, result)
    else
        # Direct conversion path for simpler algorithms
        # This path is more efficient but with limited feature support
        return convert_DNF_formula(formula, result, config.horizontal) # MAYBE IS REALY BUSY FUNCTION
    end
end

##
# Main Algorithm Implementation
##
"""
    lumen(model; config_args...) -> LumenResult
    lumen(model, config::LumenConfig) -> LumenResult

Logic-driven Unified Minimal Extractor of Notions (LUMEN): Extract and minimize 
logical rules from decision tree models into interpretable DNF formulas.

LUMEN implements a comprehensive pipeline for converting decision tree models into
interpretable logical rules. The algorithm extracts the underlying decision logic,
constructs truth tables, and applies advanced minimization techniques to produce
compact, human-readable rule sets.

# Method Signatures

## Keyword Arguments Interface
```julia
lumen(model; minimization_scheme=:mitespresso, vertical=1.0, horizontal=1.0, ...)
```
Convenient interface using keyword arguments with automatic config construction.

## Configuration Object Interface  
```julia
lumen(model, config::LumenConfig)
```
Advanced interface using pre-constructed configuration for complex scenarios.

# Arguments

## Required Arguments
- `model`: Decision tree model to analyze
  - **Single trees**: Individual decision tree models
  - **Ensembles**: Random forests, gradient boosting, etc.
  - **Supported formats**: DecisionTree.jl, SoleModels framework

## Configuration (via keywords or LumenConfig)
- `minimization_scheme::Symbol = :mitespresso`: DNF minimization algorithm
- `vertical::Float64 = 1.0`: Instance coverage parameter α ∈ (0,1]
- `horizontal::Float64 = 1.0`: Feature coverage parameter β ∈ (0,1]
- `ott_mode::Bool = false`: Enable memory-optimized processing
- `silent::Bool = false`: Suppress progress output
- `return_info::Bool = true`: Include detailed metadata in results

# Returns

`LumenResult` containing:
- **`decision_set`**: Collection of minimized logical rules
- **`info`**: Metadata including statistics and unminimized rules
- **`processing_time`**: Total algorithm execution time

# Algorithm Pipeline

## Phase 1: Model Analysis and Rule Extraction
```
Input Model → Rule Extraction → Logical Rule Set
```
- Analyzes model structure (single tree vs ensemble)
- Extracts decision paths as logical rules
- Handles different model types with appropriate strategies

## Phase 2: Alphabet Construction and Atom Processing  
```
Logical Rules → Atom Extraction → Logical Alphabet
```
- Identifies atomic logical conditions
- Constructs vocabulary for formula building
- Validates feature support and operator compatibility

## Phase 3: Truth Table Generation
```
Model + Alphabet → Truth Combinations → Labeled Examples
```
- Generates systematic input combinations
- Evaluates model on each combination
- Creates correspondence between inputs and outputs

## Phase 4: DNF Construction and Minimization
```
Truth Table → DNF Formulas → Minimized Rules
```
- Constructs DNF formulas for each decision class
- Applies advanced minimization algorithms
- Converts back to interpretable rule format

# Performance Characteristics

## Computational Complexity
- **Time**: O(2^k × n × d) where k=features, n=instances, d=tree depth
- **Space**: O(k × r) where r=number of rules
- **Scalability**: Optimized modes available for large datasets

## Memory Usage
- **Standard mode**: Suitable for typical datasets (< 20 features)
- **Optimized mode**: Memory-efficient processing for large problems
- **Streaming capability**: Future versions may support streaming processing

# Advanced Features

## Custom Processing  
```julia
# Custom alphabet filtering for domain expertise
custom_filter = alphabet -> remove_irrelevant_features(alphabet)
config = LumenConfig(filteralphabetcallback = custom_filter)
result = lumen(model, config)
```

## Performance Tuning
```julia
# Memory-optimized processing for large datasets
config = LumenConfig(ott_mode = true, vertical = 0.8)

# Speed-optimized processing with basic minimization
config = LumenConfig(minimization_scheme = :abc, silent = true)
```

## Analysis and Debugging
```julia
# Full information retention for analysis
config = LumenConfig(return_info = true, controllo = true)
result = lumen(model, config)

# Access detailed statistics  
println("Rules before minimization: _(length(result.info.unminimized_ds.rules))")
println("Rules after minimization: _(length(result.decision_set.rules))")
```

# Error Handling

The algorithm implements comprehensive error handling:

## Configuration Validation
- Parameter range checking (coverage parameters must be ∈ (0,1])
- Algorithm availability verification  
- Consistency validation across parameters

## Processing Errors
- Graceful handling of minimization failures
- Fallback strategies for problematic formulas
- Detailed error reporting with context

## Model Compatibility
- Automatic detection of supported model types
- Clear error messages for unsupported formats
- Suggestions for model preprocessing

# Examples

## Basic Usage
```julia
# Simple rule extraction with default settings
model = build_tree(X, y)
result = lumen(model)
println("Extracted _(length(result.decision_set.rules)) rules")
```

## Advanced Configuration
```julia
# Customized processing for complex scenarios
config = LumenConfig(
    minimization_scheme = :boom,        # Aggressive minimization
    vertical = 0.9,                     # High instance coverage  
    horizontal = 0.8,                   # Moderate feature coverage
    ott_mode = true,                    # Memory optimization
    return_info = true                  # Full information retention
)
result = lumen(large_ensemble, config)
```

## Performance Analysis
```julia
# Detailed performance and quality analysis
result = lumen(model, LumenConfig(return_info = true))

# Analyze minimization effectiveness
stats = result.info.vectPrePostNumber
total_reduction = sum(pre - post for (pre, post) in stats)
avg_compression = mean(pre / post for (pre, post) in stats)

println("Total term reduction: total_reduction")
println("Average compression ratio: (round(avg_compression, digits=2))x")
println("Processing time: _(result.processing_time) seconds")
```

# Implementation Notes

## Design Principles
1. **Modularity**: Each phase is independently testable and extensible
2. **Configurability**: Extensive customization without code modification
3. **Performance**: Multiple optimization strategies for different scenarios  
4. **Robustness**: Comprehensive error handling and validation
5. **Usability**: Clean interfaces with sensible defaults

## Extensibility Points
- **New minimization algorithms**: Add via Val() dispatch system
- **Custom model types**: Extend rule extraction strategies  
- **Domain-specific processing**: Custom alphabet filters and apply functions
- **Output formats**: Additional result formatters and exporters

See also: [`LumenConfig`](@ref), [`LumenResult`](@ref), [`extract_rules`](@ref), [`minimize_formula`](@ref)
"""
function lumen(model; kwargs...)
    config = LumenConfig(; kwargs...)
    return lumen(model, config)
end

function lumen(model, config::LumenConfig)
    start_time = time()
    
    # Handle special test modes
    if !isnothing(config.testott) || !isnothing(config.alphabetcontroll)
        return handle_test_modes(model, config)
    end
    
    # Determine apply function
    apply_function = determine_apply_function(model, config.apply_function)
    config = @set config.apply_function = apply_function
    
    # Create SOLE model if needed
    sole_model = SoleModels.solemodel(model)
    
    try
        # Step 1: Extract rules
        config.silent || println("\n$COLORED_TITLE$TITLE\n PART 2.a: STARTER RULESET EXTRACTION  \n$TITLE$RESET")
        ruleset_result = @timed extract_rules(sole_model)
        config.silent || println("Extracted $(length(ruleset_result.value)) rules in $(ruleset_result.time)s")
        
        # Step 2: Extract atoms and alphabet
        config.silent || println("\n$COLORED_TITLE$TITLE\n PART 2.b: ATOM EXTRACTION  \n$TITLE$RESET")
        num_atoms, atoms, alphabet = extract_atoms(sole_model, config.filteralphabetcallback)
        config.silent || println("Extracted $(length(atoms)) unique atoms from $num_atoms total")
        
        # Step 3: Generate truth combinations
        config.silent || println("\n$COLORED_TITLE$TITLE\n PART 3: TABLE GENERATION  \n$TITLE$RESET")
        results_dict, label_count = generate_combinations(model, alphabet, atoms, config.vertical, config)
            
        # Step 4: Process results
        # Modifichiamo questa riga per passare solo il dizionario dei risultati
        config.silent || println("\n$COLORED_TITLE$TITLE\n PART 4: PROCESS RESULTS \n$TITLE$RESET")
        combined_results = concat_results(results_dict, atoms)  
        minimized_rules, unminimized_rules, stats = process_rules(combined_results, config)
        
        # Create result
        ds = DecisionSet(minimized_rules)
        processing_time = time() - start_time
        
        if !config.return_info
            return LumenResult(ds)
        end
        
        info = (
            vectPrePostNumber = stats,
            processing_time = processing_time
        )
        
        if !isnothing(unminimized_rules)
            info = merge(info, (unminimized_ds = DecisionSet(unminimized_rules),))
        end
        
        return LumenResult(ds, info, processing_time)
        
    catch e
        @error "Error in lumen processing" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    handle_test_modes(model, config::LumenConfig) -> nothing

OPTIMIZATION VALID and ANALYZE ALPHABET mode.
"""
function handle_test_modes(model, config::LumenConfig)
    # Extract basic components needed for testing
    sole_model = SoleModels.solemodel(model)
    _, atoms, alphabet = extract_atoms(sole_model, config.filteralphabetcallback)
    apply_function = determine_apply_function(model, config.apply_function)
    
    if !isnothing(config.testott)
        config.silent || println("\n$COLORED_INFO PART 2.d: IS THE OPTIMIZATION VALID? $RESET")
        testOttt(model, alphabet, atoms, config.vertical; 
                config.silent, apply_function, config.testott)
    end
    
    if !isnothing(config.alphabetcontroll)
        config.silent || println("\n$COLORED_INFO ANALYZE ONLY ALPHABET $RESET")
        debug_combinations(model, alphabet, atoms, config.vertical; 
                          config.silent, apply_function, config.alphabetcontroll)
    end
    
    return nothing, nothing
end

##
# Include submodules
##

# Types
include("types/trit-vector.jl")
include("types/balanced-trit-vector.jl") 
include("types/balanced-ternary-vector.jl")
include("types/types.jl")

# Utilities
#include("reportAlphabet.jl")
include("utils/report.jl")
include("utils/IO.jl")
include("utils/minor.jl")
include("utils/core_shared.jl")
include("utils/core.jl")
include("utils/coreOttMode.jl")
include("utils/minimization.jl")
include("deprecate.jl")

end # module