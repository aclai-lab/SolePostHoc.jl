##
# Rule Extraction Strategies
##

"""
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
