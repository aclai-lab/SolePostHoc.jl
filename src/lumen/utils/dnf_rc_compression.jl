using Random

# =============================================================================
# DNF FORMULA SAMPLING UTILITIES - FIXED VERSION
# =============================================================================
# This module provides utilities for sampling Disjunctive Normal Form (DNF) 
# formulas by reducing both the number of disjuncts (horizontal sampling) and 
# the number of variables based on importance (vertical sampling).
# =============================================================================

"""
    get_variable_id_from_literal(literal) -> Int

Extract the variable ID from a SoleLogics literal.

# Arguments
- `literal`: A SoleLogics Atom literal object

# Returns
- Variable ID as an integer
"""
function get_variable_id_from_literal(literal)
    # Fix: literal is already an Atom, no need for .atom
    return literal.value.metacond.feature.i_variable
end

"""
    get_literal_info(literal) -> Tuple{Int, Any, Any}

Extract comprehensive information from a SoleLogics literal.

# Arguments
- `literal`: A SoleLogics Atom literal object

# Returns
- `Tuple{var_id, operator, threshold}`: Complete literal information
"""
function get_literal_info(literal)
    # Fix: literal is already an Atom, access .value directly
    atom_value = literal.value
    var_id = atom_value.metacond.feature.i_variable
    operator = atom_value.metacond.test_operator
    threshold = atom_value.threshold
    return (var_id, operator, threshold)
end

"""
    safe_get_variable_id_from_literal(literal) -> Union{Int, Nothing}

Safely extract the variable ID from a SoleLogics literal with error handling.

# Arguments
- `literal`: A SoleLogics literal object

# Returns
- Variable ID as an integer, or Nothing if extraction fails
"""
function safe_get_variable_id_from_literal(literal)
    try
        # Try direct access first (literal is Atom)
        return literal.value.metacond.feature.i_variable
    catch e1
        try
            # Fallback: try nested access (literal has .atom field)
            return literal.atom.value.metacond.feature.i_variable
        catch e2
            @warn "Could not extract variable ID from literal: $literal"
            @warn "Error 1 (direct access): $e1"
            @warn "Error 2 (nested access): $e2"
            return nothing
        end
    end
end

"""
    dnf_rc_compression(formula, r, c, importance_vector) -> Formula

Sample a DNF formula using two-dimensional reduction:
- **Horizontal sampling**: Reduces the number of disjuncts (OR clauses)
- **Vertical sampling**: Removes literals based on variable importance

# Arguments
- `formula`: Input DNF formula
- `r::Float64`: Fraction of disjuncts to keep ∈ (0,1]
- `c::Float64`: Fraction of variables to keep ∈ (0,1] 
- `importance_vector::Vector{Int}`: Variable importance ranks (1=most important, n=least important)

# Returns
- Sampled DNF formula or `nothing` if no valid formula remains

# Example
```julia
# Variable importance: V2 > V3 > V1 > V4
importance = [3, 1, 2, 4]  # Ranks for V1, V2, V3, V4
sampled = dnf_rc_compression(formula, 0.7, 0.5, importance)
```
"""
function dnf_rc_compression(
    formula, 
    c::Float64, 
    r::Float64, 
    importance_vector::Vector{Int};
    silent = true
)
    # Input validation
    @assert 0 < r <= 1 "r must be in (0,1]"
    @assert 0 < c <= 1 "c must be in (0,1]"
    
    original_disjuncts = formula.grandchildren
    n_original_disjuncts = length(original_disjuncts)
    
    silent || println("Original formula: $(n_original_disjuncts) disjuncts")
    
    # ==========================================================================
    # PHASE 1: HORIZONTAL SAMPLING - Reduce number of disjuncts
    # ==========================================================================
    n_keep_disjuncts = max(1, round(Int, n_original_disjuncts * r))
    
    # Random sampling without external dependencies
    disjunct_indices = collect(1:n_original_disjuncts)
    shuffle!(disjunct_indices)
    selected_indices = disjunct_indices[1:n_keep_disjuncts]
    horizontally_sampled = original_disjuncts[selected_indices]
    
    silent || println("After horizontal sampling ($(r)): $(n_keep_disjuncts) disjuncts")
    
    # ==========================================================================
    # PHASE 2: VERTICAL SAMPLING - Remove less important variables
    # ==========================================================================
    n_variables = length(importance_vector)
    n_keep_variables = max(1, round(Int, n_variables * c))
    
    # Sort variable indices by importance (ascending rank = descending importance)
    importance_order = importance_vector
    variables_to_keep = Set(importance_vector[1:n_keep_variables])

    @show variables_to_keep
    
    silent || println("Variables to keep: $(sort(collect(variables_to_keep))) (top $(n_keep_variables) most important)")
    
    # ==========================================================================
    # PHASE 3: FILTER LITERALS - Keep only literals with important variables
    # ==========================================================================
    filtered_disjuncts = []
    skipped_literals = 0
    
    for disjunct in horizontally_sampled
        filtered_literals = []
        
        # Filter literals containing only important variables
        for literal in disjunct.grandchildren
            variable_id = safe_get_variable_id_from_literal(literal)
            
            if variable_id === nothing
                skipped_literals += 1
                @warn "Skipping literal due to variable ID extraction failure"
                continue
            end
            
            if variable_id in variables_to_keep
                push!(filtered_literals, literal)
            end
        end
        
        # Keep disjuncts that have at least one literal after filtering
        if !isempty(filtered_literals)
            new_disjunct = typeof(disjunct)(filtered_literals)
            push!(filtered_disjuncts, new_disjunct)
        end
    end
    
    if skipped_literals > 0
        @warn "Skipped $skipped_literals literals due to extraction errors"
    end
    
    silent || println("Final disjuncts after vertical filtering: $(length(filtered_disjuncts))")
    
    # ==========================================================================
    # PHASE 4: CONSTRUCT NEW FORMULA
    # ==========================================================================
    if isempty(filtered_disjuncts)
        @warn " No disjuncts remaining after filtering!"
        return nothing
    end
    
    # Create new formula with same structure as original
    silent || println("pre unique",filtered_disjuncts)
    unique!(filtered_disjuncts)
    silent || println("post unique",filtered_disjuncts)
    new_formula = typeof(formula)(filtered_disjuncts)
    
    return new_formula
end
