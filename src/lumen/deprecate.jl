
"""
Generate an AND formula based on the provided binary `combination` representing a combination of atoms. Calculate conditions for each feature based on the thresholds and atom properties. Construct atoms with corresponding conditions and return the resulting AND formula.
"""
function generate_disjunct_old(
    combination::BigInt,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    comb, _ =
        process_combination(combination, num_atoms, thresholds_by_feature, atoms_by_feature)
    atoms = Vector{Atom}()

    for (feat, values) in comb
        if haskey(atoms_by_feature, feat)
            for (threshold, _) in atoms_by_feature[feat]
                if values[1] < threshold
                    mc = ScalarMetaCondition(VariableValue(feat), <)
                else
                    mc = ScalarMetaCondition(VariableValue(feat), ≥)
                end
                condition = ScalarCondition(mc, threshold)
                push!(atoms, Atom(condition))
            end
        end
    end
    return isempty(atoms) ? ⊤ : ∧(atoms...)
end


#=
    TODO VALUTARE CORRETTEZZA
=#
function generate_disjunct_all(
    combination::BitVector,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    comb, _ =
        process_combination(combination, num_atoms, thresholds_by_feature, atoms_by_feature)
    atoms = Vector{Atom}()

    for (feat, values) in comb
        if haskey(atoms_by_feature, feat)
            for (threshold, _) in atoms_by_feature[feat]
                # Creiamo un atomo per ogni condizione, senza filtrare
                if values[1] < threshold
                    # Condizione "<"
                    mc = ScalarMetaCondition(VariableValue(feat), <)
                    condition = ScalarCondition(mc, threshold)
                    push!(atoms, Atom(condition))
                else
                    # Condizione "≥"
                    mc = ScalarMetaCondition(VariableValue(feat), ≥)
                    condition = ScalarCondition(mc, threshold)
                    push!(atoms, Atom(condition))
                end
            end
        end
    end

    return isempty(atoms) ? ⊤ : ∧(atoms...)
end

# #= TODO COSTATARE CORRETTEZZA\UTILITA =#
# #= TODO SPEZZARE IN DUE E RINOMINARE  =#
# """
#     IO_print_custom_or_formula_from_mask(formula::TwoLevelDNFFormula, horizontal::Float64) !!DEPRECATE!!

# Generate a custom OR formula from a mask.

# This function takes a `TwoLevelDNFFormula` and a horizontal threshold value, and generates a set of rows representing the custom OR formula. The function applies a masking process to the formula's atoms, keeping only the most stringent atoms for each feature based on the horizontal threshold.

# Parameters:
# - `formula::TwoLevelDNFFormula`: The custom OR formula to be processed.
# - `horizontal::Float64`: The horizontal threshold value, which determines the number of features to include in the generated rows.

# Returns:
# - A vector of strings, where each string represents a row in the custom OR formula.
# """
# function IO_print_custom_or_formula_from_mask(
#     formula::TwoLevelDNFFormula,
#     horizontal::Float64,
# ) # ex generate_custom_or_formula_from_mask
#     all_rows = String[]
#     current_row = String[]

#     num_orizontal = (floor(formula.thresholds_by_feature.count * horizontal))

#     # Script per controllare e impostare i bit non stringenti a -1
#     for (j, mask) in enumerate(formula.prime_mask)
#         i = 1
#         for (feature, atoms) in formula.atoms_by_feature
#             if i <= length(atoms)
#                 # Troviamo tutti gli indici degli atomi per questa feature
#                 feature_indices = []
#                 num_atoms = 0
#                 for (value, is_greater_equal) in atoms
#                     if mask[i] != -1
#                         push!(feature_indices, (i, value, mask[i]))
#                         num_atoms += 1
#                     end
#                     i += 1
#                 end

#                 # Se abbiamo più di un atomo per questa feature
#                 if length(feature_indices) > 1
#                     # Ordiniamo gli indici per valore
#                     sort!(feature_indices, by = x -> x[2])

#                     # Per i '<', manteniamo solo il più piccolo
#                     min_less = nothing
#                     for (idx, val, bit) in feature_indices
#                         if bit == 1  # 1 rappresenta '<'
#                             if min_less === nothing || val < min_less[2]
#                                 min_less = (idx, val, bit)
#                             end
#                         end
#                     end

#                     # Per i '≥', manteniamo solo il più grande
#                     max_greater = nothing
#                     for (idx, val, bit) in feature_indices
#                         if bit == 0  # 0 rappresenta '≥'
#                             if max_greater === nothing || val > max_greater[2]
#                                 max_greater = (idx, val, bit)
#                             end
#                         end
#                     end

#                     # Impostiamo a -1 tutti tranne i più stringenti
#                     for (idx, val, bit) in feature_indices
#                         if bit == 1 && (min_less === nothing || idx != min_less[1])
#                             mask[idx] = -1
#                         elseif bit == 0 &&
#                                (max_greater === nothing || idx != max_greater[1])
#                             mask[idx] = -1
#                         end
#                     end
#                 end
#             end
#         end
#     end

#     for (j, mask) in enumerate(formula.prime_mask)
#         current_row = String[]
#         i = 1

#         for (feature, atoms) in formula.atoms_by_feature
#             if i <= length(atoms)
#                 for (value, is_greater_equal) in atoms
#                     if mask[i] != -1
#                         if feature <= num_orizontal #TODO HA SENSO ? FIXME NON HA SENSO ? 
#                             symbol = Bool(mask[i]) ? "<" : "≥"
#                             formatted_value = string(value)
#                             atom_str = "V$(feature) $(symbol) $(formatted_value)"
#                             push!(current_row, atom_str)
#                         end
#                     end
#                     i += 1
#                 end
#             end
#         end
#         if !isempty(current_row)
#             push!(all_rows, join(current_row, " ∧ "))
#         end
#     end
#     println(join(all_rows, "\n"))
# end

# DONT WORKED print_filtered_dnf
function print_filtered_dnf(f::TwoLevelDNFFormula)
    # Get all conditions
    all_conditions = conditions(f)
    
    # Initialize result string
    result = ""
    
    # Process each combination with its corresponding mask
    for (i, (combination, mask)) in enumerate(zip(eachcombination(f), f.prime_mask))
        # Skip if this is an empty combination
        if isempty(combination)
            continue
        end
        
        # Add OR separator if this isn't the first term
        if i > 1 && !isempty(result)
            result *= " ∨ "
        end
        
        # Start this conjunction term
        result *= "("
        
        # Track if we've added any atoms to this conjunction
        added_atoms = false
        
        # Process each atom position
        for (j, (bit, mask_val)) in enumerate(zip(combination, mask))
            # Only process if mask value is not -1
            if mask_val != -1
                # Add AND separator if this isn't the first atom in this conjunction
                if added_atoms
                    result *= " ∧ "
                end
                
                # Get the corresponding condition
                condition = all_conditions[j]
                
                if condition isa ScalarCondition
                    feature = condition.metacond.feature.i_variable
                    threshold = condition.threshold
                    
                    # Use bit value to determine operator (0 -> <, 1 -> ≥)
                    op_str = bit == 0 ? "<" : "≥"
                    
                    # Add the formatted condition
                    result *= "x$feature $op_str $threshold"
                    added_atoms = true
                end
            end
        end
        
        # Close this conjunction term
        result *= ")"
    end
    
    # Handle empty result
    if isempty(result)
        return "⊤"  # Return top symbol for empty formula
    end
    
    return result
end