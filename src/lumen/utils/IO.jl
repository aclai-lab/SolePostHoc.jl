"""
Prints a summary of the estimated time for a computation.

Args:
    total_estimated_time::BigFloat: The total estimated time for the computation, in seconds.

This function prints a message indicating the estimated time for the computation, formatted as seconds, minutes, or hours as appropriate.
"""
function print_estimated_time(total_estimated_time::BigFloat)
    if total_estimated_time < 60
        println("Tempo stimato: circa $(round(total_estimated_time, digits=2)) secondi")
    elseif total_estimated_time < 3600
        println("Tempo stimato: circa $(round(total_estimated_time / 60, digits=2)) minuti")
    else
        println("Tempo stimato: circa $(round(total_estimated_time / 3600, digits=2)) ore")
    end
end


"""
    print_results_summary(real_time_per_combination::BigFloat, contradictions::Int64, num_combinations::BigInt, label_count::Dict{Any,Int64})

Prints a summary of the results, including the total number of combinations, the number of valid and invalid combinations, and the distribution of labels.

Args:
    real_time_per_combination::BigFloat: The precision of the time per combination.
    contradictions::Int64: The number of logical contradictions found.
    num_combinations::BigInt: The total number of combinations.
    label_count::Dict{Any,Int64}: A dictionary mapping labels to their counts.
"""
function print_results_summary(
    real_time_per_combination::BigFloat,
    contradictions::Int64,
    num_combinations::BigInt,
    label_count::Dict{Any,Int64},
)
    valid_combinations = num_combinations - BigInt(contradictions)
    valid_percentage = (Float64(valid_combinations) / Float64(num_combinations)) * 100
    invalid_percentage = (Float64(contradictions) / Float64(num_combinations)) * 100

    println()
    println("Precisione del tempo: $real_time_per_combination")
    println("Generazione delle combinazioni completata.")
    println("Totale combinazioni: $num_combinations")
    println("Totale contraddizioni logiche trovate: $contradictions su $num_combinations")
    println("Totale combinazioni valide: $valid_combinations")

    # Aggiungere la percentuale delle valide e non valide
    println("Percentuale combinazioni valide: $(round(valid_percentage, digits=2))%")
    println("Percentuale combinazioni non valide: $(round(invalid_percentage, digits=2))%")

    println("\nDistribuzione delle etichette:")
    for (label, count) in sort(collect(label_count), by = x -> x[2], rev = true)
        percentage = (Float64(count) / Float64(valid_combinations)) * 100
        println("$label: $count ($(round(percentage, digits=2))%)")
    end
end


"""
    print_detailed_results(results::Dict{Any,Vector{BigInt}})

Prints detailed results of the computation, including the number of combinations for each result and the first 10 combinations (or all combinations if there are 10 or fewer).

Args:
    results::Dict{Any,Vector{BigInt}}: A dictionary mapping results to the corresponding combinations.
"""
function print_detailed_results(results::Dict{Any,Vector{BigInt}})
    println("\nRisultati dettagliati:")
    for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
        println("[$result] ($(length(combinations)) combinazioni):")
        if length(combinations) > 10
            println("  Prime 10 combinazioni: $(combinations[1:10])")
            println("  ... ($(length(combinations)-10) altre combinazioni)")
        else
            println("  Combinazioni: $combinations")
        end
    end
end


############################################################################################
############################################################################################
function io_formula_mask(formula::TwoLevelDNFFormula, orizontal::Float64 = 1.0)
    num_orizontal = floor(Int, formula.thresholds_by_feature.count * orizontal)
    result = Set{String}()  # Usando Set invece di Array per evitare duplicati

    for mask in formula.prime_mask
        feature_conditions = Dict{Int,Dict{String,Float64}}()
        atom_idx = 1

        # First pass: collect all conditions by feature
        for (feature, atoms) in formula.atoms_by_feature
            if atom_idx <= length(atoms)
                for (value, is_greater_equal) in atoms
                    if mask[atom_idx] != -1 && feature <= num_orizontal
                        op_type = mask[atom_idx] == 1 ? "<" : "≥"

                        if !haskey(feature_conditions, feature)
                            feature_conditions[feature] = Dict{String,Float64}()
                        end

                        if op_type == "≥"
                            # Keep maximum for ≥
                            if !haskey(feature_conditions[feature], "≥") ||
                               value > feature_conditions[feature]["≥"]
                                feature_conditions[feature]["≥"] = value
                            end
                        else
                            # Keep minimum for <
                            if !haskey(feature_conditions[feature], "<") ||
                               value < feature_conditions[feature]["<"]
                                feature_conditions[feature]["<"] = value
                            end
                        end
                    end
                    atom_idx += 1
                end
            end
        end

        # Second pass: build the formula string with most restrictive conditions
        current_atoms = String[]
        for (feature, conditions) in feature_conditions
            for (op, value) in conditions
                push!(current_atoms, "V$(feature) $(op) $(value)")
            end
        end

        if !isempty(current_atoms)
            push!(result, join(sort(current_atoms), " ∧ "))
        end
    end

    return collect(result)  # Convertiamo il Set in Array per il risultato finale
end

function printIO_custom_or_formula(
    io::IO,
    formula::TwoLevelDNFFormula,
    orizontal::Float64 = 1.0,
)
    formulas = io_formula_mask(formula, orizontal)
    println(io, "OR Formula with $(length(formulas)) combinations:")
    for (i, formula_str) in enumerate(formulas)
        println(io, "  OR[$i]: $formula_str")
    end
end