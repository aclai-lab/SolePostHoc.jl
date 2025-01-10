using Random

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

function io_formula_mask(formula::TwoLevelDNFFormula, orizontal::Float64 = 1.0)
    num_orizontal = floor(Int, formula.thresholds_by_feature.count * orizontal)
    result = Set{String}()

    if !isnothing(formula.prime_mask) && !isempty(formula.prime_mask)
        # Logica per formule semplificate (con prime_mask)
        for mask in formula.prime_mask
            feature_conditions = Dict{Int,Dict{String,Float64}}()
            current_atom_index = 1

            for (feature, atoms) in formula.atoms_by_feature
                if feature <= num_orizontal
                    for (threshold, _) in atoms
                        if current_atom_index <= length(mask) && mask[current_atom_index] != -1
                            if !haskey(feature_conditions, feature)
                                feature_conditions[feature] = Dict{String,Float64}()
                            end

                            # Se mask[current_atom_index] == 1, cerchiamo x < threshold
                            op = mask[current_atom_index] == 1 ? "<" : "≥"
                            
                            if op == "≥" && (!haskey(feature_conditions[feature], "≥") || 
                                           threshold > feature_conditions[feature]["≥"])
                                feature_conditions[feature]["≥"] = threshold
                            elseif op == "<" && (!haskey(feature_conditions[feature], "<") || 
                                               threshold < feature_conditions[feature]["<"])
                                feature_conditions[feature]["<"] = threshold
                            end
                        end
                        current_atom_index += 1
                    end
                end
            end

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
    else
        # Logica per formule non semplificate
        for term in formula.combinations 
            feature_conditions = Dict{Int,Dict{String,Float64}}()
            current_atom_index = 1
            
            for (feature, atoms) in formula.atoms_by_feature
                if feature <= num_orizontal
                    for (threshold, _) in atoms
                        if current_atom_index <= length(term)
                            is_active = term[current_atom_index]
                            
                            if is_active
                                if !haskey(feature_conditions, feature)
                                    feature_conditions[feature] = Dict{String,Float64}()
                                end
                                
                                # Se is_active (truth_row[j] == 1), cerchiamo x < threshold
                                # come in process_combination
                                op = "<"
                                
                                if !haskey(feature_conditions[feature], op) || 
                                   (op == "<" && threshold < feature_conditions[feature][op]) ||
                                   (op == "≥" && threshold > feature_conditions[feature][op])
                                    feature_conditions[feature][op] = threshold
                                end
                            else
                                # Se !is_active (truth_row[j] == 0), cerchiamo x ≥ threshold
                                if !haskey(feature_conditions, feature)
                                    feature_conditions[feature] = Dict{String,Float64}()
                                end
                                
                                op = "≥"
                                
                                if !haskey(feature_conditions[feature], op) || 
                                   (op == "<" && threshold < feature_conditions[feature][op]) ||
                                   (op == "≥" && threshold > feature_conditions[feature][op])
                                    feature_conditions[feature][op] = threshold
                                end
                            end
                        end
                        current_atom_index += 1
                    end
                end
            end

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
    end

    return collect(result)
end

function convert_DNF_formula(
    formula::TwoLevelDNFFormula,
    outcome,
    orizontal::Float64 = 1.0,
)
    formulas = io_formula_mask(formula, orizontal)    
    result = join(formulas, " ∨ ")
    
    # Specificare esplicitamente featvaltype = Real per risolvere il warning 
    φ = SoleLogics.parseformula(
        result;
        atom_parser = a->Atom(
            parsecondition(
                SoleData.ScalarCondition, 
                a; 
                featuretype = SoleData.VariableValue,
                featvaltype = Real  # Specifichiamo esplicitamente il tipo
            )
        )
    )

    # Creiamo la Rule usando l'outcome passato come parametro
    return Rule(φ, outcome)
end
#=
function convert_DNF_formula(
    formula::TwoLevelDNFFormula,
    outcome,
    orizontal::Float64 = 1.0,
)
    formulas = io_formula_mask(formula, orizontal)
    
    # Gestione caso singola condizione
    if length(formulas) == 1
        φ = Atom(parsecondition(
            SoleData.ScalarCondition, 
            first(formulas); 
            featuretype = SoleData.VariableValue,
            featvaltype = Real
        ))
    else
        result = join(formulas, " ∨ ")
        φ = SoleLogics.parseformula(
            result;
            atom_parser = a->Atom(
                parsecondition(
                    SoleData.ScalarCondition, 
                    a; 
                    featuretype = SoleData.VariableValue,
                    featvaltype = Real
                )
            )
        )
    end

    return Rule(φ, outcome)
end
=#