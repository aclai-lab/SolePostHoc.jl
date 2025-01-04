#= -------------------------------------------------------------------------------------------------------------------------
############################################################################################################################
#                                               OTTIMIZZAZIONE-MOD                                                         #
############################################################################################################################
=#


"""
Determine the optimal number of threads to use for parallel processing based on the number of combinations.

The function calculates the optimal number of threads by considering the maximum number of threads available on the system and a minimum threshold of combinations per thread. It aims to maximize the parallelism while ensuring a minimum workload per thread to avoid overhead.

Args:
    num_combinations (BigInt): The total number of combinations to be processed.

Returns:
    Int: The optimal number of threads to use for parallel processing.
"""
function determine_optimal_threads(num_combinations::BigInt)
    # Numero massimo di thread disponibili sul sistema
    max_threads = Sys.CPU_THREADS

    # Definiamo una soglia minima di combinazioni per thread
    min_combinations_per_thread = BigInt(1000)

    # Calcoliamo il numero ottimale di thread
    optimal_threads = min(
        max_threads,
        max(
            1,
            Int(
                min(num_combinations ÷ min_combinations_per_thread, BigInt(typemax(Int64))),
            ),
        ),
    )

    return optimal_threads
end


"""
Calculate a dynamic chunk size for parallel processing based on the total number of combinations and the number of threads.

The function determines an optimal chunk size that balances the number of chunks with the minimum workload per thread. It ensures the chunk size is within a reasonable range, rounds it to the nearest power of 2 for efficiency, and provides a fallback in case the calculated size is too large.

Args:
    num_combinations (BigInt): The total number of combinations to be processed.
    num_threads (Int): The number of threads to be used for parallel processing.

Returns:
    BigInt: The calculated dynamic chunk size.
"""
function calculate_dynamic_chunk_size(num_combinations::BigInt, num_threads::Int)
    # Define minimum and maximum chunk sizes
    min_chunk_size = BigInt(100)
    max_chunk_size = BigInt(10_000_000)

    # Calculate an initial chunk size
    initial_chunk_size = num_combinations ÷ (num_threads * 10)

    # Ensure the chunk size is within limits
    chunk_size = clamp(initial_chunk_size, min_chunk_size, max_chunk_size)

    # Ensure there are at least as many chunks as threads
    num_chunks = num_combinations ÷ chunk_size
    if num_chunks < num_threads
        chunk_size = num_combinations ÷ num_threads
    end

    # Round the chunk size to the nearest power of 2 for efficiency
    # Use a safe exponent calculation to avoid overflow
    if chunk_size > 0
        safe_exponent = min(63, floor(Int, log2(Float64(chunk_size))))
        chunk_size = BigInt(2)^safe_exponent
    else
        chunk_size = min_chunk_size
    end

    # Final safety check: if chunk_size is still too large, use a constant fallback
    if chunk_size > num_combinations || chunk_size <= 0
        chunk_size = min(num_combinations, BigInt(1_000_000))
    end

    return chunk_size
end


"""
Generates a set of valid combinations of feature values based on the provided model, alphabet, atoms, and vertical parameter. The function uses a parallel processing approach with dynamic chunk sizes to efficiently process the large number of combinations.

Args:
    model (Any): The model to be used for evaluating the combinations.
    alphabet (MultivariateScalarAlphabet{ScalarCondition}): The alphabet containing the feature conditions.
    atoms (Vector{<:Atom{<:ScalarCondition{Float64,<:VariableValue,<:ScalarMetaCondition{<:VariableValue,typeof(<)}}}}): The atoms representing the feature conditions.
    vertical (Float64): The vertical parameter used to scale the number of combinations.
    print_progress (bool): Whether to print progress information during the computation.

Returns:
    Tuple{Dict{Any,Vector{BigInt}}, Dict{Any,Int64}}: A tuple containing the results dictionary (mapping labels to their corresponding combinations) and the label count dictionary.
"""
function truth_combinations_ott(
    model::Any,
    alphabet,
    atoms::Vector{
        <:Atom{
            <:ScalarCondition{
                Float64,
                <:VariableValue,
                <:ScalarMetaCondition{<:VariableValue,typeof(<)},
            },
        },
    },
    vertical::Float64;
    apply_function = SoleModels.apply,
    silent = false,
    print_progress = true,
)
    thresholds_by_feature = Dict(
        subalpha.featcondition[1].feature.i_variable => sort(subalpha.featcondition[2])
        for subalpha in alphabet.subalphabets
    )
    atoms_by_feature = group_atoms_by_feature(atoms)

    num_atoms = length(atoms)
    all_combinations = BigInt(2)^num_atoms
    if isone(vertical)
        num_combinations = all_combinations
    else
        num_combinations = BigInt(round(all_combinations * vertical))
    end

    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int64}()
    contradictions = Atomic{Int64}(0)

    silent || println("$COLORED_ULTRA_OTT$TITLE\n THREADS DINAMICI \n$TITLE$RESET")
    optimal_threads = determine_optimal_threads(num_combinations)

    silent || println("Numero di combinazioni: $num_combinations")
    silent || println("Numero attuale di thread: ", Threads.nthreads())
    silent || println("Numero ottimale di thread: $optimal_threads")

    silent || println("\n\n$COLORED_ULTRA_OTT$TITLE\n CHUNK DINAMICI \n$TITLE$RESET")

    chunk_size = calculate_dynamic_chunk_size(num_combinations, Threads.nthreads())
    num_chunks = (num_combinations + chunk_size - 1) ÷ chunk_size  # Ceiling division

    silent || println("Dynamic chunk size: $chunk_size")
    silent || println("Number of chunks: $num_chunks")

    silent || println("$COLORED_ULTRA_OTT$TITLE$RESET")

    silent || println("🚀 Starting combination processing...")

    progress = Progress(num_chunks, 1)  # Creiamo una barra di progresso

    Threads.@threads for chunk = 1:num_chunks
        local_results = Dict{Any,Vector{BigInt}}()
        local_label_count = Dict{Any,Int64}()
        local_contradictions = 0

        start_i = (chunk - 1) * chunk_size
        end_i = min(start_i + chunk_size - 1, num_combinations - 1)

        for i = start_i:end_i
            if isone(vertical)
                # Systematic
                combination, has_contradiction = generate_combination_ott(
                    BigInt(i),
                    num_atoms,
                    thresholds_by_feature,
                    atoms_by_feature,
                )
            else
                has_contradiction = true
                # Find the first random non-contradicting one
                while has_contradiction
                    i_rand = rand(BigInt(0):(all_combinations-1))
                    combination, has_contradiction = generate_combination_ott(
                        BigInt(i_rand),
                        num_atoms,
                        thresholds_by_feature,
                        atoms_by_feature,
                    )
                end
            end

            if has_contradiction
                local_contradictions += 1
            else
                combination_dict = SortedDict(combination)
                combination_vector = collect(values(combination_dict))
                result = apply_function(model, combination_vector)

                push!(get!(Vector{BigInt}, local_results, result), BigInt(i))
                local_label_count[result] = get(local_label_count, result, 0) + 1
            end
        end

        # Aggiorna i risultati globali in modo thread-safe
        lock(results_lock) do
            for (result, combs) in local_results
                append!(get!(Vector{BigInt}, results, result), combs)
            end
            for (result, count) in local_label_count
                label_count[result] = get(label_count, result, 0) + count
            end
            atomic_add!(contradictions, local_contradictions)
        end

        next!(progress)  # Aggiorniamo la barra di progresso
    end

    silent || println("Combination processing completed.")

    total_contradictions = contradictions[]
    silent || begin
        valid_combinations = num_combinations - total_contradictions
        println("\nTotal combinations: $num_combinations")
        println("Logical contradictions: $total_contradictions")
        println("\nLabel distribution:")
        for (label, count) in sort(collect(label_count), by = x -> x[2], rev = true)
            percentage = (count / valid_combinations) * 100
            println("$label: $count ($(round(percentage, digits=2))%)")
        end
    end

    silent || begin
        println("\nDetailed results:")
        for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
            println("[$result] ($(length(combinations)) combinations):")
            if length(combinations) > 10
                println("  First 10: $(combinations[1:10])")
                println("  ... ($(length(combinations)-10) more)")
            else
                println("  All: $combinations")
            end
        end
    end

    return results, label_count
end


"""
Generate a valid combination of feature values based on the provided thresholds and atom properties.

This function takes the integer `i` representing a combination, the number of atoms `num_atoms`, the thresholds for each feature `thresholds_by_feature`, and the atoms grouped by feature `atoms_by_feature`. It then generates a valid combination of feature values that satisfies all the atom conditions, or returns a dictionary with `false` if no valid combination can be found.

The function works by iterating through each feature, finding a valid value for that feature based on the associated atom conditions, and building up the final combination dictionary. If no valid value can be found for a feature, the function returns an empty dictionary and `true` to indicate a contradiction.
"""
function generate_combination_ott(
    i::BigInt,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    # Convertiamo direttamente in array i digit binari
    truth_values = digits(i, base = 2, pad = num_atoms)

    function find_valid_value(
        thresholds::Vector{Float64},
        atom_list::Vector{Tuple{Float64,Bool}},
        start_idx::Int,
    )::Union{Float64,Nothing}
        # Invece di mantenere una lista di valori validi, manteniamo i limiti del dominio
        min_domain = -Inf
        max_domain = Inf

        # Processiamo ogni condizione aggiornando i limiti del dominio
        for (idx, (threshold, is_less_than)) in enumerate(atom_list)
            current_truth = truth_values[start_idx+idx-1]

            if current_truth == (is_less_than ? 1 : 0)
                # Condizione del tipo x < threshold o x ≤ threshold
                max_domain = min(max_domain, threshold)
            else
                # Condizione del tipo x ≥ threshold o x > threshold
                min_domain = max(min_domain, threshold)
            end

            # Verifichiamo subito se il dominio è vuoto
            if min_domain >= max_domain
                return nothing
            end
        end

        # Troviamo il primo valore in thresholds che cade nel dominio valido
        for val in thresholds
            if min_domain <= val < max_domain
                return val
            end
        end

        return nothing
    end

    # Inizializziamo il dizionario per la combinazione risultante
    combination = zeros(maximum(keys(thresholds_by_feature)))
    current_idx = 1

    # Processiamo ogni feature
    for (feat, atom_list) in atoms_by_feature
        thresholds = thresholds_by_feature[feat]

        # Troviamo un valore valido per questa feature
        valid_value = find_valid_value(thresholds, atom_list, current_idx)

        # Se non troviamo un valore valido, abbiamo una contraddizione
        if isnothing(valid_value)
            return Dict(enumerate(combination)), true
        end

        # Salviamo il valore valido trovato
        combination[feat] = valid_value

        # Aggiorniamo l'indice per la prossima feature
        current_idx += length(atom_list)
    end

    return Dict(enumerate(combination)), false
end


"""
Group a list of atoms by their feature index. For each feature, the atoms are sorted by their threshold value and whether the threshold is a "less than" or "greater than" condition.

Args:
    atoms (Vector{Atom}): A list of atoms to group by feature.

Returns:
    Dict{Int,Vector{Tuple{Float64,Bool}}}: A dictionary mapping feature indices to a sorted list of (threshold, is_less_than) tuples for the atoms of that feature.
"""
function group_atoms_by_feature(atoms)
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        is_less_than = atom.value.metacond.test_operator === (<)
        push!(
            get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat),
            (threshold, is_less_than),
        )
    end
    for atom_list in values(atoms_by_feature)
        sort!(atom_list)
    end
    return atoms_by_feature
end


"""
Processes a combination of feature values, applies the model to the combination, and updates the results and label count accordingly.

Args:
    i (BigInt): The index of the current combination.
    num_atoms (Int): The total number of atoms.
    thresholds_by_feature (Dict{Int,Vector{Float64}}): A dictionary mapping feature indices to a list of thresholds for that feature.
    atoms_by_feature (Dict{Int,Vector{Tuple{Float64,Bool}}}): A dictionary mapping feature indices to a sorted list of (threshold, is_less_than) tuples for the atoms of that feature.
    model: The model to apply to the combination.
    results (Dict{Any,Vector{BigInt}}): A dictionary to store the results and their corresponding combinations.
    label_count (Dict{Any,Int64}): A dictionary to store the count of each label.
    contradictions (Atomic{Int64}): An atomic counter for the number of logical contradictions.
"""
function process_combination_ott(
    i::BigInt,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
    model,
    results,
    label_count,
    contradictions,
    apply_function = SoleModels.apply,
)
    combination =
        generate_combination(i, num_atoms, thresholds_by_feature, atoms_by_feature)

    if !isnothing(combination)
        features = collect(values(combination))
        result = apply_function(model, features)

        # Thread-safe update of results and label_count
        push!(get!(() -> Vector{BigInt}(), results, result), i)
        atomic_add!(get!(() -> Atomic{Int64}(0), label_count, result), 1)
    else
        atomic_add!(contradictions, 1)
    end
end


###############################################################################################################################

"""
Generates a statistical report for the execution of the optimization mode, including information about the rules, atoms, combinations, label distribution, performance, and complexity.

Args:
    nome_file (String): The name of the file to write the report to.
    all_rules (Any): The set of all rules used in the optimization.
    ntotatoms (Int): The total number of atoms before deduplication.
    nuniqatoms (Int): The number of unique atoms after deduplication.
    results (Dict{Any,Vector{BigInt}}): A dictionary containing the results and their corresponding combinations.
    label_count (Dict{Any,Int64}): A dictionary containing the count of each label.
    combined_results (Dict{Any,TwoLevelDNFFormula}): A dictionary containing the simplified formulas for each label.
    elapsed_time (Float64): The total execution time of the optimization.
    model (Any): The model used in the optimization.
"""
function genera_report_statistiche_ott(
    nome_file,
    all_rules,
    ntotatoms,
    nuniqatoms,
    results,
    label_count,
    combined_results,
    elapsed_time,
    model,
)
    open(nome_file, "w") do file
        println(file, "=== REPORT STATISTICO DELL'ESECUZIONE ===")
        println(file, "Data e ora: ", Dates.now())
        println(file, "====================================")

        # Statistiche sulle regole e gli atomi
        println(file, "\n1. REGOLE E ATOMI")
        println(file, "  Numero totale di regole: ", length(all_rules))
        println(file, "  Numero di proposizioni prima dell'uniq: ", ntotatoms)
        println(file, "  Numero di proposizioni dopo l'uniq: ", nuniqatoms)
        println(
            file,
            "  Riduzione delle proposizioni: ",
            round((1 - nuniqatoms / ntotatoms) * 100, digits = 2),
            "%",
        )
        println(file, "  Numero di atomi: ", nuniqatoms)

        # Statistiche sulle combinazioni
        println(file, "\n2. COMBINAZIONI")
        num_combinazioni_totali = sum(length(combs) for (_, combs) in results)
        println(file, "  Numero totale di combinazioni valide: ", num_combinazioni_totali)

        # Statistiche sulle etichette
        println(file, "\n3. DISTRIBUZIONE DELLE ETICHETTE")
        #sorted_label_count = sort(collect(label_count), by=x -> x[2][], rev=true)
        for (label, count) in collect(label_count)
            percentage = (count[] / num_combinazioni_totali) * 100
            println(file, "  $label: $count (", round(percentage, digits = 2), "%)")
        end

        # Statistiche sulla semplificazione
        #=
            println(file, "\n4. SEMPLIFICAZIONE DELLE FORMULE")
            for (result, formula) in combined_results
                formula_semplificata = minimizza_dnf(formula)
                riduzione = (1 - nterms(formula_semplificata) / nterms(formula)) * 100
                println(file, "  Etichetta $result:")
                println(file, "    Termini originali: ", nterms(formula))
                println(file, "    Termini dopo la semplificazione: ", nterms(formula_semplificata))
                println(file, "    Riduzione: ", round(riduzione, digits=2), "%")
            end
        =#

        # Prestazioni
        println(file, "\n5. PRESTAZIONI")
        println(
            file,
            "  Tempo totale di esecuzione: ",
            round(elapsed_time, digits = 2),
            " secondi",
        )
        println(
            file,
            "  Tempo medio per combinazione: ",
            round(elapsed_time / num_combinazioni_totali * 1000, digits = 2),
            " millisecondi",
        )

        # Complessità
        println(file, "\n6. COMPLESSITÀ")
        println(file, "  Numero di alberi nella foresta: ", length(model.trees))

        # Aggiunta della stampa delle formule semplificate
        println(file, "\n7. FORMULE SEMPLIFICATE")
        for (result, formula) in combined_results
            println(file, "  Etichetta $result:")
            formula_semplificata = minimizza_dnf(formula)
            println(file, "    Formula originale:")
            stampa_dnf(file, formula, 3)
            println(file, "    Formula semplificata:")
            stampa_dnf(file, formula_semplificata, 3)
            println(
                file,
                "    Riduzione: ",
                round(
                    (
                        1 -
                        nterms(formula_semplificata) /
                        nterms(formula)
                    ) * 100,
                    digits = 2,
                ),
                "%",
            )
            println(file)
        end

        println(file, "\n====================================")
        println(file, "Fine del report")
    end
    println("Report generato con successo: $nome_file")
end


"""
    compare_truth_combinations(model, alpha, atom_prop, vertical; kwargs...)

Compares the results of `truth_combinations_ott` and `truth_combinations` functions.

This function executes both `truth_combinations_ott` and `truth_combinations` functions, and then compares the results. It checks if the sets of combinations for each label are identical between the two results. If any differences are found, it prints the details of the differences.

Args:
    model: The model object.
    alpha: The alpha object.
    atom_prop: The atom properties.
    vertical: A boolean indicating whether to use vertical mode.
    kwargs: Additional keyword arguments passed to the `truth_combinations_ott` and `truth_combinations` functions.

Returns:
    A boolean indicating whether the results of the two functions are identical.
"""
function compare_truth_combinations(model, alpha, atom_prop, vertical; kwargs...)
    println("\nConfrontando i risultati di truth_combinations_ott e truth_combinations...")

    println("\n\n$COLORED_INFO$TITLE\n$RESET")
    # Esegui entrambe le funzioni
    @time results_ott, label_count_ott =
        truth_combinations_ott(model, alpha, atom_prop, vertical; kwargs...)
    println("$COLORED_INFO************** ⬆️ **************$RESET")
    @time results_standard, label_count_standard =
        truth_combinations(model, alpha, atom_prop, vertical; kwargs...)
    println("$COLORED_INFO************** ⬆️ **************$RESET")
    println("$COLORED_INFO$TITLE$RESET")

    # Prepara le strutture dati necessarie per creare dnf's
    num_atoms = length(atom_prop)
    thresholds_by_feature = Dict(
        subalpha.featcondition[1].feature.i_variable => sort(subalpha.featcondition[2])
        for subalpha in alpha.subalphabets
    )
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
    for atom in atom_prop
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # Crea dnf's per entrambi i risultati
    combined_results_ott =
        concat_results(results_ott, num_atoms, thresholds_by_feature, atoms_by_feature)
    combined_results_standard =
        concat_results(results_standard, num_atoms, thresholds_by_feature, atoms_by_feature)

    # Confronta i risultati
    are_equal = true
    for (label, formula_ott) in combined_results_ott
        if !haskey(combined_results_standard, label)
            println("Etichetta $label presente solo nei risultati ottimizzati")
            are_equal = false
            continue
        end

        formula_standard = combined_results_standard[label]
        if Set(eachcombination(formula_ott)) != Set(eachcombination(formula_standard))
            println("Differenze trovate per l'etichetta $label:")
            println("  Combinazioni ottimizzate: ", nterms(formula_ott))
            println("  Combinazioni standard: ", nterms(formula_standard))
            are_equal = false
        end
    end

    for label in keys(combined_results_standard)
        if !haskey(combined_results_ott, label)
            println("Etichetta $label presente solo nei risultati standard")
            are_equal = false
        end
    end

    if are_equal
        @info "I risultati di truth_combinations_ott e truth_combinations sono identici."
    else
        @warn "Sono state rilevate differenze tra i risultati di truth_combinations_ott e truth_combinations."
    end

    return are_equal
end
