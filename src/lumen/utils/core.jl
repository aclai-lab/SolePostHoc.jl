"""
Processes a combination of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

This function takes a combination of atom values, represented as either a `BitVector` or an integer, and the dictionaries of thresholds and atom properties. It then processes the combination, filtering the valid values for each feature based on the atom values. If a feature has no valid values, the function sets a flag indicating a contradiction. Finally, it returns the dictionary of valid values for each feature and the contradiction flag.

Parameters:
- `i::Union{BitVector,Int}`: The combination of atom values to be processed.
- `num_atoms::Int`: The number of atoms in the combination.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding valid values.
- `Bool`: A flag indicating if the combination has a contradiction.
"""
function truth_combinations(
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
    print_progress = true,
    silent = false,
)
    # Inizializzazione di strutture dati per gestire le soglie e gli atomi per feature
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # Popolamento delle strutture dati con le soglie e gli atomi
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end

    # Ordinamento degli atomi per ciascuna feature
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # Calcolo del numero totale di combinazioni
    num_atoms = length(atoms)
    all_combinations = BigInt(2)^num_atoms
    println("num_atoms: ", num_atoms)
    println("2^num_atoms: ", all_combinations)
    println("vertical: ", vertical)
    
    if isone(vertical)
        num_combinations = all_combinations
    else
        num_combinations = BigInt(round(all_combinations * vertical))
    end

    # Stima del tempo di esecuzione
    silent || begin
        println("Stima del tempo di esecuzione...")
        sample_size = min(1000, num_combinations)
        start_time = time()

        for i = 1:sample_size
            process_combination(
                BigInt(i - 1),
                num_atoms,
                thresholds_by_feature,
                atoms_by_feature,
            )
        end

        end_time = time()
        estimated_time_per_combination = (end_time - start_time) / sample_size
        total_estimated_time =
            BigFloat(estimated_time_per_combination) * num_combinations

        print_estimated_time(total_estimated_time)
    end

    # Inizio della generazione delle combinazioni
    silent || println("Inizio della generazione delle combinazioni...")
    silent || println("Totale combinazioni da generare: $num_combinations")

    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()

    # Processo di generazione e valutazione delle combinazioni
    start_time = time()
    contradictions = 0
    i_rand = 0
    update_interval = max(1, num_combinations ÷ 100)  # Aggiorna ogni 1% di progresso
    
    for i = 0:(num_combinations-1)
        i_v = i                                      # Indice per la progressbar
        if isone(vertical)
            combination, has_contradiction = process_combination(
                BigInt(i),
                num_atoms,
                thresholds_by_feature,
                atoms_by_feature,
            )
        else
            has_contradiction = true
                i = rand(BigInt(0):(all_combinations-1))
                combination, has_contradiction = process_combination(
                    BigInt(i),
                    num_atoms,
                    thresholds_by_feature,
                    atoms_by_feature,
                )
        end
        if has_contradiction
            contradictions += 1
        else
            combination_dict = SortedDict(combination)
            combination_vector = vcat(collect(values(combination_dict))...)
            
            if !(combination_vector in values(results))
                result = apply_function(model, combination_vector)
                push!(get!(Vector{BigInt}, results, result), BigInt(i))
                label_count[result] = get(label_count, result, 0) + 1
            end
        end

        # Aggiornamento della barra di avanzamento testuale FIXME TOTALMENTE RIMUOVIBILE ?
        print_progress && begin
            if i_v % update_interval == 0 || i_v == num_combinations - 1
                progress = BigFloat(i_v + 1) / num_combinations
                print_progress_bar(progress)
            end
        end
    end
    print_progress && println()  # Nuova linea dopo la barra di avanzamento

    end_time = time()
    real_time_per_combination = (end_time - start_time) / num_combinations

    # Stampa dei risultati
    silent || print_results_summary(
        real_time_per_combination,
        Int64(contradictions),
        num_combinations,
        label_count,
    )
    silent || print_detailed_results(results)

    return results, label_count
end


"""
Processes a combination of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

This function takes a combination of atom values, represented as either a `BitVector` or an integer, and the dictionaries of thresholds and atom properties. It then processes the combination, filtering the valid values for each feature based on the atom values. If a feature has no valid values, the function sets a flag indicating a contradiction. Finally, it returns the dictionary of valid values for each feature and the contradiction flag.

Parameters:
- `i::Union{BitVector,Int}`: The combination of atom values to be processed.
- `num_atoms::Int`: The number of atoms in the combination.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding valid values.
- `Bool`: A flag indicating if the combination has a contradiction.
"""
function process_combination(
    i,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    combination = Dict{Int,Vector{Float64}}()

    truth_row = !isa(i, BitVector) ? digits(i, base = 2, pad = num_atoms) : i

    has_contradiction = false
    j = 1
    for (feat, atom_list) in atoms_by_feature
        thresholds = thresholds_by_feature[feat]
        valid_values = copy(thresholds)

        for (_, (threshold, _)) in enumerate(atom_list)
            if truth_row[j] == 1
                filter!(x -> x < threshold, valid_values)
            else
                filter!(x -> x >= threshold, valid_values)
            end
            j += 1
        end

        if isempty(valid_values)
            has_contradiction = true
            break
        else
            combination[feat] = valid_values
        end

    end

    if !has_contradiction
        for feat in keys(thresholds_by_feature)
            if !haskey(combination, feat)
                combination[feat] = [1605.0]
            end
        end
    end

    return combination, has_contradiction
end


"""
Concatenates the results of a custom OR formula computation into a dictionary of `TwoLevelDNFFormula` objects.

This function takes the raw results of the custom OR formula computation, along with the number of atoms and the thresholds and atom properties, and constructs a dictionary of `TwoLevelDNFFormula` objects. The results are sorted by the number of combinations in descending order, and the `TwoLevelDNFFormula` objects are created with the sorted combinations.

Parameters:
- `results::Any`: The raw results of the custom OR formula computation.
- `num_atoms::Int`: The number of atoms in the custom OR formula.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Any,TwoLevelDNFFormula}`: A dictionary mapping the original results to their corresponding `TwoLevelDNFFormula` objects.
"""
function concat_results(
    results::Any,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    results = dict_to_tritvector(results, num_atoms)
    res = Dict{Any,TwoLevelDNFFormula}()
    println("\nRisultati dettagliati:")

    for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
        println("[$result] ($(length(combinations)) combinazioni)")
        res[result] = TwoLevelDNFFormula(
            Vector{TritVector}(combinations),
            num_atoms,
            thresholds_by_feature,
            atoms_by_feature
        )
    end
    return res
end

function concat_results(results::Any, my_atoms::Vector)
    num_atoms = length(my_atoms)
    results = dict_to_tritvector(results, num_atoms)
    res = Dict{Any,TwoLevelDNFFormula}()
    println("\nRisultati dettagliati:")

    for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
        println("[$result] ($(length(combinations)) combinazioni)")
        res[result] = TwoLevelDNFFormula(my_atoms, Vector{TritVector}(combinations))
    end
    return res
end

"""
    verify_simplification(original::TwoLevelDNFFormula, simplified::TwoLevelDNFFormula)

Verifies that a simplified custom OR formula is congruent with the original formula.

This function takes the original and simplified custom OR formulas, and generates a set of random assignments to evaluate and compare the results of the two formulas. If any mismatch is found, the function returns `false`, indicating that the simplification is not correct. Otherwise, it returns `true`, confirming that the simplified formula is congruent with the original.

Parameters:
- `original::TwoLevelDNFFormula`: The original custom OR formula.
- `simplified::TwoLevelDNFFormula`: The simplified custom OR formula.

Returns:
- `true` if the simplified formula is congruent with the original, `false` otherwise.
"""
function verify_simplification(original::TwoLevelDNFFormula, simplified::TwoLevelDNFFormula)
    @info "Starting verification of simplification"

    # Funzione per generare assegnazioni basate sui BitVector
    function generate_smart_assignments(formula, num_samples)
        assignments = Set{Dict{Int,Bool}}()

        # Converti ogni BitVector in un Dict{Int,Bool}
        for combination in eachcombination(formula)
            assignment = Dict{Int,Bool}()
            for (i, bit) in enumerate(combination)
                assignment[i] = bit
            end
            push!(assignments, assignment)
            if length(assignments) >= num_samples
                break
            end
        end

        # Aggiungi alcune assegnazioni casuali per aumentare la copertura
        while length(assignments) < num_samples
            random_assignment = Dict(i => rand(Bool) for i = 1:nuberofatoms(formula))
            push!(assignments, random_assignment)
        end
        return collect(assignments)
    end

    # Funzione ottimizzata per valutare la formula
    function evaluate_custom_or_formula(formula, assignment)
        for combination in eachcombination(formula)
            all_true = true
            for (feat, atom_list) in eachatomsbyfeature(formula)
                feat_value = get(assignment, feat, false)
                for (threshold, is_less_than) in atom_list
                    if is_less_than ? (feat_value >= threshold) : (feat_value < threshold)
                        all_true = false
                        break
                    end
                end
                if !all_true
                    break
                end
            end
            if all_true
                return true
            end
        end
        return false
    end

    # Genera un set di assegnazioni "intelligenti" basate sulle combinazioni esistenti
    num_samples = min(1000, 2^nuberofatoms(original))  # Limita il numero di campioni per input molto grandi
    assignments = generate_smart_assignments(original, num_samples)

    # Verifica le formule utilizzando le assegnazioni generate
    for (i, assignment) in enumerate(assignments)
        original_result = evaluate_custom_or_formula(original, assignment)
        simplified_result = evaluate_custom_or_formula(simplified, assignment)

        if original_result != simplified_result
            @warn "Mismatch found for assignment: $assignment"
            @warn "Original result: $original_result"
            @warn "Simplified result: $simplified_result"
            return false
        end

        # Stampa il progresso ogni 100 iterazioni
        if i % 100 == 0
            @info "Processed $i out of $num_samples assignments"
        end
    end

    @info "Verification complete. Simplified formula is congruent with the original."
    return true
end
