function truth_combinations_optimized(
    model::Any,
    alphabet,
    atoms::Vector{
        <:Atom{
            <:ScalarCondition{
                Float64,
                <:VariableNamedValue,
                <:ScalarMetaCondition{<:VariableNamedValue,typeof(<)},
            },
        },
    },
    vertical::Float64;
    apply_function = SoleModels.apply,
    print_progress = true,
    silent = true,
)
    # === DATA STRUCTURE INITIALIZATION ===
    # Pre-allocate dictionaries for better memory management
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # Populate feature thresholds (sorted once)
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort!(subalpha.featcondition[2])
    end

    # Group atoms by feature with pre-allocated vectors
    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        atom_list = get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat)
        push!(atom_list, (threshold, true))
    end

    # Sort atoms once per feature (stable sort for consistency)
    for atom_list in values(atoms_by_feature)
        sort!(atom_list, by = first)
    end

    # === BINARY ID MAPPING (same order as process_combination) ===
    atom_to_bit_position = Dict{Tuple{Int,Float64},Int}()
    bit_pos = 0
    for (feat, atom_list) in atoms_by_feature
        for (threshold, _) in atom_list
            atom_to_bit_position[(feat, threshold)] = bit_pos
            bit_pos += 1
        end
    end

    # === PRE-COMPUTE VALID COMBINATIONS PER FEATURE ===
    valid_combinations_by_feature = Dict{Int,Vector{Tuple{Vector{Int8},Vector{Float64}}}}()

    for (feat, atom_list) in atoms_by_feature
        thresholds = thresholds_by_feature[feat]
        n_atoms = length(atom_list)
        valid_combos = Vector{Tuple{Vector{Int8},Vector{Float64}}}()

        # Use Int8 for truth values (memory optimization)
        for combo_idx = 0:((1<<n_atoms)-1)  # Bit shift instead of power
            truth_vals = Vector{Int8}(undef, n_atoms)

            # Extract bits efficiently using bit operations
            for i = 1:n_atoms
                truth_vals[i] = (combo_idx >> (i-1)) & 1
            end

            # Apply constraints in-place for memory efficiency
            valid_values = copy(thresholds)
            constraint_satisfied = true

            for (i, (threshold, _)) in enumerate(atom_list)
                if truth_vals[i] == 1
                    filter!(x -> x < threshold, valid_values)
                else
                    filter!(x -> x >= threshold, valid_values)
                end

                # Early termination if no valid values remain
                if isempty(valid_values)
                    constraint_satisfied = false
                    break
                end
            end

            # Only store valid combinations
            if constraint_satisfied
                push!(valid_combos, (truth_vals, valid_values))
            end
        end

        valid_combinations_by_feature[feat] = valid_combos
    end

    # Add features without atoms (using default value)
    for feat in keys(thresholds_by_feature)
        if !haskey(valid_combinations_by_feature, feat)
            valid_combinations_by_feature[feat] = [(Int8[], [1605.0])]
        end
    end

    # === EFFICIENT BINARY ID RECONSTRUCTION ===
    @inline function reconstruct_binary_id(
        combination_tuple,
        feature_order,
        atoms_by_feature,
        atom_to_bit_position,
    )
        binary_id = zero(BigInt)

        for (feat_idx, feat) in enumerate(feature_order)
            if haskey(atoms_by_feature, feat)
                truth_values = combination_tuple[feat_idx][1]
                atom_list = atoms_by_feature[feat]

                for (atom_idx, (threshold, _)) in enumerate(atom_list)
                    if truth_values[atom_idx] == 1
                        bit_pos = atom_to_bit_position[(feat, threshold)]
                        binary_id |= (BigInt(1) << bit_pos)  # Bit operations for speed
                    end
                end
            end
        end

        return binary_id
    end

    # === COMBINATION GENERATION ===
    # Maintain same order as process_combination for ID consistency
    feature_order = collect(keys(atoms_by_feature))

    # Add features without atoms
    for feat in keys(thresholds_by_feature)
        if !(feat in feature_order)
            push!(feature_order, feat)
        end
    end

    combination_sets = [valid_combinations_by_feature[feat] for feat in feature_order]

    # Apply vertical sampling efficiently
    if isone(vertical)
        # Process all valid combinations
        combinations_iter = Iterators.product(combination_sets...)
    else
        # Sample random combinations without replacement
        all_combinations = collect(Iterators.product(combination_sets...))
        n_total = length(all_combinations)
        n_sample = min(n_total, Int(round(n_total * vertical)))

        # Use reservoir sampling for large datasets
        if n_sample < n_total ÷ 10
            sampled_indices = Set{Int}()
            while length(sampled_indices) < n_sample
                push!(sampled_indices, rand(1:n_total))
            end
            combinations_iter = (all_combinations[i] for i in sampled_indices)
        else
            sample_indices = randperm(n_total)[1:n_sample]
            combinations_iter = (all_combinations[i] for i in sample_indices)
        end
    end

    # === RESULT PROCESSING ===
    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()

    # Use Set for O(1) duplicate detection with vector hash
    seen_combinations = Set{UInt64}()

    for combination_tuple in combinations_iter
        # Reconstruct original binary ID
        original_id = reconstruct_binary_id(
            combination_tuple,
            feature_order,
            atoms_by_feature,
            atom_to_bit_position,
        )

        # Build combination vector efficiently
        combination_vector = Vector{Float64}()
        for (_, valid_values) in combination_tuple
            append!(combination_vector, valid_values)
        end

        # Fast duplicate detection using hash
        combo_hash = hash(combination_vector)
        if combo_hash ∉ seen_combinations
            push!(seen_combinations, combo_hash)

            # Apply model prediction
            result = if model isa AbstractModel
                # Pre-allocate DataFrame for better performance
                df = DataFrame(reshape(combination_vector, 1, :), :auto)
                apply_function(model, df)
            else
                apply_function(model, combination_vector)
            end

            # Store result with original binary ID
            result_vec = get!(Vector{BigInt}, results, result)
            push!(result_vec, original_id)
            label_count[result] = get(label_count, result, 0) + 1
        end
    end

    # Out
    return results, label_count
end

"""
[WORK BUT IS UPGRADABLE VERSION] Processes a combination ottimize mode of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

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
function truth_combinations_ott01(
    model::Any,
    alphabet,
    atoms::Vector{
        <:Atom{
            <:ScalarCondition{
                Float64,
                <:VariableNamedValue,
                <:ScalarMetaCondition{<:VariableNamedValue,typeof(<)},
            },
        },
    },
    vertical::Float64;
    apply_function = SoleModels.apply,
    print_progress = true,
    silent = true,
)
    # Initialization of data structures
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # Populating data structures with thresholds and atoms
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end

    # Sorting atoms for each feature
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # Generate all possible truth value combinations
    num_atoms = length(atoms)
    all_combinations = BigInt(2)^num_atoms

    num_combinations = if isone(vertical)
        all_combinations
    else
        BigInt(round(all_combinations * vertical))
    end

    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()

    # Pre-calculate all valid combinations to avoid redundant process_combination calls
    valid_combinations_cache = Dict{BigInt,Tuple{Dict{Int,Vector{Float64}},Bool}}()

    # Generate only the combinations we need
    combination_indices = if isone(vertical)
        0:(all_combinations-1)
    else
        # Sample random indices
        sort(rand(BigInt(0):(all_combinations-1), Int(num_combinations)))
    end

    for i in combination_indices
        # Check if we already computed this combination
        if haskey(valid_combinations_cache, i)
            combination, has_contradiction = valid_combinations_cache[i]
        else
            combination, has_contradiction =
                process_combination(i, num_atoms, thresholds_by_feature, atoms_by_feature)
            valid_combinations_cache[i] = (combination, has_contradiction)
        end

        if !has_contradiction
            combination_dict = SortedDict(combination)
            combination_vector = vcat(collect(values(combination_dict))...)

            if !(combination_vector in values(results))
                result = if model isa AbstractModel
                    apply_function(model, DataFrame(reshape(combination_vector, 1, :), :auto))
                else
                    apply_function(model, combination_vector)
                end
                push!(get!(Vector{BigInt}, results, result), i)
                label_count[result] = get(label_count, result, 0) + 1
            end
        end
    end

    return results, label_count
end


"""
[ALFA VERSION] Processes a combination of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

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
function truth_combinations_ott00(
    model::Any,
    alphabet,
    atoms::Vector{
        <:Atom{
            <:ScalarCondition{
                Float64,
                <:VariableNamedValue,
                <:ScalarMetaCondition{<:VariableNamedValue,typeof(<)},
            },
        },
    },
    vertical::Float64;
    apply_function = SoleModels.apply,
    print_progress = true,
    silent = false,
)
    # Initialization of data structures to manage thresholds and atoms for each feature
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # Populating data structures with thresholds and atoms
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end

    # Sorting atoms for each feature
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # Calculation of the total number of combinations
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

    # Estimation of execution time
    silent || begin
        println("Estimating execution time...")
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

    # Start of combination generation
    silent || println("Starting combination generation...")
    silent || println("Total combinations to generate: $num_combinations")

    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()

    # Process of generating and evaluating combinations
    start_time = time()
    contradictions = 0
    i_rand = 0
    update_interval = max(1, num_combinations ÷ 100)  # Update every 1% of progress

    println(atoms_by_feature)
    print(
        "\n----------------------------------------------------------------------------------\n",
    )
    # Trova tutte le chiavi presenti in thresholds_by_feature per determinare le feature esistenti
    all_features = Set(keys(thresholds_by_feature))

    # Pre-calcola il numero totale di combinazioni (include epsilon per feature con atoms, 1 per feature senza atoms)
    total_combinations = 1
    for k in [4, 3, 2, 1]
        if k in all_features
            if haskey(atoms_by_feature, k)
                total_combinations *= (length(atoms_by_feature[k]) + 1)  # +1 per epsilon
            else
                total_combinations *= 1  # Solo valore simbolico
            end
        end
    end

    println("Numero totale di combinazioni: $total_combinations")
    println("Inizio generazione combinazioni...\n")

    # Ordina le chiavi nell'ordine desiderato: 4, 3, 2, 1
    keys_sorted = [1, 2, 3, 4]

    # Prepara i valori per ogni feature (con atoms o simbolico)
    symbolic_value = 1605  # Valore simbolico per feature senza atoms
    epsilon = 1e-6  # Valore di epsilon

    values_sorted = []
    existing_keys = []

    for k in keys_sorted
        if k in all_features  # Feature esiste in thresholds_by_feature
            push!(existing_keys, k)

            if haskey(atoms_by_feature, k)
                # Feature ha atoms: estrai Float64 e aggiungi epsilon
                float_values = [tuple[1] for tuple in atoms_by_feature[k]]

                if !isempty(float_values)
                    min_value = minimum(float_values)
                    epsilon_value = min_value - epsilon
                    # Inserisci epsilon_value all'inizio della lista
                    float_values = [epsilon_value; float_values]
                end

                push!(values_sorted, float_values)
            else
                # Feature non ha atoms: usa solo valore simbolico
                push!(values_sorted, [symbolic_value])
            end
        end
    end

    # Calcola la lunghezza totale del vettore binario
    total_length = 0
    for k in existing_keys
        if haskey(atoms_by_feature, k)
            total_length += length(atoms_by_feature[k]) + 1  # +1 per epsilon
        else
            total_length += length(thresholds_by_feature[k])  # Usa lunghezza da thresholds
        end
    end

    # Loop ottimizzato con generazione vettore binario
    let count = 0
        for combination in Iterators.product(values_sorted...)
            count += 1
            combination_dict = Dict(zip(existing_keys, combination))

            # Crea il vettore binario
            binary_vector = zeros(Int, total_length)
            current_pos = 1
            for (key, selected_value) in zip(existing_keys, combination)
                if haskey(atoms_by_feature, key)
                    # Feature con atoms: gestisci normalmente
                    original_float_values = [tuple[1] for tuple in atoms_by_feature[key]]
                    min_value = minimum(original_float_values)
                    epsilon_value = min_value - epsilon
                    complete_values = [epsilon_value; original_float_values]

                    # Trova l'indice del valore selezionato
                    selected_index = findfirst(x -> x == selected_value, complete_values)

                    # Imposta il pattern: 0 fino all'indice selezionato, poi tutti 1
                    for i = selected_index:length(complete_values)
                        binary_vector[current_pos+i-1] = 1
                    end

                    current_pos += length(complete_values)
                else
                    # Feature senza atoms: usa tutti 1 (assume soglia minima sempre superata)
                    threshold_length = length(thresholds_by_feature[key])
                    for i = 1:threshold_length
                        binary_vector[current_pos+i-1] = 1
                    end
                    current_pos += threshold_length
                end
            end

            binary_string = join(binary_vector)
            println("[$count/$total_combinations] $combination_dict")
            println("Vettore binario: $binary_string")
            # da vettore a intero
            integer_representation = parse(BigInt, binary_string, base = 2)
            println("Rappresentazione intera: $integer_representation")

            println("---")
            combination_dict = SortedDict(combination_dict)
            combination_vector = vcat(collect(values(combination_dict))...)
            if !(combination_vector in values(results))
                result = if model isa AbstractModel
                    apply_function(model, DataFrame(reshape(combination_vector, 1, :), :auto))
                else
                    apply_function(model, combination_vector)
                end
                println("resoult", result)
            end
            println("---")

        end
        println("\nCompletato: $count combinazioni generate.")
    end
    print(
        "\n----------------------------------------------------------------------------------\n",
    )

    for i = 0:(num_combinations-1)
        i_v = i  # Index for the progress bar
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

            println(
                "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
            )
            println(combination)
            println(
                "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
            )

            combination_dict = SortedDict(combination)
            combination_vector = vcat(collect(values(combination_dict))...)
            if !(combination_vector in values(results))
                result = if model isa AbstractModel
                    apply_function(model, DataFrame(reshape(combination_vector, 1, :), :auto))
                else
                    apply_function(model, combination_vector)
                end
                push!(get!(Vector{BigInt}, results, result), BigInt(i))
                label_count[result] = get(label_count, result, 0) + 1
            end
        end

        # Update of textual progress bar FIXME CAN BE COMPLETELY REMOVED?
        print_progress && begin
            if i_v % update_interval == 0 || i_v == num_combinations - 1
                progress = BigFloat(i_v + 1) / num_combinations
                print_progress_bar(progress)
            end
        end
    end
    print_progress && println()  # New line after progress bar

    end_time = time()
    real_time_per_combination = (end_time - start_time) / num_combinations

    # Print results
    silent || print_results_summary(
        real_time_per_combination,
        Int64(contradictions),
        num_combinations,
        label_count,
    )
    silent || print_detailed_results(results)
    print(
        "\nasssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n",
    )
    print(results)
    print(
        "\nasssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n",
    )
    return results, label_count
end
