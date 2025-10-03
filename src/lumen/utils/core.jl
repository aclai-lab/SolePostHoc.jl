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

    return results, label_count
end
