#= -------------------------------------------------------------------------------------------------------------------------
############################################################################################################################
#                                               OTTIMIZZAZIONE-MOD                                                                                                                                               #
############################################################################################################################
=#
"""
[OTT VERSION] Processes a combination ottimize mode of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

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
function truth_combinations_ott(
    model::Any,
    alphabet,
    # atoms::Vector{
    #     <:Atom{
    #         <:ScalarCondition{
    #             <:AbstractFloat,
    #             <:VariableNamedValue,
    #             <:ScalarMetaCondition{<:VariableNamedValue,typeof(<)},
    #         },
    #     },
    # },
    atoms::Vector{<:Atom{<:ScalarCondition}},
    vertical::Float64;
    apply_function = SoleModels.apply,
    print_progress = true,
    silent = true,
)
    # Same initialization as the original
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # For each feature we extract the atoms and thresholds
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    # For each atom we get it's feature index and threashold
    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end

    # We then sort each atom list (Legacy code)
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # We create an atom mapping with the same order as `process_combination`
    atom_to_bit_position = Dict{Tuple{Int,Float64}, Int}()
    bit_position = 0
    
    # Crucial: we keep the same order of iteration as `process_combination`
    for (feat, atom_list) in atoms_by_feature
        for (threshold, _) in atom_list
            atom_to_bit_position[(feat, threshold)] = bit_position
            bit_position += 1
        end
    end

    # Pre-process all the valid combination for each feature (In the same order of legacy version)
    valid_combinations_by_feature = Dict{Int, Vector{Tuple{Vector{Int}, Vector{Float64}}}}()
    
    # For each feaure we find the valid combinations
    for (feat, atom_list) in atoms_by_feature
        # Extract the thresholds
        thresholds = thresholds_by_feature[feat]
        valid_combinations_by_feature[feat] = []
        
        # Extract the atoms
        atom_thresholds = [atom[1] for atom in atom_list]
        num_atoms_for_feat = length(atom_thresholds)
        
        # Generate all truth value combinations for the feature
        for binary_combo in 0:(2^num_atoms_for_feat - 1)
            # Parse the value into binary format
            truth_values = digits(binary_combo, base=2, pad=num_atoms_for_feat)
            valid_values = copy(thresholds)
            
            # Apply the atom threshold to the value
            for (atom_idx, (threshold, _)) in enumerate(atom_list)
                truth_val = truth_values[atom_idx]
                if truth_val == 1
                    filter!(x -> x < threshold, valid_values)
                else
                    filter!(x -> x >= threshold, valid_values)
                end
            end
            
            # If any valid value is found we add it to the valid combinations vector
            if !isempty(valid_values)
                push!(valid_combinations_by_feature[feat], (truth_values, valid_values))
            end
        end
    end
    
    # Add features without any atoms
    for feat in keys(thresholds_by_feature)
        if !haskey(valid_combinations_by_feature, feat)
            push!(get!(Vector{Tuple{Vector{Int}, Vector{Float64}}}, valid_combinations_by_feature, feat), (Int[], [1605.0]))
        end
    end
    
    # This function is used to parse the combination to the relative value in base 10
    function reconstruct_binary_id(combination_tuple, feature_keys, atoms_by_feature, atom_to_bit_position)
        binary_id = BigInt(0)
        
        # As before we iterate with the same order as  `process_combination`
        for (feat, atom_list) in atoms_by_feature
            if feat in feature_keys
                feat_idx = findfirst(x -> x == feat, feature_keys)
                truth_values, _ = combination_tuple[feat_idx]
                
                # Extract the atom postition based on the value of the combination instance
                for (atom_idx, (threshold, _)) in enumerate(atom_list)
                    bit_pos = atom_to_bit_position[(feat, threshold)]
                    truth_val = truth_values[atom_idx]
                    
                    if truth_val == 1
                        binary_id += BigInt(2)^bit_pos
                    end
                end
            end
        end
        
        # Return the parsed base 10 combination value
        return binary_id
    end
    
    # Extract the feature keys
    feature_keys = [feat for (feat, _) in atoms_by_feature]  
    
    # Add features without any atoms
    for feat in keys(thresholds_by_feature)
        if !(feat in feature_keys)
            push!(feature_keys, feat)
        end
    end

    # TODO: check if this is just a name change
    combination_sets = [valid_combinations_by_feature[feat] for feat in feature_keys]
    all_combinations_iter = Iterators.product(combination_sets...)
    
    # Apply vertical sampling 
    combinations_to_process = if isone(vertical)
        all_combinations_iter
    else
        combinations_vec = collect(all_combinations_iter)
        sample_size = min(length(combinations_vec), Int(round(length(combinations_vec) * vertical)))
        sample_indices = randperm(length(combinations_vec))[1:sample_size]
        [combinations_vec[i] for i in sample_indices]
    end
    
    # Initialize results
    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()
    seen_vectors = Set{Vector{Float64}}()
    
    # For each generated combinations
    for combination_tuple in combinations_to_process
        # We reconstruct the related base 10 value
        original_binary_id = reconstruct_binary_id(combination_tuple, feature_keys, atoms_by_feature, atom_to_bit_position)
        
        # Generate the combination dict
        combination_dict = Dict{Int, Vector{Float64}}()
        for (feat_idx, feat) in enumerate(feature_keys)
            _, valid_values = combination_tuple[feat_idx]
            combination_dict[feat] = valid_values
        end
        
        # Process the combination
        combination_dict = SortedDict(combination_dict)
        combination_vector = vcat(collect(values(combination_dict))...)
        
        # Skip any duplicates
        if !(combination_vector in seen_vectors)
            push!(seen_vectors, combination_vector)
            
            # Make a prediction using the passed model
            result = if model isa AbstractModel
                apply_function(model, DataFrame(reshape(combination_vector, 1, :), :auto))
            else
                apply_function(model, combination_vector)
            end
            
            push!(get!(Vector{BigInt}, results, result), original_binary_id)
            label_count[result] = get(label_count, result, 0) + 1
        end
    end

    # Return the results
    return results, label_count
end

function testOttt(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function, testott)
    # open output file 
    open("test_ott_$testott.txt", "w") do file
        
        println(file, "🚀 Benchmark Truth Combinations vs Truth Combinations OTT")
        println(file, "=" ^ 60)

        # Test with @elapsed (Single execution)
        println(file, "\n📊 Single test with @execution:")
        println(file, "-" ^ 40)

        print(file, "truth_combinations: ")
        time1 = @elapsed result1 = truth_combinations(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
        println(file, "$(round(time1*1000, digits=3)) ms")

        print(file, "truth_combinations_ott: ")
        time2 = @elapsed result2 = truth_combinations_ott(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
        println(file, "$(round(time2*1000, digits=3)) ms")

        # Check if the results are the same 
        println(file, "\n🔍 Consistency check:")
        println(file, "-" ^ 30)
        if result1 == result2
            println(file, "✅ The results are consistent !")
        else
            println(file, "❌ WARNING: the results are not consistent !")
            println(file, "Result1 type: $(typeof(result1))")
            println(file, "Result2 type: $(typeof(result2))")

            # Check details
            if isa(result1, Array) && isa(result2, Array)
                println(file, "Result1 size: $(size(result1))")
                println(file, "Result2 size: $(size(result2))")
                if size(result1) == size(result2)
                    diff_count = sum(result1 .!= result2)
                    println(file, "Different elements: $diff_count / $(length(result1))")
                    if diff_count > 0 && diff_count <= 10
                        println(file, "First 10 differences:")
                        for i in 1:min(length(result1), 10)
                            if result1[i] != result2[i]
                                println(file, "  index $i: $(result1[i]) vs $(result2[i])")
                            end
                        end
                    end
                end
            end
        end
        
        # Multiple test needed to make the mean more reliable
        println(file, "\n🔄 Multiple test (20 iterations + warm-up):")
        println(file, "-" ^ 40)

        n_tests = 20
        n_warmup = 3

        # Warm-up to stabilize the compilation JIT
        println(file, "Warm-up...")
        for i in 1:n_warmup
            truth_combinations(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
            truth_combinations_ott(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
        end

        # Test truth_combinations
        println(file, "Testing truth_combinations...")
        times1 = Float64[]
        for i in 1:n_tests
            # Force garbage collection
            GC.gc()
            t = @elapsed truth_combinations(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
            push!(times1, t)
        end

        # Test truth_combinations_ott  
        println(file, "Testing truth_combinations_ott...")
        times2 = Float64[]
        for i in 1:n_tests
            # Force garbage collection
            GC.gc()
            t = @elapsed truth_combinations_ott(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
            push!(times2, t)
        end

        # Statistics (with outlier removal)
        # Remove the more extreme outlier (top 10% & bottom 10%)
        times1_sorted = sort(times1)
        times2_sorted = sort(times2)
        
        # Take the middle 80% (Remove the 10% at the extremes)
        start_idx = max(1, Int(round(n_tests * 0.1)))
        end_idx = min(n_tests, Int(round(n_tests * 0.9)))
        
        times1_clean = times1_sorted[start_idx:end_idx]
        times2_clean = times2_sorted[start_idx:end_idx]
        
        avg1 = sum(times1_clean) / length(times1_clean)
        avg2 = sum(times2_clean) / length(times2_clean)
        min1 = minimum(times1)
        min2 = minimum(times2)
        max1 = maximum(times1)
        max2 = maximum(times2)
        median1 = times1_sorted[div(n_tests, 2)]
        median2 = times2_sorted[div(n_tests, 2)]

        # Risultati
        println(file, "\n📈 RESULTS:")
        println(file, "=" ^ 50)
        println(file, "truth_combinations:")
        println(file, "  Mean time (no outlier): $(round(avg1*1000, digits=3)) ms")
        println(file, "  Median time:            $(round(median1*1000, digits=3)) ms")
        println(file, "  Min time:                $(round(min1*1000, digits=3)) ms") 
        println(file, "  Max time:                $(round(max1*1000, digits=3)) ms")

        println(file, "\ntruth_combinations_ott:")
        println(file, "  Mean time (no outlier): $(round(avg2*1000, digits=3)) ms")
        println(file, "  Median time:             $(round(median2*1000, digits=3)) ms")
        println(file, "  Min time:                $(round(min2*1000, digits=3)) ms")
        println(file, "  Max time:                $(round(max2*1000, digits=3)) ms")

        # Confronto
        speedup = avg1 / avg2
        speedup_median = median1 / median2
        if speedup > 1.0
            println(file, "\n🏆 truth_combinations_ott is $(round(speedup, digits=2))x faster (mean)!")
            println(file, "🏆 truth_combinations_ott is $(round(speedup_median, digits=2))x faster (median)!")
        else
            println(file, "\n⚠️  truth_combinations is $(round(1/speedup, digits=2))x faster (mean)!")
            println(file, "⚠️  truth_combinations is $(round(1/speedup_median, digits=2))x faster (median)!")
        end

        println(file, "\n📊 All times (ms):")
        println(file, "truth_combinations: ", [round(t*1000, digits=2) for t in times1])
        println(file, "truth_combinations_ott: ", [round(t*1000, digits=2) for t in times2])
        
        # Identifica outlier
        if max2 > 3 * median2
            println(file, "\n⚠️  OUTLIER FOUND in truth_combinations_ott:") 
            println(file, "   Max time $(round(max2*1000, digits=2)) ms is a lot higher than the median $(round(median2*1000, digits=2)) ms")
            println(file, "   Possible causes: GC, compilation JIT, system interference")
        end
    end
    
    println("✅ Benchmark completed! Result saved in 'test_ott.txt'")
end