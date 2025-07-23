using Dates

function debug_combinations(
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
    alphabetcontroll = "default",
)
    # Extract the threshold for each feature
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    # Count the atoms of each feature
    atoms_by_feature = Dict{Int,Int}()
    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        atoms_by_feature[feat] = get(atoms_by_feature, feat, 0) + 1
    end
    
    num_atoms = length(atoms)
    total_possible_combinations = BigInt(2)^num_atoms
    valid_combinations_count = BigInt(1)
    
    # For each feature, multiply the number of thresholds
    println(thresholds_by_feature)
    for (feat, thresholds) in thresholds_by_feature
        if  haskey(atoms_by_feature, feat)
            valid_combinations_count *= length(thresholds)
        end
    end
    
    # Additional statistics
    num_features = length(thresholds_by_feature)
    avg_thresholds_per_feature = num_features > 0 ? sum(length(thresholds) for (_, thresholds) in thresholds_by_feature) / num_features : 0
    max_thresholds = num_features > 0 ? maximum(length(thresholds) for (_, thresholds) in thresholds_by_feature) : 0
    min_thresholds = num_features > 0 ? minimum(length(thresholds) for (_, thresholds) in thresholds_by_feature) : 0
    
    # Use the timestamp as identifier
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    identifier = timestamp
    
    results = []
    push!(results, "="^60)
    push!(results, "DEBUG: COMBINATION ANALYSIS - $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    push!(results, "="^60)
    push!(results, "Number of atoms: $num_atoms")
    push!(results, "Theoretically possible combinations (2^$num_atoms): $total_possible_combinations")
    push!(results, "")
    push!(results, "FEATURE ANALYSIS:")
    push!(results, "-"^40)
    
    for (feat, thresholds) in sort(collect(thresholds_by_feature))
        num_atoms_for_feat = get(atoms_by_feature, feat, 0)
        num_thresholds = length(thresholds)
        
        push!(results, "Feature $feat:")
        push!(results, "  Available thresholds: $thresholds")
        push!(results, "  Number of atoms: $num_atoms_for_feat")
        push!(results, "  Valid combination for feature: $num_thresholds")
        push!(results, "")
    end
    
    push!(results, "FINAL RESULT:")
    push!(results, "-"^40) 
    push!(results, "Theoretically possible combinations: $total_possible_combinations")
    push!(results, "Valid combinations: $valid_combinations_count")
    if total_possible_combinations > 0
        ratio = Float64(valid_combinations_count) / Float64(total_possible_combinations) * 100
        push!(results, "valid/possible ratio: $(round(ratio, digits=2))%")
    end
    push!(results, "Removed combinations: $(total_possible_combinations - valid_combinations_count)")
    
    # Calculation details
    threshold_counts = [string(length(thresholds)) for (feat, thresholds) in sort(collect(thresholds_by_feature)) if haskey(atoms_by_feature, feat)]
    calculation = join(threshold_counts, " × ")
    push!(results, "Calculation: $calculation = $valid_combinations_count")
    push!(results, "="^60)
    
    # Print to console if not silent
    if !silent
        for line in results
            println(line)
        end
    end
    
    # Salva su file di testo con timestamp e alphabetcontroll
    # filename_txt = "debug_combinations_$(timestamp)_$(alphabetcontroll).txt"
    # open(filename_txt, "w") do file
    #     for line in results
    #         println(file, line)
    #     end
    # end
    
    # === CSV ===
    csv_filename = "debug_combinations_stats.csv"
    
    # Calculate ratio as numeric value
    validity_ratio = total_possible_combinations > 0 ? Float64(valid_combinations_count) / Float64(total_possible_combinations) : 0.0
    
    # Check if csv file exists
    file_exists = isfile(csv_filename)
    
    # Header
    header = "Identifier,NumAtoms,NumFeatures,TotalCombinations,ValidCombinations,DiscardedCombinations,ValidityRatio,AvgThresholds,MaxThresholds,MinThresholds,Vertical,AlphabetControl"
    
    # Format data
    row_data = [
        identifier,
        num_atoms,
        num_features,
        string(total_possible_combinations),
        string(valid_combinations_count),
        string(total_possible_combinations - valid_combinations_count),
        round(validity_ratio, digits=6),
        round(avg_thresholds_per_feature, digits=3),
        max_thresholds,
        min_thresholds,
        vertical,
        alphabetcontroll
    ]
    
    # Write to CSV
    open(csv_filename, file_exists ? "a" : "w") do file
        # If file doesn't exist print the header
        if !file_exists
            println(file, header)
        end
        
        # Write data
        println(file, join(row_data, ","))
    end
    
    #println("Statistiche salvate su: $filename_txt")
    println("CSV data saved in: $csv_filename")
    
    return row_data
end


function debug_combinations_old(
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
    alphabetcontroll = "default",
)
    # Estrai soglie per feature
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    # Conta atomi per feature
    atoms_by_feature = Dict{Int,Int}()
    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        atoms_by_feature[feat] = get(atoms_by_feature, feat, 0) + 1
    end

    # Calcolo semplice
    num_atoms = length(atoms)
    total_possible_combinations = BigInt(2)^num_atoms
    valid_combinations_count = BigInt(1)
    
    # Per ogni feature, moltiplica il numero di soglie
    for (feat, thresholds) in thresholds_by_feature
        valid_combinations_count *= length(thresholds)
    end
    
    results = []
    push!(results, "="^60)
    push!(results, "DEBUG: ANALISI COMBINAZIONI - $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    push!(results, "="^60)
    push!(results, "Numero di atomi: $num_atoms")
    push!(results, "Combinazioni teoricamente possibili (2^$num_atoms): $total_possible_combinations")
    push!(results, "")
    push!(results, "ANALISI PER FEATURE:")
    push!(results, "-"^40)
    
    for (feat, thresholds) in sort(collect(thresholds_by_feature))
        num_atoms_for_feat = get(atoms_by_feature, feat, 0)
        num_thresholds = length(thresholds)
        
        push!(results, "Feature $feat:")
        push!(results, "  Soglie disponibili: $thresholds")
        push!(results, "  Numero atomi: $num_atoms_for_feat")
        push!(results, "  Combinazioni valide per feature: $num_thresholds")
        push!(results, "")
    end
    
    push!(results, "RISULTATO FINALE:")
    push!(results, "-"^40)
    push!(results, "Combinazioni teoricamente possibili: $total_possible_combinations")
    push!(results, "Combinazioni effettivamente valide: $valid_combinations_count")
    if total_possible_combinations > 0
        ratio = Float64(valid_combinations_count) / Float64(total_possible_combinations) * 100
        push!(results, "Rapporto valide/possibili: $(round(ratio, digits=2))%")
    end
    push!(results, "Combinazioni scartate: $(total_possible_combinations - valid_combinations_count)")
    
    # Calcolo dettagliato
    threshold_counts = [string(length(thresholds)) for (_, thresholds) in sort(collect(thresholds_by_feature))]
    calculation = join(threshold_counts, " × ")
    push!(results, "Calcolo: $calculation = $valid_combinations_count")
    push!(results, "="^60)
    
    # Stampa su console se non silent
    if !silent
        for line in results
            println(line)
        end
    end
    
    # Salva su file con timestamp e alphabetcontroll
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "debug_combinations_$(timestamp)_$(alphabetcontroll).txt"
    
    open(filename, "w") do file
        for line in results
            println(file, line)
        end
    end
    
    println("Statistiche salvate su: $filename")
    
end