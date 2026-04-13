using SoleModels
using SolePostHoc
using DecisionTree
using DataFrames
using Random


function evaluate_forest(forest::DecisionEnsemble, X_matrix::Matrix{Float64}, y_true::Vector{String})
    df = DataFrame(X_matrix, :auto)
    
    predictions = string.(SoleModels.apply(forest, df))
    accuracy = sum(predictions .== y_true) / length(y_true)
    return round(accuracy * 100, digits=2)
end

function run_usecases()
    Random.seed!(42)

    n_samples = 500
    n_features = 10
    X_matrix = rand(n_samples, n_features)
    
    classlabels = ["First", "Second"]
    yc_int = rand(1:length(classlabels), n_samples)
    yc_string = [classlabels[i] for i in yc_int]

    train_idx = 1:400
    test_idx = 401:500
    
    f_val = X_matrix[test_idx, :]
    l_val = yc_string[test_idx]


    println("Training original Random Forest...")
    n_trees_original = 50 
    n_subfeatures = 3
    
    raw_forest = DecisionTree.build_forest(yc_int[train_idx], X_matrix[train_idx, :], n_subfeatures, n_trees_original)
    
    featurenames = ["Feature_$i" for i in 1:n_features]
    original_forest = solemodel(raw_forest; featurenames=featurenames, classlabels=classlabels)

    original_acc = evaluate_forest(original_forest, f_val, l_val)
    original_size = length(original_forest.models)

    output_filename = joinpath(@__DIR__, "orca_use_cases.txt")
    open(output_filename, "w") do io
        write(io, "=================================================\n")
        write(io, "      ORCA ALGORITHM - USE CASE COMPARISON       \n")
        write(io, "=================================================\n\n")
        
        write(io, "Dataset: 500 samples, 10 features\n")
        write(io, "Original forest:\n")
        write(io, " - Number of trees: $original_size\n")
        write(io, " - Baseline accuracy: $original_acc%\n\n")
        write(io, "-------------------------------------------------\n")

        
        modes_to_test = [:size, 
                         :depth, 
                         :alphabet,
                         :size_depth,
                         :size_alphabet,
                         :depth_alphabet,
                         :full_dimensional]
        
        for mode in modes_to_test
            println("Executing mode: $mode...")
            
            compressed_forest = SolePostHoc.Orca.compression(
                original_forest, 
                mode, 
                f_val, 
                l_val;
                population_size=20, 
                n_generations=15,
                penalty_weight=0.5 
            )
            
            comp_acc = evaluate_forest(compressed_forest, f_val, l_val)
            comp_size = length(compressed_forest.models)
            
            size_reduction = round((1 - (comp_size / original_size)) * 100, digits=1)
            acc_loss = round(original_acc - comp_acc, digits=2)

            write(io, "MODE: $mode\n")
            write(io, " - Number of trees: $comp_size (Reduction: $size_reduction%)\n")
            write(io, " - Accuracy: $comp_acc% (Variation: $(acc_loss > 0 ? "-" : "+")$(abs(acc_loss))%)\n")
            
            write(io, "-------------------------------------------------\n")
        end
    end
end

run_usecases()