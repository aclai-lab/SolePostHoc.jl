using Test
using SoleModels
using SolePostHoc
using DecisionTree  
using DataFrames
using Random

@testset "ORCA Genetic Compression Tests" begin

Random.seed!(42)

    n_samples = 100
    n_features = 4
    X_matrix = rand(n_samples, n_features)
    
    classlabels = ["classA", "classB", "classC"]
    
    yc_int = rand(1:length(classlabels), n_samples)
    
    yc_string = [classlabels[i] for i in yc_int]

    Xc = DataFrame(X_matrix, :auto)
    featurenames = string.(names(Xc))

    train_idx = 1:80
    val_idx   = 81:100
    
    n_trees_original = 10
    n_subfeatures = 2
    
    raw_forest = DecisionTree.build_forest(yc_int[train_idx], X_matrix[train_idx, :], n_subfeatures, n_trees_original)
    
    original_forest = solemodel(raw_forest; featurenames=featurenames, classlabels=classlabels)

    f_val = X_matrix[val_idx, :]
    l_val = yc_string[val_idx] 

    @test original_forest isa DecisionEnsemble
    @test length(original_forest.models) == n_trees_original


    @testset "Mode: :size (Tree Selection)" begin
        compressed_f = SolePostHoc.Orca.compression(
            original_forest, :size, f_val, l_val;
            population_size=10, n_generations=5,
            seed=42 
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) <= n_trees_original
    end

    @testset "Mode: :depth (Tree Pruning)" begin
        compressed_f = SolePostHoc.Orca.compression(
            original_forest, :depth, f_val, l_val;
            population_size=10, n_generations=5,
            seed=42 
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) == n_trees_original 
    end

    @testset "Mode: :full_dimensional (Full Pipeline)" begin
        compressed_f = SolePostHoc.Orca.compression(
            original_forest, :full_dimensional, f_val, l_val;
            population_size=10, n_generations=5,
            seed=42 
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) <= n_trees_original
    end

    @testset "Edge Cases and Errors" begin
        @test_throws ErrorException SolePostHoc.Orca.compression(
            original_forest, :modalita_inesistente, f_val, l_val
        )
    end

    @testset "Determinism Check (RNG Passing)" begin
        rng1 = MersenneTwister(1234)
        rng2 = MersenneTwister(1234)

        compressed_1 = SolePostHoc.Orca.compression(
            original_forest, :full_dimensional, f_val, l_val;
            population_size=15, n_generations=10, rng=rng1
        )

        compressed_2 = SolePostHoc.Orca.compression(
            original_forest, :full_dimensional, f_val, l_val;
            population_size=15, n_generations=10, rng=rng2
        )
        
        @test length(compressed_1.models) == length(compressed_2.models)
        
        if !isempty(compressed_1.models)
             @test SolePostHoc.Orca.tree_depth(compressed_1.models[1]) == SolePostHoc.Orca.tree_depth(compressed_2.models[1])
        end
    end

    @testset "Zero-Tree Penalty Check (Philosophical Bug)" begin
        extreme_penalty = 100.0 

        compressed_extreme = SolePostHoc.Orca.compression(
            original_forest, :size, f_val, l_val;
            population_size=15, 
            n_generations=10, 
            penalty_weight=extreme_penalty
        )
        
        @test length(compressed_extreme.models) > 0
    end


end