using Test
using SoleModels
using SolePostHoc
using MLJ
using DataFrames
using Random

@testset "ORCA Genetic Compression Tests" begin

    Xc, yc = @load_iris
    Xc = DataFrame(Xc)
    

    Random.seed!(42)
    ttpairs = MLJ.MLJBase.train_test_pairs(Holdout(; shuffle = true), 1:length(yc), yc)
    train_idx = ttpairs[1][1]
    val_idx   = ttpairs[1][2]
    
    n_trees_original = 10
    DTModel = MLJ.@load RandomForestClassifier pkg=DecisionTree verbosity=0
    model = DTModel(n_trees = n_trees_original)
    mach = machine(model, Xc, yc)
    MLJ.fit!(mach, rows = train_idx, verbosity = 0)
    
    featurenames = MLJ.report(mach).features
    classlabels = mach.fitresult[2][sortperm((mach).fitresult[3])]
    original_forest = solemodel(MLJ.fitted_params(mach).forest; featurenames, classlabels)

    f_val = Matrix{Float64}(Xc[val_idx, :])
    l_val = string.(yc[val_idx])

    @test original_forest isa DecisionEnsemble
    @test length(original_forest.models) == n_trees_original


    
    @testset "Mode: :size (Tree Selection)" begin
        compressed_f = SolePostHoc.compression(
            original_forest, :size, f_val, l_val;
            population_size=10, n_generations=5
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) <= n_trees_original
    end

    @testset "Mode: :depth (Tree Pruning)" begin
        compressed_f = SolePostHoc.compression(
            original_forest, :depth, f_val, l_val;
            population_size=10, n_generations=5
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) == n_trees_original 
    end

    @testset "Mode: :full_dimensional (Full Pipeline)" begin
        compressed_f = SolePostHoc.compression(
            original_forest, :full_dimensional, f_val, l_val;
            population_size=10, n_generations=5
        )
        
        @test compressed_f isa DecisionEnsemble
        @test length(compressed_f.models) <= n_trees_original
    end

    @testset "Edge Cases and Errors" begin
        @test_throws ErrorException SolePostHoc.compression(
            original_forest, :modalita_inesistente, f_val, l_val
        )
    end
end