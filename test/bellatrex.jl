using Test
using SoleModels
using SoleData
using SolePostHoc
using DataFrames

@testset "BELLATREX" begin
    # 1. Non-ensemble model (single DecisionTree)
    t1 = SoleModels.DecisionTree(LeafModel("classA"))
    t2 = SoleModels.DecisionTree(LeafModel("classB"))
    df = DataFrame(feature1 = [1.0, 2.0], feature2 = [3.0, 4.0])
    X = scalarlogiset(df, allow_propositional=true)
    Y = ["classA", "classB"]

    @test !isensemble(t1)

    # Regression test for Issue #95: bellatrex on a non-ensemble model must refuse with a clear error
    # (previously raised UndefVarError: model not defined or MethodError on trees(m) before isensemble guard)
    @test_throws ErrorException bellatrex(t1, X, Y)
    try
        bellatrex(t1, X, Y)
    catch err
        @test err isa ErrorException
        @test contains(err.msg, "bellatrex")
        @test contains(err.msg, "ensemble")
    end

    # Also test via extractrules interface with BellatrexRuleExtractor
    @test_throws ErrorException extractrules(BellatrexRuleExtractor(), t1, X, Y)

    # 2. Ensemble model (DecisionForest)
    forest = SoleModels.DecisionForest([t1, t2])
    @test isensemble(forest)
    @test trees(forest) == [t1, t2]
end
