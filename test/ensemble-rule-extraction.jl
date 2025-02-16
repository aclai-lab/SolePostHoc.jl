using Test
using DecisionTree: load_data, build_forest, apply_forest
using SolePostHoc
using SolePostHoc: Lumen
using SoleModels
using DataFrames


X, y = load_data("iris")
X = float.(X)
y = string.(y)


n_subfeatures = -1
n_trees = 4
partial_sampling = 0.7
min_samples_leaf = 5
max_depth = 3
min_samples_split = 2
min_purity_increase = 0.0
seed = 3

model = build_forest(
    y,
    X,
    n_subfeatures,
    n_trees,
    partial_sampling,
    max_depth,
    min_samples_leaf,
    min_samples_split,
    min_purity_increase;
    rng = seed,
)

@test_broken @test_nowarn Lumen.lumen(model; silent = true)
@test_logs (:info,) Lumen.lumen(model)

# ds = @test_logs (:warn,) extractrules(LumenRuleExtractor(), model, )
ds = @test_logs (:warn,) extractrules(LumenRuleExtractor(), model; silent = true)
apply_function = SoleModels.apply,

@test_broken @test_nowarn extractrules(InTreesRuleExtractor(), model, X, y)

@test_nowarn extractrules(InTreesRuleExtractor(), model, DataFrame((X), ["V$(i)" for i in 1:size(X, 2)]), y)

