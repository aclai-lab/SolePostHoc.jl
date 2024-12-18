using DecisionTree: load_data, build_forest, apply_forest
using SolePostHoc
using SolePostHoc: Lumen
using SoleModels

nome_dataset = "iris"

X, y = load_data(nome_dataset)
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

f = solemodel(model)

Lumen.lumen(f, model)