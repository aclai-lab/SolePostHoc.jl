using SolePostHoc.Lumen

using BenchmarkTools
using SoleXplorer
using RDatasets
using DataFrames
using Random

using SoleModels

# ---------------------------------------------------------------------------- #
#                                 iris test                                    #
# ---------------------------------------------------------------------------- #
iris_df = dataset("datasets", "iris")
y = String.(iris_df.Species)
Xdf = select(iris_df, Not(:Species))

dsc = setup_dataset(
    Xdf, y; model=RandomForestClassifier(n_trees=20, max_depth=5), rng=42)
modelc = solexplorer(dsc)
model = modelc.sole[1]

config = LumenShannonConfig(; minimization_scheme=Abc, M=4096)
lumen = Lumen.lumen_shannon(config, model);
@btime Lumen.lumen_shannon(config, model);
# 277.899 ms (23216 allocations: 2.68 MiB)

# starting point: 1.194 s (6772068 allocations: 358.79 MiB)