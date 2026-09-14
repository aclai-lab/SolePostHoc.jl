using SolePostHoc.Lumen

using BenchmarkTools
using SoleXplorer
using RDatasets
using DataFrames
using Random

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

# ---------------------------------------------------------------------------- #
#                                  config                                      #
# ---------------------------------------------------------------------------- #
config = LumenShannonConfig(; minimization_scheme=Abc, M=5000)

R=UInt32; T=Float32

lumen = Lumen.lumen_shannon(config, model)

@btime Lumen.lumen_shannon(config, model);
# 42.567 μs (429 allocations: 31.08 KiB)
# 36.340 μs (458 allocations: 41.44 KiB)
# 41.962 μs (631 allocations: 47.50 KiB)
# 43.760 μs (853 allocations: 56.88 KiB)
# 1.850 ms (1006 allocations: 548.60 KiB)

# ---------------------------------------------------------------------------- #
#                               minimization                                   #
# ---------------------------------------------------------------------------- #
raw, config = Lumen.lumen_shannon(config, model)

Lumen.run_minimization(Abc, config, raw[1])