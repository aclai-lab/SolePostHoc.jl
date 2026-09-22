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
# 41.249 μs (593 allocations: 42.13 KiB) only lumenensemble and thrsspace
# 61.730 μs (1264 allocations: 86.45 KiB) with Minimize
# 1.929 ms (682 allocations: 1.92 MiB) my leaf_extract, apply
# 2.567 ms (933 allocations: 1.34 MiB) new leaf_extract, apply
# 2.452 ms (981 allocations: 1.71 MiB) without ScratchPad
# 2.371 ms (941 allocations: 1.34 MiB)
# 2.432 ms (966 allocations: 1.60 MiB)
# 2.477 ms (933 allocations: 1.34 MiB) scratchpad
# 2.395 ms (933 allocations: 1.34 MiB) width no Int
# 2.640 ms (947 allocations: 1.34 MiB)
# 2.547 ms (947 allocations: 1.34 MiB)
# 282.968 ms (14768 allocations: 2.27 MiB)


# # ---------------------------------------------------------------------------- #
# #                                  config                                      #
# # ---------------------------------------------------------------------------- #
# config = LumenShannonConfig(; minimization_scheme=Abc, M=4096)
# lumen = Lumen.lumen_shannon(config, model)

# R=UInt32; T=Float32
# featurenames, feat_idxs =
#     assign(R, unique!(SoleModels.info(model, :featurenames)), false)
# classnames, class_idxs =
#     assign(R, unique!(SoleModels.info(model, :supporting_labels)), true)

# ensemble = LumenEnsemble(config, model);
# thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth);

# lo = R.([1, 1, 1, 1])
# hi = R.([12, 7, 15, 11])

# cube = Lumen._leaf_extract(config, thrs, ensemble, lo, hi)

# @code_warntype Lumen._leaf_extract(config, thrs, ensemble, lo, hi)
# @code_warntype gather_atoms(thrs, R.([1,1,1,1]))

# @btime Lumen._leaf_extract(config, thrs, ensemble, lo, hi);
# # 2.046 ms (27763 allocations: 2.85 MiB)

# config = LumenShannonConfig(;
#     minimization_scheme=Abc, M=4096, float_resolution=Float64)
# R=UInt32; T=Float32
# featurenames, feat_idxs =
#     assign(R, unique!(SoleModels.info(model, :featurenames)), false)
# classnames, class_idxs =
#     assign(R, unique!(SoleModels.info(model, :supporting_labels)), true)
# ensemble = LumenEnsemble(config, model);
# thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth);
# lo = ones(R, length(thrs.nlev)); hi = copy(thrs.nlev)
# atoms, cube = Lumen._leaf_extract(config, thrs, ensemble, lo, hi);

# # one backend resolved once, reused for every class
# m = Minimizer(minimization_scheme(config), config)
# minimized = run_minimization(m, atoms[2], cube[2])

# @btime begin
#     for c in $class_idxs
#         minimized = run_minimization($m, $atoms[c], $cube[c])
#     end
# end

# # whole pipeline: recursion + minimization + DecisionSet
# ds = lumen_shannon(config, model)

# @btime lumen_shannon($config, $model);
# # 393.690 ms (30721 allocations: 5.12 MiB)
# # 307.620 ms (13368 allocations: 1.96 MiB)
# # 314.979 ms (15810 allocations: 2.12 MiB)
# # 308.738 ms (16133 allocations: 2.14 MiB) # claude, danger



