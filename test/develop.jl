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

# LumenEnsemble

# ---------------------------------------------------------------------------- #
#                                  config                                      #
# ---------------------------------------------------------------------------- #
config = LumenShannonConfig(; minimization_scheme=Abc, M=5000)
lumen = Lumen.lumen_shannon(config, model)

R=UInt32; T=Float32
featurenames, feat_idxs =
    assign(R, unique!(SoleModels.info(model, :featurenames)), false)
classnames, class_idxs =
    assign(R, unique!(SoleModels.info(model, :supporting_labels)), true)

ensemble = LumenEnsemble(config, model);
thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth);

lo = R.([1, 1, 1, 1])
hi = R.([12, 7, 15, 11])

cube = Lumen._leaf_extract(config, thrs, ensemble, lo, hi)

@code_warntype Lumen._leaf_extract(config, thrs, ensemble, lo, hi)
@code_warntype gather_atoms(thrs, R.([1,1,1,1]))

@btime Lumen._leaf_extract(config, thrs, ensemble, lo, hi);
# 2.046 ms (27763 allocations: 2.85 MiB)
# ---------------------------------------------------------------------------- #
#                              lumen ensemble                                  #
# ---------------------------------------------------------------------------- #
@btime Lumen.lumen_shannon(config, model);
# 42.567 μs (429 allocations: 31.08 KiB)
# 36.340 μs (458 allocations: 41.44 KiB)
# 41.962 μs (631 allocations: 47.50 KiB)
# 43.760 μs (853 allocations: 56.88 KiB)
# 1.850 ms (1006 allocations: 548.60 KiB)
# 1.929 ms (917 allocations: 469.75 KiB)
# 1.955 ms (860 allocations: 470.33 KiB) # wrong results!
# 1.793 ms (836 allocations: 1.95 MiB) # single vector with offset
# 1.813 ms (823 allocations: 1.42 MiB) # get rid of cube
# 2.555 ms (857 allocations: 1.43 MiB)
# 2.442 ms (863 allocations: 1.43 MiB) # unique and sort
# 2.604 ms (743 allocations: 1.96 MiB) # output both atoms than cube

# ---------------------------------------------------------------------------- #
#                                    pla                                       #
# ---------------------------------------------------------------------------- #
config = LumenShannonConfig(;
    minimization_scheme=Abc, M=5000, float_resolution=Float64)
R=UInt32; T=Float32
featurenames, feat_idxs =
    assign(R, unique!(SoleModels.info(model, :featurenames)), false)
classnames, class_idxs =
    assign(R, unique!(SoleModels.info(model, :supporting_labels)), true)
ensemble = LumenEnsemble(config, model);
thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth);
atoms, cube = Lumen.lumen_shannon(config, model);

minimized = run_minimization(minimization_scheme(config), config, atoms[2], cube[2], feat_idxs)



@btime begin
    for c in $class_idxs
        minimized = run_minimization(
            minimization_scheme($config), $config, $atoms[c], $cube[c])
    end
end

# 22.686 μs (78 allocations: 8.53 KiB)
# 21.382 μs (84 allocations: 8.81 KiB) # iobuffer
# 21.485 μs (93 allocations: 9.19 KiB) # header
# 44.808 μs (846 allocations: 49.38 KiB) reintroduced cube


