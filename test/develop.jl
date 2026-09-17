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
R=UInt32; T=Float32

config = LumenShannonConfig(; minimization_scheme=Abc, M=5000)

featurenames, feat_idxs =
    assign(R, unique!(SoleModels.info(model, :featurenames)), false)
classnames, class_idxs =
    assign(R, unique!(SoleModels.info(model, :supporting_labels)), true)

ensemble = LumenEnsemble(config, model);
thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth);

lo = R.([1, 1, 1, 1])
hi = R.([12, 7, 15, 11])

@code_warntype Lumen._leaf_extract(config, thrs, ensemble, lo, hi)

@code_warntype gather_atoms(thrs, R.([1,1,1,1]))

@btime Lumen._leaf_extract(config, thrs, ensemble, lo, hi)
# 2.046 ms (27763 allocations: 2.85 MiB)
# ---------------------------------------------------------------------------- #
#                              lumen ensemble                                  #
# ---------------------------------------------------------------------------- #
@btime LumenEnsemble(config, model);
# 34.355 μs (395 allocations: 28.22 KiB)

lumen = Lumen.lumen_shannon(config, model)

# 42.567 μs (429 allocations: 31.08 KiB)
# 36.340 μs (458 allocations: 41.44 KiB)
# 41.962 μs (631 allocations: 47.50 KiB)
# 43.760 μs (853 allocations: 56.88 KiB)
# 1.850 ms (1006 allocations: 548.60 KiB)
# 1.929 ms (917 allocations: 469.75 KiB)
# 1.955 ms (860 allocations: 470.33 KiB) # wrong results!

# LumenEnsemble 33.713 μs (429 allocations: 31.09 KiB)
# LumenEnsemble new Atom struct 39.413 μs (525 allocations: 34.08 KiB)
# ThresholdSpace 47.400 μs (774 allocations: 49.70 KiB)

#  49.956 μs (782 allocations: 49.98 KiB)
# 2.122 ms (28549 allocations: 2.90 MiB)

# ---------------------------------------------------------------------------- #
#                                    pla                                       #
# ---------------------------------------------------------------------------- #
cache, thrs, raw = Lumen.lumen_shannon(config, model)

Lumen.write_pla(cache.nodes, thrs.feat_idxs, Univariate)

length.(cache.nodes)

Lumen.get_atoms(lumen)

get_feats(a::LumenAtom{R,T}, id::R) where {R<:Unsigned,T<:AbstractFloat} =
    get_thresholds(filter(a.op==id, a))

# countmap(predictions) = 3) => 2239, 2) => 9515, 1) => 2106)