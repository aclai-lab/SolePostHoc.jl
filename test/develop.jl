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
@btime LumenEnsemble(config, model);
# 34.355 μs (395 allocations: 28.22 KiB)

lumen = Lumen.lumen_shannon(config, model)

@btime Lumen.lumen_shannon(config, model);

# 42.567 μs (429 allocations: 31.08 KiB)
# 36.340 μs (458 allocations: 41.44 KiB)
# 41.962 μs (631 allocations: 47.50 KiB)
# 43.760 μs (853 allocations: 56.88 KiB)
# 1.850 ms (1006 allocations: 548.60 KiB)
# 1.929 ms (917 allocations: 469.75 KiB)
# 1.955 ms (860 allocations: 470.33 KiB) # wrong results!
# 1.731 ms (836 allocations: 1.95 MiB) # ok

# LumenEnsemble 33.713 μs (429 allocations: 31.09 KiB)
# LumenEnsemble new Atom struct 39.413 μs (525 allocations: 34.08 KiB)
# ThresholdSpace 47.400 μs (774 allocations: 49.70 KiB)

#  49.956 μs (782 allocations: 49.98 KiB)
# 2.122 ms (28549 allocations: 2.90 MiB)
# 1.806 ms (836 allocations: 1.95 MiB)

# ---------------------------------------------------------------------------- #
#                                    pla                                       #
# ---------------------------------------------------------------------------- #
conditions = Lumen.lumen_shannon(config, model)

flat_conditions = unique(collect(Iterators.flatten(Iterators.flatten(conditions))))
@btime unique(collect(Iterators.flatten(Iterators.flatten(conditions))))
# 1.100 ms (39 allocations: 2.52 MiB)

nfeats = length(thrs.feat_idxs)

# sort!(conditions; by=SD._scalarcondition_sortby)
# sort!(feat_idxs; by=syntaxstring)

# conditions = SD.removeduals(conditions)

# feat_condindxss = Vector{Vector{Int}}(undef, nfnames)
# feat_condnames = Vector{Vector{String}}(undef, nfnames)

# @inbounds for (i, feat) in enumerate(feat_idxs)
#     feat_condindxs = findall(c->SD.feature(c) == feat, conditions)
#     conds = filter(c->SD.feature(c) == feat, conditions)
#     condname = SoleLogics.syntaxstring.(conds; removewhitespaces, pretty_op)

#     feat_condindxss[i] = feat_condindxs
#     feat_condnames[i] = condname
# end

# feat_nconds = length.(feat_condindxss)

# includes, excludes = Vector{BitMatrix}(undef, nfnames),
# Vector{BitMatrix}(undef, nfnames)
# @inbounds for (i, feat_condindxs) in enumerate(feat_condindxss)
#     includes[i] = BitMatrix([
#         SD.includes(conditions[cond_i], conditions[cond_j]) for
#         cond_i in feat_condindxs, cond_j in feat_condindxs
#     ])
#     excludes[i] = BitMatrix([
#         SD.excludes(conditions[cond_j], conditions[cond_i]) for
#         cond_i in feat_condindxs, cond_j in feat_condindxs
#     ])
# end

# pla_header = if encoding == :multivariate
#     _header(feat_nconds, feat_condnames; has_offset)
# else
#     _header(conditions, feat_condnames; has_offset)
# end

# # --- ON-set rows -----------------------------------------------------
# conjuncts = _get_conjuncts(atoms)
# pla_onset_rows = Vector{String}(undef, length(conjuncts))
# Threads.@threads for i in eachindex(conjuncts)
#     row = _encode_disjunct(
#         conjuncts[i], feat_idxs, conditions, includes, excludes, feat_condindxss
#     )
#     pla_onset_rows[i] =
#         encoding == :multivariate ? _onset_rows(feat_nconds, row) : _onset_rows(row)
# end

# # --- OFF-set rows (NEW, only if `offset` given) -----------------------
# # Same conditions/includes/excludes/feat_condindxss as the ON-set: this
# # reuse is exactly what guarantees column alignment between the two
# # halves of the 
# pla_offset_rows = String[]
# if has_offset
#     offset_conjuncts = _get_conjuncts(offset)
#     pla_offset_rows = Vector{String}(undef, length(offset_conjuncts))
#     Threads.@threads for i in eachindex(offset_conjuncts)
#         row = _encode_disjunct(
#             offset_conjuncts[i], feat_idxs, conditions, includes, excludes, feat_condindxss
#         )
#         pla_offset_rows[i] =
#             encoding == :multivariate ? _offset_rows(feat_nconds, row) : _offset_rows(row)
#     end
# end

# all_rows = vcat(pla_onset_rows, pla_offset_rows)

# pla_content = join(
#     [
#         join(pla_header, "\n"),
#         ".p $(length(all_rows))",
#         join(all_rows, "\n"),
#         ".e",
#     ],
#     "\n",
# )

# return pla_content, feat_idxs