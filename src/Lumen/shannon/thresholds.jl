# ---------------------------------------------------------------------------- #
#                     grab thresholds from sole models                         #
# ---------------------------------------------------------------------------- #
"""
    ThresholdSpace{R,T}

The finite grid the model is explored on. The thresholds a forest uses on a
feature cut its axis into ordinal regions ("levels"), and the forest's
prediction is constant inside each box of regions, so classifying one
sample point per box classifies the whole input space.

For feature position `j` (`1:length(feat_idxs)`) with `n` thresholds there
are `nlev[j] = n + 1` levels. Its thresholds are stored DESCENDING, followed
by one boundary value below the smallest one (see `_thrs_boundary`), in
`thrs[thrs_offset[j] : thrs_offset[j] + n]`. Level `t` is sampled at
`thrs[thrs_offset[j] + t - 1]` and is the region `[thr_t, thr_{t-1})`, with
no upper bound for `t = 1` and no lower bound for `t = nlev[j]`.

Storage is flat on purpose: all features' thresholds live in one contiguous
`Vector{T}` addressed through an offset table, so a whole space is two
allocations rather than `2 * nfeat` of them, and the hot loops read from cache
lines instead of chasing pointers.

# Fields
- `thrs::Vector{T}`: every feature's thresholds plus its boundary value.
- `thrs_offset::Vector{R}`: where feature `j`'s block starts in `thrs`.
- `op_families::Vector{UInt8}`: the operator family of each feature; always
  `evalop(<)` after the normalization done by the constructor.
- `feat_idxs::Vector{R}`: the model's feature index of position `j`.
- `class_idxs::Vector{R}`: one entry per class; only its length is used here.
- `nlev::Vector{R}`: number of levels per feature, the grid's shape.
"""
struct ThresholdSpace{R<:Unsigned,T<:AbstractFloat}
    thrs::Vector{T}
    thrs_offset::Vector{R}
    op_families::Vector{UInt8}
    feat_idxs::Vector{R}
    class_idxs::Vector{R}
    nlev::Vector{R}
end

"""
    ThresholdSpace(ensemble, feat_idxs, class_idxs, depth) -> ThresholdSpace

Collect the split thresholds of `ensemble`, all of them when `depth ≥ 1`,
otherwise only the first `depth` fraction of each tree in breadth-first
order (`bfs_atoms`). Every split is normalized to the `<` family first,
since `x > t` and `x ≤ t` cut ℝ exactly like `x < nextfloat(t)`: after that
one threshold value is one cut, duplicates collapse, and every feature
shares the same descending layout.

A feature the model never splits on gets one level and a `NaN` sample: no
tree ever reads it, so its value is irrelevant.
"""
function ThresholdSpace(
    ensemble::LumenEnsemble{R,T},
    feat_idxs::Vector{R},
    class_idxs::Vector{R},
    depth::T
) where {R<:Unsigned,T<:AbstractFloat}
    atoms = depth < 1 ? bfs_atoms(ensemble, depth) : get_atoms(ensemble)

    # Normalize every split to the '<' partition, as classic Lumen's
    # `_normalize_atom` does:  x ≥ t  splits ℝ exactly like  x < t;
    #   x > t  ⟺  x ≥ nextfloat(t)   and   x ≤ t  ⟺  x < nextfloat(t).
    # After this there is a single family: thresholds descending, boundary
    # sample `prevfloat(min)`, and `gather_atoms` emits '<' / '≥' only.
    op_families = fill(evalop(<), length(feat_idxs))
    normthr(a) = (a.op == evalop(>) ||
        a.op == evalop(≤)) ? nextfloat(a.thr) : a.thr
    per_feat = [
        unique!(sort!([normthr(a) for a in atoms if a.feat == f], rev=true))
        for f in feat_idxs]
    nper_feat = length.(per_feat)
    # a feature the model never splits on has one level and no atoms; no tree
    # reads it, so its sample value is irrelevant and NaN stands in for it.
    thrs_boundary = [
        isempty(t) ? T(NaN) : _thrs_boundary(t, 0x01) for t in per_feat]
    thrs = reduce(
        vcat, ([t; b] for (t, b) in zip(per_feat, thrs_boundary)); init=T[])
    # each block is the feature's thresholds plus its boundary value, hence
    # the `+ one(R)` in the block length
    thrs_offset = R.(cumsum([1; nper_feat[1:end-1] .+ one(R)]))
    nlev = R.(nper_feat) .+ one(R)

    return ThresholdSpace{R,T}(
        thrs, thrs_offset, op_families, feat_idxs, class_idxs, nlev)
end

# compact form on one line, and a per-feature threshold count for the REPL
function Base.show(io::IO, ts::ThresholdSpace{R,T}) where {R,T}
    print(io, "ThresholdSpace{", R, ",", T, "}(",
        length(ts.feat_idxs), " features, ",
        length(ts.class_idxs), " classes)")
end

function Base.show(
    io::IO, ::MIME"text/plain", ts::ThresholdSpace{R,T}
) where {R,T}
    show(io, ts)
    isempty(ts.feat_idxs) && return
    println(io)
    for (i, f) in enumerate(ts.feat_idxs)
        println(io, "  V", Int(f), ": ", Int(ts.nlev[i]) - 1, " thresholds")
    end
    print(io, "  classes: ", length(ts.class_idxs))
end

# ---------------------------------------------------------------------------- #
#                            thresholds boundaries                             #
# ---------------------------------------------------------------------------- #
"""
    _thrs_boundary(thresholds, family) -> T

The extra sample value that, together with `thresholds` themselves, gives
one sample point to each of the `n + 1` ordinal regions `n` thresholds
induce. It does not modify `thresholds`; the caller appends it.

- `<` / `≤` family, thresholds sorted DESCENDING: `prevfloat(last(thresholds))`,
  for the region below the smallest threshold.
- `>` / `≥` family, thresholds sorted ASCENDING: `nextfloat(last(thresholds))`,
  for the region above the largest threshold.

`ThresholdSpace` normalizes everything to `<`, so only the first case occurs
there. The vector form maps the scalar one over features and families.
"""
function _thrs_boundary(
    thresholds::Vector{T},
    family::UInt8
) where {T<:AbstractFloat}
    # evalop(<) (descending) → boundary point is BELOW the minimum threshold
    #                     prevfloat(last) because last is the smallest value
    # evalop(>) (ascending)  → boundary point is ABOVE the maximum threshold
    #                     nextfloat(last) because last is the largest value
    return family in (evalop(<), evalop(≤)) ?
        prevfloat(last(thresholds)) :
        nextfloat(last(thresholds))
end

@inline _thrs_boundary(
    thresholds::Vector{Vector{T}},
    op_families::Vector{UInt8}
) where {T<:AbstractFloat} = _thrs_boundary.(thresholds, op_families)