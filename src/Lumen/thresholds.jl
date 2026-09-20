# ---------------------------------------------------------------------------- #
#                     grab thresholds from sole models                         #
# ---------------------------------------------------------------------------- #
# struct ThresholdSpace{R<:Unsigned,T<:AbstractFloat}
#     thresholds::Vector{T}
#     op_families::Vector{UInt8}
#     feat_idxs::Vector{R}
#     class_idxs::Vector{R}
# end

# Storage is flat on purpose: all features' thresholds live in one contiguous
# `Vector{T}` addressed through an offset table, so a whole space is two
# allocations rather than `2 * nfeat` of them, and the hot loops read from cache
# lines instead of chasing pointers.
# Thresholds are always stored DESCENDING.
struct ThresholdSpace{R<:Unsigned,T<:AbstractFloat}
    thrs::Vector{T}
    thrs_offset::Vector{R}
    op_families::Vector{UInt8}
    feat_idxs::Vector{R}
    class_idxs::Vector{R}
    nlev::Vector{R}
end

function ThresholdSpace(
    ensemble::LumenEnsemble{R,T},
    feat_idxs::Vector{R},
    class_idxs::Vector{R},
    depth::T
) where {R<:Unsigned,T<:AbstractFloat}
    atoms = get_atoms(ensemble)
    depth < 1.0 && (atoms = _take_first_percentage(atoms, depth)) # TODO check it!

    # Normalize every split to the '<' partition, as classic Lumen's
    # `_normalize_atom` does:  x ≥ t  splits ℝ exactly like  x < t;
    #   x > t  ⟺  x ≥ nextfloat(t)   and   x ≤ t  ⟺  x < nextfloat(t).
    # After this there is a single family: thresholds descending, boundary
    # sample `prevfloat(min)`, and `gather_atoms` emits '<' / '≥' only.
    op_families = fill(evalop(<), length(feat_idxs))
    normthr(a) = (a.op == evalop(>) || a.op == evalop(≤)) ? nextfloat(a.thr) : a.thr
    per_feat = [
        unique!(sort!([normthr(a) for a in atoms if a.feat == f], rev=true))
        for f in feat_idxs]
    nper_feat = length.(per_feat)
    # a feature the model never splits on has one level and no atoms; its
    # sample value is irrelevant, so any finite number will do.
    thrs_boundary = [isempty(t) ? zero(T) : _thrs_boundary(t, 0x01) for t in per_feat]
    thrs = reduce(
        vcat, ([t; b] for (t, b) in zip(per_feat, thrs_boundary)); init=T[])
    # added onr(R) to take into account the last added boundaty value
    thrs_offset = R.(cumsum([1; nper_feat[1:end-1] .+ one(R)]))
    nlev = R.(nper_feat) .+ one(R)

    return ThresholdSpace{R,T}(
        thrs, thrs_offset, op_families, feat_idxs, class_idxs, nlev)
end

function Base.show(io::IO, ts::ThresholdSpace{R,T}) where {R,T}
    print(io, "ThresholdSpace{", R, ",", T, "}(",
        length(ts.feat_idxs), " features, ",
        length(ts.class_idxs), " classes)")
end

function Base.show(io::IO, ::MIME"text/plain", ts::ThresholdSpace{R,T}) where {R,T}
    show(io, ts)
    isempty(ts.feat_idxs) && return
    println(io)
    for (i, f) in enumerate(ts.feat_idxs)
        println(io, "  V", Int(f), ": ",
            Int(ts.thrs_offset[i+1] - ts.thrs_offset[i]), " thresholds")
    end
    print(io, "  classes: ", length(ts.class_idxs))
end

# ---------------------------------------------------------------------------- #
#                            thresholds boundaries                             #
# ---------------------------------------------------------------------------- #
# Append the appropriate boundary point to `thresholds` depending on the
# operator family,
# ensuring that all `n + 1` ordinal regions induced by `n` thresholds
# are sampled.

# - `:lt` family (`<`/`≤`): thresholds sorted **descending** → appends
#   `prevfloat(last(thresholds))` 
#     to cover the region **below the smallest threshold**
#   (i.e. values smaller than every condition).

# - `:gt` family (`>`/`≥`): thresholds sorted **ascending** → appends
#   `nextfloat(last(thresholds))` 
#     to cover the region **above the largest threshold**
#   (i.e. values larger than every condition).
function _thrs_boundary(
    thresholds::Vector{T},
    family::UInt8
) where {T<:AbstractFloat}
    # 0x01 (descending) → boundary point is BELOW the minimum threshold
    #                     prevfloat(last) because last is the smallest value
    # 0x02 (ascending)  → boundary point is ABOVE the maximum threshold
    #                     nextfloat(last) because last is the largest value
    return family === 0x01 ?
                  prevfloat(last(thresholds)) :
                  nextfloat(last(thresholds))
end

@inline _thrs_boundary(
    thresholds::Vector{Vector{T}},
    op_families::Vector{UInt8}
) where {T<:AbstractFloat} = _thrs_boundary.(thresholds, op_families)
