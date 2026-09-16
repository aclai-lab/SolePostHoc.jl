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
    # nlev::Vector{R} # TODO boh
end

function ThresholdSpace(
    ensemble::LumenEnsemble{R,T},
    feat_idxs::Vector{R},
    class_idxs::Vector{R},
    depth::T
) where {R<:Unsigned,T<:AbstractFloat}
    atoms = get_atoms(ensemble)
    depth < 1.0 && (atoms = _take_first_percentage(atoms, depth)) # TODO check it!

    op_families = zeros(UInt8, length(feat_idxs))
    @inbounds for a in atoms
        fam = a.op ≤ 0x02 ? 0x01 : 0x02
        prev = op_families[a.feat]
        prev === 0x00 ? (op_families[a.feat] = fam) :
        prev === fam || throw(ArgumentError(
            "Feature $(a.feat) mixes '<'/'≤' with '>'/'≥' operators; " *
            "a single op family is required."))
    end

    per_feat = [sort!(
        get_thresholds(atoms, f), rev=op_families[f] === 0x01 ? true : false
    ) for f in feat_idxs]
    thrs_boundary = _thrs_boundary(per_feat, op_families)
    thrs = reduce(
        vcat, ([t; b] for (t, b) in zip(per_feat, thrs_boundary)); init=T[])
    thrs_offset = R.(cumsum([1; length.(per_feat)]))

    return ThresholdSpace{R,T}(
        thrs, thrs_offset, op_families, feat_idxs, class_idxs)
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
