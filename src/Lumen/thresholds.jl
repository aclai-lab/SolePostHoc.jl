struct ThresholdSpace{R<:Unsigned,T<:AbstractFloat}
    thresholds::Vector{Vector{T}}
    boundaries::Vector{T}
    op_families::Vector{UInt8}
    feat_idxs::Vector{R}
    class_idxs::Vector{R}
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

# ---------------------------------------------------------------------------- #
#                         prepare sequential context                           #
# ---------------------------------------------------------------------------- #
function ThresholdSpace(
    config::LumenShannonConfig{R,T},
    ensemble::LumenEnsemble{R,T},
    featurenames::Vector{Symbol},
    feat_idxs::Vector{R},
    classnames::Vector{String},
    class_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    nfeats = length(feat_idxs)
    atoms = get_atoms(ensemble)

    depth = config.depth
    depth < 1.0 && (atoms = _take_first_percentage(atoms, depth)) # TODO check it!

    thresholds = get_thresholds(atoms, nfeats)
    
    op_families = zeros(UInt8, nfeats)
    @inbounds for a in atoms
        fam = a.op ≤ 0x02 ? 0x01 : 0x02
        prev = op_families[a.feat]
        prev == 0x00 ? (op_families[a.feat] = fam) :
        prev == fam || throw(ArgumentError(
            "Feature $(a.feat) mixes '<'/'≤' with '>'/'≥' operators; " *
            "a single op family is required."))
    end

    thrs_boundary = _thrs_boundary(thresholds, op_families)

    return ThresholdSpace(
        thresholds, thrs_boundary, op_families, feat_idxs, class_idxs)
end
