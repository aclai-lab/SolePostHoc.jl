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
function _thrs_with_boundary(
    config::LumenConfig{R,T},
    ensemble::LumenEnsemble{R,T},
    featurenames::Vector{Symbol},
    classnames::Vector{String},
    class_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    atoms = get_atoms(ensemble)

    depth = config.depth
    depth < 1.0 && (atms = _take_first_percentage(atoms, depth)) # TODO check it!

    nfeats = R(length(featurenames))
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

    return _thrs_with_boundary(thresholds, op_families)
end

function _thrs_with_boundary(
    thresholds::Vector{T},
    family::UInt8
) where {T<:AbstractFloat}
    nthrs = length(thresholds)
    result = Vector{T}(undef, nthrs + 1)
    result[1:nthrs] .= thresholds
    # copyto!(result, thresholds)

    # 0x01 (descending) → boundary point is BELOW the minimum threshold
    #                     prevfloat(last) because last is the smallest value
    # 0x02 (ascending)  → boundary point is ABOVE the maximum threshold
    #                     nextfloat(last) because last is the largest value
    result[end] = family === 0x01 ?
                  prevfloat(last(thresholds)) :
                  nextfloat(last(thresholds))

    return result
end

@inline _thrs_with_boundary(
    thresholds::Vector{Vector{T}},
    op_families::Vector{UInt8}
) where {T<:AbstractFloat} = _thrs_with_boundary.(thresholds, op_families)
