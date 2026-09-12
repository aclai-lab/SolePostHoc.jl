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
    thresholds::Vector{T},
    family::UInt8
) where {T<:AbstractFloat}
    nthrs = length(thresholds)
    thrs_with_boundary = Vector{T}(undef, nthrs + 1)
    thrs_with_boundary[1:nthrs] .= thresholds
    # copyto!(thrs_with_boundary, thresholds)

    # 0x01 (descending) → boundary point is BELOW the minimum threshold
    #                     prevfloat(last) because last is the smallest value
    # 0x02 (ascending)  → boundary point is ABOVE the maximum threshold
    #                     nextfloat(last) because last is the largest value
    thrs_with_boundary[end] = family === 0x01 ?
                  prevfloat(last(thresholds)) :
                  nextfloat(last(thresholds))

    return thrs_with_boundary
end

@inline _thrs_with_boundary(
    thresholds::Vector{Vector{T}},
    op_families::Vector{UInt8}
) where {T<:AbstractFloat} = _thrs_with_boundary.(thresholds, op_families)

# ---------------------------------------------------------------------------- #
#                         prepare sequential context                           #
# ---------------------------------------------------------------------------- #
function _prepare_sequential_context(
    config::LumenConfig{R,T},
    ensemble::LumenEnsemble{R,T},
    nfeats::R,
    classnames::Vector{String},
    class_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    atoms = get_atoms(ensemble)
    @show length(atoms)

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

    thrs_with_boundary = _thrs_with_boundary(thresholds, op_families)

    return thresholds, thrs_with_boundary, op_families
end
