# ---------------------------------------------------------------------------- #
#                           collect atoms for rule                             #
# ---------------------------------------------------------------------------- #
function gather_atoms(
    thrs::ThresholdSpace{R,T},
    levels::AbstractVector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    out = Vector{LumenAtom{R,T}}()
    sizehint!(out, 2 * length(levels))

    @inbounds for j in eachindex(levels)
        t = levels[j]
        off = thrs.thrs_offset[j] - one(R)
        nlev = thrs.nlev[j]
        feat = thrs.feat_idxs[j]

        opin, opout = thrs.op_families[j] === 0x01 ?
            (evalop(<), evalop(≥)) :
            (evalop(≤), evalop(>))

        # region t is bounded by threshold t (inclusive side) ...
        t < nlev && push!(out, LumenAtom{R,T}(feat, thrs.thrs[off + t], opin))
        # ... and by threshold t-1 (exclusive side)
        t > 1 && push!(out, LumenAtom{R,T}(feat, thrs.thrs[off + t], opout))
    end

    return out
end

using SoleModels, DataFrames
using CategoricalArrays

# ---------------------------------------------------------------------------- #
#                        leaf extractor (Lumen legacy)                         #
# ---------------------------------------------------------------------------- #
function _leaf_extract(
    config::LumenShannonConfig{R,T,MS},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
    model
) where {R<:Unsigned,T<:AbstractFloat,MS}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    # scratch = LumenAtom{R,T}[]
    counts_c = Vector{R}(undef, nclasses)
    cursors = Vector{R}(undef, nclasses)
    raw = [Vector{Vector{LumenAtom{R,T}}}() for _ in 1:nclasses]

    widths = hi .- lo .+ one(R)
    total = prod(widths)

    tbl  = Matrix{T}(undef, total, nfeats)
    idxm = Matrix{R}(undef, total, nfeats)

    i0 = 1
    while i0 ≤ total
        this_chunk = min(config.M, total - i0 + 1)

        @inbounds for k in 1:this_chunk
            r = i0 + k - 2
            for f in 1:nfeats
                off = r % widths[f]
                r = r ÷ widths[f]
                t = R(lo[f] + off)
                tbl[k, f] = thrs.thrs[thrs.thrs_offset[f] + t - one(R)]
                idxm[k, f] = t
            end
        end

        preds = apply(ensemble, view(tbl, 1:this_chunk, :), nclasses)

        # d = PropositionalLogiset(DataFrame(tbl[1:this_chunk, :], :auto))
        # preds = levelcode.(SoleModels.apply(
        #     model, d;
        #     use_multithreads=false,
        #     suppress_parity_warning=true
        # ))

        # pass 1: count cubes per class in this chunk
        fill!(counts_c, 0)
        @inbounds for k in 1:this_chunk
            counts_c[preds[k]] += 1
        end

        # grow each class vector once
        @inbounds for c in 1:nclasses
            cursors[c] = length(raw[c])
            resize!(raw[c], cursors[c] + counts_c[c])
        end

        # pass 2: fill by cursor
        @inbounds for k in 1:this_chunk
            c = preds[k]
            raw[c][cursors[c] += 1] = gather_atoms(thrs, @view idxm[k, :])
        end

        i0 += this_chunk
    end

    raw
    # terms = Vector{Vector{LumenAtom}}(undef, nclasses)
    # classes are independent; `run_minimization` shells out to an external
    # binary, so this is both thread-safe and mostly I/O-bound.
    # Threads.@threads for c in 1:nclasses
    #     rc = raw[c]
    #     terms[c] = if isempty(rc)
    #         LumenAtom[]
    #     elseif length(rc) == 1
    #         # single cube: already minimal, skip the subprocess round-trip
    #         LumenAtom[SL.LeftmostConjunctiveForm(rc[1])]
    #     else
    #         run_minimization(config.minimization_scheme, config, rc)
    #     end
    # end

    # develop
    # rc = raw[1]
    # run_minimization(MS, config, rc)

    # return raw
end

# ---------------------------------------------------------------------------- #
#                    combining two cofactors: PLAIN UNION ONLY                 #
# ---------------------------------------------------------------------------- #
# Combine the per-class results of two sibling cofactors by plain set union,
# class by class. No re-minimization, ever -- see the module-level comment at
# the top of this file for why: it's the only combination step that's
# correct with no extra bookkeeping, given that `terms_low`/`terms_high` are
# each already fully correct for their own (disjoint) sub-rectangle, and the
# two sub-rectangles exactly partition the parent's rectangle.

# This is intentionally the ONLY combination strategy available in this file
# -- there is no flag to opt back into re-minimization here, precisely so the
# sensitivity regression that motivated this file can't silently come back.
# function _combine_cofactors(
#     terms_low::Vector{Vector{TERM}},
#     terms_high::Vector{Vector{TERM}}
# )
#     # `append!` in place instead of `vcat`: halves the copying at every
#     # internal node (and there are O(leaves) of them). Semantics unchanged —
#     # `terms_low` is a freshly built, non-aliased value from the recursion.
#     @inbounds for c in eachindex(terms_low)
#         append!(terms_low[c], terms_high[c])
#     end
#     return terms_low
# end

# ---------------------------------------------------------------------------- #
#                           extractor, shannon lumen                           #
# ---------------------------------------------------------------------------- #
function _extract(
    config::LumenShannonConfig{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{I},
    hi::Vector{I},
) where {R<:Unsigned,T<:AbstractFloat,I}
    rect_size = prod(hi .- lo .+ one(R))
    rect_size ≤ config.M &&
        return _leaf_extract(config, thrs, ensemble, lo, hi)

    jstar = argmax(hi .- lo .+ 1)

    t = lo[jstar] + (hi[jstar] - lo[jstar]) >> 1
    old_hi = hi[jstar]; hi[jstar] = t
    terms_low = _shannon_extract(config, thrs, ensemble, lo, hi)
    hi[jstar] = old_hi

    old_lo = lo[jstar]; lo[jstar] = t + one(R)
    terms_high = _shannon_extract(config, thrs, ensemble, lo, hi)
    lo[jstar] = old_lo

    # return _combine_cofactors(terms_low, terms_high)
end

# function _extract(
#     thrs::ThresholdSpace{R,T},
#     st::ShannonState
# ) where {R<:Unsigned,T<:AbstractFloat}
#     nfeat = length(thrs.feat_idxs)
#     stack = R[]
#     # empty!(stack)

#     lo = ones(R, nfeat)
#     hi = copy(st.sp.nlev)
#     append!(stack, lo)
#     append!(stack, hi)

#     while !isempty(stack)
#         base = length(stack) - 2 * nfeat
#         @inbounds for j in 1:nfeat
#             lo[j] = stack[base+j]
#             hi[j] = stack[base+nfeat+j]
#         end
#         resize!(stack, base)

#         class, ns = classify_box!(st.cf, lo, hi, st.votes, st.treeclasses, st.straddles)

#         if ns == 0
#             _emit!(st, class, lo, hi)
#             continue
#         end

#         j, c = best_split(st.straddles, ns)
#         st.nsplits += 1

#         # The stack is LIFO, so the upper half goes on first and the lower half
#         # comes off first. Emitting boxes in increasing region order is what
#         # lets `trymerge_last!` absorb siblings as they arrive.
#         let saved = lo[j]
#             lo[j] = c                       # upper half [c, hi[j]]
#             append!(stack, lo); append!(stack, hi)
#             lo[j] = saved
#         end
#         let saved = hi[j]
#             hi[j] = c - R(1)            # lower half [lo[j], c-1]
#             append!(stack, lo); append!(stack, hi)
#             hi[j] = saved
#         end
#     end
#     return nothing
# end
