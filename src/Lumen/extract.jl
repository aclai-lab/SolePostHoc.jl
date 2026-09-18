# ---------------------------------------------------------------------------- #
#                           collect atoms for rule                             #
# ---------------------------------------------------------------------------- #
# Write the atoms of the cube at `levels` into `out[pos+1:pos+n]` and return
# the new fill position. `out` must already be long enough (`count_atoms`);
# no growth happens here, so the caller can hand out views into `out`.
function gather_atoms(
    out::Vector{LumenAtom{R,T}},
    pos::Int,
    thrs::ThresholdSpace{R,T},
    levels::AbstractVector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    @inbounds for j in eachindex(levels)
        t = levels[j]
        off = thrs.thrs_offset[j] - one(R)
        nlev = thrs.nlev[j]
        feat = thrs.feat_idxs[j]
        thr = thrs.thrs[off + t]

        opin, opout = thrs.op_families[j] === 0x01 ?
            (evalop(<), evalop(≥)) :
            (evalop(≤), evalop(>))

        # region t is bounded by threshold t (inclusive side) ...
        t < nlev && (out[pos += 1] = LumenAtom{R,T}(feat, thr, opin))
        # ... and by threshold t-1 (exclusive side)
        t > 1 && (out[pos += 1] = LumenAtom{R,T}(feat, thr, opout))
    end

    return pos
end

# ---------------------------------------------------------------------------- #
#                        leaf extractor (Lumen legacy)                         #
# ---------------------------------------------------------------------------- #
function _leaf_extract(
    config::LumenShannonConfig{R,T,MS},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat,MS}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    widths = hi .- lo .+ one(R)
    total = prod(widths)
    nrows = min(Int(config.M), total)

    tbl = Matrix{T}(undef, nrows, nfeats) # apply input, one chunk at a time
    nat = Vector{Int}(undef, nrows)       # atoms produced by row k of the chunk
    preds = Vector{R}(undef, total)       # class of every row
    ncubes_c = zeros(Int, nclasses)
    natoms_c = zeros(Int, nclasses)

    # pass 1: classify every row chunk by chunk and size the per-class output
    i0 = 1
    while i0 ≤ total
        this_chunk = min(nrows, total - i0 + 1)

        @inbounds for k in 1:this_chunk
            r = i0 + k - 2
            n = 0
            for f in 1:nfeats
                off = r % widths[f]
                r = r ÷ widths[f]
                t = R(lo[f] + off)
                tbl[k, f] = thrs.thrs[thrs.thrs_offset[f] + t - one(R)]
                n += (t < thrs.nlev[f]) + (t > one(R))
            end
            nat[k] = n
        end

        chunk_preds = apply(ensemble, view(tbl, 1:this_chunk, :), nclasses)

        @inbounds for k in 1:this_chunk
            c = chunk_preds[k]
            preds[i0 + k - 1] = c
            ncubes_c[c] += 1
            natoms_c[c] += nat[k]
        end

        i0 += this_chunk
    end

    # every buffer is allocated once, at its final size: no regrowth copies
    atoms = [Vector{LumenAtom{R,T}}(undef, natoms_c[c]) for c in 1:nclasses]
    raw = [Vector{LumenCube{R,T}}(undef, ncubes_c[c]) for c in 1:nclasses]
    ccur = zeros(Int, nclasses)   # fill cursor into raw[c]
    acur = zeros(Int, nclasses)   # fill cursor into atoms[c]

    # pass 2: rewind the odometer and fill by cursor, straight from its digits
    cur = copy(lo)
    @inbounds for i in 1:total
        c = preds[i]
        buf = atoms[c]
        a0 = acur[c]
        a1 = gather_atoms(buf, a0, thrs, cur)
        acur[c] = a1
        raw[c][ccur[c] += 1] = view(buf, a0+1:a1)

        f = 1
        while f ≤ nfeats
            t = cur[f]
            if t < hi[f]
                cur[f] = t + one(R)
                break
            end
            cur[f] = lo[f]
            f += 1
        end
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
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
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
