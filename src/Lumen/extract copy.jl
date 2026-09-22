# ---------------------------------------------------------------------------- #
#                              leaf scratch space                              #
# ---------------------------------------------------------------------------- #
"""
    LeafScratch{R,T}

Every buffer a leaf needs, allocated once per run and reused by every leaf:
the apply table and per-row bookkeeping (sized by `M`), the per-class atom
pools, cube views and alphabets, and one PLA byte buffer per class. Peak
memory is therefore one leaf's worth, whatever the number of leaves, and a
leaf after the first allocates only what it grows past.
"""
struct LeafScratch{R<:Unsigned,T<:AbstractFloat}
    tbl::Matrix{T} # apply input, one chunk at a time
    nat::Vector{Int} # atoms produced by row k of the chunk
    preds::Vector{R} # class of every row of the leaf
    natoms_c::Vector{Int}
    ncubes_c::Vector{Int}
    acur::Vector{Int}
    ccur::Vector{Int}
    cube::Vector{LumenDNF{R,T}} # per class: its cubes, flat
    alphabet::Vector{Vector{LumenAtom{R,T}}} # per class: distinct atoms
    cur::Vector{R}
    widths::Vector{R}
end

function LeafScratch(
    config::LumenShannonConfig{R,T},
    thrs::ThresholdSpace{R,T}
) where {R<:Unsigned,T<:AbstractFloat}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    return LeafScratch{R,T}(
        Matrix{T}(undef, config.M, nfeats),
        Vector{Int}(undef, config.M),
        Vector{R}(undef, config.M),
        zeros(Int, nclasses), zeros(Int, nclasses),
        zeros(Int, nclasses), zeros(Int, nclasses),
        [LumenDNF{R,T}() for _ in 1:nclasses],
        [LumenAtom{R,T}[] for _ in 1:nclasses],
        Vector{R}(undef, nfeats),
        Vector{R}(undef, nfeats),
    )
end

# ---------------------------------------------------------------------------- #
#                           collect atoms for rule                             #
# ---------------------------------------------------------------------------- #
# Write the atoms of the cube at `levels` into `out[pos+1:pos+n]` and return
# the new fill position. `out` must already be long enough (`count_atoms`);
# no growth happens here, so the caller can hand out views into `out`.
function gather_atoms(
    out::LumenSlice{R,T},
    pos::Int,
    thrs::ThresholdSpace{R,T},
    levels::AbstractVector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    @inbounds for j in eachindex(levels)
        t = levels[j]
        off = thrs.thrs_offset[j] - one(R)
        nlev = thrs.nlev[j]
        feat = thrs.feat_idxs[j]
        # level t is sampled at threshold t (`_leaf_extract` feeds
        # `thrs[off + t]` to `apply`), so its region is closed on threshold t
        # and open on threshold t-1: the two bounds use DIFFERENT thresholds.
        thr_t = thrs.thrs[off + t]
        thr_p = t > 1 ? thrs.thrs[off + t - one(R)] : thr_t

        if thrs.op_families[j] === evalop(<)
            # '<' family, thresholds descending: region t = [thr_t, thr_{t-1})
            t > 1    && (out[pos += 1] = LumenAtom{R,T}(feat, thr_p, evalop(<)))
            t < nlev && (out[pos += 1] = LumenAtom{R,T}(feat, thr_t, evalop(≥)))
        else
            # '≤' family, thresholds ascending: region t = (thr_{t-1}, thr_t]
            t < nlev && (out[pos += 1] = LumenAtom{R,T}(feat, thr_t, evalop(≤)))
            t > 1    && (out[pos += 1] = LumenAtom{R,T}(feat, thr_p, evalop(>)))
        end
    end

    return pos
end

# TODO PROVA A VEDERE CHE SUCCEDE SU ELIMINI ENTRAMBI I DUALI,NON UNO SOLO
# drop atoms whose dual (same feat/thr, complementary op) was already kept.
function removeduals!(atoms::Vector{LumenAtom{R,T}}) where {R,T}
    seen = Set{LumenAtom{R,T}}()
    filter!(atoms) do a
        dual(a) in seen && return false
        push!(seen, a)
        return true
    end
end

# ---------------------------------------------------------------------------- #
#                        leaf extractor (Lumen legacy)                         #
# ---------------------------------------------------------------------------- #
"""
    _leaf_extract!(ws, config, thrs, ensemble, lo, hi) -> (alphabet, cube)

Materialise and classify every point of the sub-rectangle `[lo, hi]` and
sort its cubes by class. The result lives in `ws` (`ws.alphabet`, `ws.cube`,
backed by `ws.pool`) and is overwritten by the next leaf, so consume it
before calling again.
"""
function _leaf_extract!(
    config::LumenShannonConfig{R,T,MS},
    ws::LeafScratch{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat,MS}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    # widths = hi .- lo .+ one(R)
    # total = prod(widths)
    # nrows = min(Int(config.M), total)

    widths = ws.widths
    widths .= hi .- lo .+ one(R)
    total = prod(widths)

    # tbl = Matrix{T}(undef, nrows, nfeats) # apply input, one chunk at a time
    # nat = Vector{Int}(undef, nrows) # atoms produced by row k of the chunk
    # preds = Vector{R}(undef, total) # class of every row
    # natoms_c = zeros(Int, nclasses)
    # ncubes_c = zeros(Int, nclasses)

    tbl = ws.tbl; nat = ws.nat; preds = ws.preds
    natoms_c = fill!(ws.natoms_c, 0)
    ncubes_c = fill!(ws.ncubes_c, 0)

    # pass 1: classify every row and size the per-class output
    @inbounds for k in 1:total
        r = k - 1
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

    apply!(preds, ensemble, view(tbl, 1:total, :), nclasses)

    @inbounds for k in 1:total
        c = preds[k]
        natoms_c[c] += nat[k]
        ncubes_c[c] += 1
    end

    # per-class buffers keep their capacity across leaves: a resize! here
    # allocates only when a leaf outgrows every previous one. A cube is an
    # offset pair into the class' flat atom pool, not a 40-byte view.
    @inbounds for c in 1:nclasses
        d = ws.cube[c]
        resize!(d.atoms, natoms_c[c])
        resize!(d.offsets, ncubes_c[c] + 1)
        d.offsets[1] = 0
    end
    acur = fill!(ws.acur, 0)
    ccur = fill!(ws.ccur, 0)

    # pass 2: rewind the odometer and fill by cursor, straight from its digits
    cur = copyto!(ws.cur, lo)
    @inbounds for i in 1:total
        c = preds[i]
        d = ws.cube[c]
        a0 = acur[c]
        a1 = gather_atoms(view(d.atoms, :), a0, thrs, cur)
        acur[c] = a1
        d.offsets[(ccur[c] += 1) + 1] = a1

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

    # the cubes index into the pool, so the alphabet is gathered into its own
    # vector; the pool itself is never reordered or shrunk
    @inbounds for c in 1:nclasses
        _alphabet!(ws.alphabet[c], ws.cube[c].atoms)
    end

    return ws.alphabet, ws.cube
end

# distinct atoms of `pool`, sorted, duals removed, written into `out`
function _alphabet!(out::Vector{LumenAtom{R,T}}, pool::Vector{LumenAtom{R,T}}) where {R,T}
    empty!(out)
    seen = Set{LumenAtom{R,T}}()
    for a in pool
        a in seen && continue
        push!(seen, a); push!(out, a)
    end
    return removeduals!(sort!(out))
end

"""
    _leaf_extract(config, thrs, ensemble, lo, hi) -> (alphabet, cube)

One-off form of `_leaf_extract!` with its own scratch space.
"""
_leaf_extract(
    config::LumenShannonConfig{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat} =
    _leaf_extract!(LeafScratch(config, thrs), config, thrs, ensemble, lo, hi)

# ---------------------------------------------------------------------------- #
#                            leaf: minimise per class                          #
# ---------------------------------------------------------------------------- #
# Materialise the sub-rectangle, minimise each class' cubes in isolation and
# append the resulting terms to `out[c]`. Classes are independent and the
# backend is an external process, so they run in parallel; each thread only
# ever touches its own `out[c]`.
function _leaf_minimize!(
    config::LumenShannonConfig{R,T},
    out::Vector{LumenDNF{R,T}},
    ws::LeafScratch{R,T},
    # m::Minimizer,
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    alphabet, cube = _leaf_extract!(config, ws, thrs, ensemble, lo, hi)
    # Threads.@threads for c in eachindex(out)
    #     # append!(out[c], run_minimization(m, alphabet[c], cube[c]))
    #     append!(out[c], run_minimization(alphabet[c], cube[c]))
    # end
    return out
end

# product of the widths, saturating at `M + 1` so it can never overflow `R`
# function _rect_size(lo::Vector{R}, hi::Vector{R}, M::R) where {R<:Unsigned}
#     s = one(UInt64)
#     @inbounds for j in eachindex(lo)
#         s *= UInt64(hi[j] - lo[j]) + one(UInt64)
#         s > M && return UInt64(M) + one(UInt64)
#     end
#     return s
# end

# ---------------------------------------------------------------------------- #
#                    combining two cofactors: PLAIN UNION ONLY                 #
# ---------------------------------------------------------------------------- #
# Sibling cofactors are combined by plain set union of their terms, class by
# class, and never re-minimised: a term that dropped the literal of the split
# feature is a don't-care only INSIDE the sub-rectangle it came from, and
# feeding it back to the minimiser together with the sibling's terms would let
# other terms be discarded as "redundant" when they are not. Each cofactor is
# already exact for its own (disjoint) sub-rectangle, and the two exactly
# partition the parent's, so the union is exact too. This is intentionally the
# ONLY combination strategy in this file.

# ---------------------------------------------------------------------------- #
#                           extractor, shannon lumen                           #
# ---------------------------------------------------------------------------- #
function _extract!(
    # m::Minimizer,
    config::LumenShannonConfig{R,T},
    out::Vector{LumenDNF{R,T}},
    ws::LeafScratch{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    prod(hi .- lo .+ one(R)) ≤ config.M &&
    # _rect_size(lo, hi, config.M) ≤ config.M &&
        # return _leaf_minimize!(out, ws, m, config, thrs, ensemble, lo, hi)
        return _leaf_minimize!(config, out, ws, thrs, ensemble, lo, hi)

    # widest axis, without a temporary
    jstar = 1
    @inbounds for j in 2:length(lo)
        hi[j] - lo[j] > hi[jstar] - lo[jstar] && (jstar = j)
    end
    t = lo[jstar] + (hi[jstar] - lo[jstar]) >> 1

    old_hi = hi[jstar]; hi[jstar] = t
    # _extract!(out, ws, m, config, thrs, ensemble, lo, hi)
    _extract!(config, out, ws, thrs, ensemble, lo, hi)
    hi[jstar] = old_hi

    old_lo = lo[jstar]; lo[jstar] = t + one(R)
    # _extract!(out, ws, m, config, thrs, ensemble, lo, hi)
    _extract!(config, out, ws, thrs, ensemble, lo, hi)
    lo[jstar] = old_lo

    return out
end

"""
    _extract(m, config, thrs, ensemble, lo, hi) -> Vector{LumenDNF}

Recursive Shannon-style cofactor decomposition of the rectangle `[lo, hi]`
of the threshold-level space: halve the widest axis until a piece holds at
most `config.M` points, minimise every piece per class, union the results.
Returns one `LumenDNF` per class.
"""
function _extract(
    # m::Minimizer,
    config::LumenShannonConfig{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    out = [LumenDNF{R,T}() for _ in 1:length(thrs.class_idxs)]
    ws = LeafScratch(config, thrs)
    # return _extract!(out, ws, m, config, thrs, ensemble, copy(lo), copy(hi))
    return _extract!(config, out, ws, thrs, ensemble, copy(lo), copy(hi))
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
