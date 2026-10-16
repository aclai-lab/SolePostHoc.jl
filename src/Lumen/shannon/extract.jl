# ---------------------------------------------------------------------------- #
#                              leaf scratch space                              #
# ---------------------------------------------------------------------------- #
"""
    LeafScratch{R,T}

Every buffer a leaf needs, allocated once per run and reused by every leaf:
the apply table and per-row bookkeeping (sized by the largest possible
leaf), and per class its cubes and its alphabet. Peak memory is therefore
one leaf's worth, whatever the number of leaves, and a leaf after the first
allocates only what it grows past. The PLA text itself never lives here: it
is streamed to a file by `write_pla`.
"""
struct LeafScratch{R<:Unsigned,T<:AbstractFloat}
    tbl::Matrix{T} # apply input: one sample point per row, one feature per column
    nat::Vector{Int} # number of atoms the cube of row k will have
    preds::Vector{R} # class of every row of the leaf
    natoms_c::Vector{Int} # per class: total atoms of its cubes
    ncubes_c::Vector{Int} # per class: number of its cubes
    acur::Vector{Int} # per class: fill cursor into `cube[c].atoms`
    ccur::Vector{Int} # per class: fill cursor into `cube[c].offsets`
    cube::Vector{LumenDNF{R,T}} # per class: its cubes, flat
    alphabet::Vector{Vector{LumenAtom{R,T}}} # per class: distinct atoms
    cur::Vector{R} # odometer: the level of each feature at the current row
    widths::Vector{R} # number of levels of each feature inside the leaf
    total::Int # row capacity of the buffers, NOT the size of the current leaf
end

# scratch space sized for the largest leaf `_extract!` can produce on `thrs`
function LeafScratch(
    config::LumenConfig{R,T},
    thrs::ThresholdSpace{R,T}
) where {R<:Unsigned,T<:AbstractFloat}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    # a leaf never holds more points than the whole space nor more than M:
    # the row buffers hold min(M, prod(nlev)) rows
    total = prod(thrs.nlev)
    nrows = Int(min(config.M, total))

    return LeafScratch{R,T}(
        Matrix{T}(undef, nrows, nfeats),
        Vector{Int}(undef, nrows),
        Vector{R}(undef, nrows),
        Vector{Int}(undef, nclasses),
        Vector{Int}(undef, nclasses),
        Vector{Int}(undef, nclasses),
        Vector{Int}(undef, nclasses),
        [LumenDNF{R,T}() for _ in 1:nclasses],
        [LumenAtom{R,T}[] for _ in 1:nclasses],
        Vector{R}(undef, nfeats),
        Vector{R}(undef, nfeats),
        nrows
    )
end

# ---------------------------------------------------------------------------- #
#                           collect atoms for rule                             #
# ---------------------------------------------------------------------------- #
"""
    gather_atoms(out, pos, thrs, levels) -> Int

Write the atoms of the cube at grid point `levels` (one level per feature)
into `out[pos+1:pos+n]` and return the new fill position `pos + n`. A level
is bounded on both sides by adjacent thresholds, except the first (no upper
bound) and the last (no lower bound), so a cube has at most two atoms per
feature, and none for a feature with a single level.

`out` must already be long enough: `_leaf_extract!` sizes it in its first
pass by counting the same atoms. No growth happens here, so the cubes'
offsets into `out` stay valid.
"""
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
        # level t is sampled at threshold t (`_leaf_extract!` feeds
        # `thrs[off + t]` to `apply!`), so its region is closed on threshold t
        # and open on threshold t-1: the two bounds use DIFFERENT thresholds.
        thr_t = thrs.thrs[off + t]
        thr_p = t > 1 ? thrs.thrs[off + t - one(R)] : thr_t

        if thrs.op_families[j] === evalop(<)
            # '<' family, thresholds descending: region t = [thr_t, thr_{t-1})
            t > 1 && (out[pos+=1] = LumenAtom{R,T}(feat, thr_p, evalop(<)))
            t < nlev && (out[pos+=1] = LumenAtom{R,T}(feat, thr_t, evalop(≥)))
        else
            # '≤' family, thresholds ascending: region t = (thr_{t-1}, thr_t].
            # Unreachable today: `ThresholdSpace` normalizes every feature to
            # the '<' family; kept for a space built without normalization.
            t < nlev && (out[pos+=1] = LumenAtom{R,T}(feat, thr_t, evalop(≤)))
            t > 1 && (out[pos+=1] = LumenAtom{R,T}(feat, thr_p, evalop(>)))
        end
    end

    return pos
end

# Drop, in place, every atom whose dual (same feature and threshold,
# complementary operator) appears earlier in `atoms`: the pair stands for a
# single PLA column, so one of the two is enough to build the layout.
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
    _leaf_extract!(config, ws, thrs, ensemble, lo, hi) -> (alphabet, cube)

Materialise and classify every point of the sub-rectangle `[lo, hi]` of the
level grid (inclusive, one range per feature) and sort its cubes by class.
Returns, per class `c`, `cube[c]`: one tight cube per point classified `c`,
and `alphabet[c]`: the distinct atoms of those cubes, duals removed, which
are the columns of that class' PLA.

Two passes over the leaf's rows, row `k` being point `k` of the leaf in
odometer order (first feature fastest):
1. decode each row's levels, write its sample values into `ws.tbl` and count
   its atoms; classify all rows with one `apply!` call; size each class'
   buffers from the counts.
2. walk the same points again with an odometer and write each cube's atoms
   straight into its class' buffer, so nothing is pushed or grown.

The result lives in `ws` (`ws.alphabet`, `ws.cube`) and is overwritten by
the next leaf, so consume it before calling again.
"""
function _leaf_extract!(
    ::LumenConfig{R,T,MS},
    ws::LeafScratch{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat,MS}
    nfeats = length(thrs.feat_idxs)
    nclasses = length(thrs.class_idxs)

    ws.widths .= hi .- lo .+ one(R)
    # this leaf's own size: `ws.total` is only the buffers' capacity, and
    # looping over it would revisit the leaf's points as duplicates
    total = prod(Int, ws.widths)

    # pass 0: reset the counters
    fill!(ws.natoms_c, 0)
    fill!(ws.ncubes_c, 0)
    fill!(ws.acur, 0)
    fill!(ws.ccur, 0)

    # pass 1: classify every row and size the per-class output
    @inbounds for k in 1:total
        r = k - 1
        n = 0
        for f in 1:nfeats
            off = r % ws.widths[f]
            r = r ÷ ws.widths[f]
            t = R(lo[f] + off)
            ws.tbl[k, f] = thrs.thrs[thrs.thrs_offset[f] + t - one(R)]
            n += (t < thrs.nlev[f]) + (t > one(R))
        end
        ws.nat[k] = n
    end

    apply!(ws.preds, ensemble, view(ws.tbl, 1:total, :), nclasses)

    @inbounds for k in 1:total
        c = ws.preds[k]
        ws.natoms_c[c] += ws.nat[k]
        ws.ncubes_c[c] += 1
    end

    # per-class buffers keep their capacity across leaves: a resize! here
    # allocates only when a leaf outgrows every previous one. A cube is an
    # offset pair into the class' flat atom vector, not a 40-byte view.
    @inbounds for c in 1:nclasses
        d = ws.cube[c]
        resize!(d.atoms, ws.natoms_c[c])
        resize!(d.offsets, ws.ncubes_c[c] + 1)
        d.offsets[1] = 0
    end
    # pass 2: rewind the odometer and fill by cursor, straight from its digits;
    # it visits the points in the same order as pass 1, so `preds[i]` is the
    # class of the point `cur` holds at step i
    cur = copyto!(ws.cur, lo)
    @inbounds for i in 1:total
        c = ws.preds[i]
        d = ws.cube[c]
        a0 = ws.acur[c]
        a1 = gather_atoms(view(d.atoms, :), a0, thrs, cur)
        ws.acur[c] = a1
        d.offsets[(ws.ccur[c] += 1) + 1] = a1

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

    # the cubes index into the atom vector, so the alphabet is gathered into
    # its own vector; the atom vector itself is never reordered or shrunk
    @inbounds for c in 1:nclasses
        _alphabet!(ws.alphabet[c], ws.cube[c].atoms)
    end

    return ws.alphabet, ws.cube
end

# distinct atoms of `pool`, sorted, duals removed, written into `out`: the
# columns of one class' PLA (see `PlaLayout`)
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

One-off form of `_leaf_extract!` with its own scratch space, sized for the
whole space (or `M`) rather than reused across leaves.
"""
_leaf_extract(
    config::LumenConfig{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat} =
    _leaf_extract!(LeafScratch(config, thrs), config, thrs, ensemble, lo, hi)

# ---------------------------------------------------------------------------- #
#                            leaf: minimise per class                          #
# ---------------------------------------------------------------------------- #
# Materialise the sub-rectangle `[lo, hi]`, minimise each class' cubes in
# isolation (backend `MS`, see `run_minimization`) and append the resulting
# terms to `out[c]`. Classes are independent and the backend is an external
# process, so they run in parallel; each thread only ever touches its own
# `out[c]`.
function _leaf_minimize!(
    config::LumenConfig{R,T,MS},
    out::Vector{LumenDNF{R,T}},
    ws::LeafScratch{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat,MS}
    alphabet, cube = _leaf_extract!(config, ws, thrs, ensemble, lo, hi)
    Threads.@threads for c in eachindex(out)
        append!(out[c], run_minimization(MS, config, alphabet[c], cube[c]))
    end
    return out
end

# ---------------------------------------------------------------------------- #
#                           extractor, shannon lumen                           #
# ---------------------------------------------------------------------------- #
# In-place recursion behind `_extract`. A rectangle of at most `config.M`
# points is a leaf and goes to `_leaf_minimize!`; a larger one is cut at the
# midpoint of its widest axis and both halves recurse. `lo` and `hi` are
# mutated to describe each half and restored before returning, so the whole
# recursion allocates no rectangles. Both halves append to the same `out`:
# their regions are disjoint, so the plain union is exact.
function _extract!(
    config::LumenConfig{R,T},
    out::Vector{LumenDNF{R,T}},
    ws::LeafScratch{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    prod(hi .- lo .+ one(R)) ≤ config.M &&
        return _leaf_minimize!(config, out, ws, thrs, ensemble, lo, hi)

    # widest axis, without a temporary
    jstar = 1
    @inbounds for j in 2:length(lo)
        hi[j] - lo[j] > hi[jstar] - lo[jstar] && (jstar = j)
    end
    t = lo[jstar] + (hi[jstar] - lo[jstar]) >> 1

    # lower half [lo, t] on the cut axis
    old_hi = hi[jstar]; hi[jstar] = t
    _extract!(config, out, ws, thrs, ensemble, lo, hi)
    hi[jstar] = old_hi

    # upper half [t + 1, hi] on the cut axis
    old_lo = lo[jstar]; lo[jstar] = t + one(R)
    _extract!(config, out, ws, thrs, ensemble, lo, hi)
    lo[jstar] = old_lo

    return out
end

"""
    _extract(config, thrs, ensemble, lo, hi) -> Vector{LumenDNF}

Recursive Shannon-style cofactor decomposition of the rectangle `[lo, hi]`
of the threshold-level space: halve the widest axis until a piece holds at
most `config.M` points, minimise every piece per class, union the results.
Returns one `LumenDNF` per class, indexed by class levelcode. `lo` and `hi`
are copied, never mutated; pass `ones` and `thrs.nlev` for the whole space.

Peak memory is one leaf's scratch space (`LeafScratch`) plus the output.
The terms of different leaves are never minimised together, so the result
is exact but not globally minimal; `compact!` recovers part of the gap.
"""
function _extract(
    config::LumenConfig{R,T},
    thrs::ThresholdSpace{R,T},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    out = [LumenDNF{R,T}() for _ in 1:length(thrs.class_idxs)]
    ws = LeafScratch(config, thrs)
    return _extract!(config, out, ws, thrs, ensemble, copy(lo), copy(hi))
end
