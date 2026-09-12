# ---------------------------------------------------------------------------- #
#                                 Lumen Atom                                   #
# ---------------------------------------------------------------------------- #
struct LumenAtom{R<:Unsigned,T<:AbstractFloat}
    feat::R
    thr::T
    op::UInt8
end

# ---------------------------------------------------------------------------- #
#                              Lumen Atom Cache                                #
# ---------------------------------------------------------------------------- #
struct AtomCache
    # parts[j][r] : atoms implied by feature j being in ordinal region r
    parts::Vector{Vector{Vector{LumenAtom}}}
    # regidx[j][t] : region index for the t-th value of ctx.thrs_with_p[j]
    # regidx::Vector{Vector{Int}}
    # plen[j][r] : length(parts[j][r]), to size the cube exactly
    # plen::Vector{Vector{Int}}
end

@inline _mknode(
    feat::R,
    op::UInt8,
    thr::T
) where {R<:Unsigned,T<:AbstractFloat} = LumenAtom(feat, thr, op)

"""
    AtomCache(ctx) -> AtomCache

Materialize, once, every atom `generate_disjunct` could ever emit.

`generate_disjunct` is a pure function of `(feature, truth-row)`, and the
truth-row is itself a pure function of the ordinal region index `r` (it is
`false` on `1:r-1`, `true` on `r:n` — see `_truths_row`). So there are only
`n+1` distinct outputs per feature. Building them up-front removes *all*
`Atom`/`ScalarCondition` construction, all `BitVector` allocation and all
`findall`/`findfirst` work from the per-row hot loop.

`regidx` reproduces `_truths_by_thresholds`' `findfirst(==(value), thr)`
lookup exactly (including its duplicate-threshold behaviour and its
"not found ⇒ row n+1" boundary case), so the result is bit-identical to the
previous implementation.
"""
function AtomCache(
    # ctx::Ctx{R,T}
    thresholds::Vector{Vector{T}},
    thrs_with_boundary::Vector{Vector{T}},
    op_families::Vector{UInt8},
    nfeats::R
) where {R<:Unsigned,T<:AbstractFloat}
    parts  = Vector{Vector{Vector{LumenAtom}}}(undef, nfeats)
    # regidx = Vector{Vector{Int}}(undef, nfeats)
    # plen   = Vector{Vector{Int}}(undef, nfeats)

    @inbounds for feat in one(R):nfeats
        thr = thresholds[feat]
        n = length(thr)
        regs = Vector{Vector{LumenAtom}}(undef, n + 1)

        if n == 0
            regs[1] = LumenAtom[]
        elseif op_families[feat] === evalop(<)
            # descending thresholds: idx0 = 1:r-1 -> `< thr[r-1]`
            #                        idx1 = r:n   -> `≥ thr[r]`
            lt = [_mknode(feat, evalop(<), thr[k]) for k in 1:n]
            ge = [_mknode(feat, evalop(≥), thr[k]) for k in 1:n]
            for r in 1:(n+1)
                a = LumenAtom[]
                r > 1 && push!(a, lt[r-1])
                r ≤ n && push!(a, ge[r])
                regs[r] = a
            end
        else
            # ascending thresholds: minimum(idx0) ≡ 1, maximum(idx1) ≡ n,
            # so both atoms are region-independent.
            le = _mknode(feat, evalop(≤), thr[1])
            gt = _mknode(feat, evalop(>), thr[n])
            for r in 1:(n+1)
                a = LumenAtom[]
                r > 1 && push!(a, le)
                r ≤ n && push!(a, gt)
                regs[r] = a
            end
        end

        parts[feat]  = regs
        # plen[feat]   = Int[length(x) for x in regs]
        # regidx[feat] = Int[
        #     (k = findfirst(==(v), thr); isnothing(k) ? n + 1 : k)
        #     for v in thrs_with_boundary[feat]
        # ]
    end

    # return AtomCache(parts, regidx)
    return AtomCache(parts)
end