# ---------------------------------------------------------------------------- #
#                                 Lumen Atom                                   #
# ---------------------------------------------------------------------------- #
# struct LumenAtom{R<:Unsigned,T<:AbstractFloat}
#     feat::R
#     thr::T
#     op::UInt8
# end

# LumenAtom{R,T}() where {R<:Unsigned,T<:AbstractFloat} =
#     LumenAtom{R,T}(zero(R), T(NaN), 0xff)
# LumenAtom() = LumenAtom{UInt32,Float64}()

# isempty_atom(a::) = a.op == 0xff && iszero(a.feat)

# ---------------------------------------------------------------------------- #
#                              Lumen Atom Cache                                #
# ---------------------------------------------------------------------------- #
struct AtomCache{R<:Unsigned,T<:AbstractFloat}
    # nodes[j][r] : atoms implied by feature j being in ordinal region r
    nodes::Vector{Vector{LumenAtom}}
    # regidx[j][t] : region index for the t-th value of ctx.thrs_with_p[j]
    regidx::Vector{Vector{R}}
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
    thrs::ThresholdSpace{R,T},
) where {R<:Unsigned,T<:AbstractFloat}
    nfeats = length(thrs.feat_idxs)
    nodes = Vector{Vector{LumenAtom}}(undef, nfeats)
    regidx = Vector{Vector{R}}(undef, nfeats)

    @inbounds for f in thrs.feat_idxs
        thr = thrs.thresholds[f]
        n = length(thr)
        regs = LumenAtom[]

        if n === 0
            push!(regs, LumenAtom())
        elseif thrs.op_families[f] === evalop(<)
            # descending thresholds: idx0 = 1:r-1 -> `< thr[r-1]`
            #                        idx1 = r:n   -> `≥ thr[r]`
            lt = [_mknode(f, evalop(<), thr[k]) for k in 1:n]
            ge = [_mknode(f, evalop(≥), thr[k]) for k in 1:n]
            for r in 1:(n+1)
                r > 1 && push!(regs, lt[r-1])
                r ≤ n && push!(regs, ge[r])
            end
        else
            # ascending thresholds: minimum(idx0) ≡ 1, maximum(idx1) ≡ n,
            # so both atoms are region-independent.
            le = _mknode(f, evalop(≤), thr[1])
            gt = _mknode(f, evalop(>), thr[n])
            for r in 1:(n+1)
                r > 1 && push!(regs, le)
                r ≤ n && push!(regs, gt)
            end
        end

        nodes[f] = regs
        regidx[f] = R[
            (k = findfirst(==(v), thr); isnothing(k) ? n + 1 : k)
            for v in vcat(thrs.thresholds[f], last(thrs.thresholds[f]) + one(R))
        ]
    end

    return AtomCache{R,T}(nodes, regidx)
end