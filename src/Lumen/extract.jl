# ---------------------------------------------------------------------------- #
#                        leaf extractor (Lumen legacy)                         #
# ---------------------------------------------------------------------------- #
function _leaf_extract(
    config::LumenConfig{R,T},
    atomcache::AtomCache{R,T},
    ensemble::LumenEnsemble{R,T},
    thrs_with_boundary::Vector{Vector{T}},
    nfeats::R,
    nclasses::R,
    lo::Vector{R},
    hi::Vector{R},
) where {R<:Unsigned,T<:AbstractFloat}
    @show Int.(hi)
    @show Int.(lo)
    raw = [Vector{Vector{LumenAtom}}() for _ in one(R):nclasses]
    widths = [hi[j] - lo[j] + one(R) for j in 1:nfeats]
    total = prod(R, widths)
    @show Int.(widths)
    @show Int(total)
    batch = min(config.max_apply_batch, total)

    tbl  = Matrix{T}(undef, batch, nfeats)
    idxm = Matrix{R}(undef, batch, nfeats)

    i0 = one(R)
    while i0 ≤ total
        this_chunk = min(batch, total - i0 + one(R))

        @inbounds for k in 1:this_chunk
            r = i0 + k - 2
            for j in 1:nfeats
                off = r % widths[j]
                r   = r ÷ widths[j]
                t   = lo[j] + off
                tbl[k, j]  = thrs_with_boundary[j][t]
                idxm[k, j] = atomcache.regidx[j][t]
            end
        end

        preds = apply(ensemble, view(tbl, 1:this_chunk, :), nclasses)

        @inbounds for k in 1:this_chunk
            p  = preds[k]

            len = 0
            for j in 1:nfeats
                len += length(atomcache.nodes[j][idxm[k, j]])
            end

            cube = Vector{LumenAtom}(undef, len)   # the only per-row allocation
            q = 0
            for j in 1:nfeats
                pj = atomcache.nodes[j][idxm[k, j]]
                for a in pj
                    cube[q += 1] = a             # pointer copy, no construction
                end
            end
            push!(raw[p], cube)
        end

        i0 += this_chunk
    end

    @show length(raw[1])
    @show length(raw[2])
    @show length(raw[3])

    # raw

    # terms = Vector{Vector{TERM}}(undef, nclasses)
    # # classes are independent; `run_minimization` shells out to an external
    # # binary, so this is both thread-safe and mostly I/O-bound.
    # Threads.@threads for c in 1:nclasses
    #     rc = raw[c]
    #     terms[c] = if isempty(rc)
    #         TERM[]
    #     elseif length(rc) == 1
    #         # single cube: already minimal, skip the subprocess round-trip
    #         TERM[SL.LeftmostConjunctiveForm(rc[1])]
    #     else
    #         run_minimization(config.minimization_scheme, config, rc)
    #     end
    # end

    # return terms
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
    config::LumenConfig{R,T},
    # ctx::Ctx{R,T},
    thrs_with_p::Vector{Vector{T}},
    ensemble::LumenEnsemble{R,T},
    lo::Vector{R},
    hi::Vector{R},
    # cache::RegionCache
) where {R<:Unsigned,T<:AbstractFloat}
    rect_size = prod(hi .- lo .+ one(R))
    rect_size ≤ config.M &&
        return _leaf_extract(config, thrs_with_p, ensemble, lo, hi)

    jstar = argmax(hi .- lo .+ 1)

    t = lo[jstar] + (hi[jstar] - lo[jstar]) >> 1
    old_hi = hi[jstar]; hi[jstar] = t
    terms_low = _shannon_extract(config, thrs_with_p, ensemble, lo, hi)
    hi[jstar] = old_hi

    old_lo = lo[jstar]; lo[jstar] = t + one(R)
    terms_high = _shannon_extract(config, thrs_with_p, ensemble, lo, hi)
    lo[jstar] = old_lo

    # return _combine_cofactors(terms_low, terms_high)
end