# ---------------------------------------------------------------------------- #
#                       LUMEN SHANNON: the extraction driver                   #
#                                                                              #
#  The shape of the algorithm is unchanged: decompose the combination space    #
#  into pieces, describe each piece exactly, and recombine by plain union --   #
#  never by feeding already-minimised terms back into the minimiser, which is  #
#  what used to silently lose coverage (a term that dropped a split feature's  #
#  literal is a don't-care only INSIDE the sub-rectangle it came from).        #
#                                                                              #
#  Two things changed underneath it.                                           #
#                                                                              #
#  WHERE IT CUTS. The old recursion halved the widest axis until a piece fit   #
#  in `M` points, so its cost was set by the size of the grid -- a quantity    #
#  that is routinely astronomical and has nothing to do with the model. The    #
#  new one asks the compiled forest where it actually disagrees and cuts       #
#  there, stopping the moment a box provably gets one answer. Cost is now set  #
#  by how many regions the model genuinely distinguishes.                      #
#                                                                              #
#  WHAT IT CARRIES. A piece used to be materialised point by point: a row      #
#  tuple, a `NamedTuple` column table, a `PropositionalLogiset`, a truth-table #
#  rebuild per row per feature, and a freshly allocated `Vector{Atom}` cube.   #
#  A piece is now a pair of `Int32` bounds per feature, and a uniform region   #
#  of any size is ONE such pair rather than one cube per point inside it.      #
#                                                                              #
#  Coverage is untouched by all of this: every point of the space still ends   #
#  up in exactly one box, of exactly the class the model gives it.            #
# ---------------------------------------------------------------------------- #

"""
    ShannonState

Scratch and accumulators for one `lumen_shannon` run, allocated once up front.

Every buffer the hot loop needs lives here: the vote tally, the per-tree class
scratch, the straddle list used to choose a split, the depth-first stack, and
the per-class pending/sealed box sets. Nothing inside the decomposition
allocates.

`pending[c]` holds boxes of class `c` not yet given to the minimiser;
`sealed[c]` holds terms that already came back from it. Once a group is sealed
it is never re-minimised — that is the seal-then-combine discipline that keeps
the union step correct.
"""
mutable struct ShannonState
    sp::ThresholdSpace
    cf::CompiledForest
    cfg::LumenConfig
    session::MinimizerSession
    pending::Vector{BoxSet}
    sealed::Vector{BoxSet}
    votes::Vector{Float64}
    treeclasses::Vector{Int32}
    straddles::Vector{Int64}
    stack::Vector{Int32}
    nboxes::Int
    nsplits::Int
end

function ShannonState(sp::ThresholdSpace, cf::CompiledForest, cfg::LumenConfig,
                      session::MinimizerSession)
    nfeat = nfeatures(sp)
    nc = cf.nclass
    return ShannonState(
        sp, cf, cfg, session,
        [BoxSet(nfeat) for _ in 1:nc],
        [BoxSet(nfeat) for _ in 1:nc],
        zeros(Float64, nc),
        zeros(Int32, ntrees(cf)),
        Vector{Int64}(undef, ntrees(cf)),
        Int32[],
        0, 0,
    )
end

"""
    _emit!(st::ShannonState, class::Int32, lo, hi) -> Nothing

Record that the model answers `class` everywhere in the box `[lo, hi]`.

The box is first offered to the previous box of the same class, which absorbs it
whenever the two are adjacent along a single axis. When the pending group
reaches the configured `M`, it is sealed through the minimiser.
"""
function _emit!(st::ShannonState, class::Int32, lo::Vector{Int32}, hi::Vector{Int32})
    st.nboxes += 1
    p = st.pending[class]
    trymerge_last!(p, lo, hi) && return nothing
    pushbox!(p, lo, hi)
    (st.cfg.M > 0 && length(p) >= st.cfg.M && !st.cfg.explicit_offset) && _seal!(st, Int(class))
    return nothing
end

"""
    _seal!(st::ShannonState, c::Integer; offset_bs=nothing) -> Nothing

Minimise class `c`'s pending boxes and move the result into its sealed group.

The pending group is compacted first: merging adjacent boxes is exact and cheap,
and every box it removes is one product term the backend does not have to read,
encode and cover.
"""
function _seal!(st::ShannonState, c::Integer; offset_bs::Union{Nothing,BoxSet}=nothing)
    p = st.pending[c]
    p.n == 0 && return nothing

    st.cfg.compact && compact!(p; maxpasses=2, containment_cap=0)

    if st.cfg.minimize
        minimize_boxes!(st.sealed[c], st.session, st.sp, p;
                        offset_bs, encoding=st.cfg.encoding)
    else
        appendboxes!(st.sealed[c], p)
    end
    empty!(p)
    return nothing
end

"""
    decompose!(st::ShannonState) -> Nothing

Cut the whole combination space into boxes the model answers uniformly.

Depth-first over an explicit stack rather than the call stack: the decomposition
can descend once per distinct threshold along a path, which is thousands of
levels on a real forest, and an explicit stack neither overflows nor pays for a
stack frame per cut.

Each step asks [`classify_box!`](@ref) whether the ensemble is constant over the
current box. If it is, the box is emitted. If it is not, the box is split at the
cut most of the disagreeing trees test, which resolves that test for both halves
and so guarantees the descent terminates.
"""
function decompose!(st::ShannonState)
    nfeat = nfeatures(st.sp)
    stack = st.stack
    empty!(stack)

    lo = ones(Int32, nfeat)
    hi = copy(st.sp.nlev)
    append!(stack, lo); append!(stack, hi)

    while !isempty(stack)
        base = length(stack) - 2 * nfeat
        @inbounds for j in 1:nfeat
            lo[j] = stack[base+j]
            hi[j] = stack[base+nfeat+j]
        end
        resize!(stack, base)

        class, ns = classify_box!(st.cf, lo, hi, st.votes, st.treeclasses, st.straddles)

        if ns == 0
            _emit!(st, class, lo, hi)
            continue
        end

        j, c = best_split(st.straddles, ns)
        st.nsplits += 1

        # The stack is LIFO, so the upper half goes on first and the lower half
        # comes off first. Emitting boxes in increasing region order is what
        # lets `trymerge_last!` absorb siblings as they arrive.
        let saved = lo[j]
            lo[j] = c                       # upper half [c, hi[j]]
            append!(stack, lo); append!(stack, hi)
            lo[j] = saved
        end
        let saved = hi[j]
            hi[j] = c - Int32(1)            # lower half [lo[j], c-1]
            append!(stack, lo); append!(stack, hi)
            hi[j] = saved
        end
    end
    return nothing
end