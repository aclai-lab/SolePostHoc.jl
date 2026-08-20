# ============================================================================ #
#   LUMEN SEQUENTIAL: fused generate-and-minimize, single-pass rule extraction #
#                                                                              #
#   Overview
#   --------
#   Classic Lumen builds the WHOLE Cartesian product of feature-threshold
#   combinations, applies the model to all of it, and only then minimizes
#   each class' collected cubes. That materializes an array whose size is
#   `prod(length.(thrs_with_p))`, which for real models can be astronomically
#   large and simply won't fit in memory.
#
#   This file replaces "generate everything, then minimize" with a single
#   fused loop:
#     1. Enumerate combinations one (random) index at a time, using an O(1)
#        lazy pseudo-random permutation (`LazyPermutation`) so we NEVER
#        materialize the full combination space, and different restarts see
#        genuinely different orders.
#     2. Route each row's derived cube into the buffer of its predicted
#        class.
#     3. The moment a class buffer reaches the memory budget `M`, minimize
#        immediately: call the minimizer over (that class' CURRENT
#        minimized formula) + (everything collected since the last such
#        call). The result REPLACES the current formula. Because
#        minimization usually collapses many raw rows into fewer
#        don't-care terms, this typically pulls the buffer back well under
#        `M`, so enumeration continues; if `M` is hit again, we fold again.
#     4. If a class' formula still can't be gotten under `M` even after
#        folding, the whole restart is abandoned (too unlucky an order to
#        fit the budget).
#     5. Repeat the whole thing `N` times with independent random orders
#        ("random restarts") and keep whichever restart produced the
#        smallest total rule set. There's no guarantee any given restart
#        finds even a local optimum -- but with a good-enough shuffle across
#        `N` tries, results tend to land close to one.
#
#   Why re-feeding an already-minimized term back into the minimizer is
#   legal: a DNF term omitting the literal for some feature `f` IS exactly
#   what a don't-care on `f` means in cube / PLA form. So a previously
#   minimized term is already a perfectly valid "cube" and can be handed
#   back to the minimizer, mixed with brand-new raw cubes, in the next call.
#   This is what lets `run_minimization` itself double as the "compaction"
#   step -- no separate compaction routine is needed.
#
#   A note on minimization backends: this scheme is well-suited to Espresso
#   (`:mitespresso`), which is fine being invoked repeatedly on a formula
#   that keeps growing/shrinking incrementally (same results as calling it
#   once on the full table, just marginally slower). ABC (`:abc`) is built
#   around being handed the *complete* PLA once and doing global synthesis
#   on it, so repeatedly re-feeding it "old minimized formula + a new
#   window of rows" is off-label for it, even though the code path is
#   shared. In practice: use `:mitespresso` for `lumen_sequential`.
# ============================================================================ #


# ---------------------------------------------------------------------------- #
#                     LAZY RANDOM PERMUTATION (no materialization)             #
# ---------------------------------------------------------------------------- #
"""
    LazyPermutation

Bijective pseudo-random permutation of `{0, ..., n-1}` computed on the fly,
without ever materializing an array of size `n`.

Built on a small balanced Feistel network over the smallest power-of-two
domain `M ≥ n` (NOTE: this `M` is a Feistel-domain bound, unrelated to the
memory-budget `M` used elsewhere in this file), plus cycle-walking to fold
`[n, M)` back into `[0, n)`. This gives O(1) amortized cost per lookup
regardless of how large `n` is (the whole point: `n` here is
`prod(length.(thrs_with_p))`, which can be astronomically large and must
never be enumerated into memory).

Different `seed` values give different permutations -- this is what powers
the "random restart" loop: each restart just picks a fresh seed.

# Fields
- `n`: size of the permuted domain `{0, ..., n-1}`.
- `half_bits`: number of bits in each Feistel half (`L` and `R`).
- `mask`: bitmask selecting the low `half_bits` bits (`(1 << half_bits) - 1`).
- `rounds`: number of Feistel rounds to apply per permutation call.
- `key`: seed value mixed into every round's hash, making the permutation
  depend on `seed`.

# Notes
This is NOT cryptographically secure and does not need to be: we only need
"shuffled enough" coverage across the (possibly huge) combination space, not
an adversarial guarantee.
"""
struct LazyPermutation
    n::Int
    half_bits::Int
    mask::UInt64
    rounds::Int
    key::UInt64
end

"""
    LazyPermutation(n; seed=rand(UInt64), rounds=4)

Construct a `LazyPermutation` over `{0, ..., n-1}`.

Chooses the smallest EVEN bit-width `total_bits` with `2^total_bits >= n`
(rounded up from `log2(n)`, bumped to even so it splits cleanly into two
equal Feistel halves), then splits it into two `half_bits`-wide halves.
`seed` becomes the Feistel network's key, so different seeds yield
different (but equally well-mixed) permutations of the same domain.
"""
function LazyPermutation(n::Integer; seed::Integer=rand(UInt64), rounds::Int=4)
    n <= 0 && throw(ArgumentError("n must be positive, got $n"))
    total_bits = max(2, ceil(Int, log2(n)))
    total_bits = iseven(total_bits) ? total_bits : total_bits + 1
    half_bits = total_bits ÷ 2
    mask = (UInt64(1) << half_bits) - UInt64(1)
    return LazyPermutation(Int(n), half_bits, mask, rounds, UInt64(seed))
end

"""
    _feistel_round(l, r, round, p) -> (newl, newr)

One round of a balanced Feistel network: the new left half becomes the old
right half, and the new right half is the old left half XORed with a hash
of `(r, round, key)`, masked down to `half_bits` bits. Standard Feistel
construction -- reversibility (needed for bijectivity) doesn't actually
matter here since we only ever go in the forward direction and rely on
cycle-walking (`permute`) rather than needing to invert.
"""
@inline function _feistel_round(l::UInt64, r::UInt64, round::Int, p::LazyPermutation)
    h = (hash((r, round, p.key)) & p.mask)
    newl = r
    newr = (l ⊻ h) & p.mask
    return newl, newr
end

"""
    _feistel_permute(x, p) -> UInt64

Split `x` into its two `half_bits`-wide halves, run `p.rounds` Feistel
rounds over them, and recombine into a single value in `[0, 2^total_bits)`.
This alone is a bijection on `[0, 2^total_bits)` but NOT yet restricted to
`[0, p.n)` -- that restriction is what `permute`'s cycle-walking loop does.
"""
@inline function _feistel_permute(x::UInt64, p::LazyPermutation)
    l = (x >> p.half_bits) & p.mask
    r = x & p.mask
    for round in 1:p.rounds
        l, r = _feistel_round(l, r, round, p)
    end
    return (l << p.half_bits) | r
end

"""
    permute(p::LazyPermutation, i::Integer) -> Int

Return the `i`-th element (0-based `i`) of the pseudo-random permutation of
`{0, ..., p.n-1}` described by `p`.

Implemented via "cycle-walking": repeatedly apply the Feistel bijection
(which is only guaranteed bijective on the full power-of-two domain, not
on `[0, p.n)`) until landing back inside `[0, p.n)`. Because the Feistel
map is itself a bijection on the power-of-two domain, iterating it always
returns to `[0, p.n)` in a bounded number of steps -- in practice O(1)
amortized, since the power-of-two domain is at most ~2x the size of `p.n`.
"""
function permute(p::LazyPermutation, i::Integer)
    (0 <= i < p.n) || throw(BoundsError(p, i))
    x = UInt64(i)
    while true
        x = _feistel_permute(x, p)
        x < UInt64(p.n) && return Int(x)
    end
end


# ---------------------------------------------------------------------------- #
#                    O(1) linear-index -> combo tuple decoding                 #
# ---------------------------------------------------------------------------- #
"""
    _strides(lens::Vector{Int}) -> Vector{Int}

Compute mixed-radix strides for a product space with per-dimension sizes
`lens`, with the FIRST dimension varying fastest (consistent with
`_ProductColumn` in Lumen.jl's batch pipeline, so decoded rows line up the
same way the batch code would have produced them).

`strides[1] == 1` always; each subsequent stride is the product of all
previous dimension sizes, exactly like row-major/column-major indexing
arithmetic for an n-dimensional array flattened to 1D.
"""
function _strides(lens::Vector{Int})
    n = length(lens)
    strides = ones(Int, n)
    @inbounds for j in 2:n
        strides[j] = strides[j-1] * lens[j-1]
    end
    return strides
end

"""
    _combo_at(thrs_with_p, lens, strides, linidx0) -> Tuple

Decode a single 0-based linear index `linidx0` into the corresponding tuple
of threshold values (one per feature), without ever materializing the full
product space.

For each feature `j`, extracts that feature's mixed-radix digit via
integer division by `strides[j]` and modulo `lens[j]`, then converts the
resulting 0-based digit to a 1-based index into `thrs_with_p[j]`. This is
the "un-flatten a single flat index" counterpart to `_strides`.
"""
@inline function _combo_at(
    thrs_with_p::Vector{<:AbstractVector},
    lens::Vector{Int},
    strides::Vector{Int},
    linidx0::Int
)
    return ntuple(length(thrs_with_p)) do j
        idx = (linidx0 ÷ strides[j]) % lens[j] + 1
        @inbounds thrs_with_p[j][idx]
    end
end


# ---------------------------------------------------------------------------- #
#                       shared threshold/family derivation                     #
# ---------------------------------------------------------------------------- #
"""
    _prepare_sequential_context(config::LumenConfig, model::SM.AbstractModel)

Standalone replica of steps 1-7 of `ExtractRulesData`'s inner constructor,
stopping *before* the full Cartesian-product / model-apply step (step 8-9),
which is exactly the part this file refuses to materialize.

Kept intentionally separate from `ExtractRulesData` (rather than refactored
into a shared helper) to avoid touching the existing, working batch
pipeline.

# Steps performed
1. Extract atoms from the model (either BFS-limited by `depth`, or the full
   alphabet if `depth >= 1.0`).
2. Validate that only supported comparison operators (`<`, `≥`, `>`, `≤`)
   appear.
3. Resolve feature names from the model's stored `:featurenames` info.
4. Resolve class names from the model's stored `:supporting_labels` info.
5. For each feature, determine its operator "family" (`:lt`-like vs
   `:ge`-like) and collect + sort its thresholds accordingly.
6. Append the "plus-boundary" sentinel value per feature (`_thrs_with_boundary`),
   producing `thrs_with_p`, whose per-feature lengths define the full
   product-space shape.
7. Precompute `lens`, `strides`, and `n_total = prod(lens)` -- everything a
   caller needs to decode arbitrary linear indices later, in O(1), via
   `_combo_at` (this is `_prepare_sequential_context`'s entire reason for
   existing: enough shared setup to make random-access decoding possible,
   without doing the expensive step 8-9 Cartesian materialization at all).

# Returns
A `NamedTuple` with:
- `thresholds::Vector{Vector{<:Float}}`   -- per-feature thresholds, no boundary
- `thrs_with_p::Vector{Vector{<:Float}}`  -- per-feature thresholds + boundary
- `featurenames::Vector{Symbol}`
- `classnames::AbstractVector`
- `op_families::Vector{Symbol}`
- `lens::Vector{Int}`, `strides::Vector{Int}`, `n_total::Int`
"""
function _prepare_sequential_context(config::LumenConfig, model::SM.AbstractModel)
    depth = get_depth(config)

    atoms = unique!(_normalize_atom.(if depth < 1.0
        # Same fix as in main.jl's `ExtractRulesData`: `init` must be
        # concretely typed `Atom{ScalarCondition}[]`. `_extract_atoms_bfs_order`
        # already returns a concretely-typed vector (fixed in main.jl), but an
        # abstractly-typed `init` here would still widen the vcat'd result via
        # typejoin back to `Vector{Atom{AbstractCondition}}`, breaking the
        # `_normalize_atom.(...)` broadcast right below (invariance of
        # parametric types -- this has to be fixed at every accumulator, not
        # just at the source).
        mapreduce(vcat, SM.models(model); init=SL.Atom{SD.ScalarCondition}[]) do t
            _take_first_percentage(_extract_atoms_bfs_order(t), depth)
        end
    else
        SL.atoms(SM.alphabet(model, false))
    end))

    # Fail fast if the model contains any comparison operator we don't know
    # how to turn into a threshold-cube decision (e.g. `==`, set-membership,
    # etc.) -- better to error here than silently mis-handle it downstream.
    let unsupported = unique(op for op in get_operator.(atoms) if op ∉ _supported_operators)
        isempty(unsupported) || throw(ArgumentError(
            "Only '<', '≥', '>', '≤' operators are currently supported. " *
            "Found unsupported operators: $(unsupported)."
        ))
    end

    features = SM.featurename.(unique!(get_feature.(atoms)))
    featurenames = SM.info(model, :featurenames)
    classnames = unique!(SM.info(model, :supporting_labels))

    type = get_float_type(config)
    thresholds = Vector{Vector{type}}(undef, length(featurenames))
    op_families = Vector{Symbol}(undef, length(featurenames))

    # Per feature: figure out which operator "family" it uses (this affects
    # sort direction, since `<`-style and `≥`-style thresholds are ordered
    # oppositely for the boundary-cube construction downstream), and collect
    # + sort its thresholds. Features that never appear in any atom get an
    # empty threshold list and default to the `:lt` family.
    @inbounds for i in eachindex(featurenames)
        idx = findfirst(f -> f == featurenames[i], features)
        if isnothing(idx)
            thresholds[i] = type[]
            op_families[i] = :lt
        else
            family = _feature_op_family(atoms, features[idx])
            op_families[i] = family
            thresholds[i] = sort!(
                get_threshold.(_atoms_for_feature(atoms, features[idx]));
                rev=(family === :lt)
            )
        end
    end

    thrs_with_p = _thrs_with_boundary(thresholds, op_families)
    lens = length.(thrs_with_p)
    strides = _strides(lens)
    n_total = prod(lens)

    return (;
        thresholds, thrs_with_p, featurenames=Symbol.(featurenames),
        classnames, op_families, lens, strides, n_total
    )
end


# ---------------------------------------------------------------------------- #
#                             per-class raw buffer                             #
# ---------------------------------------------------------------------------- #
"""
    _SeqBuffer

Per-class bookkeeping used during a single restart pass.

# Fields
- `raw::Vector{Vector{SL.Atom}}` -- freshly generated, not-yet-folded-in
  cubes for this class. Each element is one row's cube (a conjunction of
  atoms/literals). Grows as rows are routed here; emptied every time we
  fold.
- `terms::Vector{Any}` -- the class' CURRENT best minimized DNF ("φ" in the
  design discussion). This is REPLACED (not appended to) every time we
  fold: there is no separate "raw results pool that keeps growing" anymore
  -- the running formula itself is the only thing carried forward, which is
  exactly what keeps total memory bounded by `M` instead of growing
  unboundedly across many folds.
- `n_rows_seen::Int` -- total rows ever routed to this class across the
  whole pass (diagnostic only, not used for any control-flow decision).
"""
mutable struct _SeqBuffer
    raw::Vector{Vector{SL.Atom}}
    terms::Vector{Any}
    n_rows_seen::Int
end

"""
    _SeqBuffer()

Construct an empty buffer: no raw cubes collected yet, no terms in the
running formula yet, zero rows seen.
"""
_SeqBuffer() = _SeqBuffer(Vector{Vector{SL.Atom}}(), Any[], 0)


# ---------------------------------------------------------------------------- #
#                      generate-and-minimize fold-in logic                     #
# ---------------------------------------------------------------------------- #
"""
    _term_to_cube(term) -> Vector{SL.Atom}

Turn a previously-minimized DNF term back into the plain "cube" shape (a
vector of atoms/literals) that `run_minimization` accepts as raw input.

This conversion is what makes folding old φ back into the minimizer legal:
a minimized term that had some literals dropped by a previous Espresso/ABC
pass IS ALREADY a valid (partial) cube -- omitting the literal for feature
`f` is exactly what a don't-care on `f` means in PLA/cube form. So handing
an old term back in as "just another row" loses no information.

Two methods:
- If `term` is already a bare `Vector{<:SL.Atom}` (e.g. it slipped through
  unminimized, or the minimizer's own output happens to already be in this
  shape), it's already a cube -- return as-is.
- Otherwise (the general case: `term` is some `SyntaxStructure`, e.g. a
  `LeftmostConjunctiveForm` or similar, as `run_minimization` typically
  returns), extract its atoms via `SL.atoms(term)` -- mirroring how
  `SL.atoms(...)` is already used elsewhere in this file to pull atoms out
  of a syntax object (see `_prepare_sequential_context`) -- and collect them
  into a concretely-typed `Vector{SL.Atom}`.
"""
_term_to_cube(term::Vector{<:SL.Atom}) = term
_term_to_cube(term) = collect(SL.Atom, SL.atoms(term))

"""
    _fold_in!(buf::_SeqBuffer, scheme::Symbol, config::LumenConfig, M::Int) -> Bool

Generate-and-minimize FUSED into a single step for one class' buffer.

If `buf.raw` is non-empty:
1. Build `combined = (buf.terms converted back to cubes) ∪ buf.raw` -- i.e.
   the class' current formula PLUS everything collected since the last
   fold, all expressed as cubes.
2. Run the minimizer (`run_minimization`, dispatching on `scheme`, so this
   is the exact same code path used for both `:abc` and `:mitespresso` --
   no backend-specific logic needed here) over `combined`.
3. The result REPLACES `buf.terms` outright (not appended -- this is the
   key difference from a naive "just keep appending minimized chunks"
   approach, and is what keeps the running formula from growing without
   bound across many folds over a huge combination space).
4. Empty `buf.raw` -- it has now been fully absorbed into `buf.terms`.

If `buf.raw` is empty, nothing is done except re-checking the existing
size (this happens for the harmless "trivial" call, e.g. when a class
received zero rows since the last fold, or zero rows overall).

# Returns
`true` if `length(buf.terms) < M` after folding (this class' memory budget
is satisfied), `false` otherwise. The caller (`_sequential_pass` /
`lumen_sequential`) treats any `false` as "this restart should be
discarded" -- there is no further attempted repair, matching the design
goal of "provaci di nuovo con un altro ordine random" rather than trying to
rescue a bad ordering mid-pass.
"""
function _fold_in!(buf::_SeqBuffer, scheme::Symbol, config::LumenConfig, M::Int)
    if !isempty(buf.raw)
        # Old formula (converted back to cube form) + everything new since
        # the last fold, fed to the minimizer as ONE combined table.
        #
        # NOTE: must convert to Vector{Vector{SL.Atom}}, NOT leave it as
        # whatever `vcat` infers. Julia's parametric types are invariant, so
        # `vcat(_term_to_cube.(buf.terms), buf.raw)` ends up a `Vector{Any}`
        # the moment `buf.terms` is empty (its first fold) or broadcast
        # can't narrow the element type away from `buf.terms`'s own
        # `Vector{Any}` container type -- even though every element already
        # IS a concrete `Vector{SL.Atom}` at runtime. This is the exact
        # invariance pitfall already flagged elsewhere in this file (see the
        # original `_flush!`'s note about `Vector{SyntaxStructure}`), and
        # it's what makes `run_minimization` -- which only has methods for
        # `Vector{Vector{Atom}}`, not `Vector{Any}` -- fail with a
        # MethodError. `convert` here does an element-wise copy-conversion
        # (each element is already the right type, so it's just a container
        # re-type), fixing the eltype back to what `run_minimization` expects.
        combined = convert(
            Vector{Vector{SL.Atom}},
            vcat(_term_to_cube.(buf.terms), buf.raw)
        )
        minimized = run_minimization(Val(scheme), config, combined)

        # REPLACE, don't append: buf.terms now IS the new (hopefully
        # smaller) formula, not "old formula + new minimized chunk".
        buf.terms = collect(Any, minimized)
        empty!(buf.raw)
    end
    return length(buf.terms) < M
end


# ---------------------------------------------------------------------------- #
#                        one sequential pass (one restart)                     #
# ---------------------------------------------------------------------------- #
"""
    _sequential_pass(config, model, ctx, scheme, seed; M, max_apply_batch)

Run exactly ONE restart: enumerate the whole combination space once, in the
pseudo-random order given by `LazyPermutation(ctx.n_total; seed)`, routing
each row's derived cube into the buffer of its predicted class, and folding
(generate-and-minimize fused, see `_fold_in!`) a class buffer the instant
it reaches `M` raw cubes.

# Batch sizing (decode/apply efficiency, decoupled from the memory budget)
Decoding one row and calling `apply` on it individually would be far too
slow, so rows are still decoded and passed to the model in batches. BUT
unlike a fixed-size chunking scheme, the batch size after the first one is
computed DYNAMICALLY every iteration as:

    room = M - (size of the currently fullest class' raw buffer)
    this_chunk = clamp(min(room, max_apply_batch, rows_left), 1, rows_left)

i.e. "decode as many rows as still comfortably fit before SOME class buffer
would need folding, but never more than `max_apply_batch` (a pure
performance cap, unrelated to when folding actually happens) and never more
than what's left to enumerate". This is what avoids "perforare sempre della
stessa misura" (always poking a fixed-size hole): the batch shrinks as
buffers fill up and grows again right after a fold frees up room.

Note: `room` is only an upper bound on how much can be decoded before
POSSIBLY needing a fold -- the actual per-row check
`length(buf.raw) >= M` inside the inner loop is still what triggers
`_fold_in!`, exactly as before. The dynamic batch sizing is purely about
not over-decoding relative to `M`; it never changes WHEN a fold happens,
only how many rows get decoded/applied together before that point.

# Arguments
- `config`, `model`: as in `lumen_sequential`.
- `ctx`: the `NamedTuple` returned by `_prepare_sequential_context`.
- `scheme`: minimization backend (`:abc` or `:mitespresso`), forwarded to
  `run_minimization` via `_fold_in!`.
- `seed`: seed for this restart's `LazyPermutation` -- different restarts
  pass different seeds here to get genuinely different enumeration orders.
- `M`: memory budget (max raw cubes per class buffer before folding).
- `max_apply_batch`: performance cap on decode/apply batch size (see
  above) -- NOT a memory-budget parameter.

# Returns
`(per_class_terms, ok)`:
- `per_class_terms[c]` is the final φ (`buf.terms`) for `ctx.classnames[c]`
  after the final fold.
- `ok` is `false` if ANY class ever failed to get back under `M` even after
  folding (i.e. this whole restart should be discarded by the caller),
  `true` otherwise.
"""
function _sequential_pass(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    scheme::Symbol,
    seed::Integer;
    M::Int,
    max_apply_batch::Int
)
    nclasses = length(ctx.classnames)
    # Map each class label to its buffer index, so predictions (which come
    # back as class labels from `apply`) can be routed to the right buffer.
    class_index = Dict(c => i for (i, c) in enumerate(ctx.classnames))
    buffers = [_SeqBuffer() for _ in 1:nclasses]

    perm = LazyPermutation(ctx.n_total; seed=seed)
    ok = true

    i0 = 0
    while i0 < ctx.n_total
        # --- Dynamic chunk sizing (decode/apply efficiency only) ---------
        # How full is the fullest class buffer right now? `room` is how
        # many more rows could be decoded before SOME buffer would reach M
        # (a per-row check inside the loop below still does the actual
        # folding trigger -- this is only about not over-decoding).
        room = M - maximum(length(b.raw) for b in buffers; init=0)
        this_chunk = clamp(min(room, max_apply_batch, ctx.n_total - i0), 1, ctx.n_total - i0)

        # --- Decode `this_chunk` pseudo-random rows, O(1) each ------------
        # `permute(perm, i0+k-1)` gives the (0-based) linear index of the
        # (i0+k)-th row in this restart's random order; `_combo_at` turns
        # that linear index into the actual per-feature threshold tuple,
        # without ever materializing more than `this_chunk` rows at once.
        rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
        @inbounds for k in 1:this_chunk
            linidx = permute(perm, i0 + k - 1)  # 0-based
            rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, linidx)
        end
        i0 += this_chunk

        # --- Batch-apply the model once per chunk --------------------------
        # Much cheaper than applying the model row-by-row: build a small
        # column-table for just this chunk and run the model's batch-apply
        # function on it once.
        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        # --- Route each row's cube into its predicted class's buffer -------
        @inbounds for k in 1:this_chunk
            ci = get(class_index, preds[k], nothing)
            isnothing(ci) && continue  # defensive: unseen label, skip

            # Turn this row's threshold tuple into the boolean truth values
            # against the (boundary-free) per-feature thresholds, then into
            # an actual cube (conjunction of atoms) via `generate_disjunct`.
            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            buf = buffers[ci]
            push!(buf.raw, cube)
            buf.n_rows_seen += 1

            # This class' buffer just hit the memory budget -- fold NOW
            # (generate-and-minimize fused), before routing any more rows.
            if length(buf.raw) >= M
                ok &= _fold_in!(buf, scheme, config, M)
            end
        end
    end

    # Whole combination space enumerated; fold whatever is left in every
    # class' buffer (there may be leftover raw rows below M that were never
    # triggered by the `>= M` check above).
    for buf in buffers
        ok &= _fold_in!(buf, scheme, config, M)
    end

    per_class_terms = [buf.terms for buf in buffers]
    return per_class_terms, ok
end


# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    lumen_sequential(config::LumenConfig, model::SM.AbstractModel;
                      M::Int=100_000,
                      N::Int=10,
                      max_apply_batch::Int=min(M, 4096),
                      score::Symbol=:n_terms) -> SM.DecisionSet

Generate-and-minimize FUSED into a single online pass over the (possibly
astronomically large) combination space, run `N` times with independent
pseudo-random enumeration orders ("random restarts"), keeping whichever
restart produced the smallest total rule set.

# How it works, end to end
1. `_prepare_sequential_context` derives everything needed to decode
   arbitrary linear indices into threshold-combination rows in O(1)
   (thresholds, per-feature lengths/strides, total combination count),
   WITHOUT materializing the product space -- this is shared setup, done
   once, independent of `M`/`N`/restarts.
2. For each of `N` restarts:
   a. Pick a fresh random seed -> a fresh `LazyPermutation` order over the
      whole combination space.
   b. Run `_sequential_pass`: enumerate in that order, routing each row's
      cube into its predicted class' buffer, and the moment any class
      buffer reaches `M` raw cubes, immediately minimize (fold in) --
      over that class' CURRENT formula plus everything collected since the
      last fold -- replacing the formula with the (usually smaller) result.
      This is exactly what makes it "generate-and-minimize fused" rather
      than "generate everything, then minimize": the minimizer never sees
      more than `M` cubes worth of "new" material at once, and the running
      formula itself never exceeds `M` terms either.
   c. If any class' formula still can't fit under `M` even after folding,
      the whole restart is discarded (`ok = false`) and doesn't count as a
      candidate.
   d. Otherwise, score the restart (currently: total term count across all
      classes) and keep it if it's the best seen so far.
3. If every single restart got discarded, raise an informative error
   (suggesting to raise `M` and/or `N`) rather than silently returning
   nothing.
4. Build the final `SM.DecisionSet` from the best restart's per-class
   terms, dropping any class that ended up with zero terms.

# Arguments
- `config::LumenConfig`: same config used by `lumen`; `minimization_scheme`
  selects `:abc` or `:mitespresso`. See the module-level notes at the top
  of this file for why `:mitespresso` is the scheme this entry point is
  actually designed around (repeated incremental re-minimization is
  off-label for `:abc`, which wants the whole PLA once).
- `model::SM.AbstractModel`: the model to extract rules from.
- `M::Int`: the ONE memory-budget knob. Maximum number of raw cubes any
  single class' buffer may accumulate before it MUST be folded
  (generate-and-minimize fused) back down. Replaces what used to be TWO
  separate knobs (a flush threshold and a results-growth threshold),
  because folding now replaces the running formula in place instead of
  letting a separate "results pool" grow across many flushes.
- `N::Int`: number of independent random-order restarts to try; the
  smallest-total-terms restart (among those that stayed under `M`) is kept.
- `max_apply_batch::Int`: PURELY a performance cap on how many rows are
  decoded and `apply`'d together in one batch inside `_sequential_pass`.
  It has NO effect on when folding happens -- that is governed only by
  `M`. It exists only so a single `apply` call doesn't try to materialize
  an unreasonably large number of rows at once when `M` itself is huge.
  Defaults to `min(M, 4096)`.
- `score::Symbol`: currently only `:n_terms` (total term count summed
  across all classes) is implemented. A literal-count-based scoring scheme
  would need a reliable, cross-representation "count atoms in a term"
  accessor that isn't assumed to exist here, since minimized terms can come
  back as arbitrary `SyntaxStructure`, not necessarily flat conjunctions.

# Returns
`SM.DecisionSet` built from the best restart's per-class terms, in the same
shape as `lumen`'s (batch pipeline's) output.

# Throws
- `ArgumentError` if `score` is anything other than `:n_terms`, or if `M`
  or `N` is not positive.
- An `error` (not `ArgumentError`, since this is a runtime/parameter-tuning
  failure rather than a malformed call) if all `N` restarts were discarded
  for exceeding the memory budget even after folding.
"""
function lumen_sequential(
    config::LumenConfig,
    model::SM.AbstractModel;
    M::Int=100_000,
    N::Int=10,
    max_apply_batch::Int=min(M, 4096),
    score::Symbol=:n_terms
)
    score === :n_terms || throw(ArgumentError(
        "Only score=:n_terms is currently supported."
    ))
    M >= 1 || throw(ArgumentError("M must be positive, got $M"))
    N >= 1 || throw(ArgumentError("N must be positive, got $N"))

    scheme = get_minimization_scheme(config)
    # Shared, restart-independent setup: thresholds, feature/class names,
    # and everything needed to O(1)-decode arbitrary linear indices later.
    ctx = _prepare_sequential_context(config, model)

    best_terms = nothing
    best_score = typemax(Int)

    for r in 1:N
        # Fresh seed per restart -> fresh LazyPermutation order. This is the
        # entirety of what makes each restart "independent": same ctx, same
        # M/scheme, different random enumeration order.
        seed = rand(UInt64)
        per_class_terms, ok = _sequential_pass(
            config, model, ctx, scheme, seed;
            M, max_apply_batch
        )
        # Restart blew the memory budget in at least one class even after
        # folding -- discard it, it's not a valid candidate.
        ok || continue

        this_score = sum(length, per_class_terms)
        if this_score < best_score
            best_score = this_score
            best_terms = per_class_terms
        end
    end

    isnothing(best_terms) && error(
        "lumen_sequential: all $N restarts exceeded the memory budget " *
        "(M=$M) even after folding-in. Raise M, or increase N to give " *
        "more random orderings a chance to fold down under budget."
    )

    # Drop classes that ended up with zero terms (i.e. the model never
    # predicted that class for any enumerated row in the winning restart).
    valid_mask = .!isempty.(best_terms)
    formulas = best_terms[valid_mask]
    classes = ctx.classnames[valid_mask]

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end