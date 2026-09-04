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
#     3. The moment a class' TOTAL occupancy (its CURRENT minimized formula
#        `terms` PLUS everything collected since the last fold, `raw`)
#        reaches the memory budget `M`, minimize immediately: call the
#        minimizer over `terms ∪ raw`. The result REPLACES the current
#        formula. Because minimization usually collapses many raw rows into
#        fewer don't-care terms, this typically pulls the buffer back well
#        under `M`, so enumeration continues; if `M` is hit again, we fold
#        again.
#        NOTE (patch): the trigger is on `length(terms) + length(raw) >= M`,
#        NOT on `length(raw) >= M` alone. This guarantees the minimizer
#        never receives more than `M` cubes in one call (previously,
#        checking `raw` alone let `combined = terms + raw` grow up to
#        `2M - 1` before folding).
#     4. If a class' formula still can't be gotten under `M` even after
#        folding, the whole restart is abandoned (too unlucky an order to
#        fit the budget).
#     5. Repeat the whole thing `N` times with independent random orders
#        ("random restarts") and keep whichever restart produced the
#        smallest total rule set. There's no guarantee any given restart
#        finds even a local optimum -- but with a good-enough shuffle across
#        `N` tries, results tend to land close to one.
#     6. NEW: once the winning restart is picked, run ONE full-domain,
#        linear (non-random) verification sweep (`_verify_and_repair_sensitivity`)
#        that checks every row against its OWN predicted class' current φ,
#        and repairs any class whose φ fails to cover rows the model
#        assigns to it (false negatives -- i.e. lost sensitivity), caused
#        by earlier folds generalizing before they had seen the whole
#        domain. See the dedicated section below for the full rationale.
#
#   Why re-feeding an already-minimized term back into the minimizer is
#   legal: a DNF term omitting the literal for some feature `f` IS exactly
#   what a don't-care on `f` means in cube / PLA form. So a previously
#   minimized term is already a perfectly valid "cube" and can be handed
#   back to the minimizer, mixed with brand-new raw cubes, in the next call.
#   This is what lets `run_minimization` itself double as the "compaction"
#   step -- no separate compaction routine is needed. It is also what makes
#   the repair pass' targeted re-folds legal for exactly the same reason.
#
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
function _prepare_sequential_context(
    config::LumenConfig{T},
    trees::Vector{SM.Branch{S}}
) where {T<:AbstractFloat,S<:SM.Label}
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
        mapreduce(vcat, trees; init=SL.Atom{SD.ScalarCondition}[]) do t
            _take_first_percentage(_extract_atoms_bfs_order(t), depth)
        end
    else
        collect(Iterators.flatten(SM.alphabet.(trees, false)))
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

"""
    _occupancy(buf::_SeqBuffer) -> Int

TOTAL occupancy of a class buffer: `terms` (current minimized formula) PLUS
`raw` (not-yet-folded cubes collected since the last fold).

This is the quantity that is now compared against `M` for deciding when to
fold (PATCH: previously only `length(buf.raw)` was checked, which let
`combined = terms + raw`, the actual input handed to the minimizer inside
`_fold_in!`, grow up to `2M - 1` before a fold was triggered). Centralizing
it here also keeps `_sequential_pass`'s trigger check and its `room`
computation for dynamic batch sizing consistent with each other.
"""
@inline _occupancy(buf::_SeqBuffer) = length(buf.terms) + length(buf.raw)


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
   PATCH: thanks to the new trigger in `_sequential_pass` (which now fires
   on `_occupancy(buf) >= M`, i.e. `length(terms) + length(raw) >= M`,
   instead of `length(raw) >= M` alone), `combined` is now guaranteed to
   have at most `M` elements at the moment it reaches this function --
   previously it could reach up to `2M - 1`.
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
its TOTAL occupancy (`terms` + `raw`, see `_occupancy`) reaches `M`.

# PATCH: trigger now on total occupancy, not on `raw` alone
Previously the fold trigger checked `length(buf.raw) >= M`, ignoring
`buf.terms`. Since `_fold_in!` builds `combined = terms ∪ raw` and hands
that whole thing to the minimizer, that meant `combined` could reach up to
`2M - 1` elements before folding -- i.e. the minimizer could be asked to
process up to roughly double the intended budget.

Now the trigger is `_occupancy(buf) >= M`, i.e.
`length(buf.terms) + length(buf.raw) >= M`. This guarantees `combined`
never exceeds `M` elements at the moment `_fold_in!` is invoked, matching
the "you stop and call the minimizer once you get close to your memory
limit" behavior (checking the TOTAL table size, not just the newly
generated part).

# Batch sizing (decode/apply efficiency, decoupled from the memory budget)
Decoding one row and calling `apply` on it individually would be far too
slow, so rows are still decoded and passed to the model in batches. BUT
unlike a fixed-size chunking scheme, the batch size after the first one is
computed DYNAMICALLY every iteration as:

    room = M - (TOTAL occupancy of the currently fullest class buffer)
    this_chunk = clamp(min(room, max_apply_batch, rows_left), 1, rows_left)

i.e. "decode as many rows as still comfortably fit before SOME class
buffer's total occupancy would need folding, but never more than
`max_apply_batch` (a pure performance cap, unrelated to when folding
actually happens) and never more than what's left to enumerate". This is
what avoids "perforare sempre della stessa misura" (always poking a
fixed-size hole): the batch shrinks as buffers fill up and grows again
right after a fold frees up room.

Note: `room` is only an upper bound on how much can be decoded before
POSSIBLY needing a fold -- the actual per-row check
`_occupancy(buf) >= M` inside the inner loop is still what triggers
`_fold_in!`. The dynamic batch sizing is purely about not over-decoding
relative to `M`; it never changes WHEN a fold happens, only how many rows
get decoded/applied together before that point.

# Arguments
- `config`, `model`: as in `lumen_sequential`.
- `ctx`: the `NamedTuple` returned by `_prepare_sequential_context`.
- `scheme`: minimization backend (`:abc` or `:mitespresso`), forwarded to
  `run_minimization` via `_fold_in!`.
- `seed`: seed for this restart's `LazyPermutation` -- different restarts
  pass different seeds here to get genuinely different enumeration orders.
- `M`: memory budget (max TOTAL occupancy -- `terms` + `raw` -- per class
  buffer before folding).
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
        # How full is the fullest class buffer right now, in TOTAL
        # occupancy (terms + raw)? `room` is how many more rows could be
        # decoded before SOME buffer would reach M (a per-row check inside
        # the loop below still does the actual folding trigger -- this is
        # only about not over-decoding).
        room = M - maximum(_occupancy(b) for b in buffers; init=0)
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

            # This class' TOTAL occupancy (terms + raw) just hit the memory
            # budget -- fold NOW (generate-and-minimize fused), before
            # routing any more rows. PATCH: was `length(buf.raw) >= M`.
            if _occupancy(buf) >= M
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
#                    sensitivity verification / repair pass                    #
# ---------------------------------------------------------------------------- #
#
#   WHY THIS EXISTS
#   ----------------
#   Each intermediate `_fold_in!` during `_sequential_pass` minimizes a
#   class' buffer having seen only the rows enumerated SO FAR for that
#   class. At that moment the minimizer treats "not yet seen" as "don't
#   care" -- but part of that unseen space will later turn out to belong to
#   OTHER classes. A term folded early can therefore end up more general
#   than it should be, and once frozen into `buf.terms` it is only ever
#   re-minimized together with new raw cubes, never re-checked against the
#   rest of the domain. The over-generalized term can end up overlapping a
#   region that truly belongs to another class -- and when the final
#   `DecisionSet` is evaluated, THAT other class loses sensitivity (its
#   instances get contended/misclassified by the too-broad rule).
#
#   This pass runs ONCE, after the winning restart has been picked, over
#   the FULL combination space (linear order this time -- no permutation
#   needed, we want total coverage, not a random sample), and asks, for
#   every row: "does this row's OWN predicted class' current φ actually
#   cover it?" Any row where the answer is no is a false negative for that
#   class -- exactly what erodes sensitivity -- and gets fed back into a
#   targeted repair fold for that class only.
#
#   SCOPE: this only targets false negatives (a class' φ failing to cover
#   rows the model assigns to it), which is what "losing sensitivity"
#   means. It does NOT attempt to fix precision/specificity issues (a
#   class' φ wrongly covering rows that belong to another class) -- that is
#   a different failure mode and was not what was asked for here.
#
#   COST: one extra O(n_total) sweep with batched model-apply, same order
#   of cost as a single restart's enumeration (actually cheaper, since it
#   uses a plain linear index instead of `LazyPermutation`'s cycle-walking).
#   Against `N` restarts already paid for, this is a ~1/N overhead -- and it
#   is still fully batched/streamed, so it does not reintroduce the
#   "materialize the whole product space" problem this file exists to
#   avoid.
#

"""
    _cube_to_formula(cube::Vector{<:SL.Atom}) -> SyntaxStructure

Wrap a bare cube (vector of atoms/literals) into the conjunctive
`SyntaxStructure` shape needed for `SL.check`. This is the same shape
`run_minimization` normally hands back for a single term; used here only
to check individual raw/repair cubes without invoking the minimizer.
"""
_cube_to_formula(cube::Vector{<:SL.Atom}) = SL.LeftmostConjunctiveForm(cube)

"""
    _class_formula(terms::Vector{Any}, float_type::Type) -> Union{Nothing,SyntaxStructure}

Build the disjunctive formula for one class' current φ (`buf.terms`, a
mix of bare cubes and minimizer-output elements -- see `_fold_in!`) so it
can be checked against a logiset with `SL.check`, exactly the way the
final `DecisionSet` will be checked downstream.

Returns `nothing` for a class with zero terms (no row was ever routed to
it in the winning restart) -- callers must treat "no formula" as "covers
nothing", not as an error.

# PATCH: forced wide element type (fixes a `check` dispatch ambiguity)
If `conjunctions` is left to whatever type Julia infers from the list
comprehension below, every element ends up concretely typed as e.g.
`LeftmostConjunctiveForm{Atom{float_type}}` (they all came out of the same
minimizer / `_cube_to_formula` path). `SL.LeftmostDisjunctiveForm(...)`
then infers a NARROW type parameter from those concrete runtime elements,
producing a `LeftmostDisjunctiveForm{LeftmostConjunctiveForm{Atom}}` that
matches SoleLogics' stricter `DNF` type alias.

That narrow match is exactly the trap already documented in
`lumen_sequential`'s own docstring for the FINAL `DecisionSet`: with a
`DNF`-matching formula, `check(::DefaultCheckAlgorithm, ::DNF, ::LogicalInstance, ...)`
and `check(::DefaultCheckAlgorithm, ::LeftmostDisjunctiveForm, ::LogicalInstance, ...)`
become EQUALLY specific for the same call, and SoleLogics raises a
`MethodError: ... is ambiguous` instead of picking one. Previously this was
only guarded against at the very end of `lumen_sequential` (when building
the returned `DecisionSet`) -- but the repair pass calls `SL.check`
*during* `_verify_and_repair_sensitivity`, on formulas built here, so the
same guard is needed at this call site too.

Fix: force `conjunctions` into the exact same WIDE element type
(`Union{LeftmostConjunctiveForm{Atom{float_type}}, SyntaxStructure}`) that
`lumen()` and the final `formulas` in `lumen_sequential` already use. A
`Union`-typed container can never match the narrower `DNF` alias, so
`check` unambiguously resolves to the `LeftmostDisjunctiveForm` method.
"""
function _class_formula(terms::Vector{Any}, float_type::Type)
    isempty(terms) && return nothing
    conjunctions = Vector{Union{
        SL.LeftmostConjunctiveForm{SL.Atom{float_type}},
        SyntaxStructure
    }}([t isa Vector{<:SL.Atom} ? _cube_to_formula(t) : t for t in terms])
    return SL.LeftmostDisjunctiveForm(conjunctions)
end

"""
    _repair_fold(terms, raw, scheme, config) -> Vector{Any}

Same merge-and-reminimize step as `_fold_in!` (old formula, converted back
to cube form via `_term_to_cube`, unioned with new raw cubes, then
re-minimized in a single call) -- but as a pure function returning the new
`terms` instead of mutating a `_SeqBuffer`, and with NO `M`-budget check.

The repair pass' job is correctness (sensitivity), not staying under the
memory budget: it always folds, regardless of how large the result ends
up. If a repair fold blows past what feels reasonable, that's a signal `M`
was too tight for the winning restart, not something this pass should
silently paper over.
"""
function _repair_fold(
    terms::Vector{Any},
    raw::Vector{Vector{SL.Atom}},
    scheme::Symbol,
    config::LumenConfig
)
    combined = convert(
        Vector{Vector{SL.Atom}},
        vcat(_term_to_cube.(terms), raw)
    )
    minimized = run_minimization(Val(scheme), config, combined)
    return collect(Any, minimized)
end

"""
    _verify_and_repair_sensitivity(
        config, model, ctx, scheme, per_class_terms, classnames;
        max_apply_batch, max_repair_rounds=3
    ) -> Vector{Vector{Any}}

Full-domain, linear (non-random) sweep that checks every row against its
OWN predicted class' current φ, collects false negatives per class, and
repairs each affected class with one targeted `_repair_fold`. Repeats up
to `max_repair_rounds` times (a repair fold can, in principle, still leave
-- or introduce -- a handful of uncovered rows; this is a fixed-point
iteration, not a one-shot guarantee), emitting a `@warn` if it never fully
converges.

# Arguments
- `config`, `model`, `scheme`: as elsewhere in this file.
- `ctx`: from `_prepare_sequential_context` (same one used for the winning
  restart).
- `per_class_terms::Vector{Vector{Any}}`: the winning restart's per-class
  φ, indexed the same way as `classnames` (i.e. `per_class_terms[i]`
  belongs to `classnames[i]`).
- `classnames`: `ctx.classnames`.
- `max_apply_batch::Int`: batch size for the verification sweep's
  decode/apply loop -- reuse the same value passed to `lumen_sequential`,
  purely a performance knob, unrelated to correctness.
- `max_repair_rounds::Int`: safety cap on repair iterations.

# Returns
The (possibly repaired) `per_class_terms`, same shape as the input.

# Note on `SL.check`
This assumes `SL.check(formula, d)` returns a per-instance `Bool`
(or `BitVector`)-like collection when given a whole logiset `d`. If your
SoleLogics version only exposes a per-instance signature
(`SL.check(formula, d, i)`), replace the batched
`covered[ci] = SL.check(formulas[ci], d)` line with a loop over `k`.
"""
function _verify_and_repair_sensitivity(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    scheme::Symbol,
    per_class_terms::Vector{Vector{Any}},
    classnames::AbstractVector;
    max_apply_batch::Int,
    max_repair_rounds::Int=3
)
    nclasses = length(classnames)
    class_index = Dict(c => i for (i, c) in enumerate(classnames))
    terms = deepcopy(per_class_terms)
    # Needed by `_class_formula` to force the wide, non-`DNF`-matching
    # element type -- see the dispatch-ambiguity note in its docstring.
    float_type = get_float_type(config)

    for round in 1:max_repair_rounds
        formulas = [_class_formula(terms[i], float_type) for i in 1:nclasses]
        repair_raw = [Vector{Vector{SL.Atom}}() for _ in 1:nclasses]
        n_false_negatives = 0

        i0 = 0
        while i0 < ctx.n_total
            this_chunk = clamp(min(max_apply_batch, ctx.n_total - i0), 1, ctx.n_total - i0)

            # Plain linear decoding this time (no LazyPermutation): the
            # repair pass needs to see EVERY row exactly once, not a random
            # sample, so there is no need for (and no benefit from) the
            # permutation's shuffling -- a linear index is cheaper too.
            rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
            @inbounds for k in 1:this_chunk
                linidx = i0 + k - 1  # 0-based, sequential
                rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, linidx)
            end
            i0 += this_chunk

            tbl = NamedTuple{Tuple(ctx.featurenames)}(
                ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
            )
            d = PropositionalLogiset(tbl)
            preds = get_apply_function(config)(
                model, d;
                use_multithreads=get_use_multithreads(config),
                suppress_parity_warning=true
            )

            # Precompute, once per batch, which instances each class'
            # CURRENT φ covers -- reusing the exact same `check` mechanism
            # the final DecisionSet is evaluated with downstream, so
            # "covered" here means exactly what "covered" will mean once
            # these rules ship.
            covered = Vector{Any}(undef, nclasses)
            @inbounds for ci in 1:nclasses
                covered[ci] = isnothing(formulas[ci]) ? nothing : SL.check(formulas[ci], d)
            end

            @inbounds for k in 1:this_chunk
                ci = get(class_index, preds[k], nothing)
                isnothing(ci) && continue  # defensive: unseen label, skip

                is_covered = !isnothing(covered[ci]) && covered[ci][k]
                is_covered && continue

                # False negative: the model predicts `classnames[ci]` for
                # this row, but that class' current φ doesn't cover it --
                # almost always because an earlier fold, with an
                # incomplete view of the domain, over-generalized a term
                # (belonging to this class or another one) that ended up
                # stealing the region. Feed the row back in as a fresh raw
                # cube for a targeted repair fold.
                n_false_negatives += 1
                truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
                cube = generate_disjunct(
                    truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
                )
                push!(repair_raw[ci], cube)
            end
        end

        n_false_negatives == 0 && break  # sensitivity fully recovered

        # One targeted, per-class repair fold: merge the class' current
        # (already minimized) φ with exactly the rows it was missing, and
        # re-minimize. Classes with zero false negatives this round are
        # left untouched.
        @inbounds for ci in 1:nclasses
            isempty(repair_raw[ci]) && continue
            terms[ci] = _repair_fold(terms[ci], repair_raw[ci], scheme, config)
        end

        if round == max_repair_rounds
            @warn(
                "lumen_sequential: sensitivity repair did not fully " *
                "converge after $max_repair_rounds round(s); " *
                "$n_false_negatives row(s) still uncovered by their " *
                "predicted class' φ. Consider raising max_repair_rounds."
            )
        end
    end

    return terms
end


# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    lumen_sequential(config::LumenConfig, model::SM.AbstractModel;
                      M::Int=100_000,
                      N::Int=10,
                      max_apply_batch::Int=min(M, 4096),
                      score::Symbol=:n_terms,
                      repair_sensitivity::Bool=true,
                      max_repair_rounds::Int=3) -> SM.DecisionSet

Generate-and-minimize FUSED into a single online pass over the (possibly
astronomically large) combination space, run `N` times with independent
pseudo-random enumeration orders ("random restarts"), keeping whichever
restart produced the smallest total rule set, then running one full-domain
sensitivity repair sweep on the winner.

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
      buffer's TOTAL occupancy (`terms` + `raw`) reaches `M`, immediately
      minimize (fold in) -- over that class' CURRENT formula plus
      everything collected since the last fold -- replacing the formula
      with the (usually smaller) result. This is exactly what makes it
      "generate-and-minimize fused" rather than "generate everything, then
      minimize": the minimizer never sees more than `M` cubes worth of
      material at once (PATCH: guaranteed now, was up to `2M-1` before),
      and the running formula itself never exceeds `M` terms either.
   c. If any class' formula still can't fit under `M` even after folding,
      the whole restart is discarded (`ok = false`) and doesn't count as a
      candidate.
   d. Otherwise, score the restart (currently: total term count across all
      classes) and keep it if it's the best seen so far.
3. If every single restart got discarded, raise an informative error
   (suggesting to raise `M` and/or `N`) rather than silently returning
   nothing.
4. NEW -- `_verify_and_repair_sensitivity`: on the WINNING restart's
   per-class φ only, run one full-domain linear sweep checking every row
   against its own predicted class' current φ. Any row the model assigns
   to class `C` but that `C`'s current φ fails to cover is a false
   negative (lost sensitivity) -- almost always caused by an earlier fold
   generalizing a term before it had seen the whole domain. False
   negatives are fed back into a targeted per-class repair fold
   (`_repair_fold`), and the sweep repeats (up to `max_repair_rounds`
   times) until no false negatives remain or the cap is hit (in which case
   a `@warn` is emitted). This step only fixes coverage
   (sensitivity/false-negatives); it does not address a class' φ being too
   broad and wrongly covering another class' rows (a precision/specificity
   concern, out of scope here). Can be disabled via
   `repair_sensitivity=false` if you want the raw, unrepaired restart
   output (e.g. for debugging or comparing against the old behavior).
5. Build the final `SM.DecisionSet` from the (possibly repaired) winning
   restart's per-class terms, dropping any class that ended up with zero
   terms.

# Arguments
- `config::LumenConfig`: same config used by `lumen`; `minimization_scheme`
  selects `:abc` or `:mitespresso`. See the module-level notes at the top
  of this file for why `:mitespresso` is the scheme this entry point is
  actually designed around (repeated incremental re-minimization is
  off-label for `:abc`, which wants the whole PLA once).
- `model::SM.AbstractModel`: the model to extract rules from.
- `M::Int`: the ONE memory-budget knob. Maximum TOTAL occupancy
  (`length(terms) + length(raw)`, see `_occupancy`) any single class'
  buffer may reach before it MUST be folded (generate-and-minimize fused)
  back down. Replaces what used to be TWO separate knobs (a flush
  threshold and a results-growth threshold), because folding now replaces
  the running formula in place instead of letting a separate "results
  pool" grow across many flushes.
- `N::Int`: number of independent random-order restarts to try; the
  smallest-total-terms restart (among those that stayed under `M`) is kept.
- `max_apply_batch::Int`: PURELY a performance cap on how many rows are
  decoded and `apply`'d together in one batch inside `_sequential_pass`
  (and reused by the repair sweep). It has NO effect on when folding
  happens -- that is governed only by `M`. It exists only so a single
  `apply` call doesn't try to materialize an unreasonably large number of
  rows at once when `M` itself is huge. Defaults to `min(M, 4096)`.
- `score::Symbol`: currently only `:n_terms` (total term count summed
  across all classes) is implemented. A literal-count-based scoring scheme
  would need a reliable, cross-representation "count atoms in a term"
  accessor that isn't assumed to exist here, since minimized terms can come
  back as arbitrary `SyntaxStructure`, not necessarily flat conjunctions.
- `repair_sensitivity::Bool`: whether to run the full-domain sensitivity
  repair sweep (step 4 above) on the winning restart before building the
  final `DecisionSet`. Default `true`. Set to `false` to skip it (e.g. to
  reproduce the pre-repair behavior, or when `M` is large enough relative
  to the model that early folds are rare/negligible and the extra sweep's
  cost isn't worth it).
- `max_repair_rounds::Int`: cap on repair-sweep iterations (see step 4
  above). Default `3`.

# Returns
`SM.DecisionSet` built from the best restart's (possibly repaired)
per-class terms, in the same shape as `lumen`'s (batch pipeline's) output.

**Type compatibility with `lumen()`**: the returned `DecisionSet`'s rule
antecedents are constructed with the exact same wide element type
(`Union{SL.LeftmostConjunctiveForm{SL.Atom{float_type}}, SyntaxStructure}`)
that `lumen()` itself uses for its `formulas` container, before broadcasting
`SL.LeftmostDisjunctiveForm.(...)` over it. This is not cosmetic: without it,
`_SeqBuffer.terms` (a `Vector{Any}` internally, needed so a class' running
formula can interchangeably hold raw cubes or minimizer output across folds)
would leak its narrow, runtime-inferred element type into the final
`LeftmostDisjunctiveForm`, causing it to match the stricter `DNF` type alias
that `lumen()`'s wider-typed result never matches. That mismatch makes
`check(::DefaultCheckAlgorithm, ::DNF, ...)` and
`check(::DefaultCheckAlgorithm, ::LeftmostDisjunctiveForm, ...)` equally
specific for the same call, and SoleLogics refuses to pick one
(`MethodError: ... is ambiguous`) -- surfacing downstream, for instance, in
`SoleModels.evaluaterule`/`calculate_decisionset_accuracy` calls on a
`lumen_sequential`-produced `DecisionSet`, even though the exact same calls
work fine on a `lumen()`-produced one. Forcing the wide type here keeps the
two extractors' outputs structurally interchangeable everywhere downstream.

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
    score::Symbol=:n_terms,
    repair_sensitivity::Bool=true,
    max_repair_rounds::Int=3
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

    # NEW: one full-domain, linear verification sweep on the winning
    # restart only (see the dedicated section above for the full
    # rationale). Fixes false negatives (sensitivity loss) caused by early
    # folds that generalized a term before having seen the whole domain.
    if repair_sensitivity
        best_terms = _verify_and_repair_sensitivity(
            config, model, ctx, scheme, best_terms, ctx.classnames;
            max_apply_batch, max_repair_rounds
        )
    end

    # Drop classes that ended up with zero terms (i.e. the model never
    # predicted that class for any enumerated row in the winning restart).
    valid_mask = .!isempty.(best_terms)
    classes = ctx.classnames[valid_mask]

    # ------------------------------------------------------------------ #
    # IMPORTANT: force the same "wide" element type that `lumen()` uses for
    # its own `formulas` container
    # (`Union{LeftmostConjunctiveForm{Atom{float_type}}, SyntaxStructure}`),
    # instead of letting `SL.LeftmostDisjunctiveForm.(...)` infer a
    # narrower type from `best_terms`'s runtime elements.
    #
    # `_fold_in!` stores each class' running formula in `buf.terms::Vector{Any}`
    # (needed internally so it can hold either bare cubes or minimizer output
    # interchangeably across folds). By the time we get here every element is
    # concretely a `LeftmostConjunctiveForm{Atom}` (or similar), but the
    # container itself is still `Vector{Any}`.
    #
    # If we broadcast `SL.LeftmostDisjunctiveForm.(...)` directly over that,
    # SoleLogics infers the resulting `LeftmostLinearForm`'s type parameter
    # from the RUNTIME elements, not from the container's static type --
    # producing a `LeftmostDisjunctiveForm{LeftmostConjunctiveForm{Atom}}`
    # that (unlike `lumen()`'s Union-typed result) matches the narrower `DNF`
    # alias. See the docstring above and `check`-dispatch ambiguity note for
    # the full explanation of why that alone breaks
    # `evaluaterule`/`calculate_decisionset_accuracy` downstream.
    #
    # Explicitly widening the container's eltype here (mirroring `lumen()`'s
    # declared `formulas` type) keeps `lumen()` and `lumen_sequential()`
    # producing bit-for-bit structurally compatible `DecisionSet`s, so
    # anything downstream behaves identically regardless of which extractor
    # produced the rules.
    # ------------------------------------------------------------------ #
    float_type = get_float_type(config)
    formulas = Vector{Vector{Union{
        SL.LeftmostConjunctiveForm{SL.Atom{float_type}},
        SyntaxStructure
    }}}(best_terms[valid_mask])

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end