# ---------------------------------------------------------------------------- #
#                     LAZY RANDOM PERMUTATION (no materialization)             #
# ---------------------------------------------------------------------------- #
"""
    LazyPermutation

Bijective pseudo-random permutation of `{0, ..., n-1}` computed on the fly,
without ever materializing an array of size `n`.

Built on a small balanced Feistel network over the smallest power-of-two
domain `M ≥ n`, plus cycle-walking to fold `[n, M)` back into `[0, n)`.
This gives O(1) amortized cost per lookup regardless of how large `n` is
(the whole point: `n` here is `prod(length.(thrs_with_p))`, which can be
astronomically large and must never be enumerated into memory).

Different `seed` values give different permutations — this is what powers
the "random restart" loop: each restart just picks a fresh seed.

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

function LazyPermutation(n::Integer; seed::Integer=rand(UInt64), rounds::Int=4)
    n <= 0 && throw(ArgumentError("n must be positive, got $n"))
    total_bits = max(2, ceil(Int, log2(n)))
    total_bits = iseven(total_bits) ? total_bits : total_bits + 1
    half_bits = total_bits ÷ 2
    mask = (UInt64(1) << half_bits) - UInt64(1)
    return LazyPermutation(Int(n), half_bits, mask, rounds, UInt64(seed))
end

@inline function _feistel_round(l::UInt64, r::UInt64, round::Int, p::LazyPermutation)
    h = (hash((r, round, p.key)) & p.mask)
    newl = r
    newr = (l ⊻ h) & p.mask
    return newl, newr
end

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
`{0, ..., p.n-1}` described by `p`. O(1) amortized via cycle-walking.
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

Mixed-radix strides for a product space with per-dimension sizes `lens`,
first dimension varying fastest (consistent with `_ProductColumn` in Lumen.jl).
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

Decode a 0-based linear index into the corresponding tuple of threshold
values, one per feature, without ever materializing the product.
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

Standalone replica of steps 1–7 of `ExtractRulesData`'s inner constructor,
stopping *before* the full Cartesian-product / model-apply step (step 8–9),
which is exactly the part we refuse to materialize here.

Kept intentionally separate from `ExtractRulesData` (rather than refactored
into a shared helper) to avoid touching the existing, working batch pipeline.

# Returns
A `NamedTuple` with:
- `thresholds::Vector{Vector{<:Float}}`   — per-feature thresholds, no boundary
- `thrs_with_p::Vector{Vector{<:Float}}`  — per-feature thresholds + boundary
- `featurenames::Vector{Symbol}`
- `classnames::AbstractVector`
- `op_families::Vector{Symbol}`
- `lens::Vector{Int}`, `strides::Vector{Int}`, `n_total::Int`
"""
function _prepare_sequential_context(config::LumenConfig, model::SM.AbstractModel)
    depth = get_depth(config)

    atoms = unique!(_normalize_atom.(if depth < 1.0
        mapreduce(vcat, SM.models(model); init=SL.Atom{SD.AbstractCondition}[]) do t
            _take_first_percentage(_extract_atoms_bfs_order(t), depth)
        end
    else
        SL.atoms(SM.alphabet(model, false))
    end))

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

Bookkeeping for one class during a single restart pass:
- `raw::Vector{Vector{SL.Atom}}` — unminimized cubes waiting to be flushed.
- `results::Vector`              — minimized terms collected from past flushes.
- `n_rows_seen::Int`             — total rows ever routed to this class
  (diagnostic only).
"""
mutable struct _SeqBuffer
    raw::Vector{Vector{SL.Atom}}
    results::Vector{Any}
    n_rows_seen::Int
end
_SeqBuffer() = _SeqBuffer(Vector{Vector{SL.Atom}}(), Any[], 0)

# ---------------------------------------------------------------------------- #
#                          flush / compaction logic                            #
# ---------------------------------------------------------------------------- #
"""
    _flush!(buf, scheme::Symbol, config, max_results_terms)

If `buf.raw` is non-empty, minimize it via `run_minimization` (dispatches on
`scheme`, so this is the SAME code path for `:abc` and `:mitespresso`) and
append the resulting terms to `buf.results`; then empty `buf.raw`.

If, after appending, `length(buf.results)` exceeds `max_results_terms`,
compact `buf.results` in place via `_refine_dnf` (cheap, no external binary,
works directly on terms/hyper-rectangles). If it's *still* over the hard
limit after compaction, return `false` (caller should treat this restart as
failed / give up on this class); otherwise `true`.
"""
function _flush!(
    buf::_SeqBuffer,
    scheme::Symbol,
    config::LumenConfig,
    max_results_terms::Int
)
    if !isempty(buf.raw)
        minimized = run_minimization(Val(scheme), config, buf.raw)
        append!(buf.results, minimized)
        empty!(buf.raw)
    end

    if length(buf.results) > max_results_terms
        # NOTE: must convert to Vector{SyntaxStructure}, NOT Vector{Any}.
        # Julia's parametric types are invariant: Vector{Any} does not match
        # the `Vector{<:SyntaxStructure}` signature of `_refine_dnf`, even
        # though every element already is a SyntaxStructure.
        refined = _refine_dnf(convert(Vector{SyntaxStructure}, buf.results))
        buf.results = convert(Vector{Any}, refined)
        if length(buf.results) > max_results_terms
            return false
        end
    end

    return true
end

# ---------------------------------------------------------------------------- #
#                        one sequential pass (one restart)                     #
# ---------------------------------------------------------------------------- #
"""
    _sequential_pass(config, model, ctx, scheme, seed;
                      chunk_size, max_buffer_rows, max_results_terms)

Enumerate the whole combination space exactly once, in the pseudo-random
order given by `LazyPermutation(ctx.n_total; seed)`, routing each row's
derived atom-cube into the buffer of its predicted class. Buffers are
flushed (minimized) whenever they cross `max_buffer_rows`.

Returns `(per_class_terms, ok)` where `per_class_terms[c]` is the (possibly
further-compacted) list of minimized terms for `ctx.classnames[c]`, and `ok`
is `false` if any class's result pool blew past `max_results_terms` even
after compaction (i.e. this restart should be discarded).
"""
function _sequential_pass(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    scheme::Symbol,
    seed::Integer;
    chunk_size::Int,
    max_buffer_rows::Int,
    max_results_terms::Int
)
    nclasses = length(ctx.classnames)
    class_index = Dict(c => i for (i, c) in enumerate(ctx.classnames))
    buffers = [_SeqBuffer() for _ in 1:nclasses]

    perm = LazyPermutation(ctx.n_total; seed=seed)
    ok = true

    i0 = 0
    while i0 < ctx.n_total
        this_chunk = min(chunk_size, ctx.n_total - i0)

        # decode `this_chunk` pseudo-random rows without materializing more
        rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
        @inbounds for k in 1:this_chunk
            linidx = permute(perm, i0 + k - 1)  # 0-based
            rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, linidx)
        end
        i0 += this_chunk

        # batch-apply the model once per chunk (cheap relative to per-row apply)
        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        # route each row's cube into its predicted class's raw buffer
        @inbounds for k in 1:this_chunk
            ci = get(class_index, preds[k], nothing)
            isnothing(ci) && continue  # defensive: unseen label, skip

            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            buf = buffers[ci]
            push!(buf.raw, cube)
            buf.n_rows_seen += 1

            if length(buf.raw) >= max_buffer_rows
                ok &= _flush!(buf, scheme, config, max_results_terms)
            end
        end
    end

    # final flush of whatever is left in each buffer
    for buf in buffers
        ok &= _flush!(buf, scheme, config, max_results_terms)
    end

    per_class_terms = [buf.results for buf in buffers]
    return per_class_terms, ok
end

# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    lumen_sequential(config::LumenConfig, model::SM.AbstractModel;
                      restarts::Int=10,
                      chunk_size::Int=4096,
                      max_buffer_rows::Int=100_000,
                      max_results_terms::Int=1_000_000,
                      score::Symbol=:n_terms) -> SM.DecisionSet

Generate-and-minimize fused into a single online pass, run `restarts` times
with independent random enumeration orders, keeping the smallest result.

Works identically for `get_minimization_scheme(config) ∈ (:abc, :mitespresso)`
— the flush step just dispatches through the existing `run_minimization`
multi-method, so no backend-specific code is needed here. (`:abc` genuinely
never sees more than `max_buffer_rows` cubes at once, so the "must have the
whole PLA" constraint is respected — it's just never handed the *whole*
table, only successive windows of it.)

# Arguments
- `config::LumenConfig`: same config used by `lumen`; `minimization_scheme`
  selects `:abc` or `:mitespresso`.
- `model::SM.AbstractModel`: the model to extract rules from.
- `restarts`: number of independent random-order passes; the best
  (fewest-terms) one is kept. Passes that blow the memory budget even after
  compaction are discarded and don't count as a valid candidate.
- `chunk_size`: rows decoded / model-applied together per batch. Larger is
  faster (fewer `apply` calls) but each chunk is fully materialized in
  memory before routing, so keep it well below `max_buffer_rows`.
- `max_buffer_rows`: per-class raw-cube buffer size that triggers a
  minimization flush (the "sto per sforare, chiamo Espresso/ABC" threshold).
- `max_results_terms`: per-class accumulated-results size that triggers a
  `_refine_dnf` compaction pass (and, if still over, marks the restart as
  failed). This is what prevents unbounded growth of `buf.results` across
  many flushes over a huge combination space.
- `score`: currently only `:n_terms` (total term count across classes) is
  implemented — literal-count scoring would need a reliable
  cross-representation "count atoms in a term" accessor that isn't assumed
  here (minimized terms can come back as arbitrary `SyntaxStructure`, not
  just flat conjunctions).

# Returns
`SM.DecisionSet` built from the best restart's per-class terms, same shape
as `lumen`'s output.
"""
function lumen_sequential(
    config::LumenConfig,
    model::SM.AbstractModel;
    restarts::Int=10,
    chunk_size::Int=4096,
    max_buffer_rows::Int=100_000,
    max_results_terms::Int=1_000_000,
    score::Symbol=:n_terms
)
    score === :n_terms || throw(ArgumentError(
        "Only score=:n_terms is currently supported."
    ))
    chunk_size <= max_buffer_rows || throw(ArgumentError(
        "chunk_size ($chunk_size) must be <= max_buffer_rows " *
        "($max_buffer_rows), otherwise a single chunk can blow the budget " *
        "before any flush has a chance to run."
    ))

    scheme = get_minimization_scheme(config)
    ctx = _prepare_sequential_context(config, model)

    best_terms = nothing
    best_score = typemax(Int)

    for r in 1:restarts
        seed = rand(UInt64)
        per_class_terms, ok = _sequential_pass(
            config, model, ctx, scheme, seed;
            chunk_size, max_buffer_rows, max_results_terms
        )
        ok || continue

        this_score = sum(length, per_class_terms)
        if this_score < best_score
            best_score = this_score
            best_terms = per_class_terms
        end
    end

    isnothing(best_terms) && error(
        "lumen_sequential: all $restarts restarts exceeded the memory " *
        "budget (max_results_terms=$max_results_terms) even after " *
        "compaction. Raise max_results_terms / max_buffer_rows, or lower " *
        "chunk_size."
    )

    valid_mask = .!isempty.(best_terms)
    formulas = best_terms[valid_mask]
    classes = ctx.classnames[valid_mask]

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end