# ============================================================================ #
#   LUMEN SHANNON: recursive exact cofactor decomposition, SOP-optimal leaves  #
#                                                                              #
#   Motivation
#   ----------
#   `lumen_sequential` fuses generate+minimize with random-order enumeration
#   and per-class folding once a buffer reaches `M`. That keeps memory
#   bounded, but a restart can be discarded (`ok = false`) whenever some
#   class can't be folded back under `M`, and the winning restart is chosen
#   only among the `N` random orders that survive. In practice this can make
#   the extracted ruleset less complete ("sensitivity") than classic Lumen,
#   which always sees every single row.
#
#   `lumen_shannon` restores exhaustive coverage (every row of the product
#   space is visited exactly once, deterministically, no discarding, no
#   restarts) while still respecting a hard memory budget `M`, by using an
#   *exact* recursive partition of the index space instead of random
#   sampling + folding:
#
#     f(x_1, ..., x_k) = f|_{x_j <= t}  ∪  f|_{x_j > t}
#
#   which is exactly Shannon's expansion theorem, generalized from a single
#   Boolean variable to the ordered multi-valued domain each feature already
#   has in this codebase (`thrs_with_p[j]`, indices `1:lens[j]`): picking a
#   cutpoint `t` inside a feature's threshold chain is itself a legitimate
#   literal for that feature (the chain is sorted per `op_families[j]` by
#   `_prepare_sequential_context`), so splitting the *index range* in half is
#   the same thing as picking two cofactors on that literal.
#
#   The two cofactors are sub-hyperrectangles of the index space (one range
#   per feature). Recursion keeps splitting the currently-widest dimension
#   in half until a sub-rectangle's size drops to <= M, at which point it is
#   materialized *in full* (classic Lumen's generate step, but restricted to
#   that sub-rectangle) and minimized once -- this is the "leaf SOP-optimal
#   expression" the BDD hangs off of. Going back up the recursion, sibling
#   leaves/subtrees are combined by set union of their (already minimal)
#   terms; if that union is still small enough (`<= M2`) it is fed back into
#   the minimizer as a fresh batch of cubes (legal because a minimized term
#   with some literals dropped IS already a valid cube -- same argument as
#   `_term_to_cube`/`_fold_in!` in the sequential file), which is where
#   opportunistic *further* compaction happens on the way up.
#
#   What this buys you
#   -------------------
#   - Exactness / completeness: identical to classic Lumen's coverage, since
#     the recursion is an exact, lossless partition of the whole product
#     space -- no row is ever skipped and no restart is ever discarded.
#   - Peak memory: still bounded by `max(M, M2)` cubes materialized/fed to
#     the minimizer in any single call, same style of guarantee as
#     `lumen_sequential`.
#   - Minimality: each LEAF is genuinely SOP-optimal (single minimizer call
#     over its whole sub-rectangle). Each INTERNAL node is minimal only when
#     its children's combined term count fits under `M2` and therefore gets
#     re-minimized; when it doesn't fit, the node's result is just the union
#     of its children's terms (a valid but not necessarily 2-level-minimal
#     formula for that subtree) -- i.e. a genuine BDD-of-SOPs, correct
#     everywhere, minimal only where budget allowed re-minimization. This
#     matches the "leaves SOP-optimal, global minimality not guaranteed"
#     design discussed for this approach.
#   - No randomness anywhere: deterministic split order (widest remaining
#     dimension), so a given `(config, model, M, M2)` always returns the
#     same ruleset -- there's no `N` restarts knob because there's nothing
#     to restart.
#
#   Relies on the same shared primitives as `lumen_sequential.jl`:
#   `_prepare_sequential_context`, `_term_to_cube`, `generate_disjunct`,
#   `_truths_by_thresholds`, `run_minimization`, `PropositionalLogiset`,
#   `get_apply_function`, `get_use_multithreads`, `get_float_type`,
#   `get_minimization_scheme`. This file assumes those are already in scope
#   (defined elsewhere in the same module), exactly like `lumen_sequential.jl`
#   does.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
#                                  leaf routine                                #
# ---------------------------------------------------------------------------- #
"""
    _leaf_extract(config, model, ctx, lo, hi, scheme; max_apply_batch) -> Vector{Vector{Any}}

Exhaustively materialize and classify every row of the sub-rectangle
`[lo[j], hi[j]]` (inclusive, 1-based, one range per feature dimension) and
return each class' minimized formula for that sub-rectangle.

This is deliberately "classic Lumen, restricted to a sub-rectangle": every
row in range is decoded and applied to the model (no sampling, no skipping),
its cube is routed to its predicted class' raw pool, and -- once the whole
sub-rectangle has been consumed -- each class' full raw pool is minimized
with exactly ONE `run_minimization` call. There is no folding/fusing here:
the caller (`_shannon_extract`) only invokes this once a sub-rectangle's
total size is already `<= M`, so a single pool per class can never exceed
`M` cubes, and doing the generate/minimize split (rather than fusing them)
is what makes this leaf provably SOP-optimal for its sub-rectangle -- the
minimizer sees the WHOLE leaf at once, same as classic (non-sequential)
Lumen would for that region.

`max_apply_batch` is purely a decode/apply chunk-size cap for performance
(mirrors the same knob in `lumen_sequential.jl`); it has no effect on
correctness or on the final minimized result, since every row in range is
visited regardless of chunk size.

# Returns
`terms::Vector{Vector{Any}}`, one entry per `ctx.classnames`, each either
the class' minimized DNF terms for this sub-rectangle or `Any[]` if no row
in range was predicted as that class.
"""
function _leaf_extract(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    lo::Vector{Int},
    hi::Vector{Int},
    scheme::Symbol;
    max_apply_batch::Int
)
    nfeat = length(ctx.featurenames)
    nclasses = length(ctx.classnames)
    class_index = Dict(c => i for (i, c) in enumerate(ctx.classnames))
    raw = [Vector{Vector{SL.Atom}}() for _ in 1:nclasses]

    # Exhaustive, deterministic enumeration of the sub-rectangle: no
    # LazyPermutation needed here at all, since we visit every index anyway
    # (order genuinely doesn't matter for a full materialization).
    dims = ntuple(j -> lo[j]:hi[j], nfeat)
    all_idx = vec(CartesianIndices(dims))
    total = length(all_idx)

    i0 = 1
    @inbounds while i0 <= total
        this_chunk = min(max_apply_batch, total - i0 + 1)
        chunk = @view all_idx[i0:(i0+this_chunk-1)]

        rows = Vector{NTuple{nfeat,eltype(ctx.thresholds[1])}}(undef, this_chunk)
        for (k, ci) in enumerate(chunk)
            rows[k] = ntuple(j -> ctx.thrs_with_p[j][ci[j]], nfeat)
        end

        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], nfeat)
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        for k in 1:this_chunk
            ci_class = get(class_index, preds[k], nothing)
            isnothing(ci_class) && continue  # defensive: unseen label, skip

            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            push!(raw[ci_class], cube)
        end

        i0 += this_chunk
    end

    terms = Vector{Vector{Any}}(undef, nclasses)
    for c in 1:nclasses
        if isempty(raw[c])
            terms[c] = Any[]
        else
            minimized = run_minimization(Val(scheme), config, raw[c])
            terms[c] = collect(Any, minimized)
        end
    end
    return terms
end


# ---------------------------------------------------------------------------- #
#                     combining two (already minimal) cofactors                #
# ---------------------------------------------------------------------------- #
"""
    _combine_cofactors(terms_low, terms_high, scheme, config; M2, remin) -> Vector{Vector{Any}}

Combine the per-class results of two sibling cofactors (`terms_low`,
`terms_high`) into the parent's per-class result.

Because `terms_low` and `terms_high` come from an EXACT partition of the
parent's index rectangle (disjoint sub-rectangles whose union is the whole
parent rectangle), the parent's correct formula for class `c` is simply the
union `terms_low[c] ∪ terms_high[c]` -- this is Shannon's theorem applied
at the level of already-derived cube sets, no literal needs to be
re-attached to either side, since each leaf's cubes were generated against
the FULL (unrestricted) threshold set (see `_leaf_extract`), not against
some abstracted/projected view of the split feature.

If that union is small enough (`length(union) <= M2`), it is opportunistically
re-minimized: each term (whether a raw cube or a previously minimized
formula) is converted back to cube form via `_term_to_cube` -- legal for the
same reason `_fold_in!` in `lumen_sequential.jl` re-feeds old terms to the
minimizer: a minimized term with some literals dropped IS already a valid
"don't care" cube. This is what lets minimality improve on the way back up
the tree, beyond what either leaf could see alone.

If the union is too large to re-minimize under `M2`, it is returned as-is
-- correct, but only as good as the union of its two children (a BDD
branch, not a re-optimized 2-level form). `remin=false` disables
re-minimization altogether (cheaper, purely a union at every internal node)
for cases where you want speed over marginal minimality gains.

# Returns
`Vector{Vector{Any}}`, one entry per class, either re-minimized or a plain
union, per the rule above.
"""
function _combine_cofactors(
    terms_low::Vector{Vector{Any}},
    terms_high::Vector{Vector{Any}},
    scheme::Symbol,
    config::LumenConfig;
    M2::Int,
    remin::Bool
)
    nclasses = length(terms_low)
    out = Vector{Vector{Any}}(undef, nclasses)

    for c in 1:nclasses
        union_terms = vcat(terms_low[c], terms_high[c])

        if remin && !isempty(union_terms) && length(union_terms) <= M2
            combined_cubes = convert(
                Vector{Vector{SL.Atom}},
                _term_to_cube.(union_terms)
            )
            minimized = run_minimization(Val(scheme), config, combined_cubes)
            out[c] = collect(Any, minimized)
        else
            out[c] = union_terms
        end
    end

    return out
end


# ---------------------------------------------------------------------------- #
#                     recursive Shannon decomposition driver                   #
# ---------------------------------------------------------------------------- #
"""
    _shannon_extract(config, model, ctx, lo, hi, scheme; M, M2, remin, max_apply_batch)

Recursively decompose the index sub-rectangle `[lo, hi]` (one inclusive
range per feature) into two exact cofactors, until each piece is small
enough to materialize directly (`_leaf_extract`), then recombine on the way
back up (`_combine_cofactors`).

# Split rule
Picks the feature dimension `j*` with the currently widest remaining range
(`argmax(hi .- lo .+ 1)`), and cuts it at its midpoint `t`, giving:
  - low  cofactor: same rectangle, but dimension `j*` restricted to `[lo[j*], t]`
  - high cofactor: same rectangle, but dimension `j*` restricted to `[t+1, hi[j*]]`
Both are legitimate sub-rectangles of the same product space (indices into
`ctx.thrs_with_p`, exactly like `_combo_at` uses elsewhere), and together
they exactly partition the parent rectangle -- nothing is lost or double
counted. Dimensions with `hi[j] == lo[j]` (feature already pinned, e.g.
zero-threshold features) are never chosen as the split axis; if `hi != lo`
somewhere then `size > 1`, and since this function only recurses when
`size > M >= 1`, such a splittable dimension is always available.

# Termination
Base case: `prod(hi .- lo .+ 1) <= M` -> `_leaf_extract`. Since each split
strictly shrinks the widest dimension, and the rectangle's total size is
finite, this always terminates (worst case depth is
`O(log2(n_total / M))`, one halving of the currently-widest axis per level).

# Returns
`Vector{Vector{Any}}`, one entry per `ctx.classnames`: this sub-rectangle's
(possibly not fully 2-level-minimal, but always CORRECT) formula per class.
"""
function _shannon_extract(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    lo::Vector{Int},
    hi::Vector{Int},
    scheme::Symbol;
    M::Int,
    M2::Int,
    remin::Bool,
    max_apply_batch::Int
)
    rect_size = prod(hi .- lo .+ 1)

    if rect_size <= M
        return _leaf_extract(config, model, ctx, lo, hi, scheme; max_apply_batch)
    end

    jstar = argmax(hi .- lo .+ 1)
    t = lo[jstar] + (hi[jstar] - lo[jstar]) ÷ 2  # midpoint cut

    lo_low, hi_low = copy(lo), copy(hi)
    hi_low[jstar] = t

    lo_high, hi_high = copy(lo), copy(hi)
    lo_high[jstar] = t + 1

    terms_low = _shannon_extract(
        config, model, ctx, lo_low, hi_low, scheme; M, M2, remin, max_apply_batch
    )
    terms_high = _shannon_extract(
        config, model, ctx, lo_high, hi_high, scheme; M, M2, remin, max_apply_batch
    )

    return _combine_cofactors(terms_low, terms_high, scheme, config; M2, remin)
end


# ---------------------------------------------------------------------------- #
#                      shared DecisionSet finalization                        #
# ---------------------------------------------------------------------------- #
"""
    _finalize_decision_set(ctx, per_class_terms, config) -> SM.DecisionSet

Build the final `SM.DecisionSet` from per-class terms, dropping classes with
zero terms and forcing the same "wide" element type `lumen()` uses for its
`formulas` container -- see the identical note in `lumen_sequential.jl`'s
tail for why this exact widening is required to avoid a `DNF` vs
`LeftmostDisjunctiveForm` dispatch ambiguity downstream in
`evaluaterule`/`calculate_decisionset_accuracy`. Factored out here so
`lumen_shannon` doesn't have to duplicate that logic (and so
`lumen_sequential` could be pointed at this same helper later, if desired).
"""
function _finalize_decision_set(
    ctx::NamedTuple,
    per_class_terms::Vector{Vector{Any}},
    config::LumenConfig
)
    valid_mask = .!isempty.(per_class_terms)
    classes = ctx.classnames[valid_mask]

    float_type = get_float_type(config)
    formulas = Vector{Vector{Union{
        SL.LeftmostConjunctiveForm{SL.Atom{float_type}},
        SyntaxStructure
    }}}(per_class_terms[valid_mask])

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end


# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    lumen_shannon(config::LumenConfig, model::SM.AbstractModel;
                  M::Int=20_000,
                  M2::Int=M,
                  remin::Bool=true,
                  max_apply_batch::Int=min(M, 4096)) -> SM.DecisionSet

Exact, deterministic rule extraction via recursive Shannon-style cofactor
decomposition of the (possibly astronomically large) threshold-combination
space, with SOP-optimal leaves. See the module-level comment at the top of
this file for the full design rationale and how this differs from
`lumen_sequential`.

Unlike `lumen_sequential`, there is no `N` (no random restarts, nothing to
discard): the recursion visits the ENTIRE combination space exactly once,
so coverage/"sensitivity" matches classic (fully materializing) Lumen
exactly, while peak memory for any single materialize-or-minimize call is
still bounded (by `M` at the leaves, `M2` on the way back up).

# Arguments
- `config`, `model`: as in `lumen_sequential`.
- `M::Int`: leaf budget -- max sub-rectangle size materialized (and hence
  max cubes handed to the minimizer) in one `_leaf_extract` call. Smaller
  `M` -> more, smaller leaves -> lower peak memory, more (cheaper)
  minimizer calls, more BDD structure (more places where minimality is
  only local).
- `M2::Int`: recombination budget -- max combined term count
  (`length(terms_low[c]) + length(terms_high[c])`) an internal node will
  re-minimize, per class. Defaults to `M`. Larger `M2` -> more of the tree
  gets re-optimized on the way up (closer to global minimality) at the cost
  of more/bigger minimizer calls near the root.
- `remin::Bool`: if `false`, internal nodes never re-minimize -- every
  internal node is a plain union of its children (fastest, least minimal,
  still always correct).
- `max_apply_batch::Int`: pure decode/apply performance cap inside leaves,
  same role as in `lumen_sequential`. No effect on correctness or on `M`/`M2`
  semantics.

# Returns
`SM.DecisionSet`, structurally compatible with both `lumen()`'s and
`lumen_sequential()`'s output (see `_finalize_decision_set`).

# Throws
- `ArgumentError` if `M` or `M2` is not positive.
"""
function lumen_shannon(
    config::LumenConfig,
    model::SM.AbstractModel;
    M::Int=20_000,
    M2::Int=M,
    remin::Bool=true,
    max_apply_batch::Int=min(M, 4096)
)
    M >= 1 || throw(ArgumentError("M must be positive, got $M"))
    M2 >= 1 || throw(ArgumentError("M2 must be positive, got $M2"))

    scheme = get_minimization_scheme(config)
    ctx = _prepare_sequential_context(config, model)

    lo = ones(Int, length(ctx.lens))
    hi = copy(ctx.lens)

    per_class_terms = _shannon_extract(
        config, model, ctx, lo, hi, scheme;
        M, M2, remin, max_apply_batch
    )

    return _finalize_decision_set(ctx, per_class_terms, config)
end