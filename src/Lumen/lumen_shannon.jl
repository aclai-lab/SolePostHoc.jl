# ============================================================================ #
#   LUMEN SHANNON: recursive exact cofactor decomposition, SOP-optimal leaves  #
#                                                                              #
#   Same design as the earlier version of this file (see its history for the  #
#   full rationale), with ONE change: there is no `remin` knob anymore.       #
#   Internal nodes ALWAYS combine their two cofactors by plain set union,     #
#   never by feeding the children's terms back into the minimizer.           #
#                                                                              #
#   Why: re-minimizing a union of two cofactors' terms requires converting    #
#   each term back into a "cube" (`_term_to_cube`), and a term that dropped   #
#   the literal for the SPLIT feature is only a valid don't-care WITHIN the   #
#   sub-rectangle that leaf/subtree was restricted to -- not globally. Mixing #
#   such a term with the sibling cofactor's terms and re-minimizing let the   #
#   minimizer treat that dropped literal as a genuine global don't-care,      #
#   which can make it discard other terms as "redundant" when they weren't -- #
#   silently losing coverage (this showed up empirically as a sensitivity     #
#   drop concentrated on rarer classes, worse the more re-minimization        #
#   happened along the path to the root -- confirmed identical on both :abc   #
#   and :mitespresso, ruling out a backend-specific bug). See                 #
#   `lumen_shannon_boundary_fix.jl` for a version that repairs this (by       #
#   re-attaching the correct boundary literal before any re-minimization)     #
#   instead of simply avoiding it.                                           #
#                                                                              #
#   Plain union has none of that risk: each cofactor's terms are already a    #
#   fully correct description of their own (disjoint) sub-rectangle, and the  #
#   union of two correct, disjoint partitions is always correct. The only     #
#   cost is a larger (less compact) final formula -- confirmed empirically to #
#   match classic Lumen's Sensitivity/Specificity/DecisionSetAccuracy exactly #
#   on the `car` dataset, at the cost of somewhat more atoms/terms.           #
#                                                                              #
#   Relies on the same shared primitives as before: `_prepare_sequential_context`,
#   `generate_disjunct`, `_truths_by_thresholds`, `run_minimization`,
#   `PropositionalLogiset`, `get_apply_function`, `get_use_multithreads`,
#   `get_float_type`, `get_minimization_scheme`.
# ============================================================================ #
const TERM = SM.LeftmostConjunctiveForm{<:SM.Atom}

_as_terms(x::SM.LeftmostDisjunctiveForm)::Vector{TERM} = collect(TERM, SL.disjuncts(x))
_as_terms(x::AbstractVector)::Vector{TERM} = collect(TERM, x)
_as_terms(::SL.Truth)::Vector{TERM} = TERM[]   # ⊤ / empty PLA

# ---------------------------------------------------------------------------- #
#                                  leaf routine                                #
# ---------------------------------------------------------------------------- #
"""
    _leaf_extract(config, model, ctx, lo, hi, scheme; max_apply_batch) -> Vector{Vector{Any}}

Exhaustively materialize and classify every row of the sub-rectangle
`[lo[j], hi[j]]` (inclusive, 1-based, one range per feature dimension) and
return each class' minimized formula for that sub-rectangle.

Unchanged from before: this is "classic Lumen, restricted to a
sub-rectangle" -- full materialization of that sub-rectangle, one
`run_minimization` call per class over its WHOLE raw pool. Leaves are still
genuinely SOP-optimal for their own sub-rectangle; the change in this file
is only in how sibling leaves/subtrees get combined afterwards (see
`_combine_cofactors` below), not in how a single leaf is computed.
"""
# @inline function _eval_tree(node::SM.Branch, row::NTuple{N,T}, featidx::Dict{Symbol,Int}) where {N,T}
#     while node isa SM.Branch
#         cond = SL.value(SM.antecedent(node))            # ScalarCondition
#         v    = SD.i_variable(SD.feature(cond))
#         j    = v isa Integer ? Int(v) : featidx[Symbol(v)]
#         node = SD.test_operator(cond)(row[j], SD.threshold(cond)) ?
#             SM.posconsequent(node) : SM.negconsequent(node)
#     end
#     return SM.outcome(node)                             # ConstantModel leaf
# end

function checkcondition(
    atom::Atom{B},
    tbl::NamedTuple{T}
)::BitVector where {B<:ScalarCondition,T}
    cond = SL.value(atom)

    cond_threshold = SD.threshold(cond)
    cond_operator = SD.test_operator(cond)
    cond_feature = SD.feature(cond)

    col = SD.i_variable(cond_feature)

    return cond_operator.(tbl[col], cond_threshold)
end

function _leaf_extract(
    config::LumenConfig{R,T},
    ctx::Ctx{R,T},
    model::FlatForest{TT,U},
    lo::Vector{R},
    hi::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat,U,TT<:AbstractFloat}
    nfeat = length(ctx.featurenames)
    nclasses = length(ctx.class_idxs)
    raw = [Vector{Vector{SM.Atom}}() for _ in 1:nclasses]

    widths = Vector{Int}(undef, nfeat)
    @inbounds for j in 1:nfeat
        widths[j] = Int(hi[j]) - Int(lo[j]) + 1
    end
    total = prod(widths)
    batch = Int(config.max_apply_batch)

    i0 = 1
    while i0 ≤ total
        this_chunk = min(batch, total - i0 + 1)
        tbl = Matrix{T}(undef, this_chunk, nfeat)

        @inbounds for k in 1:this_chunk
            r = i0 + k - 2
            for j in 1:nfeat
                off = r % widths[j]
                r = r ÷ widths[j]
                tbl[k, j] = ctx.thrs_with_p[j][Int(lo[j]) + off]
            end
        end

        preds = apply(model, tbl)

        @inbounds for k in 1:this_chunk
            label = String(preds[k])
            ci_class = searchsortedfirst(ctx.classnames, label)

            truths_row = _truths_by_thresholds(tbl[k,:], ctx.thresholds)

            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            push!(raw[ci_class], cube)
        end

        i0 += this_chunk
    end

    terms = Vector{Vector{TERM}}(undef, nclasses)
    @inbounds for c in 1:nclasses
        terms[c] = isempty(raw[c]) ?
            TERM[] :
            run_minimization(config.minimization_scheme, config, raw[c])
    end

    return terms
end


# ---------------------------------------------------------------------------- #
#                    combining two cofactors: PLAIN UNION ONLY                 #
# ---------------------------------------------------------------------------- #
"""
    _combine_cofactors(terms_low, terms_high) -> Vector{Vector{Any}}

Combine the per-class results of two sibling cofactors by plain set union,
class by class. No re-minimization, ever -- see the module-level comment at
the top of this file for why: it's the only combination step that's
correct with no extra bookkeeping, given that `terms_low`/`terms_high` are
each already fully correct for their own (disjoint) sub-rectangle, and the
two sub-rectangles exactly partition the parent's rectangle.

This is intentionally the ONLY combination strategy available in this file
-- there is no flag to opt back into re-minimization here, precisely so the
sensitivity regression that motivated this file can't silently come back.
"""
function _combine_cofactors(
    terms_low::Vector{Vector{SM.LeftmostConjunctiveForm{<:SM.Atom}}},
    terms_high::Vector{Vector{SM.LeftmostConjunctiveForm{<:SM.Atom}}}
)
    nclasses = length(terms_low)
    return [vcat(terms_low[c], terms_high[c]) for c in 1:nclasses]
end


# ---------------------------------------------------------------------------- #
#                     recursive Shannon decomposition driver                   #
# ---------------------------------------------------------------------------- #
"""
    _shannon_extract(config, model, ctx, lo, hi, scheme; M, max_apply_batch)

Recursively decompose the index sub-rectangle `[lo, hi]` into two exact
cofactors (splitting the currently-widest feature dimension in half) until
each piece is small enough to materialize directly (`_leaf_extract`), then
recombine on the way back up by plain union (`_combine_cofactors`).

Same split rule and termination argument as before: pick the widest
remaining dimension, cut it at its midpoint, recurse; terminates because
each split strictly shrinks the widest dimension and the total size is
finite (worst-case depth `O(log2(n_total / M))`).
"""
function _shannon_extract(
    config::LumenConfig{R,T},
    ctx::Ctx{R,T},
    model::FlatForest{TT,U},
    lo::Vector{R},
    hi::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat,U,TT<:AbstractFloat}
    rect_size = prod(hi .- lo .+ one(R))

    if rect_size ≤ config.M
        return _leaf_extract(config, ctx, model, lo, hi)
    end

    jstar = argmax(hi .- lo .+ 1)

    t = lo[jstar] + (hi[jstar] - lo[jstar]) >> 1
    old_hi = hi[jstar]; hi[jstar] = t
    terms_low = _shannon_extract(config, ctx, model, lo, hi)
    hi[jstar] = old_hi

    old_lo = lo[jstar]; lo[jstar] = t + one(R)
    terms_high = _shannon_extract(config, ctx, model, lo, hi)
    lo[jstar] = old_lo

    return _combine_cofactors(terms_low, terms_high)
end


# ---------------------------------------------------------------------------- #
#                      shared DecisionSet finalization                        #
# ---------------------------------------------------------------------------- #
"""
    _finalize_decision_set(ctx, per_class_terms, config) -> SM.DecisionSet

Build the final `SM.DecisionSet` from per-class terms, dropping classes with
zero terms and forcing the same "wide" element type `lumen()` uses for its
`formulas` container (see the identical note in the sequential file's tail
for why this widening is required to avoid a `DNF` vs
`LeftmostDisjunctiveForm` dispatch ambiguity downstream).
"""
function _finalize_decision_set(
    ctx::Ctx{R,T},
    per_class_terms::Vector{Vector{SM.LeftmostConjunctiveForm{<:SM.Atom}}},
    config::LumenConfig{R,T}
) where {R,T<:AbstractFloat}
    valid_mask = .!isempty.(per_class_terms)
    classes = ctx.class_idxs[valid_mask]
    formulas = per_class_terms[valid_mask]

    rules = SM.Rule{R}[
        SM.Rule(SL.LeftmostDisjunctiveForm(f), c)
        for (f, c) in zip(formulas, classes)
    ]

    return SM.DecisionSet{R}(rules)::SM.DecisionSet{R}
end


# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    lumen_shannon(config::LumenConfig, model::SM.AbstractModel;
                  M::Int=20_000,
                  max_apply_batch::Int=min(M, 4096)) -> SM.DecisionSet

Exact, deterministic rule extraction via recursive Shannon-style cofactor
decomposition of the (possibly astronomically large) threshold-combination
space, with SOP-optimal leaves and plain-union recombination.

There is no re-minimization knob here on purpose: internal nodes always
combine their two cofactors by set union, which is the one recombination
strategy that is correct with no further conditions. This means coverage
("sensitivity") always matches classic (fully materializing) Lumen exactly
-- confirmed empirically on `car` (identical Sensitivity/Specificity/
DecisionSetAccuracy across all classes) -- at the cost of a less compact
final formula than classic Lumen or a correctly-fixed re-minimizing variant
would produce (more terms/atoms, since nothing above the leaf level is ever
compacted).

# Arguments
- `config`, `model`: as before.
- `M::Int`: leaf budget -- max sub-rectangle size materialized (and hence
  max cubes handed to the minimizer) in one `_leaf_extract` call.
- `max_apply_batch::Int`: pure decode/apply performance cap inside leaves.
  No effect on correctness or on `M`'s semantics.

# Returns
`SM.DecisionSet`, structurally compatible with `lumen()`'s output.

# Throws
- `ArgumentError` if `M` is not positive.
"""

# map categorical labels onto unsigned integer of type `R`
# sort classlabels
function assign(
    ::Type{R},
    y::Vector{S}
) where {R<:Unsigned,S<:CategoricalValue}
    classlabels = sort!(unique(y))
    dict = Dict{S,R}(v => i for (i, v) in enumerate(classlabels))
    return string.(classlabels), [dict[t] for t in y]
end

function lumen_shannon(
    config::LumenConfig{R,T},
    model::SM.DecisionEnsemble{R,SM.Branch{S}},
) where {S<:CategoricalValue,R<:Unsigned,T<:AbstractFloat}
    flatmodel = flatten(model)
    classnames, class_idxs = assign(R, unique!(SM.info(model, :supporting_labels)))
    featurenames = SM.info(model, :featurenames)

    ctx = _prepare_sequential_context(config, _extract_atoms_bfs_order(model), featurenames, classnames, class_idxs)

    lo = ones(R, length(ctx.lens))
    hi = ctx.lens

    per_class_terms = _shannon_extract(config, ctx, flatmodel, lo, hi)

    return _finalize_decision_set(ctx, per_class_terms, config)
end
