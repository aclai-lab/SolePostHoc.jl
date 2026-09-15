# # ============================================================================ #
# #   LUMEN SHANNON: recursive exact cofactor decomposition, SOP-optimal leaves  #
# #                                                                              #
# #   Same design as the earlier version of this file (see its history for the  #
# #   full rationale), with ONE change: there is no `remin` knob anymore.       #
# #   Internal nodes ALWAYS combine their two cofactors by plain set union,     #
# #   never by feeding the children's terms back into the minimizer.           #
# #                                                                              #
# #   Why: re-minimizing a union of two cofactors' terms requires converting    #
# #   each term back into a "cube" (`_term_to_cube`), and a term that dropped   #
# #   the literal for the SPLIT feature is only a valid don't-care WITHIN the   #
# #   sub-rectangle that leaf/subtree was restricted to -- not globally. Mixing #
# #   such a term with the sibling cofactor's terms and re-minimizing let the   #
# #   minimizer treat that dropped literal as a genuine global don't-care,      #
# #   which can make it discard other terms as "redundant" when they weren't -- #
# #   silently losing coverage (this showed up empirically as a sensitivity     #
# #   drop concentrated on rarer classes, worse the more re-minimization        #
# #   happened along the path to the root -- confirmed identical on both :abc   #
# #   and :mitespresso, ruling out a backend-specific bug). See                 #
# #   `lumen_shannon_boundary_fix.jl` for a version that repairs this (by       #
# #   re-attaching the correct boundary literal before any re-minimization)     #
# #   instead of simply avoiding it.                                           #
# #                                                                              #
# #   Plain union has none of that risk: each cofactor's terms are already a    #
# #   fully correct description of their own (disjoint) sub-rectangle, and the  #
# #   union of two correct, disjoint partitions is always correct. The only     #
# #   cost is a larger (less compact) final formula -- confirmed empirically to #
# #   match classic Lumen's Sensitivity/Specificity/DecisionSetAccuracy exactly #
# #   on the `car` dataset, at the cost of somewhat more atoms/terms.           #
# #                                                                              #
# #   Relies on the same shared primitives as before: `_prepare_sequential_context`,
# #   `generate_disjunct`, `_truths_by_thresholds`, `run_minimization`,
# #   `PropositionalLogiset`, `get_apply_function`, `get_use_multithreads`,
# #   `get_float_type`, `get_minimization_scheme`.
# # ============================================================================ #
# const TERM = SM.LeftmostConjunctiveForm{<:SM.Atom}

# _as_terms(x::SM.LeftmostDisjunctiveForm)::Vector{TERM} = collect(TERM, SL.disjuncts(x))
# _as_terms(x::AbstractVector)::Vector{TERM} = collect(TERM, x)
# _as_terms(::SL.Truth)::Vector{TERM} = TERM[]   # ⊤ / empty PLA

# # ---------------------------------------------------------------------------- #
# #                                  leaf routine                                #
# # ---------------------------------------------------------------------------- #
# """
#     _leaf_extract(config, model, ctx, lo, hi, scheme; max_apply_batch) -> Vector{Vector{Any}}

# Exhaustively materialize and classify every row of the sub-rectangle
# `[lo[j], hi[j]]` (inclusive, 1-based, one range per feature dimension) and
# return each class' minimized formula for that sub-rectangle.

# Unchanged from before: this is "classic Lumen, restricted to a
# sub-rectangle" -- full materialization of that sub-rectangle, one
# `run_minimization` call per class over its WHOLE raw pool. Leaves are still
# genuinely SOP-optimal for their own sub-rectangle; the change in this file
# is only in how sibling leaves/subtrees get combined afterwards (see
# `_combine_cofactors` below), not in how a single leaf is computed.
# """
# # @inline function _eval_tree(node::SM.Branch, row::NTuple{N,T}, featidx::Dict{Symbol,Int}) where {N,T}
# #     while node isa SM.Branch
# #         cond = SL.value(SM.antecedent(node))            # ScalarCondition
# #         v    = SD.i_variable(SD.feature(cond))
# #         j    = v isa Integer ? Int(v) : featidx[Symbol(v)]
# #         node = SD.test_operator(cond)(row[j], SD.threshold(cond)) ?
# #             SM.posconsequent(node) : SM.negconsequent(node)
# #     end
# #     return SM.outcome(node)                             # ConstantModel leaf
# # end

# function checkcondition(
#     atom::Atom{B},
#     tbl::NamedTuple{T}
# )::BitVector where {B<:ScalarCondition,T}
#     cond = SL.value(atom)

#     cond_threshold = SD.threshold(cond)
#     cond_operator = SD.test_operator(cond)
#     cond_feature = SD.feature(cond)

#     col = SD.i_variable(cond_feature)

#     return cond_operator.(tbl[col], cond_threshold)
# end







# # ---------------------------------------------------------------------------- #
# #                      shared DecisionSet finalization                        #
# # ---------------------------------------------------------------------------- #
# """
#     _finalize_decision_set(ctx, per_class_terms, config) -> SM.DecisionSet

# Build the final `SM.DecisionSet` from per-class terms, dropping classes with
# zero terms and forcing the same "wide" element type `lumen()` uses for its
# `formulas` container (see the identical note in the sequential file's tail
# for why this widening is required to avoid a `DNF` vs
# `LeftmostDisjunctiveForm` dispatch ambiguity downstream).
# """
# function _finalize_decision_set(
#     ctx::Ctx{R,T},
#     per_class_terms::Vector{Vector{SM.LeftmostConjunctiveForm{<:SM.Atom}}},
#     config::LumenConfig{R,T}
# ) where {R,T<:AbstractFloat}
#     valid_mask = .!isempty.(per_class_terms)
#     classes = ctx.class_idxs[valid_mask]
#     formulas = per_class_terms[valid_mask]

#     rules = SM.Rule{R}[
#         SM.Rule(SL.LeftmostDisjunctiveForm(f), c)
#         for (f, c) in zip(formulas, classes)
#     ]

#     return SM.DecisionSet{R}(rules)::SM.DecisionSet{R}
# end

# map categorical labels onto unsigned integer of type `R`
# sort classlabels
function _assign(
    ::Type{R},
    y::Vector{S},
    sorted::Bool=false
) where {R<:Unsigned,S}
    labels = sorted ? sort!(unique(y)) : unique(y)
    dict = Dict{S,R}(v => i for (i, v) in enumerate(labels))
    return labels, [dict[t] for t in y]
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
function lumen_shannon(
    config::LumenShannonConfig{R,T},
    model::SM.DecisionEnsemble,
    featurenames::Vector{Symbol},
    feat_idxs::Vector{R},
    classnames::Vector{String},
    class_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    ensemble = LumenEnsemble(config, model)

    # thrs = ThresholdSpace(
    #     config,
    #     ensemble,
    #     featurenames,
    #     feat_idxs,
    #     classnames,
    #     class_idxs
    # )

    # hi = [length(t) + 1 for t in thrs.thresholds]
    # lo = ones(Int, length(hi))

    # cache = AtomCache(thrs)
    # # # per_class_terms = _extract(config, thrs_with_boundary, ensemble, lo, hi)

    # raw = _leaf_extract(config, thrs, cache, ensemble, lo, hi)

    # return cache, thrs, raw

    # return _finalize_decision_set(ctx, per_class_terms, config)
end

function lumen_shannon(
    config::LumenShannonConfig{R,T},
    model::SM.AbstractModel,
) where {R<:Unsigned,T<:AbstractFloat}
    featurenames, feat_idxs =
        _assign(R, unique!(SM.info(model, :featurenames)), false)
    classnames, class_idxs =
        _assign(R, unique!(SM.info(model, :supporting_labels)), true)

    lumen_shannon(
        config,
        model, featurenames, feat_idxs, String.(classnames), class_idxs)
end

