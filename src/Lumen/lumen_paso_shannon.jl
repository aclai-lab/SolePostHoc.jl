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

# ---------------------------------------------------------------------------- #
#                  the ONLY conversion to Sole types: DecisionSet              #
# ---------------------------------------------------------------------------- #
@inline _opfun(op::UInt8) =
    op == 0x01 ? (<) :
    op == 0x02 ? (>) :
    op == 0x03 ? (≥) : (≤)

function _toatom(a::LumenAtom{R,T}, featurenames::AbstractVector) where {R,T}
    feature = SD.VariableValue(Int(a.feat), Symbol(featurenames[a.feat]))
    return SL.Atom(SD.ScalarCondition(feature, _opfun(a.op), a.thr))
end

# an empty term is the tautology, and it swallows the whole disjunction
# function _antecedent(d::LumenDNF{R,T}, featurenames::AbstractVector) where {R,T}
#     any(isempty, d) && return SL.⊤
#     return SL.LeftmostDisjunctiveForm([
#         SL.LeftmostConjunctiveForm([_toatom(a, featurenames) for a in term])
#         for term in d
#     ])
# end

function _antecedent(
    dnf::LumenDNF{R,T},
    featurenames::AbstractVector,
) where {R,T}
    any(isempty, dnf) && return SL.⊤

    disjuncts = Vector{SyntaxStructure}(undef, length(dnf))

    for i in eachindex(dnf)
        disjuncts[i] = length(dnf[i]) == 1 ?
            _toatom(dnf[i][1], featurenames) :
            SL.LeftmostConjunctiveForm([_toatom(a, featurenames) for a in dnf[i]])
    end

    return SL.LeftmostDisjunctiveForm(disjuncts)
end

"""
    decision_set(terms, featurenames, classnames, order=eachindex(classnames))

Turn one `LumenDNF` per class into a `DecisionSet`, one rule per class with
at least one term, rules emitted in `order`. This is the single place where
Lumen atoms become `ScalarCondition`s.
"""
function decision_set(
    terms::Vector{LumenDNF{R,T}},
    featurenames::AbstractVector,
    classnames::AbstractVector,
    order::AbstractVector{O}=eachindex(classnames),
) where {R<:Unsigned,T<:AbstractFloat,O<:Integer}
    C = eltype(classnames)
    rules = SM.Rule[
        SM.Rule{C}(_antecedent(terms[c], featurenames), classnames[c])
        for c in order if !isempty(terms[c])
    ]
    return SM.DecisionSet{C}(rules)
end

# map categorical labels onto unsigned integer of type `R`
# sort classlabels
function assign(
    ::Type{R},
    y::Vector{S},
    sorted::Bool=false
) where {R<:Unsigned,S}
    labels = sorted ? sort!(unique(y)) : unique(y)
    dict = Dict{S,R}(v => i for (i, v) in enumerate(labels))
    return labels, [dict[t] for t in y]
end

# the model's feature names, column order; V<i> when the model carries none
function _featurenames(model::SM.AbstractModel, ensemble::LumenEnsemble)
    info = SM.info(model)
    haskey(info, :featurenames) && return info.featurenames
    nfeat = maximum(a.feat for a in ensemble.atoms if !isempty_atom(a); init=0)
    return [Symbol("V", i) for i in 1:nfeat]
end

# classic Lumen emits one rule per class in the order of the model's
# `supporting_labels`; keep that order when the model carries them, and fall
# back to class order otherwise. Never mutates the model's info.
function _ruleorder(model::SM.AbstractModel, classnames::AbstractVector)
    info = SM.info(model)
    haskey(info, :supporting_labels) || return collect(eachindex(classnames))
    seen = unique(_label.(info.supporting_labels))
    pos = [something(findfirst(==(c), seen), typemax(Int)) for c in classnames]
    return sortperm(pos)          # stable: unlisted classes keep class order
end

# ---------------------------------------------------------------------------- #
#                              lumen paso shannon                              #
# ---------------------------------------------------------------------------- #
function lumen_paso_shannon(
    config::LumenConfig{R,T},
    model::SM.DecisionEnsemble,
    featurenames::Vector{Symbol},
    feat_idxs::Vector{R},
    classnames::Vector{String},
    class_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    ensemble = LumenEnsemble(config, model)

    thrs = ThresholdSpace(ensemble, feat_idxs, class_idxs, config.depth)

    hi = thrs.nlev
    lo = ones(R, length(hi))

    terms = _extract(config, thrs, ensemble, lo, hi)

    # leaves were minimised in isolation: a term of one leaf may be contained
    # in, or adjacent to, a term of a neighbouring leaf. Dropping the former
    # and merging the latter never changes coverage; they are the only
    # cross-leaf compactions that are safe without re-minimising.
    foreach(compact!, terms)

    return decision_set(
        terms, featurenames, categorical(classnames), _ruleorder(model, classnames))
end

function lumen_paso_shannon(
    config::LumenConfig{R,T},
    model::SM.AbstractModel,
) where {R<:Unsigned,T<:AbstractFloat}
    featurenames, feat_idxs =
        assign(R, unique!(SM.info(model, :featurenames)), false)
    classnames, class_idxs =
        assign(R, unique!(SM.info(model, :supporting_labels)), true)

    lumen_paso_shannon(
        config, model, featurenames, feat_idxs, String.(classnames), class_idxs)
end