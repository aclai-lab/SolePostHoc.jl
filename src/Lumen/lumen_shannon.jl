# ============================================================================ #
#   LUMEN SHANNON: recursive exact cofactor decomposition, SOP-optimal leaves  #
#                                                                              #
#   Same design as the earlier version of this file (see its history for the   #
#   full rationale), with ONE change: there is no `remin` knob anymore.        #
#   Internal nodes ALWAYS combine their two cofactors by plain set union,      #
#   never by feeding the children's terms back into the minimizer.             #
#                                                                              #
#   Why: re-minimizing a union of two cofactors' terms requires converting     #
#   each term back into a "cube" (`_term_to_cube`), and a term that dropped    #
#   the literal for the SPLIT feature is only a valid don't-care WITHIN the    #
#   sub-rectangle that leaf/subtree was restricted to, not globally. Mixing    #
#   such a term with the sibling cofactor's terms and re-minimizing let the    #
#   minimizer treat that dropped literal as a genuine global don't-care,       #
#   which can make it discard other terms as "redundant" when they weren't,    #
#   silently losing coverage (this showed up empirically as a sensitivity      #
#   drop concentrated on rarer classes, worse the more re-minimization         #
#   happened along the path to the root, confirmed identical on both :abc      #
#   and :mitespresso, ruling out a backend-specific bug). See                  #
#   `lumen_shannon_boundary_fix.jl` for a version that repairs this (by        #
#   re-attaching the correct boundary literal before any re-minimization)      #
#   instead of simply avoiding it.                                             #
#                                                                              #
#   Plain union has none of that risk: each cofactor's terms are already a     #
#   fully correct description of their own (disjoint) sub-rectangle, and the   #
#   union of two correct, disjoint partitions is always correct. The only      #
#   cost is a larger (less compact) final formula -- confirmed empirically to  #
#   match classic Lumen's Sensitivity/Specificity/DecisionSetAccuracy exactly  #
#   on the `car` dataset, at the cost of somewhat more atoms/terms.            #
#                                                                              #
#   Relies on the same shared primitives as before:                            #
#   `_prepare_sequential_context`, `generate_disjunct`,                        #
#   `_truths_by_thresholds`, `run_minimization`, `PropositionalLogiset`,       #
#   `get_apply_function`, `get_use_multithreads`,                              #
#   `get_float_type`, `get_minimization_scheme`.                               #
# ============================================================================ #

# ---------------------------------------------------------------------------- #
#                                   operators                                  #
# ---------------------------------------------------------------------------- #
#     _opfun(op::UInt8) -> Function

# The comparison function behind an operator byte, the inverse of
# `evalop(::Function)`: `0x01` → `<`, `0x02` → `>`, `0x03` → `≥`, else `≤`.
@inline _opfun(op::UInt8) =
    op == 0x01 ? (<) :
    op == 0x02 ? (>) :
    op == 0x03 ? (≥) : (≤)

#     _toatom(a::LumenAtom, featurenames) -> SoleLogics.Atom

# Convert `a` to a SoleLogics `Atom` wrapping the `ScalarCondition`
# `V[a.feat] op a.thr`, the feature named `featurenames[a.feat]`.
function _toatom(a::LumenAtom{R,T}, featurenames::AbstractVector) where {R,T}
    feature = SD.VariableValue(Int(a.feat), Symbol(featurenames[a.feat]))
    return SL.Atom(SD.ScalarCondition(feature, _opfun(a.op), a.thr))
end

#     _antecedent(dnf::LumenDNF, featurenames) -> SyntaxStructure

# The SoleLogics formula of `dnf`: a `LeftmostDisjunctiveForm` with one
# disjunct per term, each a single `Atom` or a `LeftmostConjunctiveForm` of
# atoms (see [`_toatom`](@ref)). If any term is empty (the tautology) the
# whole formula is `⊤`.
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

#     _featurenames(model, ensemble::LumenEnsemble) -> AbstractVector

# The model's feature names in column order, taken from its `info`. When the
# model carries none, `[:V1, …, :Vn]` up to the largest feature index used by
# `ensemble`.
function _featurenames(model::SM.AbstractModel, ensemble::LumenEnsemble)
    info = SM.info(model)
    haskey(info, :featurenames) && return info.featurenames
    nfeat = maximum(a.feat for a in ensemble.atoms if !isempty_atom(a); init=0)
    return [Symbol("V", i) for i in 1:nfeat]
end

#     _ruleorder(model, classnames) -> Vector{Int}

# The order in which [`decision_set`](@ref) emits the rules, as indices into
# `classnames`. Classic Lumen emits one rule per class in order of first
# appearance in the model's `supporting_labels`; that order is kept when the
# model carries them, and class order is used otherwise. Classes not listed
# there come last, in class order. Never mutates the model's info.
function _ruleorder(model::SM.AbstractModel, classnames::AbstractVector)
    info = SM.info(model)
    haskey(info, :supporting_labels) || return collect(eachindex(classnames))
    seen = unique(_label.(info.supporting_labels))
    pos = [something(findfirst(==(c), seen), typemax(Int)) for c in classnames]
    return sortperm(pos) # stable: unlisted classes keep class order
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

"""
    assign(R::Type{<:Unsigned}, y::Vector, sorted::Bool=false) -> (labels, codes)

Map the values of `y` onto unsigned integer codes of type `R`. `labels` are
the distinct values of `y`, in order of first appearance or sorted when
`sorted` is `true`; `codes[i]` is the position of `y[i]` in `labels`, so
`labels[codes] == y`.

# Examples
```julia
assign(UInt32, ["b", "a", "b"])        # (["b", "a"], UInt32[1, 2, 1])
assign(UInt32, ["b", "a", "b"], true)  # (["a", "b"], UInt32[2, 1, 2])
```
"""
function assign(
    ::Type{R},
    y::Vector{S},
    sorted::Bool=false
) where {R<:Unsigned,S}
    labels = sorted ? sort!(unique(y)) : unique(y)
    dict = Dict{S,R}(v => i for (i, v) in enumerate(labels))
    return labels, [dict[t] for t in y]
end

# ---------------------------------------------------------------------------- #
#                              lumen paso shannon                              #
# ---------------------------------------------------------------------------- #
"""
    lumen_shannon(config::LumenConfig, model::AbstractModel) -> DecisionSet
    lumen_shannon(config, model::DecisionEnsemble, featurenames, feat_idxs,
                  classnames, class_idxs) -> DecisionSet

Extract from a tree ensemble an equivalent `DecisionSet` with one DNF rule
per predicted class, by exhaustive exploration of the input space the
model's thresholds induce, split into pieces of bounded size.

The steps are:
1. compile `model` into a [`LumenEnsemble`](@ref);
2. build the [`ThresholdSpace`](@ref) grid of its thresholds (only the
   BFS-first `config.depth` fraction of each tree when `depth < 1`);
3. recursively halve the grid until a piece holds at most `config.M` points,
   classify one sample per box of each piece and minimise each class' cubes
   with the backend `MS` of `config` (`Abc` or `MitEspresso`, at the setup
   level `config.p_minimization_setup`), see [`_extract`](@ref);
4. union the pieces' terms per class and [`compact!`](@ref) them;
5. convert the result to Sole types with [`decision_set`](@ref).

Pieces are disjoint and never re-minimised together, so the rules describe
the model exactly on the grid (with `depth = 1`) but are not globally
minimal. `config.M` trades memory and backend-call size against the number
of pieces. Classes the model never predicts get no rule.

The first form reads the feature names (`:featurenames`) and class labels
(`:supporting_labels`) from the model's `info`; class labels are sorted.
The second form takes them explicitly: `feat_idxs` and `class_idxs` give the
index of each feature and class position, as returned by [`assign`](@ref).
Only ensembles of trees with categorical outcomes and scalar conditions are
supported.

# Examples
```julia
config = LumenConfig(; M=4096, p_minimization_setup=Fast)
ds = lumen_shannon(config, forest)
```

See also: [`LumenConfig`](@ref), [`lumen`](@ref).
"""
function lumen_shannon(
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

function lumen_shannon(
    config::LumenConfig{R,T},
    model::SM.AbstractModel,
) where {R<:Unsigned,T<:AbstractFloat}
    featurenames, feat_idxs =
        assign(R, unique!(SM.info(model, :featurenames)), false)
    classnames, class_idxs =
        assign(R, unique!(SM.info(model, :supporting_labels)), true)

    lumen_shannon(
        config, model, featurenames, feat_idxs, String.(classnames), class_idxs)
end