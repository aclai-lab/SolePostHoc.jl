
"""
    intrees_prunerule(rc::Rule)::Rule

Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
considering the error metric

See also
[`Rule`](@ref),
[`rulemetrics`](@ref).
"""
@inline intrees_prunerule(
    r::Rule,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label};
    kwargs...
) = _prune_rule(typeof(antecedent(r)), r, X, y; kwargs...)

@inline _prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where{O} =
    intrees_prunerule(Rule(LeftmostConjunctiveForm(
        [antecedent(r)]), consequent(r), info(r)), args...; kwargs...,)


function _prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r                 :: Rule{O},
    X                 :: AbstractInterpretationSet,
    y                 :: AbstractVector{<:SoleModels.Label};
    pruning_s         :: AbstractFloat,
    pruning_decay_thr :: AbstractFloat,
    kwargs...,
) where {O}
    nruleconjuncts = nconjuncts(r)
    e_zero         = rulemetrics(r, X, y)[:error]
    valid_idxs     = 1:nruleconjuncts
    antd, cons     = antecedent(r), consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) < 2) && break

        # Indices to be considered to evaluate the rule
        other_idxs = intersect!(vcat(1:(idx-1), (idx+1):nruleconjuncts), valid_idxs)
        rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

        # Return error of the rule without idx-th pair
        e_minus_i = rulemetrics(rule, X, y)[:error]
        decay_i   = (e_minus_i - e_zero) / max(e_zero, pruning_s)

        if decay_i ≤ pruning_decay_thr 
             # Remove the idx-th pair in the vector of decisions
            deleteat!([valid_idxs...], idx)
            e_zero = e_minus_i
        end
    end

    return Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]), cons)
end

function _prune_rule(::Type{<:MultiFormula}, r::Rule{O}, args...; kwargs...) where {O}
    @assert antecedent(r) isa MultiFormula "Cannot use this function on $(antecedent(r))"
    children = [
        MultiFormula(i_modality, modant) for (i_modality, modant) in modforms(antecedent(r))
    ]

    return  length(children) < 2 ? r :
            intrees_prunerule(
                Rule(LeftmostConjunctiveForm(children), consequent(r), info(r)),
                args...;
                kwargs...,
            )
end
