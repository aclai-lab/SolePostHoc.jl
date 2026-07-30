# ---------------------------------------------------------------------------- #
#                           Intrees pruning utility                            #
# ---------------------------------------------------------------------------- #
"""
    _prune_rule(::Type{<:Atom}, r::Rule, args...; kwargs...) -> Rule

Wrap a single-atom antecedent in a `LeftmostConjunctiveForm` and delegate to
the conjunctive-form pruning method.
"""
function _prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where {O}
    r = Rule(LeftmostConjunctiveForm([antecedent(r)]), consequent(r), info(r))
    _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
end

"""
    _prune_rule(
        ::Type{<:LeftmostConjunctiveForm},
        r::Rule,
        X::AbstractInterpretationSet,
        y::AbstractVector{<:SoleModels.Label};
        s::AbstractFloat,
        decay_threshold::AbstractFloat,
        kwargs...
    ) -> Rule

Greedily remove conjuncts from `r`'s antecedent whose removal does not
increase the rule's error by more than `decay_threshold` (normalized by
`s`), as in the original InTrees pruning procedure.
"""
function _prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r::Rule{O},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label};
    decay_threshold::Float32,
    percentage_degradation::Bool,
    s::Float32,
    kwargs...,
) where {O}
    nruleconjuncts = SoleModels.nconjuncts(r)
    e_zero = SoleModels.rulemetrics(r, X, y)[:error]
    valid_idxs = 1:nruleconjuncts
    antd, cons = SoleModels.antecedent(r), SoleModels.consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) < 2) && break

        # indices to be considered to evaluate the rule
        other_idxs = setdiff(valid_idxs, idx)
        rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

        # return error of the rule without idx-th pair
        e_minus_i = SoleModels.rulemetrics(rule, X, y)[:error]
        decay_i = percentage_degradation ?
            (e_minus_i - e_zero) / max(e_zero, s) :
            (e_minus_i - e_zero)

        if decay_i ≤ decay_threshold
            # remove the idx-th pair in the vector of decisions
            valid_idxs = setdiff(valid_idxs, idx)
            e_zero = e_minus_i
        end
    end

    return Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]), cons)
end

"""
    _prune_rule(::Type{<:MultiFormula}, r::Rule, args...; kwargs...) -> Rule

Flatten a `MultiFormula` antecedent into its per-modality children and
delegate to the conjunctive-form pruning method.
"""
function _prune_rule(::Type{<:MultiFormula}, r::Rule{O}, args...; kwargs...) where {O}
    @assert antecedent(r) isa MultiFormula "Cannot use this function on $(antecedent(r))"
    children = [
        MultiFormula(i_modality, modant) for (i_modality, modant) in modforms(antecedent(r))
    ]

    return length(children) < 2 ? r : begin
        r = Rule(LeftmostConjunctiveForm(children), consequent(r), info(r))
        _prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
    end
end

"""
    _prune_ruleset(
        ruleset::Vector{<:Rule},
        X::AbstractInterpretationSet,
        y::AbstractVector{<:SoleModels.Label},
        config::InTreesConfig
    ) -> Vector{<:Rule}

Prune every rule in `ruleset` in parallel via [`_prune_rule`](@ref), using the
pruning parameters stored in `config`. Rules with a `BooleanTruth` antecedent
(e.g. produced by XGBoost) are kept unchanged.
"""
function _prune_ruleset(
    ruleset::Vector{<:Rule},
    X::AbstractInterpretationSet,
    y::AbstractVector{<:SoleModels.Label},
    config::InTreesConfig
)
    pruned = similar(ruleset)

    @inbounds Threads.@threads for i in eachindex(ruleset)
        pruned[i] = if ruleset[i].antecedent isa SoleLogics.BooleanTruth
            ruleset[i]  # keep BooleanTruth rules as-is (e.g., from XGBoost)
        else
            _prune_rule(
                typeof(antecedent(ruleset[i])), ruleset[i], X, y;
                decay_threshold=get_decay_threshold(config),
                percentage_degradation=get_percentage_degradation(config),
                s=get_s(config)
            )
        end
    end

    return pruned
end