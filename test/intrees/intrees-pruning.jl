
    """
        intrees_prunerule(rc::Rule)::Rule

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`Rule`](@ref),
    [`rulemetrics`](@ref).
    """
    intrees_prunerule(r::Rule, X, y; kwargs...) = _prune_rule(typeof(antecedent(r)), r, X, y; kwargs...)
    
    function _prune_rule(
        ::Type{<:Atom},
        r::Rule{O},
        args...; kwargs...
    ) where {O}
      intrees_prunerule(Rule(LeftmostConjunctiveForm([antecedent(r)]),consequent(r),info(r)), args...; kwargs...)
    end
    
    function _prune_rule(
        ::Type{<:LeftmostConjunctiveForm},
        r::Rule{O},
        X, y; pruning_s, pruning_decay_threshold, kwargs...
    ) where {O}
        # @show r
        nruleconjuncts = nconjuncts(r)
        E_zero = rulemetrics(r, X, y)[:error]
        valid_idxs = collect(1:nruleconjuncts)
        antd, cons = antecedent(r), consequent(r)
        for idx in reverse(valid_idxs)
            (length(valid_idxs) < 2) && break

            # Indices to be considered to evaluate the rule
            other_idxs = intersect!(vcat(1:(idx-1),(idx+1):nruleconjuncts),valid_idxs)
            rule = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]), cons)

            # Return error of the rule without idx-th pair
            E_minus_i = rulemetrics(rule, X, y)[:error]

            decay_i = (E_minus_i - E_zero) / max(E_zero, pruning_s)

            if decay_i <= pruning_decay_threshold
                # Remove the idx-th pair in the vector of decisions
                deleteat!(valid_idxs, idx)
                E_zero = E_minus_i
            end
        end
        prunedr = Rule(LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]), cons)
        # @show prunedr
        return prunedr
    end

    function _prune_rule(
        ::Type{<:MultiFormula},
        r::Rule{O},
        args...; kwargs...
    ) where {O}
        @assert antecedent(r) isa MultiFormula "Cannot use this function on $(antecedent(r))"
        children = [
            MultiFormula(i_modality,modant)
            for (i_modality,modant) in modforms(antecedent(r))
        ]

        return length(children) < 2 ?
            r :
            intrees_prunerule(Rule(LeftmostConjunctiveForm(children),consequent(r),info(r)), args...; kwargs...)
    end
