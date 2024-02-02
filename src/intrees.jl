# using DecisionTree
using ComplexityMeasures

############################################################################################
# Rule extraction from random forest
############################################################################################

"""
    intrees(args...; kwargs...)::DecisionList

Extract rules from a model, reduces the length of each rule (number of variable value
pairs), and applies a sequential coverage approach to obtain a set of relevant and
non-redundant rules

# Arguments
- `model::Union{AbstractModel,DecisionForest}`: input model
- `X::Any`: dataset
- `Y::AbstractVector{<:Label}`: label vector

# Keywords
- `prune_rules::Bool=true`: access to prune or not
- `pruning_s::Float=nothing`: parameter that limits the denominator in the pruning metric calculation
- `pruning_decay_threshold::Float=nothing`: threshold used in pruning to remove or not a joint from the rule
- `method_rule_selection::Symbol=:CBC`: method of rule selection
- `accuracy_rule_selection::Float=nothing`: percentage of rules that rule selection must follow
- `min_frequency::Float=nothing`: minimum frequency that the rules must have at the beginning of the definition of the learner

# Returns
- `DecisionList`: decision list that represent a new learner

TODO cite paper and specify that the method is for forests, but was extended to work with any other model

See also
[`AbstractModel`](@ref),
[`DecisionForest`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
# Extract rules from a forest, with respect to a dataset
function intrees(
    model::Union{AbstractModel,DecisionForest},
    X,
    Y::AbstractVector{<:Label};
    #
    prune_rules = true,
    pruning_s = nothing,
    pruning_decay_threshold = nothing,
    #
    method_rule_selection = :CBC,
    accuracy_rule_selection = nothing,
    min_frequency = nothing,
    #
    rng::AbstractRNG = MersenneTwister(1),
    kwargs...,
)
    isnothing(pruning_s) && !isnothing(pruning_decay_threshold) && (prune_rules = false)
    isnothing(pruning_decay_threshold) && !isnothing(pruning_s) && (prune_rules = false)
    isnothing(pruning_s) && (pruning_s = 1.0e-6)
    isnothing(pruning_decay_threshold) && (pruning_decay_threshold = 0.05)
    isnothing(accuracy_rule_selection) && (accuracy_rule_selection = 0.0)
    isnothing(min_frequency) && (min_frequency = 0.01)

    @assert model isa DecisionForest || SoleModels.issymbolic(model) "Cannot extract rules for model of type $(typeof(model))."

    """
        prune_rule(rc::Rule)::Rule

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`Rule`](@ref),
    [`rulemetrics`](@ref).
    """
    function prune_rule(
        r::Rule{O,<:LeftmostConjunctiveForm} # TODO use nchildren..?
    ) where {O}
        nruleconjuncts = nconjuncts(r)
        E_zero = rulemetrics(r,X,Y)[:error]
        valid_idxs = collect(1:nruleconjuncts)

        for idx in reverse(valid_idxs)
            (length(valid_idxs) < 2) && break

            # Indices to be considered to evaluate the rule
            other_idxs = intersect!(vcat(1:(idx-1),(idx+1):nruleconjuncts),valid_idxs)
            rule = r[other_idxs]

            # Return error of the rule without idx-th pair
            E_minus_i = rulemetrics(rule,X,Y)[:error]

            decay_i = (E_minus_i - E_zero) / max(E_zero, pruning_s)

            if decay_i <= pruning_decay_threshold
                # Remove the idx-th pair in the vector of decisions
                deleteat!(valid_idxs, idx)
                E_zero = E_minus_i
            end
        end

        return r[valid_idxs]
    end

    function prune_rule(
        r::Rule{O,<:MultiFormula}
    ) where {O}
        children = [
            MultiFormula(i_modality,modant)
            for (i_modality,modant) in modforms(antecedent(r))
        ]

        return length(children) < 2 ?
            r :
            prune_rule(Rule(LeftmostConjunctiveForm(children),consequent(r),info(r)))
    end

    """
        cfs()

    Prunes redundant or irrelevant conjuncts of the antecedent of the input rule cascade
    considering the error metric

    See also
    [`Rule`](@ref),
    [`rulemetrics`](@ref).
    """
    function cfs(
        X,
        Y::Vector{<:Label},
    )
        entropyd(x) = ComplexityMeasures.entropy(probabilities(x))
        midd(x, y) = -entropyd(collect(zip(x, y)))+entropyd(x)+entropyd(y)
        information_gain(f1, f2) = entropyd(f1) - (entropyd(f1) - midd(f1, f2))
        su(f1, f2) = (2.0 * information_gain(f1, f2) / (entropyd(f1) + entropyd(f2)))
        function merit_calculation(X, Y::Vector{<:Label})
            n_samples, n_features = size(X)
            rff = 0
            rcf = 0

            for i in collect(1:n_features)
                fi = X[:, i]
                rcf += su(fi, Y)  # su is the symmetrical uncertainty of fi and y
                for j in collect(1:n_features)
                    if j > i
                        fj = X[:, j]
                        rff += su(fi, fj)
                    end
                end
            end

            rff *= 2
            merits = rcf / sqrt(n_features + rff)
        end

        n_samples, n_features = size(X)
        F = [] # vector of returned indices samples
        M = [] # vector which stores the merit values

        while true
            merit = -100000000000
            idx = -1

            for i in collect(1:n_features)
                if i ∉ F
                    append!(F,i)
                    idxs_column = F[findall(F .> 0)]
                    t = merit_calculation(X[:,idxs_column],Y)

                    if t > merit
                        merit = t
                        idx = i
                    end

                    pop!(F)
                end
            end

            append!(F,idx)
            append!(M,merit)

            (length(M) > 5) &&
                (M[length(M)] <= M[length(M)-1]) &&
                    (M[length(M)-1] <= M[length(M)-2]) &&
                        (M[length(M)-2] <= M[length(M)-3]) &&
                            (M[length(M)-3] <= M[length(M)-4]) && break
        end

        valid_idxs = findall(F .> 0)
        @show F
        @show valid_idxs
        return F[valid_idxs]
    end

    function selectRuleRRF(
        matrixrulemetrics::Matrix{Float64},
        X,
        Y::Vector{<:Label},
    )
        #coefReg = 0.95 .- (0.01*matrixrulemetrics[:,3]/max(matrixrulemetrics[:,3]...))
        #@show coefReg
        rf = DecisionTree.build_forest(Y,X,2,50,0.7,-1; rng=rng)
        imp = begin
            #importance = impurity_importance(rf, coefReg)
            importance = DecisionTree.impurity_importance(rf)
            importance/max(importance...)
        end
        @show imp

        feaSet = findall(imp .> 0.01)
        @show feaSet
        ruleSetPrunedRRF = hcat(matrixrulemetrics[feaSet,:],imp[feaSet],feaSet)

        finalmatrix = sortslices(ruleSetPrunedRRF, dims=1, by=x->(x[4],x[2],x[3]), rev=true)

        return Integer.(finalmatrix[:,5])
    end

    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    println("Rule extraction in...")
    ruleset = @time begin
        if model isa DecisionForest
            rs = unique([listrules(tree; use_shortforms=true) for tree in trees(model)])
            # TODO maybe also sort?
            rs isa Vector{<:Vector{<:Any}} ? reduce(vcat,rs) : rs
        else
            listrules(model; use_shortforms=true)
        end
    end
    ########################################################################################

    ########################################################################################
    # Prune rules with respect to a dataset
    ########################################################################################
    prune_rules ? (println("Pruning phase in...")) : (println("Skipped pruning in..."))
    ruleset = @time begin
        if prune_rules
            afterpruningruleset = Vector{Rule}(undef, length(ruleset))
            Threads.@threads for (i,r) in collect(enumerate(ruleset))
                afterpruningruleset[i] = prune_rule(r)
            end
            #ruleset = @time prune_rule.(ruleset)
            afterpruningruleset
        else
            ruleset
        end
    end

    println("# Rules to checking: $(length(ruleset))")

    ########################################################################################
    # Rule selection to obtain the best rules
    ########################################################################################
    println("Selection phase with $(method_rule_selection) in...")
    best_rules = @time begin
        if method_rule_selection == :CBC
            matrixrulemetrics = Matrix{Float64}(undef,length(ruleset),3)
            afterselectionruleset = Vector{BitVector}(undef, length(ruleset))
            Threads.@threads for (i,rule) in collect(enumerate(ruleset))
                eval_result = rulemetrics(rule, X, Y)
                afterselectionruleset[i] = eval_result[:antsat]
                matrixrulemetrics[i,1] = eval_result[:support]
                matrixrulemetrics[i,2] = eval_result[:error]
                matrixrulemetrics[i,3] = eval_result[:length]
            end
            #M = hcat([evaluaterule(rule, X, Y)[:antsat] for rule in ruleset]...)
            M = hcat(afterselectionruleset...)

            #best_idxs = findcorrelation(cor(M); threshold = accuracy_rule_selection)
            #best_idxs = cfs(M,Y)
            best_idxs = selectRuleRRF(matrixrulemetrics,M,Y)
            @show best_idxs
            ruleset[best_idxs]
        else
            error("Unexpected method specified: $(method)")
        end
    end
    ########################################################################################

    println("# Rules to checking: $(length(best_rules))")
    ########################################################################################
    # Construct a rule-based model from the set of best rules
    ########################################################################################
    println("STEL phase starting, end soon")

    D = deepcopy(X) # Copy of the original dataset
    L = deepcopy(Y)
    R = Rule[]      # Ordered rule list
    # S is the vector of rules left
    S = [deepcopy(best_rules)..., Rule(bestguess(Y; suppress_parity_warning = true))]

    # Rules with a frequency less than min_frequency
    S = begin
        rules_support = Vector{AbstractFloat}(undef,length(S))
        Threads.@threads for (i,s) in collect(enumerate(S))
            rules_support[i] = rulemetrics(s,X,Y)[:support]
        end

        idxs_undeleted = findall(rules_support .>= min_frequency)
        S[idxs_undeleted]
    end

    println("# rules in S: $(length(S))")
    println("# rules in R: $(length(R))")

    while true
        # Metrics update based on remaining instances
        rules_support = Vector{AbstractFloat}(undef,length(S))
        rules_error = Vector{AbstractFloat}(undef,length(S))
        rules_length = Vector{Integer}(undef,length(S))

        Threads.@threads for (i,s) in collect(enumerate(S))
            metrics = rulemetrics(s,D,L)
            rules_support[i] = metrics[:support]
            rules_error[i] = metrics[:error]
            rules_length[i] = metrics[:length]
        end

        println("Rules error:")
        println(rules_error)
        println("Minimo: $(min(rules_error...))")

        #metrics = [rulemetrics(s,D,Y) for s in S]
        #println(metrics[1])
        #rules_support = [metrics[i][:support] for i in eachindex(metrics)]
        #rules_error = [metrics[i][:error] for i in eachindex(metrics)]
        #rules_length = [metrics[i][:length] for i in eachindex(metrics)]

        # Best rule index
        idx_best = begin
            idx_best = nothing
            # First: find the rule with minimum error
            idx = findall(rules_error .== min(rules_error...))
            (length(idx) == 1) && (idx_best = idx)
            println("Idxs min error: $(idx)")

            # If not one, find the rule with maximum frequency
            if isnothing(idx_best)
                idx_support = findall(rules_support .== max(rules_support[idx]...))
                (length(intersect!(idx, idx_support)) == 1) && (idx_best = idx)
            end
            println("Idxs min support: $(idx)")

            # If not one, find the rule with minimum length
            if isnothing(idx_best)
                idx_length = findall(rules_length .== min(rules_length[idx]...))
                (length(intersect!(idx, idx_length)) == 1) && (idx_best = idx)
            end
            println("Idxs min length: $(idx)")

            # Final case: more than one rule with minimum length
            # Randomly choose a rule
            isnothing(idx_best) && (idx_best = rand(rng,idx))

            length(idx_best) > 1 && error("More than one best indexes")
            println("Idx best: $(idx_best)")
            idx_best[1]
        end
        println("Idx best: $(idx_best)")

        # Add at the end the best rule
        push!(R, S[idx_best])
        println("# rules in R: $(length(R))")

        # Indices of the remaining instances
        idx_remaining = begin
            sat_unsat = evaluaterule(S[idx_best], D, L)[:antsat]
            # Remain in D the rule that not satisfying the best rule'pruning_s condition
            findall(sat_unsat .== false)
        end
        D = length(idx_remaining) > 0 ?
                slicedataset(D, idx_remaining; return_view=true) : nothing
        L = L[idx_remaining]

        if idx_best == length(S)
            return DecisionList(R[1:end-1],consequent(R[end]))
        elseif isnothing(D) || ninstances(D) == 0
            return DecisionList(R,bestguess(Y; suppress_parity_warning = true))
        end

        # Delete the best rule from S
        deleteat!(S,idx_best)
        # Update of the default rule
        S[end] = Rule(bestguess(L; suppress_parity_warning = true))
    end

    return error("Unexpected error in intrees!")
end
