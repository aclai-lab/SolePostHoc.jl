import DecisionTree as DT
using ComplexityMeasures
using SoleModels

include("intrees-pruning.jl")

############################################################################################
# Rule extraction from random forest
############################################################################################

"""
Deng, Houtao. "Interpreting tree ensembles with intrees." International Journal of Data Science and Analytics 7.4 (2019): 277-287.

See also [`extractrules`](@ref), [`intrees`](@ref), [`RuleExtractor`](@ref).
"""
struct InTreesRuleExtractor <: RuleExtractor end

function extractrules(::InTreesRuleExtractor, m, args...; kwargs...)
  dl = intrees(m, args...; kwargs...)
  return listrules(dl)
end

"""
    intrees(model::Union{AbstractModel,DecisionForest}, X, y::AbstractVector{<:Label}; kwargs...)::DecisionList

Return a decision list which approximates the behavior of the input `model` on the specified supervised dataset.
The set of relevant and
non-redundant rules in the decision list are obtained by means of rule selection, rule pruning,
and sequential covering (STEL).

# Keyword Arguments
- `prune_rules::Bool=true`: access to prune or not
- `pruning_s::Float=nothing`: parameter that limits the denominator in the pruning metric calculation
- `pruning_decay_threshold::Float=nothing`: threshold used in pruning to remove or not a joint from the rule
- `rule_selection_method::Symbol=:CBC`: rule selection method. Currently only supports `:CBC`.
- `min_coverage::Float=nothing`: minimum rule coverage for STEL

TODO cite paper and specify that the method is for forests, but was extended to work with any other models.
Reference: Deng, Houtao. "Interpreting tree ensembles with intrees." International Journal of Data Science and Analytics 7.4 (2019): 277-287.

See also
[`AbstractModel`](@ref),
[`DecisionForest`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
# Extract rules from a forest, with respect to a dataset
function intrees(
    model,
    X,
    y::AbstractVector{<:Label};
    #
    prune_rules = true,
    pruning_s = nothing,
    pruning_decay_threshold = nothing,
    #
    rule_selection_method = :CBC,
    rule_length_metric = :natoms,
# - `accuracy_rule_selection::Float=nothing`: percentage of rules that rule selection must follow
    # accuracy_rule_selection = nothing,
    min_coverage = nothing,
    silent = false,
    #
    rng::AbstractRNG = MersenneTwister(1),
    return_info::Bool = false,
    # kwargs...,
)
    isnothing(pruning_s) && !isnothing(pruning_decay_threshold) && (prune_rules = false)
    isnothing(pruning_decay_threshold) && !isnothing(pruning_s) && (prune_rules = false)
    isnothing(pruning_s) && (pruning_s = 1.0e-6)
    isnothing(pruning_decay_threshold) && (pruning_decay_threshold = 0.05)
    # isnothing(accuracy_rule_selection) && (accuracy_rule_selection = 0.0)
    isnothing(min_coverage) && (min_coverage = 0.01)

    if !(X isa AbstractInterpretationSet)
        X = SoleData.scalarlogiset(X; silent, allow_propositional = true)
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
        y::AbstractVector{<:Label},
    )
        entropyd(_x) = ComplexityMeasures.entropy(probabilities(_x))
        midd(_x, _y) = -entropyd(collect(zip(_x, _y)))+entropyd(_x)+entropyd(_y)
        information_gain(f1, f2) = entropyd(f1) - (entropyd(f1) - midd(f1, f2))
        su(f1, f2) = (2.0 * information_gain(f1, f2) / (entropyd(f1) + entropyd(f2)))
        function merit_calculation(X, y::AbstractVector{<:Label})
            n_samples, n_features = size(X)
            rff = 0
            rcf = 0

            for i in collect(1:n_features)
                fi = X[:, i]
                rcf += su(fi, y)  # su is the symmetrical uncertainty of fi and y
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
                    t = merit_calculation(X[:,idxs_column],y)

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
        # @show F
        # @show valid_idxs
        return F[valid_idxs]
    end

    function starterruleset(model; kwargs...)
        unique(reduce(vcat, [listrules(subm; kwargs...) for subm in SoleModels.models(model)]))
        # TODO maybe also sort?
    end

    if !haslistrules(model)
        model = solemodel(model)
    end

    info = (;)

    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    silent || println("Extracting starting rules...")
    listrules_kwargs = (;
        use_shortforms=true,
        # flip_atoms = true,
        normalize = true,
        # normalize_kwargs = (; forced_negation_removal = true, reduce_negations = true, allow_atom_flipping = true, rotate_commutatives = false)
    )
    ruleset = isensemble(model) ? starterruleset(model; listrules_kwargs...) : listrules(model; listrules_kwargs...)

    ########################################################################################
    # Prune rules with respect to a dataset
    ########################################################################################
    if prune_rules
        silent || println("Pruning $(length(ruleset)) rules...")
        if return_info
            info = merge(info, (; unpruned_ruleset = ruleset))
        end
        ruleset = @time begin
            afterpruningruleset = Vector{Rule}(undef, length(ruleset))
            Threads.@threads for (i,r) in collect(enumerate(ruleset))
                afterpruningruleset[i] = intrees_prunerule(r, X, y; pruning_s, pruning_decay_threshold)
            end
            afterpruningruleset
        end
    end

    ########################################################################################
    # Rule selection to obtain the best rules
    ########################################################################################
    silent || println("Selecting via $(string(rule_selection_method)) from a pool of $(length(ruleset)) rules...")
    ruleset = @time begin
        if return_info
            info = merge(info, (; unselected_ruleset = ruleset))
        end
        if rule_selection_method == :CBC
            matrixrulemetrics = Matrix{Float64}(undef,length(ruleset),3)
            afterselectionruleset = Vector{BitVector}(undef, length(ruleset))
            Threads.@threads for (i,rule) in collect(enumerate(ruleset))
                eval_result = rulemetrics(rule, X, y)
                afterselectionruleset[i] = eval_result[:checkmask,]
                matrixrulemetrics[i,1] = eval_result[:coverage]
                matrixrulemetrics[i,2] = eval_result[:error]
                matrixrulemetrics[i,3] = eval_result[rule_length_metric]
            end
            #M = hcat([evaluaterule(rule, X, y)[:checkmask,] for rule in ruleset]...)
            M = hcat(afterselectionruleset...)

            #best_idxs = findcorrelation(cor(M); threshold = accuracy_rule_selection)
            #best_idxs = cfs(M,y)
            
            #coefReg = 0.95 .- (0.01*matrixrulemetrics[:,3]/max(matrixrulemetrics[:,3]...))
            #@show coefReg
            rf = DT.build_forest(y,M,2,50,0.7,-1; rng=rng)
            importances = begin
                #importance = impurity_importance(rf, coefReg)
                importance = DT.impurity_importance(rf)
                importance/max(importance...)
            end
            # @show importances
            best_idxs = begin
                selected_features = findall(importances .> 0.01)
                # @show selected_features
                ruleSetPrunedRRF = hcat(matrixrulemetrics[selected_features,:],importances[selected_features],selected_features)
                finalmatrix = sortslices(ruleSetPrunedRRF, dims=1, by=x->(x[4],x[2],x[3]), rev=true)
                best_idxs = Int.(finalmatrix[:,5])
            end
            # @show best_idxs
            ruleset[best_idxs]
        else
            error("Unexpected rule selection method specified: $(rule_selection_method)")
        end
    end
    silent || println("# rules selected: $(length(ruleset)).")
    
    ########################################################################################
    # Construct a rule-based model from the set of best rules
    ########################################################################################
    silent || println("Applying STEL...")
    
    dl = STEL(ruleset, X, y; min_coverage, rule_length_metric, rng, silent)

    if return_info
        return dl, info
    else
        return dl
    end
end





function STEL(ruleset, X, y; min_coverage, rule_length_metric = :natoms, rng::AbstractRNG = MersenneTwister(1), silent = false)
    D = deepcopy(X) # Copy of the original dataset
    L = deepcopy(y)
    R = Rule[]      # Ordered rule list
    # S is the vector of rules left
    S = [deepcopy(ruleset)..., Rule(bestguess(y; suppress_parity_warning = true))]

    # Rules with a frequency less than min_coverage
    S = begin
        rules_coverage = Vector{Float64}(undef, length(S))
        Threads.@threads for (i,s) in collect(enumerate(S))
            rules_coverage[i] = rulemetrics(s,X,y)[:coverage]
        end

        idxs_undeleted = findall(rules_coverage .>= min_coverage)
        S[idxs_undeleted]
    end

    silent || println("# rules in S: $(length(S))")
    silent || println("# rules in R: $(length(R))")

    # Metrics update based on remaining instances
    rules_coverage = Vector{Float64}(undef, length(S))
    rules_error = Vector{Float64}(undef, length(S))
    rules_length = Vector{Int}(undef, length(S))

    while true
        silent || println()
        silent || println()
        silent || println()
        
        resize!(rules_coverage, length(S))
        resize!(rules_error, length(S))
        resize!(rules_length, length(S))

        silent || println("Rules left: $(length(rules_coverage)).")

        Threads.@threads for (i,s) in collect(enumerate(S))
            metrics = rulemetrics(s,D,L)
            rules_coverage[i] = metrics[:coverage]
            rules_error[i] = metrics[:error]
            rules_length[i] = metrics[rule_length_metric]
            # @show metrics
        end
        silent || println("Rules error:")
        silent || println(rules_error)

        #metrics = [rulemetrics(s,D,y) for s in S]
        #silent || println(metrics[1])
        #rules_coverage = [metrics[i][:coverage] for i in eachindex(metrics)]
        #rules_error = [metrics[i][:error] for i in eachindex(metrics)]
        #rules_length = [metrics[i][rule_length_metric] for i in eachindex(metrics)]

        # Best rule index
        idx_best = begin
            idx_best = nothing
            # First: find the rule with minimum error
            silent || println("By error")
            m = minimum(filter(!isnan, rules_error))
            silent || println("Minimum rule error: $(m)")
            best_idxs = findall(rules_error .== m)
            (length(best_idxs) == 1) && (idx_best = best_idxs[1])
            silent || println("Idxs min error: $(best_idxs)")

            # If not one, find the rule with maximum coverage
            if isnothing(idx_best)
                silent || println("By coverage")
                # @show rules_coverage[best_idxs]
                m = maximum(filter(!isnan, rules_coverage[best_idxs]))
                silent || println("Max coverage: $(m)")
                idx_coverage = findall(rules_coverage[best_idxs] .== m)
                best_idxs = best_idxs[idx_coverage]
                silent || println("Idxs max coverage: $(best_idxs)")
                (length(best_idxs) == 1) && (idx_best = best_idxs[1])
            end

            # If not one, find the rule with minimum length
            if isnothing(idx_best)
                silent || println("By length")
                # @show rules_length[best_idxs]
                m = minimum(filter(!isnan, rules_length[best_idxs]))
                silent || println("Min length: $(m)")
                idx_length = findall(rules_length[best_idxs] .== m)
                best_idxs = best_idxs[idx_length]
                silent || println("Idxs min length: $(best_idxs)")
                (length(best_idxs) == 1) && (idx_best = best_idxs[1])
            end

            # Final case: more than one rule with minimum length
            # Randomly choose a rule
            if isnothing(idx_best)
                idx_best = rand(rng, best_idxs)
            end

            idx_best
        end
        silent || println("Idx best: $(idx_best)")

        # Add at the end the best rule
        push!(R, S[idx_best])
        silent || println("# rules in R: $(length(R))")

        # Indices of the remaining instances
        idx_remaining = begin
            sat_unsat = evaluaterule(S[idx_best], D, L)[:checkmask,]
            # Remain in D the rule that not satisfying the best rule's pruning_s condition
            findall(sat_unsat .== false)
        end

        # Exit condition
        if idx_best == length(S)
            # @show S
            # @show R
            # @show DecisionList(R[1:end-1], consequent(R[end]))
            return DecisionList(R[1:end-1], consequent(R[end]))
        elseif length(idx_remaining) == 0
            return DecisionList(R, bestguess(y; suppress_parity_warning = true))
        end
        
        D = slicedataset(D, idx_remaining; return_view=true)
        L = L[idx_remaining]
        # Delete the best rule from S
        deleteat!(S, idx_best)
        # Update of the default rule
        S[end] = Rule(bestguess(L; suppress_parity_warning = true))
    end
    error("Unexpected error.")
end
