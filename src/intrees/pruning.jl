using SoleLogics: value, grandchildren

using SoleData: AbstractCondition, RangeScalarCondition, scalartiling,
    _rangescalarcond_to_scalarconds_in_conjunction

function Rprune_rule(
    set::ClassificationRule{T},
    X::Matrix{S},
    y::Vector{<:Label},
    featurenames::Vector{F};
    max_decay::Float64=0.05,
    type_decay::Int=2
) where {T,S<:Float,F}
    pred, err = measure_rule(set, X, y, featurenames)
    conds = LeftmostConjunctiveForm(set)
    nconds = length(conds)
    nconds ≤ 1 && return set, pred, err

    for i in nconds:-1:1
        remaining_conds = copy(conds)
        deleteat!(remaining_conds, i)
        _, err_new = measure_rule(remaining_conds, X, y, featurenames; pred)

        decay = type_decay == 1 ?
            (err_new - err) / max(err, 1e-6) :
            err_new - err

        decay ≤ max_decay && begin
            deleteat!(conds, i)
            set = scalartiling(value.(grandchildren(remaining_conds)))
            err = err_new
            length(conds) ≤ 1 && break
        end
    end
    # return scalartiling(conds), pred, err
    return set, pred, err
end

# ---------------------------------------------------------------------------- #
#                           Intrees pruning utility                            #
# ---------------------------------------------------------------------------- #
function prune_rule(::Type{<:Atom}, r::Rule{O}, args...; kwargs...) where{O}
    r = Rule(LeftmostConjunctiveForm([antecedent(r)]), consequent(r), info(r))
    prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
end

function prune_rule(
    ::Type{<:LeftmostConjunctiveForm},
    r::Rule{O},
    X::Matrix{<:Float},
    y::Vector{<:Label};
    decay_threshold::Float64=0.05,
    decay_type::Int=2,
    decay_s::Float64=1.0e-6
) where {O}
    nruleconjuncts = SoleModels.nconjuncts(r)
    e_zero = SoleModels.rulemetrics(r, X, y)[:error]
    valid_idxs = 1:nruleconjuncts
    antd, cons = SoleModels.antecedent(r), SoleModels.consequent(r)

    for idx in reverse(valid_idxs)
        (length(valid_idxs) ≤ 1) && break

        # indices to be considered to evaluate the rule
        other_idxs = setdiff(valid_idxs, idx)
        rule = Rule(
            LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[other_idxs]),
            cons
        )

        # return error of the rule without idx-th pair
        e_minus_i = SoleModels.rulemetrics(rule, X, y)[:error]

        decay_i = decay_type == 1 ?
            (e_minus_i - e_zero) / max(e_zero, decay_s) :
            e_minus_i - e_zero

        if decay_i ≤ decay_threshold 
            # remove the idx-th pair in the vector of decisions
            valid_idxs = setdiff(valid_idxs, idx)
            e_zero = e_minus_i
        end
    end

    return Rule(
        LeftmostConjunctiveForm(SoleLogics.grandchildren(antd)[valid_idxs]),
        cons
    )
end

function prune_rule(
    ::Type{<:MultiFormula},
    r::Rule{O},
    args...;
    kwargs...
) where {O}
    # @assert antecedent(r) isa MultiFormula
    #     "Cannot use this function on $(antecedent(r))"
    children = [
        MultiFormula(i_modality, modant)
            for (i_modality, modant) in modforms(antecedent(r))
    ]

    return  length(children) ≤ 1 ? r : begin
        r = Rule(LeftmostConjunctiveForm(children), consequent(r), info(r))
        prune_rule(typeof(antecedent(r)), r, args...; kwargs...,)
    end
end

function prune_ruleset(
    ruleset::Vector{<:Rule},
    args...;
    kwargs...
)
    pruned = similar(ruleset)
    
    @inbounds Threads.@threads for i in eachindex(ruleset)
        pruned[i] = if ruleset[i].antecedent isa SoleLogics.BooleanTruth
            ruleset[i]  # keep BooleanTruth rules as-is (e.g., from XGBoost)
        else
            prune_rule(
                typeof(antecedent(ruleset[i])),
                ruleset[i],
                args...;
                kwargs...
            )
        end
    end
    
    return pruned
end

# function _prune_ruleset(
#     ruleset::Vector{<:Rule},
#     X::AbstractInterpretationSet,
#     y::AbstractVector{<:SoleModels.Label},
#     extractor::InTreesRuleExtractor
# )
#     pruned = similar(ruleset)
    
#     @inbounds Threads.@threads for i in eachindex(ruleset)
#         pruned[i] = if ruleset[i].antecedent isa SoleLogics.BooleanTruth
#             ruleset[i]  # keep BooleanTruth rules as-is (e.g., from XGBoost)
#         else
#             prune_rule(
#                 typeof(antecedent(ruleset[i])), ruleset[i], X, y;
#                 pruning_s=get_pruning_s(extractor),
#                 pruning_decay_thr=get_pruning_decay_threshold(extractor)
#             )
#         end
#     end
    
#     return pruned
# end
