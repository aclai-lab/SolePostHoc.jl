using SoleLogics
using SoleModels
using SoleModels: DecisionList, rulebase, defaultconsequent, info

function GENERATE_CROSSOVER_FUNCTION()
    function crossover(i1::DLIndividual, i2::DLIndividual; rng::AbstractRNG = default_rng())
        cls = classes(i1)
        alp = alphabet(i1)
        adls = alldls(i1)
        dl1, dl2 = dl(i1), dl(i2)
        dl1, dl2 = swapintervalsrule(dl1, dl2; rng = rng)

        return (DLIndividual(dl1, alp, cls, adls), DLIndividual(dl2, alp, cls, adls))
    end

    return crossover
end

############################################################################################
############################################################################################
############################################################################################

function swapintervalsrule(
    dl1::DecisionList,
    dl2::DecisionList;
    rng::AbstractRNG = default_rng(),
)
    dl1_rules = rulebase(dl1)
    dl2_rules = rulebase(dl2)

    if length(dl1_rules) == 0 && length(dl2_rules) == 0
        return (dl1, dl2)
    elseif length(dl1_rules) == 0
        from_dl2, to_dl2 = from_to(length(dl2_rules); rng = rng, one = true)

        dl1_rules = dl2_rules[from_dl2:to_dl2]
        dl2_rules = deleteat!(dl2_rules, from_dl2:to_dl2)

    elseif length(dl2_rules) == 0
        from_dl1, to_dl1 = from_to(length(dl1_rules); rng = rng, one = true)

        dl2_rules = dl1_rules[from_dl1:to_dl1]
        dl1_rules = deleteat!(dl1_rules, from_dl1:to_dl1)
    else
        from_dl1, to_dl1 = from_to(length(dl1_rules); rng = rng, one = true)
        from_dl2, to_dl2 = from_to(length(dl2_rules); rng = rng, one = true)

        tmp_rules = dl1_rules[from_dl1:to_dl1]
        dl1_rules = [
            dl1_rules[1:(from_dl1-1)]...,
            dl2_rules[from_dl2:to_dl2]...,
            dl1_rules[(to_dl1+1):end]...,
        ]
        dl2_rules =
            [dl2_rules[1:(from_dl2-1)]..., tmp_rules..., dl2_rules[(to_dl2+1):end]...]
    end

    #=
    TODO: fix ancestors
    ancestors = unique([
        i1,
        i2,
        (haskey(info(dl1), :evo_ancestors) ? info(dl1).ancestors : [])...,
        (haskey(info(dl2), :evo_ancestors) ? info(dl2).ancestors : [])...,
    ])

    i1 = Base.structdiff(info(dl1), (; tree_id = nothing))
    i1 = merge(i1, (; evo_ancestors = ancestors))
    i2 = Base.structdiff(info(dl2), (; tree_id = nothing))
    i2 = merge(i2, (; evo_ancestors = ancestors))
    =#

    return (
        DecisionList(dl1_rules, defaultconsequent(dl1), info(dl1)), #, i1
        DecisionList(dl2_rules, defaultconsequent(dl2), info(dl2)), #, i2
    )
end
