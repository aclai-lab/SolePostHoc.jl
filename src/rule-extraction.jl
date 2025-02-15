module RuleExtraction

using Reexport
using Test
using ThreadsX
using Random
using StatsBase

using SoleData
using SoleData: slicedataset
using SoleData: feature, i_variable

using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet, LeftmostConjunctiveForm
using SoleLogics: conjuncts, nconjuncts

using SoleModels
using SoleModels: AbstractModel, AbstractLogiset
using SoleModels: Rule, antecedent, consequent, rulemetrics, trees, info
using SoleModels: LeafModel, Branch, DecisionForest, DecisionList
using SoleModels: listrules, checkantecedent
using SoleModels: CLabel, RLabel, Label
using SoleModels: bestguess, evaluaterule

# using SoleFeatures: findcorrelation

# using Statistics: cor

using SoleModels: RuleExtractor
import SoleModels: isexact, extractrules

function _get_rule_extractor_docstring(ruleextractorname::String, method)
    return """Extract rules from a symbolic model using [`$(string(method))`](ref).""" *
    "\n\n" *
    """See also [`extractrules`](@ref), [`RuleExtractor`](@ref)."""
end

export InTreesRuleExtractor
export intrees
include("intrees/intrees.jl")


"""$(_get_rule_extractor_docstring("InTreesRuleExtractor", intrees))"""
@kwdef struct InTreesRuleExtractor <: RuleExtractor
    prune_rules::Bool = true
    pruning_s::Union{Float64,Nothing} = nothing
    pruning_decay_threshold::Union{Float64,Nothing} = nothing
    rule_selection_method::Symbol = :CBC
    rule_complexity_metric::Symbol = :natoms
    # accuracy_rule_selection = nothing
    min_coverage::Union{Float64,Nothing} = nothing
end

function extractrules(::InTreesRuleExtractor, m, args...; kwargs...)
  dl = intrees(m, args...; kwargs...)
  return listrules(dl)
end


export BellatrexRuleExtractor
export bellatrex
include("bellatrex.jl")


"""$(_get_rule_extractor_docstring("BellatrexRuleExtractor", bellatrex))"""
@kwdef struct BellatrexRuleExtractor <: RuleExtractor
    ntrees::Float64
    ndims::Union{Nothing,Int}
    nclusters::Int
end

function extractrules(::BellatrexRuleExtractor, m, args...; kwargs...)
  dl = bellatrex(m, args...; kwargs...)
  return listrules(dl)
end


export lumen, Lumen
export LumenRuleExtractor

export batrees, BATrees

using SoleModels: RuleExtractor
import SoleModels: isexact, extractrules

include("lumen/main.jl")
@reexport using .Lumen

"""$(_get_rule_extractor_docstring("LumenRuleExtractor", lumen))"""
struct LumenRuleExtractor <: RuleExtractor end

function extractrules(::LumenRuleExtractor, m, args...; kwargs...)
  dl = lumen(m, args...; kwargs...)
  return listrules(dl)
end

include("BA-Trees/src/main.jl")
@reexport using .BATrees


"""$(_get_rule_extractor_docstring("BATreesRuleExtractor", batrees))"""
struct BATreesRuleExtractor <: RuleExtractor end

function extractrules(::BATreesRuleExtractor, m, args...; kwargs...)
  dl = batrees(m, args...; kwargs...)
  return listrules(dl)
end

end
