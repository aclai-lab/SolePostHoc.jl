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
import SoleModels: isexact, modalextractrules

function _get_rule_extractor_docstring(ruleextractorname::String, method)
    return """Extract rules from a symbolic model using [`$(string(method))`](ref).""" *
    "\n\n" *
    """See also [`modalextractrules`](@ref), [`RuleExtractor`](@ref)."""
end

export InTreesRuleExtractor
export intrees
include("intrees/intrees.jl")
include("intrees/apiIntrees.jl")


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

function modalextractrules(::InTreesRuleExtractor, m, args...; kwargs...)
  dl = intrees(m, args...; kwargs...)
  ll = listrules(dl) # decision list to list of rules
  rules_obj = convert_classification_rules(ll)
  ds = DecisionSet(rules_obj)  
  return ds
end

#= 
  IS NOT RULE EXTRACTION BUT MAYBE IN FUTURE 
    export BellatrexRuleExtractor
    export bellatrex
    include("bellatrex.jl")


    """$(_get_rule_extractor_docstring("BellatrexRuleExtractor", bellatrex))"""
    @kwdef struct BellatrexRuleExtractor <: RuleExtractor
        ntrees::Float64
        ndims::Union{Nothing,Int}
        nclusters::Int
    end

    function modalextractrules(::BellatrexRuleExtractor, m, args...; kwargs...)
      dl = bellatrex(m, args...; kwargs...)
      return listrules(dl)
    end
=#

export lumen, Lumen
export LumenRuleExtractor

export batrees, BATrees

using SoleModels: RuleExtractor
import SoleModels: isexact, modalextractrules

include("lumen/main.jl")
@reexport using .Lumen

"""$(_get_rule_extractor_docstring("LumenRuleExtractor", lumen))"""
struct LumenRuleExtractor <: RuleExtractor end

function modalextractrules(::LumenRuleExtractor, m, args...; kwargs...)
  ds = lumen(m, args...; kwargs...)
  return ds
end

#= TODO: Include BATrees =#
include("BA-Trees/src/main.jl")
@reexport using .BATrees


"""$(_get_rule_extractor_docstring("BATreesRuleExtractor", batrees))"""
struct BATreesRuleExtractor <: RuleExtractor end

function modalextractrules(::BATreesRuleExtractor, m, args...; kwargs...)
  ds = batrees(m, dsOutput = true, args...; kwargs...)
  return ds
end

#
include("Refne/src/main.jl")
include("Refne/src/apiREFNESole.jl")
@reexport using .REFNE


"""$(_get_rule_extractor_docstring("REFNERuleExtractor", REFNE))"""
struct REFNERuleExtractor <: RuleExtractor end

function modalextractrules(::REFNERuleExtractor, m, args...; kwargs...)
  dl = refne(m, args...; kwargs...)
  ds = convertApi(dl)
  return ds
end

# TODO see and test this
#==============================================================================#
include("Trepan/src/main.jl")
include("Trepan/src/apiTREPANSole.jl")
@reexport using .TREPAN


"""$(_get_rule_extractor_docstring("TREPANRuleExtractor", TREPAN))"""
struct TREPANRuleExtractor <: RuleExtractor end

function modalextractrules(::TREPANRuleExtractor, m, args...; kwargs...)
  dl = trepan(m, args...; kwargs...)
  ds = convertApi(dl)
  return ds
end

#=
include("RuleCosiplus/src/main.jl")
#include("RuleCosiplus/src/apiRuleCosi.jl")
@reexport using .RULECOSIPLUS


"""$(_get_rule_extractor_docstring("RULECOSIPLUSRuleExtractor", RULECOSIPLUS))"""
struct RULECOSIPLUSRuleExtractor <: RuleExtractor end

function modalextractrules(::RULECOSIPLUSRuleExtractor, m, args...; kwargs...)
  ds = rulecosiplus(m, args...; kwargs...)  # for now we return directly decision set object
  return ds
end
=#

end
