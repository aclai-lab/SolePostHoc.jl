module RuleExtraction

using Reexport
# using Test     --- not needed
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

export convert_classification_rules, refne_classification_rules
include("shared_utils.jl")

#======================================================================================================================================
                                                        InTrees
======================================================================================================================================#

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
  ll = listrules(dl, use_shortforms=false) # decision list to list of rules
  rules_obj = convert_classification_rules(dl, ll)
  dsintrees = DecisionSet(rules_obj)
  return dsintrees
end

#======================================================================================================================================
                                                        Lumen
======================================================================================================================================#

export lumen, Lumen
export LumenRuleExtractor

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

#======================================================================================================================================
                                                        Batrees
======================================================================================================================================#

export batrees, BATrees
export BATreesRuleExtractor

include("BA-Trees/src/main.jl")
@reexport using .BATrees


"""$(_get_rule_extractor_docstring("BATreesRuleExtractor", batrees))"""
struct BATreesRuleExtractor <: RuleExtractor end

function modalextractrules(::BATreesRuleExtractor, m, args...; kwargs...)
  dsbatrees = batrees(m, dsOutput=true, args...; kwargs...)
  return dsbatrees
end

#======================================================================================================================================
                                                        Refne
======================================================================================================================================#
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

#======================================================================================================================================
                                                        TrePan
======================================================================================================================================#
include("Trepan/src/main.jl")
# include("Trepan/src/apiTREPANSole.jl")
@reexport using .TREPAN


"""$(_get_rule_extractor_docstring("TREPANRuleExtractor", TREPAN))"""
struct TREPANRuleExtractor <: RuleExtractor end

function modalextractrules(::TREPANRuleExtractor, m, args...; kwargs...)
  dl = trepan(m, args...; kwargs...)
  ds = convertApi(dl)
  return ds
end

#======================================================================================================================================
                                                        RULECOSIPLUS
======================================================================================================================================#
export RULECOSIPLUSRuleExtractor

include("RuleCosiplus/src/main.jl")
include("RuleCosiplus/src/apiRuleCosi.jl")
@reexport using .RULECOSIPLUS


"""$(_get_rule_extractor_docstring("RULECOSIPLUSRuleExtractor", RULECOSIPLUS))"""
struct RULECOSIPLUSRuleExtractor <: RuleExtractor end

function modalextractrules(::RULECOSIPLUSRuleExtractor, m, args...; kwargs...)
  dl = rulecosiplus(m, args...; kwargs...) # decision list   
  ll = listrules(dl, use_shortforms=false) # decision list to list of rules
  rules_obj = convert_classification_rules(dl, ll)
  dsrulecosiplus = DecisionSet(rules_obj)
  return dsrulecosiplus
end


end
