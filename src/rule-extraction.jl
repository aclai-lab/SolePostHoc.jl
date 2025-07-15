module RuleExtraction

using Reexport
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

# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #
const Optional{T}     = Union{T, Nothing}
const OptFloat64      = Optional{Float64}

# ---------------------------------------------------------------------------- #
#                                   intrees                                    #
# ---------------------------------------------------------------------------- #
export InTreesRuleExtractor
export intrees
include("intrees/intrees.jl")
include("intrees/apiIntrees.jl")


"""$(_get_rule_extractor_docstring("InTreesRuleExtractor", intrees))"""
struct InTreesRuleExtractor <: RuleExtractor
    prune_rules             :: Bool
    pruning_s               :: OptFloat64
    pruning_decay_threshold :: OptFloat64
    rule_selection_method   :: Symbol
    rule_complexity_metric  :: Symbol
    max_rules               :: Int
    min_coverage            :: OptFloat64
    rng                     :: AbstractRNG

    function InTreesRuleExtractor(;
        prune_rules             :: Bool        = true,
        pruning_s               :: OptFloat64  = nothing,
        pruning_decay_threshold :: OptFloat64  = nothing,
        rule_selection_method   :: Symbol      = :CBC,
        rule_complexity_metric  :: Symbol      = :natoms,
        max_rules               :: Int         = -1,
        min_coverage            :: OptFloat64  = nothing,
        rng                     :: AbstractRNG = TaskLocalRNG()
    )
        return new(
            prune_rules,
            pruning_s,
            pruning_decay_threshold,
            rule_selection_method,
            rule_complexity_metric,
            max_rules,
            min_coverage,
            rng
        )
    end
end

function Base.iterate(e::InTreesRuleExtractor, state=1)
    fields = fieldnames(typeof(e))
    if state > length(fields)
        return nothing
    end
    field = fields[state]
    return (field => getfield(e, field)), state + 1
end

function Base.show(io::IO, info::InTreesRuleExtractor)
    println(io, "InTreesRuleExtractor:")
    for field in fieldnames(InTreesRuleExtractor)
        value = getfield(info, field)
        println(io, "  ", rpad(String(field) * ":", 25), value)
    end
end

function modalextractrules(e::InTreesRuleExtractor, m, args...; kwargs...)
  dl = intrees(m, args...; e..., kwargs...)
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
struct LumenRuleExtractor <: RuleExtractor 
    minimization_scheme    :: Symbol
    vertical               :: Real
    horizontal             :: Real
    ott_mode               :: Bool
    controllo              :: Bool
    minimization_kwargs    :: NamedTuple
    filteralphabetcallback :: Function
    vetImportance          :: Vector

    function LumenRuleExtractor(;
        minimization_scheme    :: Symbol     = :mitespresso,
        vertical               :: Real       = 1.0,
        horizontal             :: Real       = 1.0,
        ott_mode               :: Bool       = false,
        controllo              :: Bool       = false,
        minimization_kwargs    :: NamedTuple = (;),
        filteralphabetcallback :: Function   = identity,
        vetImportance          :: Vector     = []
    )
        return new(
            minimization_scheme,
            vertical,
            horizontal,
            ott_mode,
            controllo,
            minimization_kwargs,
            filteralphabetcallback,
            vetImportance
        )
    end
end

function Base.iterate(e::LumenRuleExtractor, state=1)
    fields = fieldnames(typeof(e))
    if state > length(fields)
        return nothing
    end
    field = fields[state]
    return (field => getfield(e, field)), state + 1
end

function Base.show(io::IO, info::LumenRuleExtractor)
    println(io, "LumenRuleExtractor:")
    for field in fieldnames(LumenRuleExtractor)
        value = getfield(info, field)
        println(io, "  ", rpad(String(field) * ":", 25), value)
    end
end

function modalextractrules(e::LumenRuleExtractor, m, args...; kwargs...)
  ds = lumen(m, args...; e..., kwargs...)
  return ds
end

#======================================================================================================================================
                                                        Batrees
======================================================================================================================================#

export batrees, BATrees

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
