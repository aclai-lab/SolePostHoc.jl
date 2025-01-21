module RuleExtraction

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

function _get_rule_extractor_docstring(::Type{<:RuleExtractor}, method)
    return """Extract rules from a symbolic model using [`$(string(method))`](ref).""" *
    "\n\n" *
    """See also [`extractrules`](@ref), [`RuleExtractor`](@ref)."""
end

using SoleModels: RuleExtractor
import SoleModels: isexact, extractrules

export InTreesRuleExtractor
include("intrees/intrees.jl")

export BellatrexRuleExtractor
include("bellatrex.jl")

export lumen, Lumen
export LumenRuleExtractor

#export batrees, BATrees #TODO here

using SoleModels: RuleExtractor
import SoleModels: isexact, extractrules

"""
TODO explain
Pagliarini, Giovanni, et al. "Minimal Rules from Decision Forests: a Systematic Approach." OVERLAY@ AI* IA. 2024.

See also [`extractrules`](@ref), [`intrees`](@ref), [`RuleExtractor`](@ref).
"""
struct LumenRuleExtractor <: RuleExtractor end

extractrules(::LumenRuleExtractor, m, args...; kwargs...) = Lumen.lumen(m, args...; kwargs...)

include("lumen/main.jl")
#include("BA-Trees/src/main.jl") #TODO here

end
