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

using Statistics: cor

using SoleModels: RuleExtractor
import SoleModels: isexact, extractrules

export InTreesRuleExtractor
include("intrees.jl")

export BellatrexRuleExtractor
include("bellatrex.jl")

export lumen, Lumen
export LumenRuleExtractor
include("lumen/main.jl")

end
