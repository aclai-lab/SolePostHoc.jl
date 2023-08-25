module RuleExtraction

using SoleBase
using Test
using ThreadsX
using Random

using SoleData
using SoleData: slicedataset

using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet, LeftmostConjunctiveForm
using SoleLogics: conjuncts, nconjuncts

using SoleModels
using SoleModels: AbstractModel, AbstractLogiset
using SoleModels: Rule, antecedent, consequent, rulemetrics, trees, info
using SoleModels: LeafModel, Branch, DecisionForest, DecisionList
using SoleModels: listrules, antecedenttops
using SoleModels: CLabel, RLabel, Label
using SoleModels: bestguess, evaluaterule

using SoleFeatures: findcorrelation

using Statistics: cor

export intrees
include("intrees.jl")

export bellatrex
include("bellatrex.jl")

end
