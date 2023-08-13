module RuleExtraction

using SoleBase

using SoleData
using SoleData: slicedataset

using SoleLogics
using SoleLogics: ⊤, AbstractInterpretationSet, nconjuncts, LeftmostConjunctiveForm

using SoleModels
using SoleModels: AbstractModel, AbstractLogiset
using SoleModels: Rule, antecedent, consequent, rulemetrics, trees
using SoleModels: LeafModel, Branch, DecisionForest, DecisionList
using SoleModels: listrules
using SoleModels: CLabel, RLabel, Label
using SoleModels: bestguess, evaluaterule

using SoleFeatures: findcorrelation

using Statistics: cor

export intrees
include("intrees.jl")

export bellatrex
include("bellatrex.jl")

end
