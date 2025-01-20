# TODO


## LUMEN
- Remove Italian sentences
- 

## BA-Trees
BA-Trees must be integrated well as a `RuleExtraction` object. Some first steps:
- Separate forest learning, forest simplification, and algorithm testing logics.
- Write a function that translates decision forests to BA-tree txt format.
- Write a function that translates ANY (non-bag) symbolic model to BA-tree txt format. Remember that any non-bag model can be seen as a particular case of a forest (ensemble of trees), via a polynomial reduction
- Build object `BATreesRuleExtraction`, which outputs a SoleModel.DecisionTree. Be careful: clarify the assumpions on the symbolic model type, or the underlying logical language
- Can we extend it to the modal level?
