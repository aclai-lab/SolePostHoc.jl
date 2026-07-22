using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

using RCall
@rlibrary inTrees

Xc, yc = @load_iris
Xc = DataFrame(Xc)

@rput Xc yc

# ---------------------------------------------------------------------------- #
R"""
library(randomForest)
set.seed(1)

rf <- randomForest(as.data.frame(Xc), as.factor(yc), ntree=2)

for (i in seq_len(rf$ntree)) {
    cat("\n--- Tree", i, "---\n")
    print(getTree(rf, k = i, labelVar = TRUE))
}
"""

# --- Tree 1 ---
#    left daughter right daughter    split var split point status prediction
# 1              2              3 petal_length        2.45      1       <NA>
# 2              0              0         <NA>        0.00     -1     setosa
# 3              4              5 petal_length        4.95      1       <NA>
# 4              6              7  petal_width        1.65      1       <NA>
# 5              8              9 petal_length        5.05      1       <NA>
# 6              0              0         <NA>        0.00     -1 versicolor
# 7              0              0         <NA>        0.00     -1  virginica
# 8             10             11 sepal_length        6.50      1       <NA>
# 9              0              0         <NA>        0.00     -1  virginica
# 10             0              0         <NA>        0.00     -1  virginica
# 11             0              0         <NA>        0.00     -1 versicolor

# --- Tree 2 ---
#    left daughter right daughter    split var split point status prediction
# 1              2              3  petal_width        0.75      1       <NA>
# 2              0              0         <NA>        0.00     -1     setosa
# 3              4              5  petal_width        1.75      1       <NA>
# 4              6              7  petal_width        1.55      1       <NA>
# 5              8              9  petal_width        1.85      1       <NA>
# 6             10             11 petal_length        4.95      1       <NA>
# 7             12             13 petal_length        5.40      1       <NA>
# 8             14             15 sepal_length        5.95      1       <NA>
# 9              0              0         <NA>        0.00     -1  virginica
# 10             0              0         <NA>        0.00     -1 versicolor
# 11             0              0         <NA>        0.00     -1  virginica
# 12             0              0         <NA>        0.00     -1 versicolor
# 13             0              0         <NA>        0.00     -1  virginica
# 14            16             17 petal_length        4.95      1       <NA>
# 15             0              0         <NA>        0.00     -1  virginica
# 16             0              0         <NA>        0.00     -1 versicolor
# 17             0              0         <NA>        0.00     -1  virginica

featurenames = [
    :sepal_length,
    :sepal_width,
    :petal_length,
    :petal_width,
]

leaf(label) = ConstantModel(label, (;
    supporting_predictions = [label],
    supporting_labels = [label],
))

split(feature_id::Int, threshold::Real) = Atom(
    ScalarCondition(
        VariableValue(feature_id, featurenames[feature_id]),
        <,
        threshold,
    ),
)

node(condition, left, right) = Branch(condition, left, right, (;
    supporting_predictions = [
        left.info[:supporting_predictions]...,
        right.info[:supporting_predictions]...,],
    supporting_labels = [
        left.info[:supporting_labels]...,
        right.info[:supporting_labels]...,],
))

# --------------------------------------------------------------------------- #
# Tree 1 from randomForest::getTree(rf, k = 1, labelVar = TRUE)
#
# 1: petal_length <= 2.45
# ├─ 2: setosa
# └─ 3: petal_length <= 4.95
#    ├─ 4: petal_width <= 1.65
#    │  ├─ 6: versicolor
#    │  └─ 7: virginica
#    └─ 5: petal_length <= 5.05
#       ├─ 8: sepal_length <= 6.50
#       │  ├─ 10: virginica
#       │  └─ 11: versicolor
#       └─ 9: virginica
# --------------------------------------------------------------------------- #
tree1_root = node(
    split(3, 2.45), # petal_length <= 2.45
    leaf("setosa"),
    node(
        split(3, 4.95), # petal_length <= 4.95
        node(
            split(4, 1.65), # petal_width <= 1.65
            leaf("versicolor"),
            leaf("virginica"),
        ),
        node(
            split(3, 5.05), # petal_length <= 5.05
            node(
                split(1, 6.50), # sepal_length <= 6.50
                leaf("virginica"),
                leaf("versicolor"),
            ),
            leaf("virginica"),
        ),
    ),
)

tree1 = SoleModels.DecisionTree(tree1_root, (;
    featurenames = featurenames,
    supporting_predictions = tree1_root.info[:supporting_predictions],
    supporting_labels = tree1_root.info[:supporting_labels],
))

# --------------------------------------------------------------------------- #
# Tree 2 from randomForest::getTree(rf, k = 2, labelVar = TRUE)
#
# 1: petal_width <= 0.75
# ├─ 2: setosa
# └─ 3: petal_width <= 1.75
#    ├─ 4: petal_width <= 1.55
#    │  ├─ 6: petal_length <= 4.95
#    │  │  ├─ 10: versicolor
#    │  │  └─ 11: virginica
#    │  └─ 7: petal_length <= 5.40
#    │     ├─ 12: versicolor
#    │     └─ 13: virginica
#    └─ 5: petal_width <= 1.85
#       ├─ 8: sepal_length <= 5.95
#       │  ├─ 14: petal_length <= 4.95
#       │  │  ├─ 16: versicolor
#       │  │  └─ 17: virginica
#       │  └─ 15: virginica
#       └─ 9: virginica
# --------------------------------------------------------------------------- #
tree2_root = node(
    split(4, 0.75), # petal_width <= 0.75
    leaf("setosa"),
    node(
        split(4, 1.75), # petal_width <= 1.75
        node(
            split(4, 1.55), # petal_width <= 1.55
            node(
                split(3, 4.95), # petal_length <= 4.95
                leaf("versicolor"),
                leaf("virginica"),
            ),
            node(
                split(3, 5.40), # petal_length <= 5.40
                leaf("versicolor"),
                leaf("virginica"),
            ),
        ),
        node(
            split(4, 1.85), # petal_width <= 1.85
            node(
                split(1, 5.95), # sepal_length <= 5.95
                node(
                    split(3, 4.95), # petal_length <= 4.95
                    leaf("versicolor"),
                    leaf("virginica"),
                ),
                leaf("virginica"),
            ),
            leaf("virginica"),
        ),
    ),
)

tree2 = SoleModels.DecisionTree(tree2_root, (;
    featurenames = featurenames,
    supporting_predictions = tree2_root.info[:supporting_predictions],
    supporting_labels = tree2_root.info[:supporting_labels],
))

trees = [tree1, tree2]

ensemble_info = (;
    featurenames = featurenames,
    supporting_predictions = vcat([
        tree.info[:supporting_predictions] for tree in trees
    ]...),
    supporting_labels = vcat([
        tree.info[:supporting_labels] for tree in trees
    ]...),
)

solem_rf = SoleModels.DecisionEnsemble{eltype(yc)}(
    trees,
    ensemble_info;
    parity_func = x -> argmax(x),
)

# ---------------------------------------------------------------------------- #
R"""
library(inTrees)

treeList <- RF2List(rf)
exec <- extractRules(treeList,Xc)
exec[1:2,]

presentRules(exec, c(
    "sepal_length",
    "sepal_width",
    "petal_length",
    "petal_width"
))
"""

#       condition                                                                                       
#  [1,] "petal_length<=2.45"                                                                               
#  [2,] "petal_length>2.45 & petal_length<=4.95 & petal_width<=1.65"                                       
#  [3,] "petal_length>2.45 & petal_length<=4.95 & petal_width>1.65"                                        
#  [4,] "sepal_length<=6.5 & petal_length>2.45 & petal_length>4.95 & petal_length<=5.05"                   
#  [5,] "sepal_length>6.5 & petal_length>2.45 & petal_length>4.95 & petal_length<=5.05"                    
#  [6,] "petal_length>2.45 & petal_length>4.95 & petal_length>5.05"                                        
#  [7,] "petal_width<=0.75"                                                                                
#  [8,] "petal_length<=4.95 & petal_width>0.75 & petal_width<=1.75 & petal_width<=1.55"                    
#  [9,] "petal_length>4.95 & petal_width>0.75 & petal_width<=1.75 & petal_width<=1.55"                     
# [10,] "petal_length<=5.4 & petal_width>0.75 & petal_width<=1.75 & petal_width>1.55"                      
# [11,] "petal_length>5.4 & petal_width>0.75 & petal_width<=1.75 & petal_width>1.55"                       
# [12,] "sepal_length<=5.95 & petal_length<=4.95 & petal_width>0.75 & petal_width>1.75 & petal_width<=1.85"
# [13,] "sepal_length<=5.95 & petal_length>4.95 & petal_width>0.75 & petal_width>1.75 & petal_width<=1.85" 
# [14,] "sepal_length>5.95 & petal_width>0.75 & petal_width>1.75 & petal_width<=1.85"                      
# [15,] "petal_width>0.75 & petal_width>1.75 & petal_width>1.85"  

extractor = InTreesRuleExtractor()
extracted_rules =
    RuleExtraction.extractrules(extractor, solem_rf, Xc, yc)

# set = ClassificationRule{String}
#  [▣ ([petal_length] ∈ [-Inf,2.45))  ↣  setosa
# , ▣ (([petal_length] ∈ [2.45,4.95))) ∧ (([petal_width] ∈ [-Inf,1.65)))  ↣  versicolor
# , ▣ (([petal_length] ∈ [2.45,4.95))) ∧ (([petal_width] ∈ [1.65,Inf]))  ↣  virginica
# , ▣ (([petal_length] ∈ [4.95,5.05))) ∧ (([sepal_length] ∈ [-Inf,6.5)))  ↣  virginica
# , ▣ (([petal_length] ∈ [4.95,5.05))) ∧ (([sepal_length] ∈ [6.5,Inf]))  ↣  versicolor
# , ▣ ([petal_length] ∈ [5.05,Inf])  ↣  virginica
# , ▣ ([petal_width] ∈ [-Inf,0.75))  ↣  setosa
# , ▣ (([petal_width] ∈ [0.75,1.55))) ∧ (([petal_length] ∈ [-Inf,4.95)))  ↣  versicolor
# , ▣ (([petal_width] ∈ [0.75,1.55))) ∧ (([petal_length] ∈ [4.95,Inf]))  ↣  virginica
# , ▣ (([petal_width] ∈ [1.55,1.75))) ∧ (([petal_length] ∈ [-Inf,5.4)))  ↣  versicolor
# , ▣ (([petal_width] ∈ [1.55,1.75))) ∧ (([petal_length] ∈ [5.4,Inf]))  ↣  virginica
# , ▣ (([petal_width] ∈ [1.75,1.85))) ∧ (([sepal_length] ∈ [-Inf,5.95))) ∧ (([petal_length] ∈ [-Inf,4.95)))  ↣  versicolor
# , ▣ (([petal_width] ∈ [1.75,1.85))) ∧ (([sepal_length] ∈ [-Inf,5.95))) ∧ (([petal_length] ∈ [4.95,Inf]))  ↣  virginica
# , ▣ (([petal_width] ∈ [1.75,1.85))) ∧ (([sepal_length] ∈ [5.95,Inf]))  ↣  virginica
# , ▣ ([petal_width] ∈ [1.85,Inf])  ↣  virginica
# ]

R"""
ruleMetric <- getRuleMetric(exec,Xc,yc)
ruleMetric[1:2,]
# pruneRule <- pruneRule(ruleMetric,Xc,yc)

presentRules(ruleMetric, c(
    "sepal_length",
    "sepal_width",
    "petal_length",
    "petal_width"
))
"""

#       pred        
#  [1,] "setosa"    
#  [2,] "versicolor"
#  [3,] "virginica" 
#  [4,] "virginica" 
#  [5,] "versicolor"
#  [6,] "virginica" 
#  [7,] "setosa"    
#  [8,] "versicolor"
#  [9,] "virginica" 
# [10,] "versicolor"
# [11,] "virginica" 
# [12,] "versicolor"
# [13,] "virginica" 
# [14,] "virginica" 
# [15,] "virginica"

# ---------------------------------------------------------------------------- #
# function listrules is validated against R implementation,
# giving the same results in the same order
# ---------------------------------------------------------------------------- #
