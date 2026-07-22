using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

using RCall
@rlibrary inTrees

Xc, yc = @load_iris
Xc = DataFrame(Xc)

@rput Xc yc

R"""
library(randomForest)

rf <- randomForest(as.data.frame(Xc), as.factor(yc), ntree=2)

# dtree <- rpart(target ~ ., data = df, parms = list(split = "entropy"))

# dtree <- rpart(target ~ ., data = df, method = "class", parms = list(split = "information"), control = rpart.control(
#         cp = 0,
#         minsplit = 2,
#         minbucket = 1,
#         maxdepth = 30,
#         xval = 0,
#     )
# )
for (i in seq_len(rf$ntree)) {
    cat("\n--- Tree", i, "---\n")
    print(getTree(rf, k = i, labelVar = TRUE))
}
"""

# --- Tree 1 ---
#   left daughter right daughter    split var split point status prediction
# 1             2              3 petal_length        2.60      1       <NA>
# 2             0              0         <NA>        0.00     -1     setosa
# 3             4              5  petal_width        1.75      1       <NA>
# 4             0              0         <NA>        0.00     -1 versicolor
# 5             0              0         <NA>        0.00     -1  virginica

# --- Tree 2 ---
#    left daughter right daughter    split var split point status prediction
# 1              2              3  petal_width        0.80      1       <NA>
# 2              0              0         <NA>        0.00     -1     setosa
# 3              4              5  petal_width        1.75      1       <NA>
# 4              6              7  sepal_width        2.25      1       <NA>
# 5              8              9 sepal_length        6.00      1       <NA>
# 6             10             11  petal_width        1.25      1       <NA>
# 7             12             13  petal_width        1.65      1       <NA>
# 8             14             15 petal_length        4.85      1       <NA>
# 9              0              0         <NA>        0.00     -1  virginica
# 10             0              0         <NA>        0.00     -1 versicolor
# 11             0              0         <NA>        0.00     -1  virginica
# 12             0              0         <NA>        0.00     -1 versicolor
# 13            16             17 petal_length        4.75      1       <NA>
# 14             0              0         <NA>        0.00     -1 versicolor
# 15             0              0         <NA>        0.00     -1  virginica
# 16             0              0         <NA>        0.00     -1  virginica
# 17             0              0         <NA>        0.00     -1 versicolor

featurenames = [
    :sepal_length,
    :sepal_width,
    :petal_length,
    :petal_width,
]

setosa     = yc[findfirst(==("setosa"), yc)]
versicolor = yc[findfirst(==("versicolor"), yc)]
virginica  = yc[findfirst(==("virginica"), yc)]

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
# 1: petal_length < 2.60
# ├─ 2: setosa
# └─ 3: petal_width < 1.75
#    ├─ 4: versicolor
#    └─ 5: virginica
# --------------------------------------------------------------------------- #
tree1_root = node(
    split(3, 2.60), # petal_length < 2.60
    leaf(setosa),
    node(
        split(4, 1.75), # petal_width < 1.75
        leaf(versicolor),
        leaf(virginica),
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
# 1: petal_width < 0.80
# ├─ 2: setosa
# └─ 3: petal_width < 1.75
#    ├─ 4: sepal_width < 2.25
#    │  ├─ 6: petal_width < 1.25
#    │  │  ├─ 10: versicolor
#    │  │  └─ 11: virginica
#    │  └─ 7: petal_width < 1.65
#    │     ├─ 12: versicolor
#    │     └─ 13: petal_length < 4.75
#    │        ├─ 16: virginica
#    │        └─ 17: versicolor
#    └─ 5: sepal_length < 6.00
#       ├─ 8: petal_length < 4.85
#       │  ├─ 14: versicolor
#       │  └─ 15: virginica
#       └─ 9: virginica
# --------------------------------------------------------------------------- #
tree2_root = node(
    split(4, 0.80), # petal_width < 0.80
    leaf(setosa),
    node(
        split(4, 1.75), # petal_width < 1.75
        node(
            split(2, 2.25), # sepal_width < 2.25
            node(
                split(4, 1.25), # petal_width < 1.25
                leaf(versicolor),
                leaf(virginica),
            ),
            node(
                split(4, 1.65), # petal_width < 1.65
                leaf(versicolor),
                node(
                    split(3, 4.75), # petal_length < 4.75
                    leaf(virginica),
                    leaf(versicolor),
                ),
            ),
        ),
        node(
            split(1, 6.00), # sepal_length < 6.00
            node(
                split(3, 4.85), # petal_length < 4.85
                leaf(versicolor),
                leaf(virginica),
            ),
            leaf(virginica),
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
#                          intrees rules extraction                            #
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
#  [1,] "petal_length<=2.6"                                                                             
#  [2,] "petal_length>2.6 & petal_width<=1.75"                                                          
#  [3,] "petal_length>2.6 & petal_width>1.75"                                                           
#  [4,] "petal_width<=0.8"                                                                              
#  [5,] "sepal_width<=2.25 & petal_width>0.8 & petal_width<=1.75 & petal_width<=1.25"                   
#  [6,] "sepal_width<=2.25 & petal_width>0.8 & petal_width<=1.75 & petal_width>1.25"                    
#  [7,] "sepal_width>2.25 & petal_width>0.8 & petal_width<=1.75 & petal_width<=1.65"                    
#  [8,] "sepal_width>2.25 & petal_length<=4.75 & petal_width>0.8 & petal_width<=1.75 & petal_width>1.65"
#  [9,] "sepal_width>2.25 & petal_length>4.75 & petal_width>0.8 & petal_width<=1.75 & petal_width>1.65" 
# [10,] "sepal_length<=6 & petal_length<=4.85 & petal_width>0.8 & petal_width>1.75"                     
# [11,] "sepal_length<=6 & petal_length>4.85 & petal_width>0.8 & petal_width>1.75"                      
# [12,] "sepal_length>6 & petal_width>0.8 & petal_width>1.75"  

extractor = InTreesRuleExtractor()
extracted_rules =
    RuleExtraction.extractrules(extractor, solem_rf, Xc, yc)
