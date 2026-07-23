using Test
using SoleModels, SolePostHoc, SoleLogics, SoleData
using MLJ
using DataFrames, Random

using SoleData: AbstractCondition, RangeScalarCondition,
    honors_minval, honors_maxval

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

featurenames = [:sepal_length, :sepal_width, :petal_length, :petal_width]

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

presentRules(exec, colnames(Xc))
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

presentRules(ruleMetric, colnames(Xc))
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
R"""
set.seed(1)

pruneRules <- pruneRule(ruleMetric,Xc,yc)
pruneRules[1:2,]
selectRules <- selectRuleRRF(pruneRules,Xc,yc)

# presentRules(ruleMetric, colnames(Xc))
"""

# ---------------------------------------------------------------------------- #
# using DataFrames
# using Statistics

# _empty_metric() = (
#     len = -1,
#     freq = -1.0,
#     err = -1.0,
#     condition = "",
#     pred = "",
# )

# # R equivalent of names(which.max(table(ys))) for character-like labels.
# function _modal_label(ys)
#     labels = string.(ys)
#     candidates = sort!(unique(labels))
#     counts = [count(label -> label == candidate, labels) for candidate in candidates]
#     return candidates[argmax(counts)]
# end

# # R's is.numeric(target), except Bool is treated as categorical.
# _is_numeric_target(target) =
#     all(value -> value isa Real && !(value isa Bool), target)

# function _matching_indices(rule_exec, X, target, rule2table::Function)
#     membership = rule2table(rule_exec, X, target)

#     membership isa AbstractVector ||
#         throw(ArgumentError("rule2table must return a vector"))

#     length(membership) == size(X, 1) ||
#         throw(DimensionMismatch(
#             "rule2table returned $(length(membership)) entries, " *
#             "but X has $(size(X, 1)) rows",
#         ))

#     # Equivalent to R's which(...). Missing values are not matches.
#     return findall(value -> !ismissing(value) && !iszero(value), membership)
# end

# ---------------------------------------------------------------------------- #
R"""
pred <- NULL
regMethod <- "mean"
X <- Xc
target <- yc
ruleExec <- ruleMetric[2,]["condition"]

print(ruleExec)
len <- length(unlist(strsplit(ruleExec, split=" & ")))
print(len)
origRule <- ruleExec
ruleExec <- paste("which(", ruleExec, ")")
print(ruleExec)
ixMatch <- eval(parse(text=ruleExec))
print(ixMatch)

# if(length(ixMatch)==0){
# v <- c("-1","-1", "-1", "", "")
# names(v) <- c("len","freq","err","condition","pred")
# return(v)
# }
# ys <- target[ixMatch]
# freq <- round(length(ys)/nrow(X),digits=3)
"""

#   if(is.numeric(target))
#   {
#       if(regMethod == "median"){
#         ysMost = median(ys)
#       }else{
#         ysMost <- mean(ys)
#       }
#       err <- sum((ysMost - ys)^2)/length(ys)   
#   }else{ 
#     if(length(pred)>0){ #if pred of the rule is provided
#       ysMost = as.character(pred)
#     }else{
#       ysMost <- names(which.max(  table(ys))) # get back the first max
#     }
#     ly <- sum(as.character(ys)==ysMost)
#     conf <- round(ly/length(ys),digits=3)    
#     err <- 1 - conf
#   }
#   rule <- origRule

#   v <- c(len, freq, err, rule, ysMost)
#   names(v) <- c("len","freq","err","condition","pred")
#   return(v)

# measureRule <-
# function(ruleExec,X,target,pred=NULL,regMethod="mean"){
#   len <- length(unlist(strsplit(ruleExec, split=" & ")))
#   origRule <- ruleExec
#   ruleExec <- paste("which(", ruleExec, ")")
#   ixMatch <- eval(parse(text=ruleExec)) 
#   if(length(ixMatch)==0){
#     v <- c("-1","-1", "-1", "", "")
#     names(v) <- c("len","freq","err","condition","pred")
#     return(v)
#   }
#   ys <- target[ixMatch]
#   freq <- round(length(ys)/nrow(X),digits=3)

#   if(is.numeric(target))
#   {
#       if(regMethod == "median"){
#         ysMost = median(ys)
#       }else{
#         ysMost <- mean(ys)
#       }
#       err <- sum((ysMost - ys)^2)/length(ys)   
#   }else{ 
#     if(length(pred)>0){ #if pred of the rule is provided
#       ysMost = as.character(pred)
#     }else{
#       ysMost <- names(which.max(  table(ys))) # get back the first max
#     }
#     ly <- sum(as.character(ys)==ysMost)
#     conf <- round(ly/length(ys),digits=3)    
#     err <- 1 - conf
#   }
#   rule <- origRule

#   v <- c(len, freq, err, rule, ysMost)
#   names(v) <- c("len","freq","err","condition","pred")
#   return(v)
# }

# ---------------------------------------------------------------------------- #
function measure_rule(
    rule_exec::AbstractString,
    X,
    target::AbstractVector;
    pred=nothing,
    reg_method=:mean,
    rule2table::Function,
)
    clause_count = isempty(rule_exec) ? 0 : length(split(rule_exec, " & "))
    ix_match = _matching_indices(rule_exec, X, target, rule2table)

    isempty(ix_match) && return _empty_metric()

    ys = target[ix_match]
    freq = round(length(ys) / size(X, 1); digits=3)

    if _is_numeric_target(target)
        ys_most =
            (reg_method == :median || reg_method == "median") ?
            median(ys) :
            mean(ys)

        err = sum(abs2, ys .- ys_most) / length(ys)
    else
        # Matches R: use the supplied rule prediction when present;
        # otherwise use the modal target value.
        ys_most = pred === nothing ? _modal_label(ys) : string(pred)

        correct = count(y -> string(y) == ys_most, ys)
        confidence = round(correct / length(ys); digits=3)
        err = 1 - confidence
    end

    # A NamedTuple is the Julia equivalent of R's named vector.
    return (
        len = clause_count,
        freq = freq,
        err = err,
        condition = String(rule_exec),
        pred = ys_most,
    )
end

# function prune_single_rule(
#     rule,
#     X,
#     target::AbstractVector,
#     max_decay::Real,
#     type_decay::Integer;
#     rule2table::Function,
# )
#     condition = String(rule[:condition])
#     pred = rule[:pred]

#     # Deliberately does not pass `pred`, matching the R implementation.
#     new_metric = measure_rule(
#         condition,
#         X,
#         target;
#         rule2table=rule2table,
#     )

#     err_orig = new_metric.err
#     rule_parts = split(condition, " & ")

#     length(rule_parts) <= 1 && return new_metric

#     # Keep `err_orig` fixed throughout, as in the R implementation.
#     for i in length(rule_parts):-1:1
#         remaining_parts = copy(rule_parts)
#         deleteat!(remaining_parts, i)

#         rest_rule = join(remaining_parts, " & ")
#         metric_tmp = measure_rule(
#             rest_rule,
#             X,
#             target;
#             pred=pred,
#             rule2table=rule2table,
#         )

#         err_new = metric_tmp.err

#         decay = if type_decay == 1
#             (err_new - err_orig) / max(err_orig, 1e-6)
#         else
#             err_new - err_orig
#         end

#         if decay <= max_decay
#             rule_parts = remaining_parts
#             new_metric = metric_tmp

#             length(rule_parts) <= 1 && break
#         end
#     end

#     return new_metric
# end

# function _metrics_dataframe(metrics)
#     lengths = Int[]
#     frequencies = Float64[]
#     errors = Float64[]
#     conditions = String[]
#     predictions = Any[]

#     for metric in metrics
#         push!(lengths, metric.len)
#         push!(frequencies, metric.freq)
#         push!(errors, metric.err)
#         push!(conditions, metric.condition)
#         push!(predictions, metric.pred)
#     end

#     return DataFrame(
#         len=lengths,
#         freq=frequencies,
#         err=errors,
#         condition=conditions,
#         pred=predictions,
#     )
# end

# """
#     prune_rule(rules, X, target, max_decay=0.05, type_decay=2;
#                rule2table)

# Julia equivalent of R's `pruneRule`.
# """
# function prune_rule(
#     rules::AbstractDataFrame,
#     X,
#     target::AbstractVector,
#     max_decay::Real=0.05,
#     type_decay::Integer=2;
#     rule2table::Function,
# )
#     all(name -> name in propertynames(rules), (:condition, :pred)) ||
#         throw(ArgumentError("rules must contain :condition and :pred columns"))

#     metrics = [
#         prune_single_rule(
#             rule,
#             X,
#             target,
#             max_decay,
#             type_decay;
#             rule2table=rule2table,
#         )
#         for rule in eachrow(rules)
#     ]

#     # This replaces repeatedly calling R's rbind().
#     return _metrics_dataframe(metrics)
# end

# _select_rules_cbc(set, X, y, extractor)
# ruleset = ClassificationRule{String}[
#   ▣ ([petal_length] ∈ [-Inf,2.45))  ↣  setosa
# , ▣ ([petal_width] ∈ [-Inf,0.75))  ↣  setosa
# , ▣ ([petal_length] ∈ [2.45,4.95))  ↣  versicolor
# , ▣ ([petal_width] ∈ [0.75,1.55))  ↣  versicolor
# , ▣ ([petal_width] ∈ [1.55,1.75))  ↣  versicolor
# , ▣ ([petal_width] ∈ [1.75,1.85))  ↣  versicolor
# , ▣ (([petal_length] ∈ [4.95,5.05))) ∧ (([sepal_length] ∈ [-Inf,6.5)))  ↣  virginica
# ]

R"""
library("RRF")
set.seed(1)

ruleI = sapply(pruneRules[,"condition"],rule2Table,Xc,yc)
print(exec[,"condition"])
coefReg <- 0.95 - 0.01*as.numeric(pruneRules[,"len"])/max(as.numeric(pruneRules[,"len"]))
rf <- RRF(ruleI,as.factor(yc), flagReg = 1, coefReg=coefReg, mtry = (ncol(ruleI)*1/2) , ntree=50, maxnodes= 10,replace=FALSE) 
imp <- rf$importance/max(rf$importance)
feaSet <- which(imp > 0.01)
ruleSetPrunedRRF <- cbind(pruneRules[feaSet,,drop=FALSE],impRRF=imp[feaSet])
ix = order(as.numeric(ruleSetPrunedRRF[,"impRRF"]),
            - as.numeric(ruleSetPrunedRRF[,"err"]),
            - as.numeric(ruleSetPrunedRRF[,"len"]),
            decreasing=TRUE)
ruleSelect <- ruleSetPrunedRRF[ix,,drop=FALSE]
"""

# Extract rules from model
model = solem_rf
listrules_kwargs = (use_shortforms=true, normalize=true)
set = SolePostHoc.RuleExtraction._starterruleset(model; listrules_kwargs...)


# function select_rule_rrf(
#     prune_rules::ClassificationRule{T},
#     X,
#     y;
#     rule2table::Function,
#     rrf_fit::Function,
#     rng::AbstractRNG = MersenneTwister(1),
# ) where T
    if nrow(prune_rules) == 0
        result = copy(prune_rules)
        result[!, :impRRF] = Float64[]
        return result
    end

    rule_columns = rule2table.(
        prune_rules[!, :condition],
        Ref(X),
        Ref(y),
    )
    ruleI = Float64.(hcat(rule_columns...))

    # R: 0.95 - 0.01 * len / max(len)
    lengths = Float64.(prune_rules[!, :len])
    coef_reg = 0.95 .- 0.01 .* lengths ./ maximum(lengths)

    # `y` does not need explicit factor conversion if rrf_fit handles labels.
    # R's ncol(ruleI) * 1 / 2 is converted to an integer by RRF internally.
    rf = rrf_fit(
        ruleI,
        y;
        flagReg = 1,
        coefReg = coef_reg,
        mtry = size(ruleI, 2) ÷ 2,
        ntree = 50,
        maxnodes = 10,
        replace = false,
        rng = rng,
    )

    raw_importance = Float64.(vec(rf.importance))
    length(raw_importance) == size(ruleI, 2) ||
        throw(ArgumentError("RRF importance must contain one value per rule"))

    max_importance = maximum(raw_importance)
    imp = iszero(max_importance) ?
          zeros(Float64, length(raw_importance)) :
          raw_importance ./ max_importance

    # R: which(imp > 0.01)
    feature_set = findall(>(0.01), imp)

    # R: cbind(pruneRules[feaSet, , drop = FALSE], impRRF = imp[feaSet])
    rule_set_pruned_rrf = copy(prune_rules[feature_set, :])
    selected_importance = imp[feature_set]
    rule_set_pruned_rrf[!, :impRRF] = selected_importance

    # Equivalent to:
    # order(impRRF, -err, -len, decreasing = TRUE)
    # i.e. importance descending, error ascending, length ascending.
    err = Float64.(rule_set_pruned_rrf[!, :err])
    len = Float64.(rule_set_pruned_rrf[!, :len])
    ix = sortperm(
        eachindex(selected_importance);
        by = i -> (-selected_importance[i], err[i], len[i], i),
    )

    return rule_set_pruned_rrf[ix, :]
# end

# function check(rules::Rule{T}, X::Vector, featurenames::vector{T}) where T
rules = set[2]
rules isa Rule{T} where T
Xm = Matrix(Xc)
X = Xm[2,:]
featurenames = names(Xc)

for a in atoms(rules)
    fidx = findfirst(==(i_name(a)), Symbol.(featurenames))
    checkcondition(a, X[fidx])
end

a = atoms(rules)[1]

checkcondition(set, X, featurenames)
checkcondition(set[5], Xm, featurenames)
# ---------------------------------------------------------------------------- #
value(a::Atom) = a.value
i_name(a::Atom) = a.value.feature.i_name

checkcondition(r::AbstractCondition, featval::Real) =
    error("Please, provide method for $(typeof(r)), $(typeof(featval))")

@inline function checkcondition(r::RangeScalarCondition, featval::Real)
    honors_minval(r, featval) && honors_maxval(r, featval)
end

@inline checkcondition(a::Atom, featval::Real) =
    checkcondition(value(a), featval)

@inline function checkcondition(
    rule::ClassificationRule{T},
    x::AbstractVector, 
    featurenames::Vector{<:Union{String,Symbol}}
) where T
    return all(atoms(rule)) do a
        fidx = findfirst(==(T.(i_name(a))), T.(featurenames))
        checkcondition(value(a), x[fidx])
    end
end

@inline function checkcondition(
    rule::ClassificationRule{T},
    X::AbstractArray{S}, 
    args...
) where {T,S}
    return [checkcondition(rule, x, args...) for x in eachrow(X)]
end

@inline function checkcondition(
    set::Vector{ClassificationRule{T}},
    args...
) where T
    return [checkcondition(r, args...) for r in set]
end

# ---------------------------------------------------------------------------- #

function rule_matches(rule, x::AbstractVector, featurenames)
    feature_symbols = Symbol.(featurenames)

    return all(atoms(rule)) do atom
        feature_index = findfirst(==(i_name(atom)), feature_symbols)

        feature_index === nothing && throw(ArgumentError(
            "Feature $(repr(i_name(atom))) does not occur in featurenames",
        ))

        checkcondition(value(atom), x[feature_index])
    end
end

matches = rule_matches(rules, X, featurenames)