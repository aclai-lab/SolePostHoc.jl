using Test

using SolePostHoc

using MLJ
using MLJBase
using DataFrames

using MLJDecisionTreeInterface
using SoleModels

import DecisionTree as DT

X, y = @load_iris
X = DataFrame(X)

train_ratio = 0.8

train, test = partition(eachindex(y), train_ratio, shuffle=true)
X_train, y_train = X[train, :], y[train]
X_test, y_test = X[test, :], y[test]

println("Training set size: ", size(X_train), " - ", size(y_train))
println("Test set size: ", size(X_test), " - ", size(y_test))
println("Training set type: ", typeof(X_train), " - ", typeof(y_train))
println("Test set type: ", typeof(X_test), " - ", typeof(y_test))

Forest = MLJ.@load RandomForestClassifier pkg=DecisionTree

model = Forest(
  max_depth=3,
  min_samples_leaf=1,
  min_samples_split=2,
  n_trees = 10,
)

# Bind the model and data into a machine
mach = machine(model, X_train, y_train)
# Fit the model
fit!(mach)


classlabels = (mach).fitresult[2]
classlabels = classlabels[sortperm((mach).fitresult[3])]
featurenames = report(mach).features
solem = solemodel(fitted_params(mach).forest; classlabels, featurenames)
solem = solemodel(fitted_params(mach).forest; classlabels, featurenames, keep_condensed = false)


SolePostHoc.RuleExtraction.intrees(solem)

bellatrex()

