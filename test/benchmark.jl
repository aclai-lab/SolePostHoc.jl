using Test

using SolePostHoc
const SP = SolePostHoc

using SoleXplorer
const SX = SoleXplorer

using BenchmarkTools
using RDatasets
using DataFrames
using Random

# ---------------------------------------------------------------------------- #
#                                 iris test                                    #
# ---------------------------------------------------------------------------- #
iris_df = dataset("datasets", "iris")
y = String.(iris_df.Species)
Xdf = select(iris_df, Not(:Species))

dsc = setup_dataset(Xdf, y; model=SX.RandomForestClassifier(n_trees=20), rng=42)
modelc = solexplorer(dsc)
solem = modelc.sole[1]

config = SP.Lumen.LumenConfig(minimization_scheme=:mitespresso)
SP.Lumen.super_lumen(config, solem, M=10000, N=10)

# # ---------------------------------------------------------------------------- #
# #                                caravan test                                  #
# # ---------------------------------------------------------------------------- #
# df = dataset("ISLR", "Caravan")
# y = String.(df.Purchase)
# Xdf = select(df, Not(:Purchase))
# X = Matrix(Xdf)
# vnames = String.(names(Xdf))

# dt = SX.setup_dataset(
#     X, y;
#     model=SX.DecisionTreeClassifier(),
#     vnames,
#     rng=Xoshiro(42)
# )
# classlabels = unique(y)
# X_train = SX.get_X(dt, :train)
# X_train = Matrix(hcat(X_train...))
# y_train = Vector(get_y(dt, :train)[1])
# X_test = SX.get_X(dt, :test)
# X_test = Matrix(hcat(X_test...))
# y_test = Vector(get_y(dt, :test)[1])

# # ---------------------------------------------------------------------------- #
# #                                  pima test                                   #
# # ---------------------------------------------------------------------------- #
# pima_df = dataset("MASS", "Pima.tr")
# y = String.(pima_df.Type)
# Xdf = select(pima_df, Not(:Type))
# X = Matrix(Xdf)
# vnames = String.(names(Xdf))

# dt = SX.setup_dataset(
#     X, y;
#     model=SX.DecisionTreeClassifier(),
#     vnames,
#     rng=Xoshiro(42)
# )
# classlabels = unique(y)
# X_train = SX.get_X(dt, :train)
# X_train = Matrix(hcat(X_train...))
# y_train = Vector(get_y(dt, :train)[1])
# X_test = SX.get_X(dt, :test)
# X_test = Matrix(hcat(X_test...))
# y_test = Vector(get_y(dt, :test)[1])

# # ---------------------------------------------------------------------------- #
# #                                 boston test                                  #
# # ---------------------------------------------------------------------------- #
# boston = dataset("MASS", "Boston")
# yr = boston.MedV
# Xr = select(boston, Not(:MedV))
# X = Matrix{Float32}(Xr)
# y = Float32.(yr)
# vnames = String.(names(Xr))

# dt = SX.setup_dataset(
#     X, y;
#     model=SX.DecisionTreeRegressor(),
#     vnames,
#     rng=Xoshiro(42)
# )
# classlabels = unique(y)
# X_train = SX.get_X(dt, :train)
# X_train = Matrix(hcat(X_train...))
# y_train = Vector(get_y(dt, :train)[1])
# X_test = SX.get_X(dt, :test)
# X_test = Matrix(hcat(X_test...))
# y_test = Vector(get_y(dt, :test)[1])
