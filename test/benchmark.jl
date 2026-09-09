using Test

using SolePostHoc
const SP = SolePostHoc

using SoleXplorer
const SX = SoleXplorer

using SoleModels
const SM = SoleModels

using BenchmarkTools
using RDatasets
using DataFrames
using Random

using InteractiveUtils
using JET

using CategoricalArrays

# ---------------------------------------------------------------------------- #
#                                 iris test                                    #
# ---------------------------------------------------------------------------- #
iris_df = dataset("datasets", "iris")
y = String.(iris_df.Species)
Xdf = select(iris_df, Not(:Species))

dsc = setup_dataset(Xdf, y; model=SX.RandomForestClassifier(n_trees=20, max_depth=5), rng=42)
modelc = solexplorer(dsc)
solem = modelc.sole[1]

# config = SP.Lumen.LumenConfig(; minimization_scheme=:mitespresso)
# SP.Lumen.super_lumen(config, solem, M=100000, N=10);

config = SP.Lumen.LumenConfig(; minimization_scheme=Abc, M=5000)
lumen = SP.Lumen.lumen_shannon(config, solem);
@btime SP.Lumen.lumen_shannon(config, solem);
# 3.915 s (15192644 allocations: 801.34 MiB)
# 2.043 s (15223660 allocations: 755.99 MiB)
# 2.087 s (14735030 allocations: 745.08 MiB)
# 2.059 s (14733296 allocations: 745.04 MiB)
# 2.054 s (14772318 allocations: 744.97 MiB)
# 2.064 s (14759051 allocations: 744.59 MiB)
# 2.057 s (14759046 allocations: 744.59 MiB)

@code_warntype SP.Lumen.lumen_shannon(config, solem)
result = @report_opt SP.Lumen.lumen_shannon(config, solem)
open("jet_report.txt", "w") do io
    show(io, MIME"text/plain"(), result)
end

ctx, per_class_terms, config = SP.Lumen.lumen_shannon(config, solem);

@code_warntype SP.Lumen._finalize_decision_set(ctx, per_class_terms, config)
result = @report_opt SP.Lumen._extract_atoms_bfs_order(model)
open("jet_report.txt", "w") do io
    show(io, MIME"text/plain"(), result)
end


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
