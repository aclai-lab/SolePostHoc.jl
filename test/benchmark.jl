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

config = SP.Lumen.LumenConfig(; minimization_scheme=:abc, M=5000)
lumen = SP.Lumen.lumen_shannon(config, solem)
@btime SP.Lumen.lumen_shannon(config, solem);
# 3.915 s (15192644 allocations: 801.34 MiB)
# typed float_type in config
# 2.416 s (15192236 allocations: 801.31 MiB)
# avoid get method for internal use
# 2.395 s (15192236 allocations: 801.31 MiB)
# config Float32 by default (but propably it doesn't pass)
# 2.452 s (15192236 allocations: 801.31 MiB)
# M and max_apply_batch in config
# 2.410 s (15208042 allocations: 798.60 MiB)
# passing branch instead AbstractModel
# 3.316 s (10240912 allocations: 411.73 MiB) # laptop on batteries
@code_warntype SP.Lumen.lumen_shannon(config, solem)
result = @report_opt SP.Lumen.lumen_shannon(config, solem)
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

# ---------------------------------------------------------------------------- #
#                             apply graph walking                              #
# ---------------------------------------------------------------------------- #
function apply(
    tree::MetaGraph{R,G},
    v::R,
    instance::SubArray{T, 1}
) where {R<:Unsigned,T<:AbstractFloat,G<:SimpleGraph}
    while length(tree.graph.fadjlist[v]) > 1 
        dir = instance[tree[v].feature] < tree[v].threshold ?
            1 : 2
        v = tree.graph.fadjlist[v][dir+1]
    end

    return tree[v].outcome
end

function apply(
    ::Type{L},
    tree::MetaGraph{R,G},
    X::Matrix{T}
) where {L<:ClassificationLoss,R<:Unsigned,T<:AbstractFloat,G<:SimpleGraph}
    predictions = Vector{R}(undef, size(X,1))

    for i in axes(X, 1)
        predictions[i] = R(apply(tree, one(R), @view X[i, :]))
    end

    return predictions
end

function apply(
    ::Type{L},
    tree::MetaGraph{R,G},
    X::Matrix{T}
) where {L<:RegressionLoss,R<:Unsigned,T<:AbstractFloat,G<:SimpleGraph}
    predictions = Vector{T}(undef, size(X,1))

    for i in axes(X, 1)
        predictions[i] = apply(tree,  one(R), @view X[i,:])
    end

    return predictions
end

function apply(trees::Vector{SM.Branch{S}}, args...) where{S<:SM.Labels}
    ntrees = length(trees)

    for t in 1:ntrees

    end
end

config = SP.Lumen.LumenConfig(; minimization_scheme=:abc, M=5000)
model = solem
trees = SM.models(solem)
featurenames = String.(SM.info(solem, :featurenames))
classnames = String.(unique!(SM.info(solem, :supporting_labels)))
depth = 1.0

atoms = unique!(SP.Lumen._normalize_atom.(if depth < 1.0
    mapreduce(vcat, SM.models(model); init=SL.Atom{SD.ScalarCondition}[]) do t
        _take_first_percentage(_extract_atoms_bfs_order(t), depth)
    end
else
    SM.atoms(SM.alphabet(model, false))
end))

btoms = unique!(SP.Lumen._normalize_atom.(if depth < 1.0
    mapreduce(vcat, trees; init=SL.Atom{SD.ScalarCondition}[]) do t
        _take_first_percentage(_extract_atoms_bfs_order(t), depth)
    end
else
    collect(Iterators.flatten(SM.alphabet.(trees, false)))
end))

features = SM.featurename.(unique!(SP.Lumen.get_feature.(atoms)))
beatures = String.(SM.featurename.(unique!(SP.Lumen.get_feature.(btoms))))

family = SP.Lumen._feature_op_family(atoms, features[1])
family = SP.Lumen._feature_op_family(btoms, beatures[1])

feat_atoms = SP.Lumen._atoms_for_feature(atoms, features[1])
ops = unique(get_operator.(feat_atoms))

beat_atoms = SP.Lumen._atoms_for_feature(btoms, beatures[1])
ops = unique(get_operator.(feat_atoms))

ctx = SP.Lumen._prepare_sequential_context(config, trees, featurenames, classnames);
ctxTest = SP.Lumen._prepare_sequential_context(config, solem);

