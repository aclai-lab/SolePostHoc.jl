using SoleModels
using MultivariateStats
using Clustering
using NearestNeighbors
using LinearAlgebra: norm
using Random

# hyperparameter: ntrees, nclusters

function concordingtrees(
    m::Vector{<:DecisionTree},
    X::AbstractLogiset,
    y::Label;
    ntrees::Union{Nothing,Integer}=nothing,
    kwargs...,
)
    y_preds = [apply(t,X) for t in m]

    treesloss = begin
        if y isa CLabel
            # Classification problem: find an alternative of euclidean norm
            y_preds
        else
            # Regression problem: euclidean norm
            [norm(y-y_preds[k]) for k in eachindex(y_preds)]
        end
    end

    preselectionidxs = sortperm(treesloss)[1:ntrees]

    #return m[y_preds .== y]
    return (m[preselectionidxs], preselectionidxs)
end

# To see: tree_splits_to_vector() in Bellatrex/code_scripts/TreeRepresentation_utils.py
function tree2vector(
    m::DecisionTree
)
end

function projectionvectors(
    matrixtrees::Array{Number, 2};
    projection_method::Symbol = :MDS,
    kwargs...,
)
    Xtrees = hcat(vtrees...)
    maxoutdim = min(size(Xtrees)...)

    space = begin
        if projection_method == :MDS
            # Xtr is adjoint training dataset because each column is an observation
            # Maybe maxoutdims is the number of classes in Y
            M = fit(PCA, Xtrees; maxoutdim=maxoutdim)
            predict(M, X)
        elseif projection_method == :PCA
            M = fit(MDS, Xtrees; maxoutdim=maxoutdim, distances=false)
            predict(M)
        else
            error("$(projection_method) is not in possible methods")
        end
    end

    return space
end

function clusteringvectors(
    matrixtrees::Array{Number, 2};
    rng::AbstractRNG = MersenneTwister(1),
    nclusters::Integer = 2,
    kwargs...,
)
    R = kmeans(matrixtrees, nclusters; maxiter=5000, display=:iter, rng=rng)

    # .centers gets the cluster centers
    return R.centers
end

# Repository web path: https://github.com/Klest94/Bellatrex
# Bellatrex is a a
function bellatrex(
    m::DecisionForest,
    X::AbstractLogiset, # testing dataset
    Y::AbstractVector{<:Label};
    rng::AbstractRNG = MersenneTwister(1),
    kwargs...,
)
    # 1) For each instance in X, we would extract the most adaptive rules for it
    for idx in 1:ninstances(X)

        # 2) Only trees that correctly predict the considered instance are considered
        ctrees, idxsctrees =
            concordingtrees(trees(m), slicedataset(X,[idx]; return_view=true), y[idx])

        # 3) Representing trees as a vector
        vtrees = tree2vector.(ctrees)

        # https://juliastats.org/MultivariateStats.jl/dev/pca/
        # 4) Vector projection (one of them):
        #     - Principal Component Analysis (PCA)
        #     - Multidimensional Scaling (MDS)
        ptrees = projectionvectors(vtrees; kwargs...)

        # https://juliastats.org/Clustering.jl/dev/kmeans.html
        # 5) Clustering using k-means++ method
        ktrees = clusteringvectors(ptrees; rng=rng, kwargs...)

        # 6) Surrogate model: combination between rules and instance
        # KDTree construction: https://github.com/KristofferC/NearestNeighbors.jl
        kdtrees = KDTree(ptrees, Minkowski(3.5); leafsize = 40)
        knn_idxs, _ = knn(kdtrees, ktrees, 1, sortres = true)
        flatten_knn_idxs = vcat(knn_idxs)
        final_trees_idx = [idxsctrees[k] for k in flatten_knn_idxs]

        # https://github.com/Klest94/Bellatrex/blob/main/code_scripts/LocalMethod_class.py row 241
    end
end
