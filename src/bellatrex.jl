using SoleModels
using MultivariateStats
using Clustering
using NearestNeighbors
using LinearAlgebra: norm
using Random
using ModalDecisionTrees: DForest

# hyperparameter: ntrees, nclusters

function concordingtrees(
    m::Vector{<:DecisionTree},
    X::AbstractLogiset,
    y::Label;
    nrules::Integer=20,
    kwargs...,
)
    y_preds = [apply(t,X) for t in m]

    treesloss = begin
        if y isa CLabel
            # Classification problem: find an alternative of euclidean norm
            y_preds .== y
        else
            # Regression problem: euclidean norm
            [norm(y-y_preds[k]) for k in eachindex(y_preds)]
        end
    end

    preselectionidxs = sortperm(treesloss)[1:nrules]

    #return m[y_preds .== y]
    return (m[preselectionidxs], preselectionidxs)
end

# To see: tree_splits_to_vector() in Bellatrex/code_scripts/TreeRepresentation_utils.py
function rule2vector(
    m::Rule{O,<:LeftmostConjunctiveForm},
    X::AbstractLogiset;
    maxconjuncts,
    kwargs...,
) where {O}
    infos = info(m)
    sumweighted = 0
    cons = consequent(m)
    nodes = conjuncts(m)
    nnodes = nconjuncts(m)
    totinstances = ninstances(X)
    fallinstances = collect(1:totinstances)

    transvector = Vector{AbstractFloat}(undef,maxconjuncts)
    for inode in 1:nnodes
        rnode = inode == 1 ?
                    Rule(nodes[inode],cons,infos) :
                    Rule(LeftmostConjunctiveForm(nodes[1:inode]),cons,infos)

        truthvalues = antecedenttops(rnode,slicedataset(X,fallinstances; return_view=true))
        idxstrues = findall(truthvalues .== true)

        setdiff!(fallinstances,idxstrues)
        sumweighted += (length(idxstrues)/totinstances)
        transvector[inode] = sumweighted * inode
    end
    transvector[
        filter(idx-> isassigned(transvector,idx) == false,1:length(transvector))
    ] .= 0

    return transvector
end

function projectionvectors(
    matrixrules::Array{<:Number, 2},
    X::AbstractLogiset,
    Y::AbstractVector{<:Label};
    projection_method::Symbol = :MDS,
    kwargs...,
)
    maxoutdim = min(size(matrixrules)...,length(unique(Y)))

    space = begin
        if projection_method == :MDS
            # Xtr is adjoint training dataset because each column is an observation
            # Maybe maxoutdims is the number of classes in Y
            M = fit(PCA, matrixrules; maxoutdim=maxoutdim)
            predict(M, X)
        elseif projection_method == :PCA
            M = fit(MDS, matrixrules; maxoutdim=maxoutdim, distances=false)
            predict(M)
        else
            error("$(projection_method) is not in possible methods")
        end
    end

    return space
end

function clusteringvectors(
    matrixtrees::Array{<:Number, 2};
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
    m::Union{AbstractModel,DecisionForest,DForest},
    X::AbstractLogiset, # testing dataset
    Y::AbstractVector{<:Label};
    rng::AbstractRNG = MersenneTwister(1),
    kwargs...,
)
    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    ftrees = trees(m)
    ruleset = begin
        if m isa DecisionForest
            unique([listrules(tree; use_shortforms=false, use_leftmostlinearform=true) for tree in ftrees])
        else
            listrules(model)
        end
    end

    # 1) For each instance in X, we would extract the most adaptive rules for it
    for idx in 1:ninstances(X)

        println("Concording phase")
        # 2) Only trees that correctly predict the considered instance are considered
        ctrees, idxsctrees =
            concordingtrees(ftrees, slicedataset(X,[idx]; return_view=true), Y[idx])

        # Translate trees into rules
        rtrees = begin
            rs = ruleset[idxsctrees]
            rs isa Vector{<:Vector{<:Any}} ? reduce(vcat,rs) : rs
        end

        println("Rule to vector phase")
        # 3) Representing trees as a vector
        maxconjuncts = max(nconjuncts.(rtrees)...)
        vrules = [rule2vector(t,X; maxconjuncts=maxconjuncts) for t in rtrees]

        # https://juliastats.org/MultivariateStats.jl/dev/pca/
        # 4) Vector projection (one of them):
        #     - Principal Component Analysis (PCA)
        #     - Multidimensional Scaling (MDS)
        println("Projection phase")
        prules = projectionvectors(hcat(vrules...),X,Y; kwargs...)

        # https://juliastats.org/Clustering.jl/dev/kmeans.html
        # 5) Clustering using k-means++ method
        println("Clustering phase")
        krules = clusteringvectors(prules; rng=rng, kwargs...)

        # 6) Surrogate model: combination between rules and instance
        # KDTree construction: https://github.com/KristofferC/NearestNeighbors.jl
        println("Final phase")
        kdrules = KDTree(prules, Minkowski(3.5); leafsize = 40)
        knn_idxs, _ = knn(kdrules, krules, 1, sortres = true)
        final_rules_idx = vcat(knn_idxs)

        # https://github.com/Klest94/Bellatrex/blob/main/code_scripts/LocalMethod_class.py row 241
        @show final_rules_idx
    end
end
