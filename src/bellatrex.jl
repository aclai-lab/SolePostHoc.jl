using SoleModels
using MultivariateStats
using Clustering
using NearestNeighbors
using LinearAlgebra: norm
using Random
using ModalDecisionTrees: DForest

# hyperparameter: ntrees, nclusters

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
    ndatasetinstances = ninstances(X)
    ruleset = begin
        if m isa DecisionForest
            unique([listrules(tree; use_shortforms=true)for tree in ftrees])
        else
            listrules(model)
        end
    end

    # 1) For each instance in X, we would extract the most adaptive rules for it
    for idx in 1:ndatasetinstances

        println("Running instances $(idx)/$(ndatasetinstances)")
        currentinstance = slicedataset(X,[idx]; return_view=true)

        # 2) Only trees that correctly predict the considered instance are considered
        println("Concording phase in...")
        ctrees, idxsctrees = @time
            concordingtrees(ftrees, currentinstance, Y[idx])

        # 3) Representing rules as a vector
        println("Rule to vector phase in...")
        vrules = @time begin
            rtrees = begin
                rs = ruleset[idxsctrees]
                rs = rs isa Vector{<:Vector{<:Any}} ? reduce(vcat,rs) : rs
                rs[antecedenttops(antnode,currentinstance)]
            end

            [rule2vector(t,X) for t in rtrees]
        end

        # https://juliastats.org/MultivariateStats.jl/dev/pca/
        # 4) Vector projection (one of them):
        #     - Principal Component Analysis (PCA)
        #     - Multidimensional Scaling (MDS)
        println("Projection phase in...")
        prules = @time projectionvectors(hcat(vrules...),X,Y; kwargs...)

        # https://juliastats.org/Clustering.jl/dev/kmeans.html
        # 5) Clustering using k-means++ method
        println("Clustering phase in...")
        krules, cluster_sizes = @time clusteringvectors(prules; rng=rng, kwargs...)

        # 6) Surrogate model: combination between rules and instance
        # KDTree construction: https://github.com/KristofferC/NearestNeighbors.jl
        println("Final phase: k-nearest neighbors in...")
        final_rules_idx = @time begin
            kdrules = KDTree(prules, Minkowski(3.5); leafsize = 40)
            knn_idxs, _ = knn(kdrules, krules, 1, sortres = true)

            vcat(knn_idxs)
        end

        # https://github.com/Klest94/Bellatrex/blob/main/code_scripts/LocalMethod_class.py row 241
        final_trees_idx = [idxsctrees[k] for k in final_rules_idx]
        #localprediction(final_rules_idx,cluster_sizes)
    end
end

############################################################################################
########################### Supporting functions of BELLATREX ##############################
############################################################################################

function concordingtrees(
    m::Vector{<:DecisionTree},
    X::AbstractLogiset,
    y::Label;
    ntrees::Integer=20,
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

    preselectionidxs = sortperm(treesloss)[1:ntrees]

    #return m[y_preds .== y]
    return (m[preselectionidxs], preselectionidxs)
end

                    #############################################
                    #############################################

# Features extraction: returns all the features in the SyntaxTree
featuresextract(m::MultiFormula) =
    reduce(vcat,[featuresextract(modant) for (_,modant) in modforms(m)])
featuresextract(m::SyntaxTree) = i_variable.(feature.(atom.(propositions(m))))

function rule2vector(m::Rule,X::AbstractLogiset)
    rulevector = zeros(nfeatures(X))
    datanextnodes = rule2vector(antecedent(m),X)

    icovariate = countmap(
                    datanextnodes[:featuresnodes],
                    datanextnodes[:nfallinstances]/ninstances(X)
                )

    return rulevector[1:length(icovariate)] += icovariate
end

function rule2vector(
    m::M,
    X::AbstractLogiset,
) where {M<:Union{MultiFormula,LeftmostConjunctiveForm}}
    featuresnodes = []
    nfallinstances = []
    idxsremaininstances = collect(1:ninstances(X))

    conjs = begin
        if m isa MultiFormula
            [modant for (_,modant) in modforms(m)]
        elseif m isa LeftmostConjunctiveForm
            children(m)
        end
    end

    for conj in conjs
        truthvalues = antecedenttops(conj,X)

        # features and ninstances fall in the current node
        push!(featuresnodes,featuresextract(conj))
        idxsfallinstances = intersect(findall(truthvalues .== true),idxsremaininstances)
        push!(nfallinstances,length(idxsfallinstances))

        # updating remain instances
        idxsremaininstances = idxsfallinstances
    end

    return (; featuresnodes=featuresnodes, nfallinstances=nfallinstances)
end

                    #############################################
                    #############################################

function projectionvectors(
    matrixrules::Array{<:Number, 2},
    X::AbstractLogiset,
    Y::AbstractVector{<:Label};
    projection_method::Symbol = :MDS,
    kwargs...,
)
    maxoutdim = min(size(matrixrules)...)

    space = begin
        if projection_method == :MDS
            M = fit(MDS, matrixrules; maxoutdim=maxoutdim, distances=false)
            predict(M)
        elseif projection_method == :PCA
            # Xtr is adjoint training dataset because each column is an observation
            # Maybe maxoutdims is the number of classes in Y
            M = fit(PCA, matrixrules; maxoutdim=maxoutdim)
            predict(M, X)
        else
            error("$(projection_method) is not in possible methods")
        end
    end

    return space
end

                    #############################################
                    #############################################

function clusteringvectors(
    matrixtrees::Array{<:Number, 2};
    rng::AbstractRNG = MersenneTwister(1),
    nclusters::Integer = 2,
    kwargs...,
)
    R = kmeans(matrixtrees, nclusters; maxiter=5000, display=:iter, rng=rng)

    # .centers gets the cluster centers
    return R.centers, counts(R)
end

                    #############################################
                    #############################################
#=
function localprediction(
    final_trees_idx::Vector{Integer},
    cluster_sizes::Vector{Integer},
)
    @assert length(final_trees_idx) == length(cluster_sizes) "$(final_trees_idx) !=" *
        " $(cluster_sizes)"

    bestguess
end
=#
