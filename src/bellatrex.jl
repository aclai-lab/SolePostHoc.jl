using SoleModels
using MultivariateStats
using Clustering
using NearestNeighbors
using LinearAlgebra: norm
using Random
using Hyperopt
using ModalDecisionTrees: DForest

# hyperparameter: ntrees, nclusters

"""
Dedja, Klest, et al. "BELLATREX: Building explanations through a locally accurate rule extractor." Ieee Access 11 (2023): 41348-41367.

See also [`extractrules`](@ref), [`bellatrex`](@ref), [`RuleExtractor`](@ref).
"""
struct BellatrexRuleExtractor <: RuleExtractor end

extractrules(::BellatrexRuleExtractor, m, args...; kwargs...) = bellatrex(m, args...; kwargs...)

# Repository web path: https://github.com/Klest94/Bellatrex
# Bellatrex is a a
function bellatrex(
    m::Union{AbstractModel,DecisionForest,DForest},
    X::AbstractLogiset, # testing dataset
    Y::AbstractVector{<:Label};
    exec_ntrees::Vector{Float64}=[0.2,0.5,0.8],
    exec_ndims::Vector=[2,5,nothing],
    exec_nclusters::Vector{Int64}=[1,2,3],
    rng::AbstractRNG = MersenneTwister(1),
    kwargs...,
)
    ########################################################################################
    # Extract rules from each tree, obtain full ruleset
    ########################################################################################
    ftrees = trees(m)
    ndatasetinstances = ninstances(X)
    final_predictions = []
    ruleset = begin
        if isensemble(m)
            unique([listrules(tree; use_shortforms=true) for tree in ftrees])
        else
            listrules(model)
        end
    end

    # 1) For each instance in X, we would extract the most adaptive rules for it
    for idx in 1:ndatasetinstances
        dict_predictions = Dict()

        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " *
                "Running instances $(idx)/$(ndatasetinstances) " *
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        )
        currentinstance = slicedataset(X,[idx]; return_view=true)

        global ho
        ho = @hyperopt for i = (
                            length(exec_ntrees)*length(exec_ndims)*length(exec_nclusters)
                        ),
                        ntrees = exec_ntrees,
                        ndims = exec_ndims,
                        nclusters = exec_nclusters

            #println("Iteration number: $(i)")
            #println("Parameters:")
            #println(" - # Trees: $(ntrees)")
            #println(" - # Dimensions: $(ndims)")
            #println(" - # Clusters: $(nclusters)")
            # 2) Only trees that correctly predict the considered instance are considered
            #println("Concording phase in...")
            ctrees, idxsctrees = concordingtrees(
                ftrees, currentinstance, Y[idx]; ntrees=ntrees, kwargs...,
            )

            # 3) Representing rules as a vector
            #println("Rule to vector phase in...")
            rtrees, vrules = begin
                rtrees = begin
                    rs = ruleset[idxsctrees]
                    rs = rs isa Vector{<:Vector{<:Any}} ? reduce(vcat,rs) : rs

                    brs = Vector{Bool}(undef,length(rs))
                    Threads.@threads for (i,r) in collect(enumerate(rs))
                        brs[i] = first(checkantecedent(r,currentinstance))
                    end

                    rs[findall(brs .== true)]
                end

                (rtrees,[rule2vector(t,X) for t in rtrees])
            end

            # https://juliastats.org/MultivariateStats.jl/dev/pca/
            # 4) Vector projection (one of them):
            #     - Principal Component Analysis (PCA)
            #     - Multidimensional Scaling (MDS)
            #println("Projection phase in...")
            prules = projectionvectors(
                hcat(vrules...),X,Y; maxoutdim=ndims, kwargs...,
            )

            # https://juliastats.org/Clustering.jl/dev/kmeans.html
            # 5) Clustering using k-means++ method
            #println("Clustering phase in...")
            krules, cluster_sizes = clusteringvectors(
                prules; rng=rng, nclusters=nclusters, kwargs...,
            )

            # 6) Surrogate model: combination between rules and instance
            # KDTree construction: https://github.com/KristofferC/NearestNeighbors.jl
            #println("Final phase: k-nearest neighbors in...")
            final_rules_idx = begin
                kdrules = KDTree(prules, Minkowski(3.5); leafsize = 40)
                knn_idxs, _ = knn(kdrules, krules, 1, true)

                vcat(knn_idxs)
            end

            # https://github.com/Klest94/Bellatrex/blob/main/code_scripts/LocalMethod_class.py row 241
            #final_trees_idx = [idxsctrees[k] for k in final_rules_idx]
            #localprediction(final_rules_idx,cluster_sizes)
            local_predictions = [
                first(apply(r,currentinstance))
                for r in rtrees[reduce(vcat,final_rules_idx)]
            ]
            final_prediction =
                bestguess(local_predictions,cluster_sizes; suppress_parity_warning=true)

            dict_predictions[[ntrees,ndims,nclusters]] = final_prediction

            final_prediction == Y[idx] ? 0.0 : 1.0
        end

        best_hyperparam = ho.minimizer
        @show best_hyperparam
        best_prediction = dict_predictions[best_hyperparam]
        @show best_prediction
        push!(final_predictions,best_prediction)
    end

    return final_predictions
end

############################################################################################
########################### Supporting functions of BELLATREX ##############################
############################################################################################

function concordingtrees(
    m::Vector{<:SoleModels.DecisionTree},
    X::AbstractLogiset,
    y::Label;
    ntrees::AbstractFloat=0.2,
    kwargs...,
)
    #y_preds = [apply(t,X) for t in m]
    ntrees = Integer(ntrees * length(m))
    y_preds = Vector{Label}(undef,length(m))
    Threads.@threads for (i,t) in collect(enumerate(m))
        y_preds[i] = first(apply(t,X))
    end

    preselectionidxs = begin
        if y isa CLabel
            # Classification problem: find an alternative of euclidean norm
            treesloss = y_preds .== y
            sortperm(treesloss; rev=true)[1:ntrees]
        else
            # Regression problem: euclidean norm
            treesloss = [norm(y-y_preds[k]) for k in eachindex(y_preds)]
            sortperm(treesloss)[1:ntrees]
        end
    end

    #return m[y_preds .== y]
    return (m[preselectionidxs], preselectionidxs)
end

                    #############################################
                    #############################################

# Features extraction: returns all the features in the SyntaxTree
featuresextract(m::MultiFormula) =
    reduce(vcat,[featuresextract(modant) for (_,modant) in modforms(m)])
featuresextract(m::SyntaxTree) = i_variable.(feature.(value.(unique(atoms(m)))))

function rule2vector(m::Rule,X::AbstractLogiset)
    rulevector = zeros(nfeatures(X))
    datanextnodes = rule2vector(antecedent(m),X)

    icovariate = countmap(
                    datanextnodes[:featuresnodes],
                    datanextnodes[:nfallinstances]/ninstances(X)
                )

    return rulevector[1:length(icovariate)] .+= values(icovariate)
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
            SoleLogics.grandchildren(m)
        end
    end

    for conj in conjs
        truthvalues = check(conj,X)

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
    maxoutdim::Union{Integer,Nothing} = nothing,
    kwargs...,
)
    maxoutdim = isnothing(maxoutdim) ?
                    min(size(matrixrules)...) :
                    min(maxoutdim,size(matrixrules)...)

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
