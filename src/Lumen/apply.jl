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

@inline apply(d::DecisionGraph{L,R,T}, args...) where{L,R,T} =
    apply(L, d.graph, args...)
