struct FlatForest{T<:AbstractFloat,R<:Integer}
    feat::Vector{Int32}      # >0: split variable; -1: leaf
    thr::Vector{T}
    left::Vector{Int32}      # taken when antecedent holds
    right::Vector{Int32}
    leaf::Vector{R}          # level code, valid where feat < 0
    roots::Vector{Int32}
    pool::CategoricalPool{String,R}
end

function flatten(m::DecisionEnsemble{U,Branch{S}}) where {U,S<:CategoricalValue}
    T = Float64
    feat, thr = Int32[], T[]
    left, right, leaf = Int32[], Int32[], UInt32[]
    roots = Int32[]

    function push!!(node)::Int32
        push!(feat, 0); push!(thr, zero(T))
        push!(left, 0); push!(right, 0); push!(leaf, 0)
        id = Int32(length(feat))

        if node isa ConstantModel
            feat[id] = -1
            leaf[id] = UInt32(levelcode(outcome(node)))
        else
            cond = SL.value(antecedent(node))
            feat[id] = Int32(SD.i_variable(SD.feature(cond)))
            thr[id]  = T(SD.threshold(cond))
            left[id]  = push!!(posconsequent(node))
            right[id] = push!!(negconsequent(node))
        end
        return id
    end

    for sm in models(m)
        push!(roots, push!!(sm))
    end
    pool = CategoricalArrays.pool(first(info(m).supporting_predictions))::CategoricalPool{String,UInt32}
    return FlatForest{T,UInt32}(feat, thr, left, right, leaf, roots, pool)
end

function apply(
    f::FlatForest{TT,R},
    d::Matrix{T}
) where {T<:AbstractFloat,R,TT<:AbstractFloat}
    n = size(d, 1)
    pool = f.pool
    preds  = Vector{CategoricalValue{String,R}}(undef, n)
    counts = zeros(Int, length(levels(pool)))

    @inbounds for i in 1:n
        fill!(counts, 0)
        for r in f.roots
            node = r
            while f.feat[node] > 0
                node = d[i, f.feat[node]] < f.thr[node] ? f.left[node] : f.right[node]
            end
            counts[f.leaf[node]] += 1
        end
        preds[i] = pool[argmax(counts)]
    end
    return preds
end

function apply(
    f::FlatForest{TT,R},
    d::AbstractMatrix{T}
) where {T<:AbstractFloat,R,TT<:AbstractFloat}
    n = size(d, 1)
    pool = f.pool
    preds  = Vector{CategoricalValue{String,R}}(undef, n)
    counts = zeros(Int, length(levels(pool)))

    @inbounds for i in 1:n
        fill!(counts, 0)
        for r in f.roots
            node = r
            while f.feat[node] > 0
                node = d[i, f.feat[node]] < f.thr[node] ? f.left[node] : f.right[node]
            end
            counts[f.leaf[node]] += 1
        end
        preds[i] = pool[argmax(counts)]
    end
    return preds
end
