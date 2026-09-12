@inline evalop(::typeof(<)) = 0x01
@inline evalop(::typeof(≤)) = 0x02
@inline evalop(::typeof(>)) = 0x03
@inline evalop(::typeof(≥)) = 0x04

@inline function evalop(op::UInt8, x::T, thr::T)::Bool where {T<:AbstractFloat}
    op == 0x01 ? (x < thr) :
    op == 0x02 ? (x ≤ thr) :
    op == 0x03 ? (x > thr) : (x ≥ thr)
end

# ---------------------------------------------------------------------------- #
#                                 Lumen Node                                   #
# ---------------------------------------------------------------------------- #
struct LumenNode{R<:Unsigned,T<:AbstractFloat}
    feat::R
    thr::T
    op::UInt8
    left::R
    right::R
    leaf::R
end

Base.:(==)(a::LumenNode, b::LumenNode) =
    a.feat == b.feat && a.op == b.op && a.thr == b.thr
Base.isequal(a::LumenNode, b::LumenNode) = a == b
Base.hash(n::LumenNode, h::UInt) =
    hash(n.thr, hash(n.op, hash(n.feat, h)))
@inline isleaf(n::LumenNode{R}) where {R} = n.leaf != zero(R)

# ---------------------------------------------------------------------------- #
#                               Lumen Ensemble                                 #
# ---------------------------------------------------------------------------- #
struct LumenEnsemble{R<:Unsigned,T<:AbstractFloat}
    nodes::Vector{LumenNode{R,T}}
    roots::Vector{R}
end

function LumenEnsemble(
    ::LumenConfig{R,T},
    model::DecisionEnsemble{U,SM.Branch{S}}
) where {R<:Unsigned,T<:AbstractFloat,U,S<:CategoricalValue}
    nodes = LumenNode{R,T}[]
    roots = R[]

    function fillensemble(node)::R
        # reserve slot first so children get higher ids
        push!(nodes, LumenNode{R,T}(
            zero(R), zero(T), zero(R), zero(R), zero(R), zero(R)))
        id = R(length(nodes))

        if node isa SM.ConstantModel
            nodes[id] = LumenNode{R,T}(
                zero(R), zero(T), zero(R), zero(R), zero(R),
                R(levelcode(outcome(node))))
        else
            cond = SL.value(antecedent(node))
            left = fillensemble(posconsequent(node))
            right = fillensemble(negconsequent(node))
            nodes[id] = LumenNode{R,T}(
                R(SD.i_variable(SD.feature(cond))),
                T(SD.threshold(cond)),
                evalop(SD.test_operator(cond)),
                left, right, zero(R)
            )
        end

        return id
    end

    for m in models(model)
        push!(roots, fillensemble(m))
    end

    return LumenEnsemble{R,T}(nodes, roots)
end

Base.length(e::LumenEnsemble) = length(e.nodes)
Base.getindex(e::LumenEnsemble, i::Integer) = e.nodes[i]
Base.getindex(e::LumenEnsemble, idxs::AbstractVector{<:Integer}) = e.nodes[idxs]
Base.iterate(e::LumenEnsemble, s...) = iterate(e.nodes, s...)
Base.eltype(::LumenEnsemble{R,T}) where {R,T} = LumenNode{R,T}

@inline get_atoms(e::LumenEnsemble) = unique!(filter(!isleaf, e.nodes))
function get_thresholds(
    atoms::Vector{LumenNode{R,T}},
    nfeats::R
) where {R<:Unsigned,T<:AbstractFloat}
    thresholds = [T[] for _ in 1:nfeats]
    @inbounds for a in atoms
        push!(thresholds[a.feat], a.thr)
    end
    return thresholds
end
@inline get_thresholds(atoms::Vector{LumenNode{R,T}}) where {R,T} =
    [a.thr for a in atoms]
@inline get_op(atoms::Vector{LumenNode{R,T}}) where {R,T} =
    [a.op for a in atoms]
@inline featidxs(atoms::Vector{LumenNode{R,T}}, feat::R) where {R,T} =
    findall(a -> a.feat == feat, atoms)

# ---------------------------------------------------------------------------- #
#                                   apply                                      #
# ---------------------------------------------------------------------------- #
function apply(
    f::LumenEnsemble{R,T},
    d::SubArray{T},
    nclasses::R
) where {R<:Unsigned,T<:AbstractFloat}
    n = size(d, 1)
    preds = Vector{R}(undef, n)
    counts = Vector{R}(undef, nclasses)

    @inbounds for i in 1:n
        fill!(counts, zero(R))
        for r in f.roots
            node = f.nodes[r]
            while !isleaf(node)
                r = evalop(node.op, d[i, node.feat], node.thr) ?
                    node.left : node.right
                node = f.nodes[r]
            end
            counts[node.leaf] += one(R)
        end
        preds[i] = argmax(counts)
    end

    return preds
end
