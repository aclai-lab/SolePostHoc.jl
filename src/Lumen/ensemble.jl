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
#                                 Lumen Atom                                   #
# ---------------------------------------------------------------------------- #
struct LumenAtom{R<:Unsigned,T<:AbstractFloat}
    feat::R
    thr::T
    op::UInt8
end

Base.:(==)(a::LumenAtom, b::LumenAtom) =
    a.feat == b.feat && a.op == b.op && a.thr == b.thr
Base.isequal(a::LumenAtom, b::LumenAtom) = a == b
Base.hash(n::LumenAtom, h::UInt) =
    hash(n.thr, hash(n.op, hash(n.feat, h)))

LumenAtom{R,T}() where {R<:Unsigned,T<:AbstractFloat} =
    LumenAtom{R,T}(zero(R), T(NaN), 0xff)
LumenAtom() = LumenAtom{UInt32,Float64}()

isempty_atom(a::LumenAtom) = a.op == 0xff && iszero(a.feat)

@inline get_thresholds(atoms::Vector{LumenAtom{R,T}}) where {R,T} =
    [a.thr for a in atoms]
@inline get_thresholds(atoms::Vector{LumenAtom{R,T}}, id::R) where {R,T} =
    get_thresholds(filter(a -> a.feat == id, atoms))
@inline get_op(atoms::Vector{LumenAtom{R,T}}) where {R,T} =
    [a.op for a in atoms]
@inline featidxs(atoms::Vector{LumenAtom{R,T}}, feat::R) where {R,T} =
    findall(a -> a.feat == feat, atoms)

const LumenCube{R,T} = 
    SubArray{LumenAtom{R,T},1,Vector{LumenAtom{R,T}},Tuple{UnitRange{Int}},true}

# ---------------------------------------------------------------------------- #
#                                 Lumen Node                                   #
# ---------------------------------------------------------------------------- #
struct LumenNode{R<:Unsigned}
    left::R
    right::R
    leaf::R
end

Base.:(==)(a::LumenNode, b::LumenNode) =
    a.left == b.left && a.right == b.right && a.leaf == b.leaf
Base.isequal(a::LumenNode, b::LumenNode) = a == b
Base.hash(n::LumenNode, h::UInt) =
    hash(n.leaf, hash(n.right, hash(n.left, h)))
@inline isleaf(n::LumenNode{R}) where {R} = n.leaf != zero(R)

# ---------------------------------------------------------------------------- #
#                               Lumen Ensemble                                 #
# ---------------------------------------------------------------------------- #
struct LumenEnsemble{R<:Unsigned,T<:AbstractFloat}
    atoms::Vector{LumenAtom{R,T}}
    # nodes::Vector{LumenNode{R,T}}
    nodes::Vector{LumenNode{R}}
    roots::Vector{R}
end

function LumenEnsemble(
    ::LumenShannonConfig{R,T},
    model::AbstractModel
) where {R<:Unsigned,T<:AbstractFloat}
    LumenEnsemble(model, R, T)
end

function LumenEnsemble( 
    model::DecisionEnsemble{U,SM.Branch{S}},
    R::Type{<:Unsigned}=UInt32,
    T::Type{<:AbstractFloat}=Float32
) where {U,S<:CategoricalValue}
    atoms = LumenAtom{R,T}[]
    nodes = LumenNode{R}[]
    roots = R[]

    function fillensemble(node)::R
        # reserve slot first so children get higher ids
        push!(nodes, LumenNode{R}(zero(R), zero(R), zero(R)))
        push!(atoms, LumenAtom{R,T}())
        id = R(length(nodes))

        if node isa SM.ConstantModel
            atoms[id] = LumenAtom{R,T}()
            nodes[id] = LumenNode{R}(
                zero(R), zero(R),
                R(levelcode(outcome(node)))
            )
        else
            cond = SL.value(antecedent(node))
            left = fillensemble(posconsequent(node))
            right = fillensemble(negconsequent(node))
            atoms[id] = LumenAtom{R,T}(
                R(SD.i_variable(SD.feature(cond))),
                T(SD.threshold(cond)),
                evalop(SD.test_operator(cond)),                
            )
            nodes[id] = LumenNode{R}(left, right, zero(R))
        end

        return id
    end

    for m in models(model)
        push!(roots, fillensemble(m))
    end

    return LumenEnsemble{R,T}(atoms, nodes, roots)
end

Base.length(e::LumenEnsemble) = length(e.nodes)
Base.getindex(e::LumenEnsemble, i::Integer) = e.nodes[i]
Base.getindex(e::LumenEnsemble, idxs::AbstractVector{<:Integer}) = e.nodes[idxs]
Base.iterate(e::LumenEnsemble, s...) = iterate(e.nodes, s...)
Base.eltype(::LumenEnsemble{R,T}) where {R,T} = LumenAtom{R,T}

@inline get_thresholds(e::LumenEnsemble) = get_thresholds(e.atoms)
@inline get_thresholds(e::LumenEnsemble, id::R) where R =
    get_thresholds(e.atoms, id)
@inline get_atoms(e::LumenEnsemble) = unique!(filter(!isempty_atom, e.atoms))

function Base.show(io::IO, e::LumenEnsemble{R,T}) where {R,T}
    nleaves = count(isleaf, e.nodes)
    print(io, "LumenEnsemble{", R, ",", T, "}(",
        length(e.roots), " trees, ",
        length(e.nodes), " nodes, ",
        nleaves, " leaves)")
end

function Base.show(io::IO, ::MIME"text/plain", e::LumenEnsemble{R,T}) where {R,T}
    show(io, e)
    isempty(e.nodes) && return
    nsplits = length(e.nodes) - count(isleaf, e.nodes)
    feats = unique(a.feat for (a, n) in zip(e.atoms, e.nodes) if !isleaf(n))
    println(io)
    println(io, "  splits:   ", nsplits)
    println(io, "  features: ", length(feats))
    print(io,   "  classes:  ", length(unique(n.leaf for n in e.nodes if isleaf(n))))
end

# ---------------------------------------------------------------------------- #
#                                   apply                                      #
# ---------------------------------------------------------------------------- #
function apply(
    f::LumenEnsemble{R,T},
    d::SubArray{T},
    nclasses::Integer
) where {R<:Unsigned,T<:AbstractFloat}
    n = size(d, 1)
    preds = Vector{R}(undef, n)
    counts = Vector{R}(undef, nclasses)

    @inbounds for i in 1:n
        fill!(counts, zero(R))

        for r in f.roots
            node = f.nodes[r]
            while !isleaf(node)
                a = f.atoms[r]
                r = evalop(a.op, d[i, a.feat], a.thr) ?
                    node.left : node.right
                node = f.nodes[r]
            end
            counts[node.leaf] += one(R)
        end

        best = 1
        tie = false
        @inbounds for j in 2:nclasses
            if counts[j] > counts[best]
                best, tie = j, false
            elseif counts[j] == counts[best]
                tie = true
            end
        end
        preds[i] = best
    end

    return preds
end
