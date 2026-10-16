# ---------------------------------------------------------------------------- #
#                                 operators                                    #
# ---------------------------------------------------------------------------- #
@inline evalop(::typeof(<)) = 0x01
@inline evalop(::typeof(>)) = 0x02
@inline evalop(::typeof(≥)) = 0x03
@inline evalop(::typeof(≤)) = 0x04

function evalop(op::UInt8)
    op == 0x01 ? '<' :
    op == 0x02 ? '>' :
    op == 0x03 ? '≥' : '≤'
end

function evalop(op::UInt8, x::T, thr::T)::Bool where {T<:AbstractFloat}
    op == 0x01 ? (x < thr) :
    op == 0x02 ? (x > thr) :
    op == 0x03 ? (x ≥ thr) : (x ≤ thr)
end

# dual operator: < ↔ ≥ , > ↔ ≤
@inline dualop(op::UInt8) =
    op == 0x01 ? 0x03 :
    op == 0x03 ? 0x01 :
    op == 0x02 ? 0x04 : 0x02

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

# sort key mirroring SoleData._scalarcondition_sortby:
# (feature, operator, threshold) — operators ordered so that the
# "inclusive" side (≥, ≤) sorts consistently relative to (<, >).
@inline function _lumenatom_sortby(a::LumenAtom)
    (a.feat, a.op, a.thr)
end

Base.isless(a::LumenAtom, b::LumenAtom) =
    isless(_lumenatom_sortby(a), _lumenatom_sortby(b))

LumenAtom{R,T}() where {R<:Unsigned,T<:AbstractFloat} =
    LumenAtom{R,T}(zero(R), T(NaN), 0xff)
LumenAtom() = LumenAtom{UInt32,Float64}()

isempty_atom(a::LumenAtom) = a.op == 0xff && iszero(a.feat)

@inline get_features(atoms::Vector{LumenAtom{R,T}}) where {R,T} =
    [a.feat for a in atoms]
@inline get_thresholds(atoms::Vector{LumenAtom{R,T}}) where {R,T} =
    [a.thr for a in atoms]
@inline get_thresholds(atoms::Vector{LumenAtom{R,T}}, id::R) where {R,T} =
    get_thresholds(filter(a -> a.feat == id, atoms))
@inline get_op(atoms::Vector{LumenAtom{R,T}}) where {R,T} =
    [a.op for a in atoms]
@inline featidxs(atoms::Vector{LumenAtom{R,T}}, feat::R) where {R,T} =
    findall(a -> a.feat == feat, atoms)

@inline dual(a::LumenAtom{R,T}) where {R,T} =
    LumenAtom{R,T}(a.feat, a.thr, dualop(a.op))

const LumenSlice{R,T} = SubArray{
    LumenAtom{R,T},1,
    Vector{LumenAtom{R,T}},Tuple{Base.Slice{Base.OneTo{Int64}}},true
}
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
struct LumenEnsemble{R<:Unsigned,T<:AbstractFloat,S,A<:Base.Callable}
    atoms::Vector{LumenAtom{R,T}}
    nodes::Vector{LumenNode{R}}
    roots::Vector{R}
    classes::Vector{S}
    aggregate::A
end

# the plain label behind an outcome (a categorical value unwraps to its level)
_label(x::CategoricalValue) = CategoricalArrays.unwrap(x)
_label(x) = x

function LumenEnsemble(
    ::LumenConfig{R,T},
    model::AbstractModel
) where {R<:Unsigned,T<:AbstractFloat}
    LumenEnsemble(model, R, T)
end

function LumenEnsemble( 
    model::DecisionEnsemble{U,SM.Branch{S}},
    R::Type{<:Unsigned},
    T::Type{<:AbstractFloat}
) where {U,S<:CategoricalValue}
    atoms = LumenAtom{R,T}[]
    nodes = LumenNode{R}[]
    roots = R[]
    classes = S[]

    function fillensemble(node)::R
        # reserve slot first so children get higher ids
        push!(nodes, LumenNode{R}(zero(R), zero(R), zero(R)))
        push!(atoms, LumenAtom{R,T}())
        id = R(length(nodes))

        if node isa SM.ConstantModel
            o = outcome(node)
            c = levelcode(o)
            if c > length(classes)
                resize!(classes, length(CategoricalArrays.levels(o)))
            end
            classes[c] = o
            atoms[id] = LumenAtom{R,T}()
            nodes[id] = LumenNode{R}(zero(R), zero(R), R(c))
        elseif node isa SM.Branch
            cond = SL.value(antecedent(node))
            left = fillensemble(posconsequent(node))
            right = fillensemble(negconsequent(node))
            atoms[id] = LumenAtom{R,T}(
                R(SD.i_variable(SD.feature(cond))),
                T(SD.threshold(cond)),
                evalop(SD.test_operator(cond)),                
            )
            nodes[id] = LumenNode{R}(left, right, zero(R))
        else
            throw(ArgumentError("unsupported model node $(typeof(node))"))
        end

        return id
    end

    for m in models(model)
        push!(roots, fillensemble(m))
    end

    aggregate = SM.weighted_aggregation(model)
    return LumenEnsemble{R,T,S,typeof(aggregate)}(atoms, nodes, roots, classes, aggregate)
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

"""
    bfs_atoms(e::LumenEnsemble, depth) -> Vector{LumenAtom}

The atoms of each tree in breadth-first order (positive child before
negative), truncated to the first `ceil(n * depth)` per tree, then
concatenated over the trees and made unique. This is classic Lumen's
`depth` semantics (`_extract_atoms_bfs_order` + `_take_first_percentage`),
reproduced on the compiled forest.
"""
function bfs_atoms(e::LumenEnsemble{R,T}, depth::Real) where {R,T}
    out = LumenAtom{R,T}[]
    queue = R[]
    order = LumenAtom{R,T}[]
    for r in e.roots
        empty!(order); empty!(queue)
        push!(queue, r)
        while !isempty(queue)
            id = popfirst!(queue)
            n = e.nodes[id]
            isleaf(n) && continue
            push!(order, e.atoms[id])
            push!(queue, n.left); push!(queue, n.right)
        end
        ntake = min(Int(ceil(length(order) * depth)), length(order))
        append!(out, view(order, 1:ntake))
    end
    return unique!(out)
end

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
    apply!(Vector{R}(undef, size(d, 1)), f, d, nclasses)
end

# in-place form: `preds` must hold at least `size(d, 1)` entries
function apply!(
    preds::Vector{R},
    f::LumenEnsemble{R,T},
    d::SubArray{T},
    nclasses::Integer
) where {R<:Unsigned,T<:AbstractFloat}
    n = size(d, 1)
    counts = Vector{Float64}(undef, nclasses)
    roots = f.roots
    votes = Vector{R}(undef, length(roots)) # leaf class of every tree
    tied = similar(f.classes, length(roots)) # touched only on a tie

    @inbounds for i in 1:n
        fill!(counts, 0.0)

        for k in eachindex(roots)
            r = roots[k]
            node = f.nodes[r]
            while !isleaf(node)
                a = f.atoms[r]
                r = evalop(a.op, d[i, a.feat], a.thr) ?
                    node.left : node.right
                node = f.nodes[r]
            end
            # counts[node.leaf] += weights[k]
            counts[node.leaf] += 1
            votes[k] = node.leaf
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
        # A tie is broken by the model's own aggregation, fed the leaf
        # outcomes in tree order exactly as `SM.apply` does: its parity rule
        # depends on the model (default `bestguess` takes `argmax` of a
        # countmap, the DecisionTreeExt `:alphanumeric` tiebreaker the lowest
        # voted class), so no fixed rule here can match every model.
        if tie
            @inbounds for k in eachindex(votes)
                tied[k] = f.classes[votes[k]]
            end
            best = levelcode(f.aggregate(tied; suppress_parity_warning=true))
        end
        preds[i] = best
    end

    return preds
end