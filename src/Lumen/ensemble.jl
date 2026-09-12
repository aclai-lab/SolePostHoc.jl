const OpCondition = Dict{UInt8,Function}(
    0x01 => (<),
    0x02 => (>),
    0x03 => (≤),
    0x04 => (≥)
)

const OpCode = Dict{Function,UInt8}(v => k for (k, v) in OpCondition)

@inline opcondition(code::Unsigned) = OpCondition[UInt8(code)]
@inline opcode(f::Function) = OpCode[f]

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
                R(levelcode(outcome(node)))
            )
        else
            cond = SL.value(antecedent(node))
            left = fillensemble(posconsequent(node))
            right = fillensemble(negconsequent(node))
            nodes[id] = LumenNode{R,T}(
                R(SD.i_variable(SD.feature(cond))),
                T(SD.threshold(cond)),
                opcode(SD.test_operator(cond)),
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

@inline atoms(e::LumenEnsemble) = filter(!isleaf, e.nodes)

# ---------------------------------------------------------------------------- #
#                                 Base.show                                    #
# ---------------------------------------------------------------------------- #
@inline Base.show(io::IO, n::LumenNode{R,T}) where {R,T} =
    isleaf(n) ?
        print(io, "LumenNode{$R,$T}(leaf=", n.leaf, ")") :
        print(io, "LumenNode{$R,$T}(V", n.feat, " ",
            opcondition(n.op), " ", n.thr, ")")

function Base.show(
    io::IO,
    ::MIME"text/plain",
    n::LumenNode{R,T}
) where {R,T}
    println(io, "LumenNode{$R,$T}")
    if isleaf(n)
        println(io, "  leaf     :", n.leaf)
    else
        println(io, "  feature  :", n.feat)
        println(io, "  operator :", opcondition(n.op))
        println(io, "  threshold:", n.thr)
        println(io, "  left     :", n.left)
        println(io, "  right    :", n.right)
    end
end

@inline Base.show(io::IO, e::LumenEnsemble{R,T}) where {R,T} =
    print(io, "LumenEnsemble{$R,$T}(",
        length(e.roots), " trees, ", length(e.nodes), " nodes)")

function Base.show(
    io::IO,
    ::MIME"text/plain",
    e::LumenEnsemble{R,T}
) where {R,T}
    println(io, "LumenEnsemble{$R,$T}")
    println(io, "  trees :", length(e.roots))
    println(io, "  nodes :", length(e.nodes))
    println(io, "  leaves:", count(isleaf, e.nodes))
end