const OpCondition{R<:Unsigned} = Dict(
    R(1) => (<),
    R(2) => (>),
    R(3) => (≤),
    R(4) => (≥)
)

# ---------------------------------------------------------------------------- #
#                               Lumen Ensemble                                 #
# ---------------------------------------------------------------------------- #
struct LumenEnsemble{R<:Unsigned,T<:AbstractFloat}
    feat::Vector{R}
    thr::Vector{T}
    op::Vector{R}
    left::Vector{R}
    right::Vector{R}
    leaf::Vector{R}
    roots::Vector{R}

    function LumenEnsemble(
        ::LumenConfig{R,T},
        model::DecisionEnsemble{U,SM.Branch{S}}
    ) where {R<:Unsigned,T<:AbstractFloat,U,S<:CategoricalValue}
        roots = R[]
        feat, thr, op = R[], T[], Symbol[]
        left, right, leaf = R[], R[], R[]

        function fillensemble(node)::Int32
            push!(feat, zero(R))
            push!(thr, zero(T))
            push!(op, zero(R))
            push!(left, zero(R))
            push!(right, zero(R))
            push!(leaf, zero(R))
            id = R(length(feat))

            if node isa SM.ConstantModel
                feat[id] = zero(R)
                leaf[id] = R(levelcode(outcome(node)))
            else
                cond = SL.value(antecedent(node))
                feat[id] = R(SD.i_variable(SD.feature(cond)))
                thr[id] = T(SD.threshold(cond))
                # op[id] = Symbol(SD.test_operator(cond))
                left[id] = fillensemble(posconsequent(node))
                right[id] = fillensemble(negconsequent(node))
            end

            return id
        end

        for m in models(model)
            push!(roots, fillensemble(m))
        end
        
        return new{R,T}(feat, thr, left, right, leaf, roots)
    end
end

Base.length(e::LumenEnsemble) = length(e.feat)
Base.getindex(e::LumenEnsemble, i::Integer) =
    (feat=e.feat[i], thr=e.thr[i], op=e.op[i])
Base.getindex(e::LumenEnsemble, idxs::AbstractVector{<:Integer}) =
    [e[i] for i in idxs]

@inline atoms(e::LumenEnsemble{R}) where {R<:Unsigned} =
    e[findall(==(zero(R)), e.leaf)]
