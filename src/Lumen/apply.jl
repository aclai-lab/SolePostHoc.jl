# ---------------------------------------------------------------------------- #
#                            apply on matrix data                              #
# ---------------------------------------------------------------------------- #
function apply(
    m::DecisionEnsemble{U,Branch{S}},
    d::Matrix{T};
    suppress_parity_warning::Bool=true
)::Vector{S} where {S,U,T<:AbstractFloat}
    preds = hcat([apply(subm, d) for subm in models(m)]...)

    preds = [
        SM.aggregation(m)(preds[i,:]; suppress_parity_warning) for i in axes(preds,1)
    ]

    return preds
end

function apply(
    m::Branch{S},
    d::Union{Matrix{T},SubArray{T}}
)::Vector{S} where {S,T<:AbstractFloat}
    checkmask = checkantecedent(m, d)
    negmask = .!checkmask
    preds = Vector{outputtype(m)}(undef, length(checkmask))
    preds[checkmask] .= apply(posconsequent(m), @view d[checkmask, :])
    preds[negmask] .= apply(negconsequent(m), @view d[negmask, :])
    return preds
end

function apply(
    m::ConstantModel{S},
    ::Union{Matrix{T},SubArray{T}}
)::S where {S,T<:AbstractFloat}
    outcome(m)
end

function checkantecedent(
    m::Branch{S},
    d::Union{Matrix{T},SubArray{T}}
)::BitVector where {S,T<:AbstractFloat}
    check(antecedent(m), d)
end

function check(
    atom::Atom{<:SM.ScalarCondition}, d::Union{Matrix{T},SubArray{T}}
)::BitVector where {T<:AbstractFloat}
    cond = SL.value(atom)
    return checkcondition(cond, d)
end

function checkcondition(
    cond::ScalarCondition,
    d::Union{Matrix{T},SubArray{T}}
)::BitVector where {T<:AbstractFloat}
    cond_threshold = SD.threshold(cond)
    cond_operator = SD.test_operator(cond)
    cond_feature = SD.feature(cond)

    col = SD.i_variable(cond_feature)
    
    return cond_operator.(@view(d[:, col]), cond_threshold)
end
