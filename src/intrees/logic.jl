using SoleData: AbstractCondition, RangeScalarCondition, ScalarCondition,
    honors_minval, honors_maxval, apply_test_operator, test_operator, threshold

using SoleModels: Label

const Float = Union{Float32,Float64}

# ---------------------------------------------------------------------------- #
#                              check condition                                 #
# ---------------------------------------------------------------------------- #
value(a::Atom) = a.value
i_name(a::Atom{<:RangeScalarCondition}) = a.value.feature.i_name
i_name(a::Atom{<:ScalarCondition}) = a.value.metacond.feature.i_name

checkcondition(r::AbstractCondition, featval::Real) =
    error("Please, provide method for $(typeof(r)), $(typeof(featval))")

@inline function checkcondition(r::RangeScalarCondition, featval::Real)
    honors_minval(r, featval) && honors_maxval(r, featval)
end
@inline function checkcondition(c::ScalarCondition, featval::Real)
    apply_test_operator(test_operator(c), featval, threshold(c))
end
@inline checkcondition(a::Atom, featval::Real) =
    checkcondition(value(a), featval)

@inline function checkcondition(
    rule::ClassificationRule{T},
    x::AbstractVector, 
    featurenames::Vector{<:Union{String,Symbol}}
) where T
    # return all(atoms(rule)) do a
    return all(atoms(antecedent(rule))) do a
        fidx = findfirst(==(T.(i_name(a))), T.(featurenames))
        checkcondition(value(a), x[fidx])
    end
end

@inline function checkcondition(
    rule::ClassificationRule{T},
    X::AbstractArray{S}, 
    args...
) where {T,S}
    return [checkcondition(rule, x, args...) for x in eachrow(X)]
end

@inline function checkcondition(
    set::Vector{ClassificationRule{T}},
    args...
) where T
    return [checkcondition(rule, args...) for rule in set]
end

# ---------------------------------------------------------------------------- #
#                                 measure_rule                                 #
# ---------------------------------------------------------------------------- #
function measure_rule(
    set::ClassificationRule{T},
    X::Matrix{S},
    y::Vector{<:Label},
    featurenames::Vector{F};
    reg_func::Base.Callable=mean
) where {T,S<:Float,F}
    idx_match = checkcondition(set, X, featurenames)
    
    y_match = y[idx_match]
    y_most = mode(sort!(y_match))
    correct = count(==(y_most), y_match)
    confidence = correct / length(y_match)
    err = 1 - confidence

    return y_most, err
end