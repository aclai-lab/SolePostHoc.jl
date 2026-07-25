using SoleData: AbstractCondition, RangeScalarCondition, ScalarCondition,
    LeftmostConjunctiveForm,
    honors_minval, honors_maxval, apply_test_operator, test_operator, threshold,
    minval, maxval, feature, threshold,
    _isgreater_test_operator, _isless_test_operator

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
    return all(atoms(antecedent(rule))) do a
        fidx = findfirst(==(T.(i_name(a))), T.(featurenames))
        checkcondition(value(a), x[fidx])
    end
end

@inline function checkcondition(
    rule::Vector{ScalarCondition},
    x::AbstractVector,
    featurenames::Vector{<:Union{String,Symbol}}
)
@show "PASO"
    return all(atoms(rule)) do a
        name = i_name(a)
        fidx = findfirst(==(name), typeof(name).(featurenames))
        checkcondition(a, x[fidx])
    end
end

@inline function checkcondition(
    rule::Union{ClassificationRule{T},Vector{ScalarCondition}},
    X::AbstractArray{S}, 
    args...
) where {T,S}
@show "PASO"
    return [checkcondition(rule, x, args...) for x in eachrow(X)]
end

@inline function checkcondition(
    set::Vector{ClassificationRule{T}},
    args...
) where T
    return [checkcondition(rule, args...) for rule in set]
end

# ---------------------------------------------------------------------------- #
#                ClassificationRule to LeftmostConjunctiveForm                 #
# ---------------------------------------------------------------------------- #
function LeftmostConjunctiveForm(set::ClassificationRule{T}) where T
    cond = antecedent(set)
    conds = Atom[]

    for c in cond
        minv = minval(c.value)
        !isnothing(minv) && !isinf(minv) && push!(
            conds,
            Atom(ScalarCondition(
                feature(c.value),
                _isgreater_test_operator(c.value),
                minv
            ))
        )
        maxv = maxval(c.value)
        !isnothing(maxv) && !isinf(maxv) && push!(
            conds,
            Atom(ScalarCondition(
                feature(c.value),
                _isless_test_operator(c.value),
                maxv
            ))
        )
    end

    return LeftmostConjunctiveForm(conds)
end

# ---------------------------------------------------------------------------- #
#                                 measure_rule                                 #
# ---------------------------------------------------------------------------- #
function measure_rule(
    set::Union{ClassificationRule{T},Vector{ScalarCondition}},
    X::Matrix{S},
    y::Vector{<:Label},
    featurenames::Vector{F};
    pred::Union{Nothing,Label}=nothing,
    reg_func::Base.Callable=mean
) where {T,S<:Float,F}
    idx_match = checkcondition(set, X, featurenames)
    
    y_match = y[idx_match]
    isnothing(pred) && (pred = mode(sort(y_match)))
    confidence = count(==(pred), y_match) / length(y_match)
    err = 1 - confidence

    return pred, err
end
