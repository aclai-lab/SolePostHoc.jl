# map categorical labels onto unsigned integer of type `R`
# sort classlabels
function assign(
    ::Type{R},
    y::Vector{S}
) where {R<:Unsigned,S<:CategoricalValue}
    classlabels = sort!(unique(y))
    dict = Dict{S,R}(v => i for (i, v) in enumerate(classlabels))
    return string.(classlabels), [dict[t] for t in y]
end