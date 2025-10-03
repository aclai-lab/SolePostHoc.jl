import SoleLogics: alphabet
import Base: show

############################################################################################
############################################################################################
############################################################################################

struct DLIndividual
    dl::DecisionList
    _alphabet::Base.RefValue{<:AbstractConditionalAlphabet}
    _classes::Base.RefValue{<:AbstractVector{<:Label}}
    _alldls::Base.RefValue{<:AbstractVector{<:DecisionList}}

    function DLIndividual(
        dl::DecisionList,
        _alphabet::Base.RefValue{<:AbstractConditionalAlphabet},
        _classes::Base.RefValue{<:AbstractVector{<:Label}},
        _alldls::Base.RefValue{<:AbstractVector{<:DecisionList}},
    )
        new(dl, _alphabet, _classes, _alldls)
    end

    #function DLIndividual()
    #    return new(DecisionList(Rule[],NaN), Ref(), Ref([NaN]))
    #end
end

dl(i::DLIndividual) = i.dl
alphabet(i::DLIndividual) = i._alphabet
classes(i::DLIndividual) = i._classes
alldls(i::DLIndividual) = i._alldls

Base.show(io::IO, indiv::DLIndividual) =
    println(io, "DLIndividual\n$(displaymodel(dl(indiv); header = false))")

meandelay(i::DLIndividual; kwargs...) = meandelay(dl(i); kwargs...)
symbolnumber(i::DLIndividual) = symbolnumber(dl(i))
kappa(i::DLIndividual; kwargs...) = kappa(dl(i); kwargs...)
accuracy(i::DLIndividual; kwargs...) = accuracy(dl(i); kwargs...)
_error(i::DLIndividual; kwargs...) = _error(dl(i); kwargs...)
nrules(i::DLIndividual; kwargs...) = nrules(dl(i); kwargs...)
absnrulescomplexity(i::DLIndividual; kwargs...) = absnrulescomplexity(dl(i); kwargs...)
nrulescomplexity(i::DLIndividual; kwargs...) = nrulescomplexity(dl(i); kwargs...)
nsymbolscomplexity(i::DLIndividual; kwargs...) = nsymbolscomplexity(dl(i); kwargs...)

#=
function finalmetrics(
    i::DLIndividual,
    finalmetricsfuns::Vector{<:Function};
    kwargs...
)
   return map(f->f(i; kwargs...), finalmetricsfuns)
end
=#

function finalmetrics(
    i::DLIndividual;
    X::AbstractLogiset,
    Y::AbstractVector{<:Label},
    memostruct,
)
    return [
        kappa(i; X = X, Y = Y, memostruct = memostruct),
        _error(i; X = X, Y = Y, memostruct = memostruct),
        nrules(i),
        symbolnumber(i),
        meandelay(i; X = X, memostruct = memostruct),
        absnrulescomplexity(i; Y = Y),
    ]
end

############################################################################################
################################ Evolutionary ##############################################
############################################################################################

Evolutionary.default_values(x::DLIndividual) = x

function Evolutionary.EvolutionaryObjective(
    f::TC,
    x::DLIndividual,
    F::Union{Real,AbstractArray{<:Real}} = zero(f(x));
    eval::Symbol = :serial,
) where {TC}
    defval = default_values(x)
    # convert function into the in-place one
    TF = typeof(F)
    fn, TN = if funargnum(f) == 2 && F isa AbstractArray
        ff = (Fv, xv) -> (Fv .= f(xv))
        ff, typeof(ff)
    else
        f, TC
    end
    EvolutionaryObjective{TN,TF,typeof(x),Val{eval}}(fn, F, defval, 0)
end

############################################################################################
############################################################################################
############################################################################################
