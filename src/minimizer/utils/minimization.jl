function minimizza_dnf(_minimization_scheme::Val, formula::TwoLevelDNFFormula, kwargs...)
    error("Unknown minimization scheme: $(_minimization_scheme)!")
end

function minimizza_dnf(
    ::Val{:mitespresso},
    formula::TwoLevelDNFFormula;
    silent = true,
    mitespresso_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)
    formula = SoleData.espresso_minimize(formula, silent; mitespresso_kwargs...)
    silent || (println(); @show formula)
    formula = convert(TwoLevelDNFFormula, formula)
    silent || (println(); @show formula)
    return formula
end
