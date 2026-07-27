module InTrees

import DecisionTree as DT
using ComplexityMeasures
using DataFrames
using Random

using SoleData

using SoleLogics
using SoleLogics: Atom, LeftmostConjunctiveForm
using SoleLogics: AbstractInterpretationSet, BooleanTruth

using SoleData
using SoleData: scalarlogiset


using SoleModels
using SoleModels: Label, AbstractModel, Rule, antecedent, consequent, info
using SoleModels: rulemetrics, bestguess, evaluaterule
using SoleModels: DecisionList, ConstantModel, isensemble, listrules
using SoleModels: RuleExtractor
using SoleModels: MultiFormula, modforms

abstract type AbstractConfig end
abstract type AbstractRuleSelection end
abstract type AbstractPostProcess end

include("config.jl")
include("pruning.jl")
include("ruleselection.jl")
include("postprocess.jl")

export intrees, InTreesConfig, PruningConfig, CBC
export get_dns

# ---------------------------------------------------------------------------- #
#                                   InTrees                                    #
# ---------------------------------------------------------------------------- #
@inline _starterruleset(model::AbstractModel; kwargs...) = unique!(
    reduce(vcat, [listrules(subm; kwargs...) for subm in SoleModels.models(model)]),
)

"""
    intrees(config::InTreesConfig, model::AbstractModel, X, y::AbstractVector{<:Label})
        -> DecisionList

Return a decision list which approximates the behavior of the input `model` on
the specified supervised dataset. The set of relevant and non-redundant rules
in the decision list is obtained by means of rule extraction, rule pruning,
CBC rule selection, and sequential covering (STEL), using the parameters
encoded in `config`.

# References
- Deng, Houtao. "Interpreting tree ensembles with intrees." International
  Journal of Data Science and Analytics 7.4 (2019): 277-287.

# Pipeline
1. Extract the starting ruleset from `model` (per-tree rules for ensembles,
   `listrules` otherwise).
2. If `get_prune_rules(config)`, prune every rule's antecedent
   ([`_prune_ruleset`](@ref)).
3. Select a compact subset of rules via CBC ([`_select_rules_cbc`](@ref)).
4. Build the final decision list via sequential covering ([`_stel`](@ref)).

Although the method was originally presented for forests it is hereby extended
to work with any symbolic model.

# Arguments
- `config::InTreesConfig`: Algorithm configuration (pruning, CBC, STEL
  parameters).
- `model::AbstractModel`: A single (possibly ensemble) symbolic model.
- `X::AbstractInterpretationSet`: The dataset used to evaluate rules.
- `y::AbstractVector{<:SoleModels.Label}`: Ground-truth labels for `X`.

# Returns
- `DecisionList`: The extracted decision list.

---

    intrees(model::AbstractModel, X, y::AbstractVector{<:Label}; kwargs...)
        -> DecisionList

Convenience method that builds an [`InTreesConfig`](@ref) internally from
`kwargs` and forwards the call to `intrees(config, model, X, y)`. Any keyword
argument accepted by `InTreesConfig` (e.g. `prune_rules`, `pruning_s`,
`pruning_decay_threshold`, `rule_selection_method`, `rule_complexity_metric`,
`min_coverage`, `max_rules`, `n_subfeatures`, `n_trees`, `partial_sampling`,
`max_depth`, `rng`, `dns`, `cbc_threshold`) can be passed here.

`X` can be an `AbstractInterpretationSet` or an `AbstractDataFrame` (in the
latter case it is converted via `SoleData.scalarlogiset`).

# Examples
```julia
# Default configuration
dl = intrees(model, X, y)

# Explicit config object
config = InTreesConfig(prune_rules=true, max_rules=20)
dl = intrees(config, model, X, y)

# Custom CBC / pruning parameters via keyword arguments
dl = intrees(model, X, y; n_trees=100, pruning_decay_threshold=0.1)
```

See also
[`InTreesConfig`](@ref),
[`AbstractModel`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
function intrees(
    config::InTreesConfig{T,S},
    model::AbstractModel,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:Label}
) where {T,S}
    # Extract rules from model
    listrules_kwargs = (use_shortforms=true, normalize=false)
    set = isensemble(model) ?
          _starterruleset(model; listrules_kwargs...) :
          listrules(model; listrules_kwargs...)

    # prune rules if enabled
    get_prune_rules(config) && (set = _prune_ruleset(set, X, y, config))

    # rule selection
    T isa Type{Nothing} || (set = ruleselection(config, set, X, y))

    # construct final decision list via sequential covering
    return S isa Type{Nothing} ?
        DecisionList(set, bestguess(y; suppress_parity_warning=true)) :
        postprocess(config, set, X, y)
end

intrees(
    config::InTreesConfig,
    X::AbstractInterpretationSet,
    y::AbstractVector{<:Label},
    model::AbstractModel
) = intrees(config, model, X, y)

intrees(config::InTreesConfig, m, X::AbstractDataFrame, y) =
    intrees(config, m, scalarlogiset(X; allow_propositional=true), y)

function intrees(
    model::AbstractModel,
    X,
    y::AbstractVector{<:Label};
    kwargs...
)
    config = InTreesConfig(; kwargs...)
    return intrees(config, model, X, y)
end

end