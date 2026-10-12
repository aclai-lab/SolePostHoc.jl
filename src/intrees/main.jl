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

export intrees, InTreesConfig, PruningConfig, CBC, STEL
export get_dns

# ---------------------------------------------------------------------------- #
#                                   InTrees                                    #
# ---------------------------------------------------------------------------- #

# Extract the starting ruleset from an ensemble, with per-tree debug logging.
function _starterruleset(model::AbstractModel; kwargs...)
    submodels = SoleModels.models(model)
    n = length(submodels)
    @info "[InTrees] Extracting rules from ensemble" n_trees=n

    per_tree = Vector{Any}(undef, n)
    t_extract = time()
    total = 0
    for (i, subm) in enumerate(submodels)
        t = @elapsed rs = listrules(subm; kwargs...)
        per_tree[i] = rs
        total += length(rs)
        @info "[InTrees]   tree $i/$n" rules_from_tree=length(rs) total_rules_so_far=total secs=round(t, digits=3)
    end
    @info "[InTrees] Extraction done" total_rules=total secs=round(time() - t_extract, digits=2)

    t_cat = @elapsed allrules = reduce(vcat, per_tree)
    @info "[InTrees] vcat done" n_rules=length(allrules) secs=round(t_cat, digits=2)

    @info "[InTrees] Running unique! on rules (can be slow: needs rule equality/hashing)..."
    t_unique = @elapsed unique!(allrules)
    @info "[InTrees] unique! done" unique_rules=length(allrules) secs=round(t_unique, digits=2)

    return allrules
end

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
    t_total = time()
    @info "[InTrees] START" config_type=typeof(config) n_instances=SoleLogics.ninstances(X) n_labels=length(y) is_ensemble=isensemble(model)

    # 1. Rule extraction
    listrules_kwargs = (use_shortforms=true, normalize=false)
    t1 = @elapsed begin
        set = isensemble(model) ?
              _starterruleset(model; listrules_kwargs...) :
              listrules(model; listrules_kwargs...)
    end
    @info "[InTrees] STEP 1/4 rule extraction finished" n_rules=length(set) secs=round(t1, digits=2)

    # 2. Pruning
    if get_prune_rules(config)
        @info "[InTrees] STEP 2/4 pruning START" n_rules=length(set)
        t2 = @elapsed (set = _prune_ruleset(set, X, y, config))
        @info "[InTrees] STEP 2/4 pruning DONE" n_rules=length(set) secs=round(t2, digits=2)
    else
        @info "[InTrees] STEP 2/4 pruning SKIPPED"
    end

    # 3. Rule selection
    if T isa Type{Nothing}
        @info "[InTrees] STEP 3/4 rule selection SKIPPED"
    else
        @info "[InTrees] STEP 3/4 rule selection START" n_rules=length(set) method=T
        t3 = @elapsed (set = ruleselection(config, set, X, y))
        @info "[InTrees] STEP 3/4 rule selection DONE" n_rules=length(set) secs=round(t3, digits=2)
    end

    # 4. Decision list construction (sequential covering)
    @info "[InTrees] STEP 4/4 decision list START" n_rules=length(set) method=S
    t4 = @elapsed begin
        dl = S isa Type{Nothing} ?
             DecisionList(set, bestguess(y; suppress_parity_warning=true)) :
             postprocess(config, set, X, y)
    end
    @info "[InTrees] STEP 4/4 decision list DONE" secs=round(t4, digits=2)

    @info "[InTrees] FINISHED" total_secs=round(time() - t_total, digits=2)
    return dl
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