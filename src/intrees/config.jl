# ---------------------------------------------------------------------------- #
#                               pruning config                                 #
# ---------------------------------------------------------------------------- #
"""
    PruningConfig

Configuration object for the pruning rule-extraction algorithm.

# Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `prune_rules` | `Bool` | `true` | Whether to prune each rule's antecedent
before selection. |
| `decay_threshold` | `Float64` | `0.05` | Maximum tolerated error-decay before
a conjunct is dropped from a rule. |
| `percentage_degradation` | `Bool` | `false` | Choose between percentage
degradation(`true`) and absolute degradation(`false`). |
| `s` | `Float64` | `1.0e-6` | Denominator floor used when computing percentage
degradation. |

# Validation

The constructor throws `ArgumentError` when:
- `pruning_s` or `pruning_decay_threshold` is negative.

See also: [`intrees`](@ref), [`InTreesConfig`](@ref)
"""
struct PruningConfig
    prune_rules::Bool
    decay_threshold::Float32
    percentage_degradation::Bool
    s::Float32

    function PruningConfig(;
        prune_rules::Bool=true,
        decay_threshold::Float64=0.05,
        percentage_degradation::Bool=false,
        s::Float64=1.0e-6,
    )
        # validate non-negative parameters
        if s < 0.0 || decay_threshold < 0.0
            throw(ArgumentError(
                "s and decay_threshold must be non-negative. Got " *
                "s=$(s), decay_threshold=$(decay_threshold)."
            ))
        end

        new(
            prune_rules,
            decay_threshold,
            percentage_degradation,
            s
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                 cbc config                                   #
# ---------------------------------------------------------------------------- #
"""
    CBC <: AbstractRuleSelection

Configuration object for the CBC rule-selection algorithm.

# Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `rule_complexity_metric` | `Symbol` | `:natoms` | Metric used to estimate rule
complexity (must be a key returned by `SoleModels.rulemetrics`). |
| `threshold` | `Float64` | `0.01` | Minimum normalized feature importance for a
rule to survive CBC selection. |

| `max_rules` | `Int64` | `-1` | Maximum number of rules in the final decision
list (excluding the default rule). `-1` means unlimited. |
| `nsubfeatures` | `Int64` | `2` | Number of candidate features considered at
each split of the CBC random forest. |
| `ntrees` | `Int64` | `50` | Number of trees in the CBC random forest. |
| `partial_sampling` | `Float64` | `0.7` | Fraction of samples used to build
each tree of the CBC random forest, ∈ (0, 1]. |
| `max_depth` | `Int64` | `5` | Maximum depth of each tree in the CBC random
forest. |

# Validation

The constructor throws `ArgumentError` when:
- `partial_sampling` is outside (0.0, 1.0].
- `ntrees` or `max_depth` is not positive.
- `nsubfeatures` is negative.
- `max_rules` is neither `-1` nor a positive integer.

See also: [`intrees`](@ref), [`InTreesConfig`](@ref)
"""

struct CBC <: AbstractRuleSelection
    threshold::Float32
    nsubfeatures::UInt32
    ntrees::UInt32
    partial_sampling::Float32
    max_depth::UInt32

    function CBC(;
        threshold::Float64=0.01,
        nsubfeatures::Int64=2,
        ntrees::Int64=50,
        partial_sampling::Float64=0.7,
        max_depth::Int64=5
    )
        # validate CBC random-forest parameters
        if partial_sampling ≤ 0.0 || partial_sampling > 1.0
            throw(ArgumentError(
                "partial_sampling must be in range (0.0, 1.0]. " *
                "Got: $(partial_sampling)."
            ))
        end

        # validate non-negative parameters
        if threshold < 0.0
            throw(ArgumentError(
                "threshold must be non-negative. Got threshold=$(threshold)."
            ))
        end

        ntrees > 0 || throw(ArgumentError(
            "ntrees must be a positive integer. Got: $(ntrees)."
        ))

        nsubfeatures ≥ 0 || throw(ArgumentError(
            "nsubfeatures must be non-negative. Got: $(nsubfeatures)."
        ))

        max_depth > 0 || throw(ArgumentError(
            "max_depth must be a positive integer. Got: $(max_depth)."
        ))

        new(
            threshold,
            nsubfeatures,
            ntrees,
            partial_sampling,
            max_depth
        )
    end
end

# ---------------------------------------------------------------------------- #
#                               InTrees config                                 #
# ---------------------------------------------------------------------------- #
"""
    InTreesConfig <: AbstractConfig

Configuration object for the InTrees rule-extraction algorithm.

Bundles every tunable parameter into a single, validated, immutable struct.
All fields are set through the keyword constructor, which performs range
and consistency validation before storing anything.

# Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `dns` | `Bool` | `true` | Whether the starting ruleset is built in
"decision-node-set" mode (one rule per branch/leaf) rather than per-path. |
| `pruning` | `PruningConfig` | `default` | Set pruning algorithm parameters. |
| `rule_selection` | `AbstractRuleSelection` | `default` | Set rule_selection
algorithm and parameters. |
| `min_coverage` | `Float64` | `0.01` | Minimum rule coverage required to enter
the STEL sequential-covering step. |
| `rng` | `AbstractRNG` | `Random.TaskLocalRNG()` | RNG used for any randomized
step (CBC forest, tie-breaking in STEL). |

# Validation

The constructor throws `ArgumentError` when:
- `min_coverage`, or `cbc_threshold` is negative.

# Examples

```julia
# Default configuration
cfg = InTreesConfig()

# Custom pruning and coverage parameters
cfg = InTreesConfig(
    pruning = PruningConfig(decay_threshold = 0.1),
    min_coverage = 0.02,
    max_rules = 20,
)

# Tune the underlying CBC random forest
cfg = InTreesConfig(
    n_trees = 100,
    max_depth = 8,
    partial_sampling = 0.8,
)
```

See also: [`intrees`](@ref), [`RuleExtractor`](@ref), [`PruningConfig`](@ref)
"""

struct InTreesConfig{T} <: AbstractConfig
    pruning::PruningConfig
    rule_selection::Union{Nothing,AbstractRuleSelection}
    complexity_metric::Symbol
    max_rules::UInt32
    min_coverage::Float32
    dns::Bool
    rng::AbstractRNG

    function InTreesConfig(;
        pruning::PruningConfig=PruningConfig(),
        rule_selection::T=nothing,
        complexity_metric::Symbol=:natoms,
        max_rules::Int64=0,
        min_coverage::Float64=0.01,
        dns::Bool=true,
        rng::AbstractRNG=Random.TaskLocalRNG()
    ) where {T<:Union{Nothing,AbstractRuleSelection}}
        # validate non-negative parameters
        if min_coverage < 0.0
            throw(ArgumentError(
                "min_coverage must be non-negative. Got " *
                "min_coverage=$(min_coverage)."
            ))
        end

        new{T}(
            pruning,
            rule_selection,
            complexity_metric,
            max_rules,
            min_coverage,
            dns,
            rng
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
"""
    get_complexity_metric(r::InTreesConfig) -> Symbol

Return the rule complexity metric identifier stored in `r`.
"""
@inline get_complexity_metric(r::InTreesConfig) =
    r.complexity_metric

"""
    get_max_rules(r::InTreesConfig) -> UInt32

Return the maximum number of rules in the final decision list stored in `r`.
`0` means unlimited.
"""
@inline get_max_rules(r::InTreesConfig) = r.max_rules

"""
    get_min_coverage(r::InTreesConfig) -> Float64

Return the minimum rule coverage required for STEL, stored in `r`.
"""
@inline get_min_coverage(r::InTreesConfig) = r.min_coverage

"""
    get_dns(r::InTreesConfig) -> Bool

Return `true` if the starting ruleset is built in decision-node-set mode.
"""
@inline get_dns(r::InTreesConfig) = r.dns

"""
    get_rng(r::InTreesConfig) -> AbstractRNG

Return the RNG stored in `r`.
"""
@inline get_rng(r::InTreesConfig) = r.rng

# ---------------------------------------------------------------------------- #
"""
    get_prune_rules(r::InTreesConfig) -> Bool

Return `true` if rule pruning is enabled in `r`.
"""
@inline get_prune_rules(r::InTreesConfig) = r.pruning.prune_rules

"""
    get_decay_threshold(r::InTreesConfig) -> Float32

Return the pruning decay threshold stored in `r`.
"""
@inline get_decay_threshold(r::InTreesConfig) = r.pruning.decay_threshold

"""
    get_percentage_degradation(r::InTreesConfig) -> Bool

Return `true` if percentage_degradation is applied,
    `false` if absolute degradation is applied.
"""
@inline get_percentage_degradation(r::InTreesConfig) =
    r.pruning.percentage_degradation

"""
    get_s(r::InTreesConfig) -> Float32

Return the denominator floor used in the pruning decay metric stored in `r`.
"""
@inline get_s(r::InTreesConfig) = r.pruning.s

# ---------------------------------------------------------------------------- #
"""
    get_cbc_threshold(r::InTreesConfig{T}) where {T<:CBC}
        -> Float32

Return the minimum normalized feature importance required by CBC selection,
stored in `r`.
"""
@inline get_threshold(r::InTreesConfig{T}) where {T<:CBC} =
    r.rule_selection.threshold

"""
    get_ntrees(r::InTreesConfig{T}) where {T<:CBC}
        -> UInt32

Return the number of trees in the CBC random forest stored in `r`.
"""
@inline get_ntrees(r::InTreesConfig{T}) where {T<:CBC} =
    r.rule_selection.ntrees

"""
    get_nsubfeatures(r::InTreesConfig{T}) where {T<:CBC}
        -> UInt32

Return the number of candidate features per split of the CBC random forest,
stored in `r`.
"""
@inline get_nsubfeatures(r::InTreesConfig{T}) where {T<:CBC} =
    r.rule_selection.nsubfeatures

"""
    get_partial_sampling(r::InTreesConfig{T}) where {T<:CBC} -> Float32

Return the per-tree sampling fraction of the CBC random forest stored in `r`.
"""
@inline get_partial_sampling(r::InTreesConfig{T}) where {T<:CBC} =
    r.rule_selection.partial_sampling

"""
    get_max_depth(r::InTreesConfig{T}) where {T<:CBC} -> UInt32

Return the maximum depth of each tree in the CBC random forest stored in `r`.
"""
@inline get_max_depth(r::InTreesConfig{T}) where {T<:CBC} =
    r.rule_selection.max_depth
