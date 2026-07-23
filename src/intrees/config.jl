# ---------------------------------------------------------------------------- #
#                               InTrees config                                 #
# ---------------------------------------------------------------------------- #
"""
    InTreesConfig <: RuleExtractor

Configuration object for the InTrees rule-extraction algorithm.

Bundles every tunable parameter into a single, validated, immutable struct.
All fields are set through the keyword constructor, which performs range
and consistency validation before storing anything.

# Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `dns` | `Bool` | `true` | Whether the starting ruleset is built in "decision-node-set" mode (one rule per branch/leaf) rather than per-path. |
| `prune_rules` | `Bool` | `true` | Whether to prune each rule's antecedent before selection. |
| `pruning_s` | `Float64` | `1.0e-6` | Denominator floor used when computing the pruning decay metric. |
| `pruning_decay_threshold` | `Float64` | `0.05` | Maximum tolerated error-decay before a conjunct is dropped from a rule. |
| `cbc_threshold` | `Float64` | `0.01` | Minimum normalized feature importance for a rule to survive CBC selection. |
| `min_coverage` | `Float64` | `0.01` | Minimum rule coverage required to enter the STEL sequential-covering step. |
| `max_rules` | `Int64` | `-1` | Maximum number of rules in the final decision list (excluding the default rule). `-1` means unlimited. |
| `rule_selection_method` | `Symbol` | `:CBC` | Rule selection method. Currently only `:CBC` is supported. |
| `rule_complexity_metric` | `Symbol` | `:natoms` | Metric used to estimate rule complexity (must be a key returned by `SoleModels.rulemetrics`). |
| `n_subfeatures` | `Int64` | `2` | Number of candidate features considered at each split of the CBC random forest. |
| `n_trees` | `Int64` | `50` | Number of trees in the CBC random forest. |
| `partial_sampling` | `Float64` | `0.7` | Fraction of samples used to build each tree of the CBC random forest, ∈ (0, 1]. |
| `max_depth` | `Int64` | `5` | Maximum depth of each tree in the CBC random forest. |
| `rng` | `AbstractRNG` | `Random.TaskLocalRNG()` | RNG used for any randomized step (CBC forest, tie-breaking in STEL). |

# Validation

The constructor throws `ArgumentError` when:
- `rule_selection_method` is not `:CBC`.
- `pruning_s`, `min_coverage`, or `cbc_threshold` is negative.
- `pruning_decay_threshold` is negative.
- `partial_sampling` is outside (0.0, 1.0].
- `n_trees` or `max_depth` is not positive.
- `n_subfeatures` is negative.
- `max_rules` is neither `-1` nor a positive integer.

# Examples

```julia
# Default configuration
cfg = InTreesConfig()

# Custom pruning and coverage parameters
cfg = InTreesConfig(
    pruning_decay_threshold = 0.1,
    min_coverage            = 0.02,
    max_rules               = 20,
)

# Tune the underlying CBC random forest
cfg = InTreesConfig(
    n_trees          = 100,
    max_depth        = 8,
    partial_sampling = 0.8,
)
```

See also: [`intrees`](@ref), [`RuleExtractor`](@ref)
"""
struct InTreesConfig <: RuleExtractor
    dns::Bool
    prune_rules::Bool
    pruning_s::Float64
    pruning_decay_threshold::Float64
    cbc_threshold::Float64
    min_coverage::Float64
    max_rules::Int64
    rule_selection_method::Symbol
    rule_complexity_metric::Symbol
    n_subfeatures::Int64
    n_trees::Int64
    partial_sampling::Float64
    max_depth::Int64
    rng::AbstractRNG

    function InTreesConfig(;
        dns::Bool=true,
        prune_rules::Bool=true,
        pruning_s::Float64=1.0e-6,
        pruning_decay_threshold::Float64=0.05,
        cbc_threshold::Float64=0.01,
        min_coverage::Float64=0.01,
        max_rules::Int64=-1,
        rule_selection_method::Symbol=:CBC,
        rule_complexity_metric::Symbol=:natoms,
        n_subfeatures::Int64=2,
        n_trees::Int64=50,
        partial_sampling::Float64=0.7,
        max_depth::Int64=5,
        rng::AbstractRNG=Random.TaskLocalRNG()
    )
        # validate rule selection method - only :CBC is currently implemented
        rule_selection_method === :CBC || throw(ArgumentError(
            "rule_selection_method must be :CBC. Got: $(rule_selection_method)."
        ))

        # validate non-negative parameters
        if pruning_s < 0.0 || min_coverage < 0.0 || cbc_threshold < 0.0 ||
           pruning_decay_threshold < 0.0
            throw(ArgumentError(
                "pruning_s, min_coverage, cbc_threshold and " *
                "pruning_decay_threshold must be non-negative. Got " *
                "pruning_s=$(pruning_s), min_coverage=$(min_coverage), " *
                "cbc_threshold=$(cbc_threshold), " *
                "pruning_decay_threshold=$(pruning_decay_threshold)."
            ))
        end

        # validate CBC random-forest parameters
        if partial_sampling ≤ 0.0 || partial_sampling > 1.0
            throw(ArgumentError(
                "partial_sampling must be in range (0.0, 1.0]. " *
                "Got: $(partial_sampling)."
            ))
        end

        n_trees > 0 || throw(ArgumentError(
            "n_trees must be a positive integer. Got: $(n_trees)."
        ))

        max_depth > 0 || throw(ArgumentError(
            "max_depth must be a positive integer. Got: $(max_depth)."
        ))

        n_subfeatures ≥ 0 || throw(ArgumentError(
            "n_subfeatures must be non-negative. Got: $(n_subfeatures)."
        ))

        # validate max_rules - either unlimited (-1) or a positive integer
        (max_rules == -1 || max_rules > 0) || throw(ArgumentError(
            "max_rules must be -1 (unlimited) or a positive integer. " *
            "Got: $(max_rules)."
        ))

        new(
            dns,
            prune_rules,
            pruning_s,
            pruning_decay_threshold,
            cbc_threshold,
            min_coverage,
            max_rules,
            rule_selection_method,
            rule_complexity_metric,
            n_subfeatures,
            n_trees,
            partial_sampling,
            max_depth,
            rng
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
"""
    get_dns(r::InTreesConfig) -> Bool

Return `true` if the starting ruleset is built in decision-node-set mode.
"""
@inline get_dns(r::InTreesConfig) = r.dns

"""
    get_prune_rules(r::InTreesConfig) -> Bool

Return `true` if rule pruning is enabled in `r`.
"""
@inline get_prune_rules(r::InTreesConfig) = r.prune_rules

"""
    get_pruning_s(r::InTreesConfig) -> Float64

Return the denominator floor used in the pruning decay metric stored in `r`.
"""
@inline get_pruning_s(r::InTreesConfig) = r.pruning_s

"""
    get_pruning_decay_threshold(r::InTreesConfig) -> Float64

Return the pruning decay threshold stored in `r`.
"""
@inline get_pruning_decay_threshold(r::InTreesConfig) = r.pruning_decay_threshold

"""
    get_cbc_threshold(r::InTreesConfig) -> Float64

Return the minimum normalized feature importance required by CBC selection,
stored in `r`.
"""
@inline get_cbc_threshold(r::InTreesConfig) = r.cbc_threshold

"""
    get_min_coverage(r::InTreesConfig) -> Float64

Return the minimum rule coverage required for STEL, stored in `r`.
"""
@inline get_min_coverage(r::InTreesConfig) = r.min_coverage

"""
    get_max_rules(r::InTreesConfig) -> Int64

Return the maximum number of rules in the final decision list stored in `r`.
`-1` means unlimited.
"""
@inline get_max_rules(r::InTreesConfig) = r.max_rules

"""
    get_rule_selection_method(r::InTreesConfig) -> Symbol

Return the rule selection method stored in `r`.
"""
@inline get_rule_selection_method(r::InTreesConfig) = r.rule_selection_method

"""
    get_rule_complexity_metric(r::InTreesConfig) -> Symbol

Return the rule complexity metric identifier stored in `r`.
"""
@inline get_rule_complexity_metric(r::InTreesConfig) = r.rule_complexity_metric

"""
    get_n_subfeatures(r::InTreesConfig) -> Int64

Return the number of candidate features per split of the CBC random forest,
stored in `r`.
"""
@inline get_n_subfeatures(r::InTreesConfig) = r.n_subfeatures

"""
    get_n_trees(r::InTreesConfig) -> Int64

Return the number of trees in the CBC random forest stored in `r`.
"""
@inline get_n_trees(r::InTreesConfig) = r.n_trees

"""
    get_partial_sampling(r::InTreesConfig) -> Float64

Return the per-tree sampling fraction of the CBC random forest stored in `r`.
"""
@inline get_partial_sampling(r::InTreesConfig) = r.partial_sampling

"""
    get_max_depth(r::InTreesConfig) -> Int64

Return the maximum depth of each tree in the CBC random forest stored in `r`.
"""
@inline get_max_depth(r::InTreesConfig) = r.max_depth

"""
    get_rng(r::InTreesConfig) -> AbstractRNG

Return the RNG stored in `r`.
"""
@inline get_rng(r::InTreesConfig) = r.rng