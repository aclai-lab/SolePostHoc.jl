"""
    orca

Entry point for the bitvector-driven Random Forest compression algorithm.

Orchestrates the pipeline stages selected via `parts_to_run`:
  1. Run the Genetic Algorithm (`parts.core`) optimising only the requested stages
  2. Apply the winning mask through the selected compression stages:
       `first_part`  – tree selection   (presence bits)
       `second_part` – tree pruning     (depth bits)
       `third_part`  – alphabet/feature selection (alphabet bits)

# Dependencies
- `types.jl`  – `GeneticForest` struct and accessors
- `parts.jl`  – `first_part`, `second_part`, `third_part`, `core`
"""
module Orca

using SoleModels
using SoleLogics
using StatsBase: mode
using Evolutionary
using DataFrames

export compression

include("types.jl")
include("utils.jl")
include("parts.jl")


"""
    compression(original_f; ..., parts_to_run=[1,2,3]) → (final_forest, best_bitvector, best_cost)

Run the compression pipeline on the inserted forest.

## Keyword arguments
| Name              | Default   | Description                                                 |
|-------------------|-----------|-------------------------------------------------------------|
| `original_f`      | required  | The random forest that will be optimized                    |
| `mode`            | required  | The stages of the compression pipeline that will be done    |
| `f_val`           | required  | The features matrix that will be used to optimize the forest|
| `l_val`           | required  | The labels vector that will be used to optimize the forest  |
| `population_size` | 50        | GA population size                                          |
| `n_generations`   | 100       | GA iterations                                               |
| `penalty_weight`  | 0.3       | Compression-vs-accuracy trade-off for GA fitness            |

## `mode` examples
- `[:size]`             – tree selection only (no pruning, no alhabet modification)
- `[:size_depth]`       – selection + pruning (no alphabet modification)
- `[:depth_alphabet]`   – pruning + alphabet  (all trees kept, no selection)
- `[:full_dimensional]` – full pipeline 
"""

const PARTS_MAP = Dict(
    :size             => [1],
    :depth            => [2],
    :alphabet         => [3],
    :size_depth       => [1, 2],
    :size_alphabet    => [1, 3],
    :depth_alphabet   => [2, 3],
    :full_dimensional => [1, 2, 3]
)

function compression(
        original_f::DecisionEnsemble,
        mode::Symbol,
        f_val::Matrix{Float64}, 
        l_val::Vector{String};
        population_size::Int      = 50,
        n_generations::Int        = 100,
        penalty_weight::Float64   = 0.3)


    # ── Validate parts_to_run ─────────────────────────────────────────────────
    if !haskey(PARTS_MAP, mode)
        error("Error: please use one of the following symbol: ", join(keys(PARTS_MAP), ", "))
    end  

    # ─────────────────────────────────────────────────────────────────────────
    # EVOLUTIONARY OPTIMISATION
    #
    # The GA searches for the BitVector of length 3*n_trees that:
    #   • minimises classification error on the validation set
    #   • minimises the number of active bits (1s) → smaller forest
    #
    # fitness = error_rate + penalty_weight * (n_ones / n_bits)
    #
    # penalty_weight (passed as keyword argument):
    #   higher → prefer compression over accuracy
    #   lower  → prefer accuracy over compression
    # ─────────────────────────────────────────────────────────────────────────
    parts_to_run = PARTS_MAP[mode]


    n_trees = length(original_f.models)

    best_bitvector, best_cost = core(
        n_trees, original_f, f_val, l_val;
        population_size = population_size,
        n_generations   = n_generations,
        penalty_weight  = penalty_weight,
        active_parts    = parts_to_run
    )


    # ─────────────────────────────────────────────────────────────────────────
    # APPLY THE BEST BITVECTOR — only for the requested stages
    # ─────────────────────────────────────────────────────────────────────────

    # ── Stage 1 — tree selection ──────────────────────────────────────────────
    if 1 in parts_to_run
        selected_forest, presence_mask = first_part(n_trees, original_f, best_bitvector, parts_to_run)
    else
        selected_forest = original_f
        presence_mask   = trues(n_trees)
    end


    # ── Stage 2 — tree pruning ────────────────────────────────────────────────
    if 2 in parts_to_run
        pruned_forest, _, _ =
            second_part(n_trees, presence_mask, original_f, best_bitvector, parts_to_run)
    else
        # Stage 2 skipped: build a forest from the kept trees (no pruning)
        kept_indices  = findall(presence_mask)
        pruned_forest = DecisionEnsemble(original_f.models[kept_indices])
    end


    # ── Stage 3 — alphabet modification ──────────────────────────────────────
    if 3 in parts_to_run
        final_forest, _, _ =
            third_part(n_trees, presence_mask, pruned_forest, best_bitvector, parts_to_run)
    else
        final_forest   = pruned_forest
    end


    return final_forest

end  # compression

end  # module Orca