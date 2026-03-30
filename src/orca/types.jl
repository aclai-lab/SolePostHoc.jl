"""
    types

Defines the core data structure used throughout the compression pipeline.

`GeneticForest` pairs a trained `DecisionEnsemble` with a `BitVector` that
encodes how the forest should be compressed. All other logic (training, pruning,
pipeline execution) lives in `utils.jl`, `parts.jl`, and `pipeline.jl`.

# Exports
- `GeneticForest`    – the main struct
- `vector`           – accessor for the bitvector
- `original_forest`  – accessor for the forest
- `with_vector`      – functional updater for the bitvector
- `with_forest`      – functional updater for the forest
"""
module types

using SoleModels

export GeneticForest, vector, original_forest, with_vector, with_forest


# ─────────────────────────────────────────────────────────────────────────────
# STRUCT
# ─────────────────────────────────────────────────────────────────────────────

"""
    GeneticForest

Immutable struct that pairs a trained `DecisionEnsemble` with a `BitVector`
encoding which trees to keep and how to prune them.

The bitvector has length `3 * n_trees` and is structured as:

  [ presence (N) | depth/pruning (N) | alphabet (N) ]

# Fields
- `vector::BitVector`                 : Bitvector of length `3 * n_trees(forest)`.
- `original_forest::DecisionEnsemble` : Full trained forest, used as the source
                                        of trees throughout the pipeline.
"""
struct GeneticForest
    vector::BitVector
    original_forest::DecisionEnsemble
end

"""
    GeneticForest(vector, forest) → GeneticForest

Outer constructor with a compatibility check. Raises an error if
`length(vector) ≠ n_trees(forest) * 3`, ensuring exactly three bits per tree.
"""
function GeneticForest(vector::BitVector, forest::DecisionEnsemble, n_trees)
    expected_len = n_trees * 3
    length(vector) == expected_len || error(
        "Incompatible vector and forest: expected $(expected_len) bits, got $(length(vector))"
    )
    return GeneticForest(vector, forest)
end


# ─────────────────────────────────────────────────────────────────────────────
# ACCESSORS
# ─────────────────────────────────────────────────────────────────────────────

"""
    vector(gf) → BitVector

Return the bitvector stored in `gf`.
"""
vector(gf::GeneticForest) = gf.vector

"""
    original_forest(gf) → DecisionEnsemble

Return the reference forest stored in `gf`.
"""
original_forest(gf::GeneticForest) = gf.original_forest


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTIONAL UPDATERS
# ─────────────────────────────────────────────────────────────────────────────

"""
    with_vector(gf, new_vector) → GeneticForest

Return a new `GeneticForest` identical to `gf` but with `vector` replaced.
The original forest reference is preserved unchanged.
"""
with_vector(gf::GeneticForest, new_vector::BitVector) =
    GeneticForest(new_vector, gf.original_forest)

"""
    with_forest(gf, new_forest) → GeneticForest

Return a new `GeneticForest` identical to `gf` but with `original_forest`
replaced. The bitvector is preserved unchanged.
"""
with_forest(gf::GeneticForest, new_forest::DecisionEnsemble) =
    GeneticForest(gf.vector, new_forest)

end  # module types