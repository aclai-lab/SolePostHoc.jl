module REFNE

# using Revise
# using Pkg
using IterTools, DataFrames
using DecisionTree: load_data, build_forest, apply_forest

include("utils.jl")
include("apiREFNESole.jl")
# set of classification parameters and respective default values

# n_subfeatures: number of features to consider randomly for each split (default: -1, sqrt(# features))
# n_trees: number of trees to train (default: 10)
# partial_sampling: fraction of samples to train each tree on (default: 0.7)
# max_depth: maximum depth of the decision trees (default: no maximum (-1))
# min_samples_leaf: minimum number of samples each leaf must have (default: 5)
# min_samples_split: minimum number of samples required for a split (default: 2)
# min_purity_increase: minimum purity required for a split (default: 0.0)
# keyword rng: the random number generator or seed to use (default Random.GLOBAL_RNG)
# multi-threaded forests must be initialized with an `Int`
export refne

"""
    refne(m, Xmin, Xmax; L=100, perc=1.0, max_depth=-1, n_subfeatures=-1, 
          partial_sampling=0.7, min_samples_leaf=5, min_samples_split=2, 
          min_purity_increase=0.0, seed=3)

Extract interpretable rules from a trained neural network ensemble using decision tree approximation.

This implementation follows the REFNE-a (Rule Extraction From Neural Network Ensemble) algorithm,
which approximates complex neural network behavior with an interpretable decision tree model.

# Arguments
- `m`: Trained neural network model to extract rules from
- `Xmin`: Minimum values for each input feature
- `Xmax`: Maximum values for each input feature
- `L`: Number of samples to generate in the synthetic dataset (default: 100)
- `perc`: Percentage of generated samples to use (default: 1.0)
- `max_depth`: Maximum depth of the decision tree (default: -1, unlimited)
- `n_subfeatures`: Number of features to consider at each split (default: -1, all)
- `partial_sampling`: Fraction of samples used for each tree (default: 0.7)
- `min_samples_leaf`: Minimum number of samples required at a leaf node (default: 5)
- `min_samples_split`: Minimum number of samples required to split a node (default: 2)
- `min_purity_increase`: Minimum purity increase required for a split (default: 0.0)
- `seed`: Random seed for reproducibility (default: 3)

# Returns
- A forest-decision trees representing the extracted rules

# Description
The algorithm works by:
1. Generating a synthetic dataset spanning the input space
2. Using the neural network to label these samples
3. Training a decision tree to approximate the neural network's behavior

# References
- Zhi-Hua, Zhou, et al. Extracting Symbolic Rules from Trained Neural Network Ensembles

# Example
```julia
model = load_decision_tree_model()
refne(model, Xmin, Xmax)

See also
[`AbstractModel`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
function refne(
    f,
    Xmin,
    Xmax;
    L = 100,
    perc = 1.0,
    max_depth = -1,
    n_subfeatures = -1,
    partial_sampling = 0.7,
    min_samples_leaf = 5,
    min_samples_split = 2,
    min_purity_increase = 0.0,
    seed = 3,
    ott_mode = false,
)

    if ott_mode == true
        ddf = generate_univers_of_combinations_ott(Xmin, Xmax, L)
    else
        ddf = generate_univers_of_combinations(Xmin, Xmax, L)
    end

    y_pred =
        apply(f, SoleData.scalarlogiset(DataFrame(ddf, :auto); allow_propositional = true))

    n_trees = 1

    model = build_forest(
        y_pred,
        ddf,
        n_subfeatures,
        n_trees,
        partial_sampling,
        max_depth,
        min_samples_leaf,
        min_samples_split,
        min_purity_increase;
        rng = seed,
    )


    println(model)

    f = solemodel(model) # f = solemodel(model; classlabels = labels)
    println(f)

    return f
end

end
