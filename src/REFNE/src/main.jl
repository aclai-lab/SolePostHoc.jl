module REFNE

using Revise
using Pkg
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

function refne(f, Xmin, Xmax; L=nothing, perc=1.0, max_depth=-1, n_subfeatures=-1, partial_sampling=0.7, min_samples_leaf=5, min_samples_split=2, min_purity_increase=0.0, seed=3)

    ddf = generate_univers_of_combinations(Xmin, Xmax, L)
    y_pred = apply(
        f,
        SoleData.scalarlogiset(
            DataFrame(ddf, :auto);
            allow_propositional=true,
        ),
    )

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
        rng=seed,
    )


    println(model)

    f = solemodel(model) # f = solemodel(model; classlabels = labels)
    println(f)

    return f
end

end
