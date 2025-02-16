using Pkg
Pkg.activate(".")
using Revise
using Random
using Logging
using Dates
using DataStructures
using ScikitLearn
using SoleModels
using DecisionTree: load_data, build_forest, apply_forest
using AbstractTrees
using SoleData
using SoleData: UnivariateScalarAlphabet

using ModalDecisionTrees
using SoleLogics
using DataFrames

using Base.Threads: Atomic, atomic_add!
using Profile
using ConcurrentCollections
using DelimitedFiles
using StatsBase

using DecisionTree

function learn_and_convert(
    numero_alberi::Int,
    nome_dataset::String,
    max_depth::Int=-1,
)
    start_time = time()
    println(
        "\n\nPART 0 DATASET CONFIGURATION \n",
    )


    features, labels = load_data(nome_dataset)
    features = float.(features)
    labels = string.(labels)
    

    @info "dataset loaded: $nome_dataset correctly... good luck!"

    println(
        "\n\n PART 1 GENERATION OF THE FOREST with decisionTree.jl \n",
    )

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

    n_subfeatures = -1
    n_trees = numero_alberi
    partial_sampling = 0.7
    #max_depth = -1              # from 6 it becomes too much... already with 1...
    min_samples_leaf = 5
    min_samples_split = 2
    min_purity_increase = 0.0
    seed = 3

    model = build_forest(
        labels,
        features,
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

    println(
        "\n\n PART 2 CONVERSION OF THE FOREST INTO SOLE with solemodel variant \n",
    )

    f = solemodel(model) # f = solemodel(model; classlabels = labels)
    println(f)

    return f, model, start_time
end