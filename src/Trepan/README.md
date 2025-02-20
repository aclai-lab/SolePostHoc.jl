# TREPAN 🧠

## Overview

TREPAN (Rule Extraction From Neural Network Ensemble) is a Julia implementation for extracting interpretable rules from trained neural network ensembles using decision tree approximation.

## Function Index

```.jl
📦 TREPAN
┣━━ MAIN FUNCTION
┃   ┗━━ TREPAN
┣━━ UTILITIES
┃   ┗━━ generate_univers_of_combinations
┗━━ API INTERFACE
    ┗━━ apiTREPANSole
```

## Dependencies

Required Julia packages:
- IterTools
- DataFrames
- DecisionTree
- Revise
- SoleData
- SoleLogics
- SoleModels

## Main Function

```julia
TREPAN(
    f,                            # trained model
    Xmin,                         # minimum feature values
    Xmax;                         # maximum feature values
    L=100,                        # number of samples
    perc=1.0,                     # sample percentage
    max_depth=-1,                 # maximum tree depth
    n_subfeatures=-1,             # features per split
    partial_sampling=0.7,         # sampling fraction
    min_samples_leaf=5,          # minimum leaf samples
    min_samples_split=2,         # minimum split samples
    min_purity_increase=0.0,     # minimum purity increase
    seed=3                       # random seed
)
```

### Parameters
- `f`: Trained neural network model
- `X`: All dataset (NOT ONLY TRAINNIG)

- `max_depth`: Maximum decision tree depth
- `n_subfeatures`: Number of features per split
- `partial_sampling`: Sample fraction per tree
- `min_samples_leaf`: Minimum samples in leaf nodes
- `min_samples_split`: Minimum samples for splitting
- `min_purity_increase`: Required purity increase
- `seed`: Random seed for reproducibility

### Example Usage

```julia
# Load or create your model
model = load_neural_network_model()

X = compleate_dataset()

# Extract rules
tree = TREPAN(model, X)
```

## Algorithm Steps

1. Generate synthetic dataset spanning input space
2. Label samples using neural network
3. Train decision tree to approximate network behavior
4. Extract interpretable rules from tree structure

## References

- Mark W. Craven, et al. "Extracting Thee-Structured Representations of Thained Networks"

> [!NOTE]
> ### Note
> TREPAN is particularly useful for understanding and interpreting complex neural network behavior through simpler decision tree rules.

> [!TIP]
> ### Optimization Tips
> - Tune `max_depth` for rule complexity vs. accuracy
> - Use `partial_sampling` to control overfitting
> - Experiment with `min_samples_leaf` for rule stability

> [!IMPORTANT]
> ### Performance Considerations
- More complex neural networks may require larger synthetic datasets