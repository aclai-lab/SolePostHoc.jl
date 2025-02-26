# BATrees 🌳 

## Overview

BATrees (Born Again Trees) is a Julia module for building and training adaptive binary decision trees. It provides a flexible framework for creating decision tree ensembles with customizable parameters.

## Function Index

```.jl
📦 BATrees
┣━━ MAIN FUNCTIONS
┃   ┗━━ batrees
┣━━ CORE COMPONENTS
┃   ┣━━ RunModule
┃   ┗━━ TradModule
┗━━ API INTERFACE
    ┗━━ apiBaSole
```

## Dependencies

The module requires the following Julia packages:
- SoleForest
- SoleData
- SoleLogics
- SoleModels

## Main Function

```julia
batrees(
    f=nothing;                    # input SoleForest (optional)
    dataset_name="iris",          # dataset name
    num_trees=10,                 # number of trees
    max_depth=10,                 # maximum tree depth
    dsOutput=true                 # output format flag
)
```

### Parameters
- `f`: Optional SoleForest input model
- `dataset_name`: Name of the dataset (default: "iris")
- `num_trees`: Number of trees to build (default: 10)
- `max_depth`: Maximum depth of each tree (default: 10)
- `dsOutput`: If true, returns DecisionSet; if false, returns SoleTree (default: true)

### Example Usage

```julia
# Basic usage with default parameters
result = batrees()

# Custom configuration with existing forest
forest = create_forest()  # Your forest creation code
ds = batrees(forest, num_trees=15, max_depth=8)
```

## Core Components

### RunModule
Handles the execution and training of binary decision trees.

### TradModule
Manages the translation and conversion of tree structures.

## Output Formats

1. **DecisionSet (dsOutput=true)**
   - Comprehensive set of decision rules
   - Suitable for ensemble analysis

2. **SoleTree (dsOutput=false)**
   - Single tree representation
   - Useful for individual tree analysis

> [!NOTE]
> ### Note
> BATrees is designed for efficient binary decision tree creation and manipulation. While it supports complex tree structures, consider memory usage when working with large datasets or deep trees.

> [!TIP]
> ### Best Practices
> - Start with default parameters for initial testing
> - Adjust `max_depth` based on your data complexity
> - Use `dsOutput=true` for ensemble analysis
> - Consider dataset size when setting `num_trees`