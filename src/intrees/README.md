# InTrees Algorithm Documentation 🌳

## Overview

InTrees (Interpreting Tree Ensembles) is a Julia implementation of an algorithm designed to extract interpretable rules from tree ensemble models. The implementation extends the original algorithm to work with any symbolic model, making it more versatile than the original paper's implementation.

## Core Components

### Main Function
```julia
intrees(model, X, y::AbstractVector{<:Label}; kwargs...)::DecisionList
```

### Key Parameters

- `model`: Input model (AbstractModel or DecisionForest)
- `X`: Feature dataset
- `y`: Label vector
- `prune_rules`: Boolean to enable/disable rule pruning (default: true)
- `pruning_s`: Denominator limit in pruning metric (default: 1.0e-6)
- `pruning_decay_threshold`: Threshold for joint removal in pruning (default: 0.05)
- `rule_selection_method`: Rule selection method (currently only :CBC supported)
- `rule_complexity_metric`: Metric for rule complexity (default: :natoms)
- `min_coverage`: Minimum rule coverage for STEL (default: 0.01)

## Algorithm Components

### 1. Rule Extraction
- Extracts initial ruleset from the model
- Handles both ensemble and single models
- Uses shortforms and normalization

### 2. Rule Pruning
- Optional pruning phase for rule optimization
- Parallel processing for improved performance
- Removes redundant or irrelevant conjuncts

### 3. Rule Selection (CBC)
Features:
- Correlation-based feature selection
- Random Forest importance calculation
- Threshold-based feature filtering (>0.01 importance)
- Multi-criteria sorting (importance, error, complexity)

### 4. STEL (Sequential Covering)
Key features:
- Iterative rule selection
- Coverage-based instance removal
- Default rule handling
- Multi-criteria optimization:
  1. Minimum error
  2. Maximum coverage
  3. Minimum length
  4. Random selection for ties

## Dependencies

- DecisionTree
- ComplexityMeasures
- SoleModels (and related SOLE framework components)

## Workflow

1. Rule Extraction
2. Optional Rule Pruning
3. Rule Selection via CBC
4. Sequential Covering (STEL)
5. Decision List Construction

## Metrics Used

- Coverage
- Error rate
- Rule complexity
- Information gain
- Symmetrical uncertainty

## Key Features

- Parallel processing support
- Extensible to any symbolic model
- Built-in optimization techniques
- Comprehensive rule evaluation
- Automatic default rule generation

## Output

Returns a DecisionList containing:
- Selected and pruned rules
- Optimized rule order
- Default rule for uncovered cases

## Reference

Deng, Houtao. "Interpreting tree ensembles with intrees." International Journal of Data Science and Analytics 7.4 (2019): 277-287.

## Usage Example

```julia
# Basic usage
model = load_your_model()
X, y = load_your_data()
rulelist = intrees(model, X, y)

# With custom parameters
rulelist = intrees(model, X, y;
    prune_rules = true,
    pruning_s = 1.0e-6,
    min_coverage = 0.02,
    rule_selection_method = :CBC
)
```

## Notes

🔬 **Research Implementation**
- Focuses on interpretability
- Handles complex ensemble models
- Supports custom symbolic models
- Includes performance optimizations

## Performance Tips

1. Enable parallel processing for large datasets
2. Adjust pruning parameters for balance
3. Tune min_coverage for rule set size
4. Consider rule complexity metric choice