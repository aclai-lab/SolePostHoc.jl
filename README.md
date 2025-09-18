<div align="center"><a href="https://github.com/aclai-lab/Sole.jl"><img src="logo.png" alt="" title="This package is part of Sole.jl" width="200"></a></div>

# SolePostHoc.jl – Post-Hoc Analysis for Symbolic Learning

[![Stable](https://img.shields.io/badge/docs-stable-9558B2.svg)](https://aclai-lab.github.io/SolePostHoc.jl/dev/)
[![CI](https://img.shields.io/badge/CI-5464F4.svg)](https://github.com/aclai-lab/SolePostHoc.jl/actions/workflows/CI.yml)
[![Last Commit](https://img.shields.io/github/last-commit/aclai-lab/SolePostHoc.jl?color=5464F4)](https://github.com/aclai-lab/SolePostHoc.jl/commits/main)
[![Coverage](https://codecov.io/gh/aclai-lab/SolePostHoc.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aclai-lab/SolePostHoc.jl) 
[![License](https://img.shields.io/github/license/aclai-lab/SolePostHoc.jl?color=389826)](https://github.com/aclai-lab/SolePostHoc.jl/blob/main/LICENSE)
[![Julia](https://img.shields.io/badge/julia-1.10%2B-389826)](https://julialang.org/)
[![Issues](https://img.shields.io/github/issues/aclai-lab/SolePostHoc.jl?color=9558B2)](https://github.com/aclai-lab/SolePostHoc.jl/issues)

## In a nutshell

*SolePostHoc.jl* provides a comprehensive suite of post-hoc analysis and optimization tools for symbolic learning models. The package enables knowledge extraction from both symbolic and non-symbolic models through a unified interface, facilitating the comparison of different interpretation methods while maintaining consistency and ease of use.

## Key Features

### Knowledge Extraction Algorithms
- **LUMEN**;
- **InTrees**;
- **TREPAN**;
- **REFNE**;
- **RuleCOSI+**;
- **BATrees**.

### Model Analysis Capabilities
- Rule extraction from decision trees, random forests, and black-box models
- Model transformation and enhancement through rule minimization
- Performance optimization while preserving interpretability
- Support for surrogate model generation and knowledge distillation

### Unified Interface
All algorithms are accessible through a consistent API via the `SolePostHoc.modalextractrules` function, enabling seamless comparison between different extraction methods.

## Usage

### Direct Algorithm Access
```julia
# Call specific algorithms directly
extracted_rules = lumen(model)
extracted_rules = intrees(model, X_test, y_test)
```

### Unified Interface (Recommended)
```julia
# Use the unified interface with rule extractors
extractor_lumen = LumenRuleExtractor()
extractor_intrees = IntreesRuleExtractor()

decision_set_lumen = modalextractrules(extractor_lumen, model)
decision_set_intrees = modalextractrules(extractor_intrees, model, X_test, y_test)
```

The unified interface converts outputs into `DecisionSet` objects—vectors of propositional logical rules in Disjunctive Normal Form (DNF), with one rule per class label.

### Example: Rule Extraction from Random Forest
```julia
using SolePostHoc

# Assume we have a trained Random Forest on the Iris dataset
rf_model = train_random_forest(X_train, y_train)

# Extract interpretable rules using LUMEN
extractor = LumenRuleExtractor()
interpretable_rules = modalextractrules(extractor, rf_model)

# The result is a DecisionSet with logical rules explaining the model's decisions
```

## Algorithm Categories

**Surrogate Trees**: Approximate complex models with interpretable decision trees.

**Knowledge Distillation**: Transfer knowledge from complex to transparent models. 

**Rule Extraction**: Derive clear logical rules from any machine learning model.

## Integration

*SolePostHoc.jl* is part of the [*Sole.jl*](https://github.com/aclai-lab/Sole.jl) ecosystem:
- [SoleLogics.jl](https://github.com/aclai-lab/SoleLogics.jl): Logical foundations
- [SoleData.jl](https://github.com/aclai-lab/SoleData.jl): Data handling  
- [SoleModels.jl](https://github.com/aclai-lab/SoleModels.jl): Model definitions
- [SoleFeatures.jl](https://github.com/aclai-lab/SoleFeatures.jl): Feature engineering

## References

For theoretical foundations, see: [*Modal Symbolic Learning: from theory to practice*, G. Pagliarini (2024)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=FRo4yrcAAAAJ&citation_for_view=FRo4yrcAAAAJ:LkGwnXOMwfcC)

## About

Developed by the [ACLAI Lab](https://aclai.unife.it/en/) @ University of Ferrara.
