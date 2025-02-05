<div align="center"><a href="https://github.com/aclai-lab/Sole.jl"><img src="logo.png" alt="" title="This package is part of Sole.jl" width="200"></a></div>

# SolePostHoc.jl – Post-Hoc Analysis for Symbolic Learning
🚧 This package is under construction. 🚧


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aclai-lab.github.io/SolePostHoc.jl/stable)
[![Build Status](https://api.cirrus-ci.com/github/aclai-lab/SolePostHoc.jl.svg?branch=main)](https://cirrus-ci.com/github/aclai-lab/SolePostHoc.jl)
[![Coverage](https://codecov.io/gh/aclai-lab/SolePostHoc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/aclai-lab/SolePostHoc.jl)

## In a nutshell

*SolePostHoc.jl* is dedicated to post-hoc analysis and optimization of symbolic learning models. It provides a comprehensive suite of algorithms for:
- Rule extraction from both symbolic and non-symbolic models
- Rule minimization and optimization
- Model transformation and enhancement
- Interpretability analysis

## Key Features

### Rule Extraction and Model Optimization
- Extraction of comprehensible rules from complex models
- Support for various source models:
  - Decision trees and random forests
  - Black-box models
- A clean Rule extraction interface (`SolePostHoc.extractrules`)
- Implementation of state-of-the-art algorithms:
  - LUMEN (L: Logic-driven U: Unified M: Minimal E: Extractor of N: Notions)
  - InTrees (Interpret Tree Ensembles)
  - BELLATRIX
- have binding with other state-of-the-art algorithms
	- ModalETEL
	- binding to intrees
	- binding to RuleCOSI(+)
  - BATrees (Born Again Trees)
### Through these we guarantee
- Rule minimization techniques
- Model simplification while preserving accuracy
- Performance enhancement through post-processing

### Integration
- Seamless integration with other Sole.jl packages


## Usage Example

```julia
# Load packages
using SolePostHoc
using SoleModels
using MLJ

# Load and prepare a model (e.g., a random forest)
🌳 = load_model("your_model.jl")

# Extract rules
🍃 = extractrules(🌳, method = :LUMEN)

# View metrics
printmetrics(🍃)
```


## Want to know more?
For the theoretical foundations of Sole framework, refer to:
[*Modal Symbolic Learning: from theory to practice*, G. Pagliarini (2024)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=FRo4yrcAAAAJ&citation_for_view=FRo4yrcAAAAJ:LkGwnXOMwfcC)

## About

The package is developed by the [ACLAI Lab](https://aclai.unife.it/en/) @ University of Ferrara.

*SolePostHoc.jl* is part of the [*Sole.jl*](https://github.com/aclai-lab/Sole.jl) ecosystem, working alongside:
- [SoleLogics.jl](https://github.com/aclai-lab/SoleLogics.jl): Logical foundations
- [SoleData.jl](https://github.com/aclai-lab/SoleData.jl): Data handling
- [SoleModels.jl](https://github.com/aclai-lab/SoleModels.jl): Model definitions
- [SoleFeatures.jl](https://github.com/aclai-lab/SoleFeatures.jl): Feature engineering