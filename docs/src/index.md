```@meta
CurrentModule = SolePostHoc
```

# SolePostHoc

Welcome!!! to the documentation for [SolePostHoc](https://github.com/aclai-lab/SolePostHoc.jl).

## Installation

To install SolePostHoc, simply launch:
```julia
using Pkg
Pkg.add("SolePostHoc")
```

## [Feature](@id feature)
SolePostHoc.jl provides knowledge extraction algorithms through a uniform interface, allowing for the comparison of different post-hoc interpretation methods while maintaining a coherent and intuitive user experience:

```
struct ALGORITHMNAME <: RuleExtractor end
modalextractrules(:ALGORITHMNAME, model, args...)
```

SolePostHoc.jl integrates a wide range of algorithms for knowledge extraction, including:

- **Surrogate Trees**, algorithms that approximate complex models such as neural networks or random forests with more interpretable decision trees;
- **Knowledge Distillation**, techniques for transferring knowledge from complex models to simpler and more transparent ones;
- **Rule Extraction**, methods for deriving clear and understandable logical rules from any machine learning model.

## About

The package is developed by the [ACLAI Lab](https://aclai.unife.it/en/) @ University of Ferrara.

*ModalAssociationRules.jl* lives in the context of [*Sole.jl*](https://github.com/aclai-lab/Sole.jl), an open-source framework for *symbolic machine learning*, originally designed for machine learning based on modal logics (see [Eduard I. Stan](https://eduardstan.github.io/)'s PhD thesis *'Foundations of Modal Symbolic Learning'* [here](https://www.repository.unipr.it/bitstream/1889/5219/5/main.pdf)).
