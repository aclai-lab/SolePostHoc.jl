# [Getting started](@id getting-started)

In this introductory section, you will learn about the main building blocks of SolePostHoc.jl.
The above introduces two important ideas for using post-hoc explanation algorithms.
Further on in the documentation, the potential of SolePostHoc.jl will become apparent: this package's primary purpose is to provide a uniform interface for knowledge extraction algorithms, enabling the comparison of different post-hoc interpretation methods while maintaining a coherent and intuitive user experience.

## Fast introduction

Consider a machine learning model trained on a generic dataset. For example, let us consider a Random Forest Classifier learned on the [Iris](https://www.kaggle.com/datasets/uciml/iris) dataset to classify 3 different species of flowers. We are interested in extracting interpretable rules that explain the model's decision process. SolePostHoc.jl offers two primary methods for accomplishing this task.

The first approach is to directly call the specific algorithm function. For example:

```julia
# Extract rules using the LUMEN algorithm directly
extracted_rules = lumen(model, X_test, y_test, args...)
```

The second approach uses the unified interface through rule extractors:

```julia
# Extract rules using the unified interface
extractor = LumenRuleExtractor()
decision_set = modalextractrules(extractor, model, X_test, y_test, args...)
```

The key advantage of the second approach is that it not only executes the original algorithm (equivalent to calling `lumen(...)` directly) but also converts the output into a `DecisionSet`. A `DecisionSet` is a vector of propositional logical rules in Disjunctive Normal Form (DNF), with one rule per class/label.

Consider a trained model that classifies hand gestures. Using SolePostHoc.jl, we might extract the following decision set:

```
Class "Iris-setosa": IF (SepalLengthCm < -0.5) AND (SepalWidthCm < 8.2) THEN predict "Iris-setosa"
Class "Iris-versicolor": IF (SepalLengthCm > 0.5) AND (SepalWidthCm < 3.25) THEN predict "Iris-versicolor"
Class "Iris-virginica": IF (PetalWidthCm > 2.0) THEN predict "Iris-virginica"
```

## Core definitions

The foundation of SolePostHoc.jl lies in providing interpretable explanations for complex machine learning models through rule extraction.

```julia
abstract type RuleExtractor
```

A `RuleExtractor` is an abstract type that defines the interface for all post-hoc explanation algorithms. Each concrete implementation represents a specific knowledge extraction method.

A `DecisionSet` represents the extracted knowledge as a collection of logical rules, where each rule corresponds to a specific class or decision outcome in Disjunctive Normal Form.

The main entry point for rule extraction is:

```julia
modalextractrules(extractor::RuleExtractor, model, args...)
```

## Algorithm Types

SolePostHoc.jl integrates a wide range of algorithms for knowledge extraction, categorized into three main types:

### Surrogate Trees
Algorithms that approximate complex models such as neural networks or random forests with more interpretable decision trees.

```julia
struct REFNERuleExtractor <: RuleExtractor end
struct BATreesRuleExtractor <: RuleExtractor end 
struct TREPANRuleExtractor <: RuleExtractor end
```

### Knowledge Distillation
Techniques for transferring knowledge from complex models (teacher) to simpler and more transparent ones (student).

```julia
struct RuleCOSIPLUSRuleExtractor <: RuleExtractor end
struct InTreesRuleExtractor <: RuleExtractor end
```

### Rule Extraction
Methods for deriving clear and understandable logical rules from any machine learning model.

```julia
struct LUMENRuleExtractor <: RuleExtractor end
```

## Direct Algorithm Access

For users who prefer to use algorithms in their original form without the unified interface, SolePostHoc.jl provides direct access to each algorithm:

```@docs
intrees
lumen
batrees
refne
trepan
rulecosiplus
```

## Rule Extraction, simplification and Optimization

One of the key features of SolePostHoc.jl is its ability to extract, simplify and optimize extracted rules while maintaining their expressiveness.

For example, consider this decision forest:

```
├[1/2]┐ (V3 < 2.45)
│     ├✔ Iris-setosa
│     └✘ (V4 < 1.75)
│       ├✔ (V3 < 4.65)
│       │ ├✔ Iris-versicolor
│       │ └✘ Iris-versicolor
│       └✘ Iris-virginica
└[2/2]┐ (V4 < 0.8)
      ├✔ Iris-setosa
      └✘ (V1 < 5.65)
        ├✔ (V4 < 1.2)
        │ ├✔ Iris-versicolor
        │ └✘ Iris-versicolor
        └✘ (V3 < 4.75)
          ├✔ Iris-versicolor
          └✘ Iris-virginica
```

SolePostHoc.jl can leverage logical reasoning to obtain a more succinct and equally expressive theory:

```
▣
├[1/3] ((V3 ≥ 2.45) ∧ (V4 ≥ 1.75)) ∨ ((V1 ≥ 5.65) ∧ (V3 ≥ 4.75) ∧ (V4 ≥ 0.8))  ↣  Iris-virginica
├[2/3] ((V3 ≥ 2.45) ∧ (V3 < 4.65) ∧ (V4 < 1.75)) ∨ ((V3 ≥ 4.65) ∧ (V3 < 4.75) ∧ (V4 < 1.75)) ∨ ((V1 < 5.65) ∧ (V3 ≥ 4.65) ∧ (V4 < 1.75)) ∨ ((V3 ≥ 2.45) ∧ (V4 < 0.8))  ↣  Iris-versicolor
└[3/3] (V3 < 2.45)  ↣  Iris-setosa

```


## Customization and Extension

Users can implement their own rule extraction algorithms by extending the `RuleExtractor` interface:

```julia

function algorithm(model, args...)
    # ordinary function
    return output  # regular generic type of output
end

struct MyCustomExtractor <: RuleExtractor
    # algorithm-specific parameters
end

function modalextractrules(extractor::MyCustomExtractor, model, args...)
    # implement your custom convert `generic type of output in decision set` logic
    # return a DecisionSet
end
```

## Integration with Sole.jl Ecosystem

SolePostHoc.jl seamlessly integrates with the broader Sole.jl ecosystem, particularly:

- **SoleLogics.jl**: For modal logic reasoning and formula manipulation
- **SoleData.jl**: For handling multivariate time series and relational data structures
- **SoleModels.jl**: For interpretable model training and symbolic learning
