# Documentazione Progetto Lumen 🕯️

## Indice delle Funzioni:

```.jl
📦 Lumen
┣━━ NUOVA STRUTTURA CREATA
┃   ┣━━ TwoLevelDNFFormula
┃   ┣━━ generate_disjunct
┃   ┣━━ generate_custom_or_formula_from_mask
┃   ┗━━ evaluate (per TwoLevelDNFFormula, Atom, ScalarCondition, Connectives, SyntaxBranch)
┣━━ FUNZIONI PRINCIPALI
┃   ┣━━ truth_combinations
┃   ┣━━ truth_combinations_ott (versione ottimizzata)
┃   ┣━━ process_combination
┃   ┣━━ get_alphabet_to_rules
┃   ┣━━ extract_rules_by_label
┃   ┣━━ rules_of_the_forest
┃   ┣━━ concat_results
┃   ┣━━ stampa_dnf
┃   ┣━━ stampa_disjunct
┃   ┗━━ minimizza_dnf
┣━━ FUNZIONI AUSILIARIE
┃   ┣━━ print_rules_of_the_forest
┃   ┣━━ create_table
┃   ┣━━ print_estimated_time
┃   ┣━━ print_results_summary
┃   ┣━━ print_detailed_results
┃   ┣━━ print_progress_bar
┃   ┗━━ generate_binary_row
┣━━ FUNZIONI DI TESTING E VERIFICA
┃   ┣━━ verify_simplification
┃   ┗━━ compare_truth_combinations
┣━━ FUNZIONI DI REPORTING
┃   ┣━━ genera_report_statistiche
┃   ┗━━ genera_report_statistiche_ott (versione ottimizzata)
┗━━ FUNZIONE PRINCIPALE
    ┗━━ experimenter
```
## Overview

Lumen🕯️  
**L**: Logic-driven
**U**: Unified
**M**: Minimal
**E**: Extractor of
**N**: Notions.

is an framework of Julia packet using the SOLE (SymbOlic LEarning) framework. It includes functionality for generating truth combinations, simplifying logical formulas, and producing detailed statistical reports.
The ultimate goal is the extraction of Minimal Rules from Decision Forests with a Systematic Approach.
it is positioned as a script in SolePostHoc.jl

> [!NOTE] 
>## Dependencies
> 
>The script uses several Julia packages, including:
>
> -AbstractTrees
>
> -Base
>
> -Base.Threads
>
> -BenchmarkTools
>
> -StatProfilerHTML
>
> -Profile
>
> -Test
>
> -ConcurrentCollections
>
> -DataFrames
>
> -DataStructures
>
> -Dates
>
> -DecisionTree
>
> -Logging
>
> -ModalDecisionTrees
>
> -PrettyTables
>
> -ProgressMeter
>
> -Random
>
> -Revise
>
> -SoleData
>
> -SoleDecisionTreeInterface
>
> -SoleLogics
>
> -SoleModels

## Main Components

### 1. Data Structures

- `MultivariateScalarAlphabet`: A type for handling multiple scalar alphabets.
- `TwoLevelDNFFormula`: A custom type representing a large OR of AND formulas.


#### Description of my TwoLevelDNFFormula Struct: In-Depth Analysis
The `TwoLevelDNFFormula` struct is a specialized data structure designed to represent and manipulate large OR formulas composed of AND terms. This structure is particularly useful for handling complex logical expressions in decision tree analysis and rule-based systems.

#### Structure Definition and Visualization

##### Code Structure
```julia
struct TwoLevelDNFFormula <: Formula
    combinations::Vector{Int}
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
    prime_mask::Vector{Vector{Int}}
end
```

##### Component Hierarchy
```ascii
TwoLevelDNFFormula
│
├── combinations
│   └── [BitVector, BitVector, BitVector, ...]
│
├── num_atoms: Int
│
├── thresholds_by_feature
│   ├── 1 => [Float64, Float64, ...]
│   ├── 2 => [Float64, Float64, ...]
│   └── ...
│
├── atoms_by_feature
│   ├── 1 => [(Float64, Bool), (Float64, Bool), ...]
│   ├── 2 => [(Float64, Bool), (Float64, Bool), ...]
│   └── ...
│
└── prime_mask
    └── [[Int, Int, ...], [Int, Int, ...], ...]
```

##### Formula Structure Visualization
```ascii
                       OR Formula
                           │
                 ┌─────────┴─────────┐
                AND                 AND
                 │                   │
        ┌────────┼────────┐     ┌────┴────┐
    Atom1     Atom2    Atom3   Atom1   Atom4
```

#### Component Details

##### 1. Combinations Vector
Binary representation of AND terms:
```ascii
Vector[BitVector[0101], BitVector[1100], BitVector[0011], ...]
                 │││└─ Atom 4 (<)   (LSB)
                 ││└── Atom 3 (≥)
                 │└─── Atom 2 (<)
                 └──── Atom 1 (≥)   (MSB)
```

Extended visualization:
```ascii
Combination Layout
┌─────┬─────┬─────┬─────┐
│ A1  │ A2  │ A3  │ A4  │
├─────┼─────┼─────┼─────┤
│  0  │  1  │  0  │  1  │
└─────┴─────┴─────┴─────┘
```

##### 2. Feature Thresholds
```ascii
Thresholds Structure
┌─Feature 1─────┐  ┌─Feature 2──┐  ┌─Feature 3──────────┐
│ 2.5, 3.7, 4.2 │  │ 1.0, 2.0   │  │ 0.5, 1.5, 2.5, 3.5 │
└───────────────┘  └────────────┘  └────────────────────┘
```

##### 3. Atomic Propositions
```ascii
Feature 1                 Feature 2
┌───────┬─────────┐      ┌───────┬─────────┐
│Value  │Operator │      │Value  │Operator │
├───────┼─────────┤      ├───────┼─────────┤
│ 2.5   │   <     │      │ 1.0   │   <     │
│ 3.7   │   ≥     │      │ 2.0   │   ≥     │
│ 4.2   │   <     │      └───────┴─────────┘
└───────┴─────────┘
```

##### 4. Prime Mask Structure
```ascii
Prime Mask Representation
┌─Combination 1─┐  ┌─Combination 2─┐
│  -1  1  -1  0 │  │  1  -1  0  -1 │
└───────────────┘  └───────────────┘
     │                    │
     └─ -1: Non-essential │
        0/1: Essential    │
        bit positions     └─ Different mask
```

##### 5. Logical Flow
```ascii
Input Features  ─────┐
                     │
                     v
┌────────────────────────────┐
│     Threshold Comparison   │
└────────────────────────────┘
            │
            v
┌────────────────────────────┐
│    Binary Combination      │
└────────────────────────────┘
            │
            v
┌────────────────────────────┐
│      Prime Analysis        │
└────────────────────────────┘
            │
            v
     Final Evaluation
```

#### Detailed Components

##### Thresholds by Feature
Example structure:
```julia
{
    1 => [2.5, 3.7, 4.2],
    2 => [1.0, 2.0],
    3 => [0.5, 1.5, 2.5, 3.5]
}
```

##### Atoms by Feature
Example structure:
```julia
{
    1 => [(2.5, true), (3.7, false), (4.2, true)],
    2 => [(1.0, true), (2.0, false)],
    3 => [(0.5, false), (1.5, true), (2.5, false)]
}
```

##### Prime Mask Structure
Example:
```julia
[
    [-1, 1, -1, 0],  # For combination 1
    [1, -1, 0, -1],  # For combination 2
    ...
]
```

#### Logical Representation
The formula follows Disjunctive Normal Form (DNF):
```ascii
(A₁ ∧ A₂ ∧ ... ∧ Aₙ) ∨ (B₁ ∧ B₂ ∧ ... ∧ Bₙ) ∨ ... ∨ (Z₁ ∧ Z₂ ∧ ... ∧ Zₙ)

Where:
A₁, B₁, Z₁ = Atomic propositions
∧ = AND operator
∨ = OR operator
```

#### Implementation Example
```ascii
Decision Boundary

     │
  Y  │    ┌──────┐
     │    │ True │
     │    │      │
     │    └──────┘
     │
     └──────────────
            X
```

#### Core Functionalities
1. **Formula Evaluation**
```ascii
Input → Threshold Check → Combination Match → Result
  │           │               │                │
  └───────────┴───────────────┴────────────────┘
```

2. **Prime Analysis**
```ascii
Combinations → McQuinss → Essential Bits → Optimization
     │            │            │              │
     └────────────┴────────────┴──────────────┘
```

#### Advantages and Limitations

##### Advantages
```ascii
┌─ Memory Efficiency
├─ Fast Evaluation
├─ Flexible Structure
├─ Built-in Semplification
└─ Analysis Capabilities
```

##### Limitations
```ascii
┌─ Fixed After Init
├─ Complex Management
└─ Careful Mapping Required
```

### 2. Core Functions

#### `experimenter(numero_alberi::Int64, nome_dataset::String ,max_depth::Int64=-1,  vertical::Float64=1.0, horizontal::Float64=1.0, ott_mode::Bool=false, check_ott_mode::Bool=false)`

The 'experimenter' function that orchestrates the entire workflow. It takes the following parameters:

- `numero_alberi`: Number of trees in the forest
- `nome_dataset`: Name of the dataset to use
- `max_depth`: Maximum depth of the decision trees (-1 for no limit)
- `ott_mode`: Boolean flag to use optimized truth table generation
- `check_ott_mode`: Boolean flag to check optimized truth table algorithm
- `horizontal`: Percentage (0.0 to 1.0) of features to consider in truth table
  - Controls feature selection during rule generation
  - 1.0 means use all features
  - 0.5 means use 50% of features
  - Example: With 8 features, horizontal=0.75 will use 6 features
- `vertical`: Percentage (0.0 to 1.0) of combinations to consider in truth table
  - Controls how many rows of the truth table to include
  - 1.0 means use all combinations
  - 0.5 means use 50% of combinations
  - Example: With 16 combinations, vertical=0.25 will use 4 combinations

###### Truth Table Selection Visualization:

```
Original Truth Table:        With horizontal=0.5, vertical=0.75:
Features →                   Features →
│ F1 F2 F3 F4              │ F1 F2 -- --   (50% of features)
├─────────────             ├───────
│ 1  1  1  0               │ 1  1          ┐
│ 1  1  0  0               │ 1  1          │
│ 1  0  1  0               │ 1  0          │ 75% of
│ 0  1  1  0               │ 0  1          │ combinations that take L1
│ 1  0  0  1               │ --            ┘
│ 0  1  0  1               │ --
```

###### Selection Process:

1. **Horizontal (Feature) Selection:**
```
8 Total Features Example with horizontal = 0.75:
F1 F2 F3 F4 F5 F6 F7 F8  →  F1 F2 F3 F4 F5 F6 -- --
                             └───── Selected Features (6) ─────┘
```

2. **Vertical (Combinations) Selection:**
```
16 Total Combinations Example with vertical = 0.25:
Combination 1    →   Combination 1
Combination 2    →   Combination 2
Combination 3    →   Combination 3
Combination 4    →   Combination 4
Combination 5    →   ---
...              →   ---
Combination 16   →   ---
```

###### Usage Examples:

```julia
# Use all features and combinations
experimenter(10, "iris", -1, 1.0, 1.0, true, false)

# Use 75% of features, all combinations
experimenter(10, "iris", -1, 1.0, 0.75, true, false)

# Use all features, 50% of combinations
experimenter(10, "iris", -1, 0.5, 1.0, true, false)

# Use 60% of features, 30% of combinations
experimenter(10, "iris", -1, 0.3, 0.6, true, false)
```

###### Implementation Effects:

1. **Feature Selection (horizontal):**
   - Reduces the dimensionality of the problem
   - Affects the complexity of the generated rules
   - Can improve computational efficiency
   - May impact rule accuracy and generalization

2. **Combination Selection (vertical):**
   - Controls the number of rules considered
   - Affects the coverage of the solution space
   - Can significantly reduce computation time
   - May impact rule completeness

###### Statistical Reporting:

The report will include:
- Number of features selected vs total features
- Number of combinations used vs total combinations
- Impact on rule coverage and accuracy
- Computation time savings

#### `truth_combinations_ott(model, alphabet, atoms)`

An optimized version of truth table generation.

#### `truth_combinations(model, alphabet, atoms)`

The standard version of truth table generation.

#### `minimizza_dnf(formula::TwoLevelDNFFormula)`

Simplifies a TwoLevelDNFFormula through a Quine-McCluskey map-based simplification algorithm.
Returns custom OR formulas along with their associated prime implicant masks (prime_mask).

#### `generate_custom_or_formula_from_mask(formula::TwoLevelDNFFormula)`
This function generates custom OR formulas from a TwoLevelDNFFormula and its prime_mask through two main steps:

Simplification Phase:

Analyzes each feature and its associated atoms
Keeps only the most stringent conditions per feature:

For '<' operators, retains only the smallest value
For '≥' operators, retains only the largest value


Eliminates redundant conditions by setting them to -1 in the mask


Generation Phase:

Converts the simplified mask into readable strings
Creates expressions like 'Vfeature < value' or 'Vfeature ≥ value'
Combines conditions using AND (∧) operators
Prints the resulting rows separated by newlines



Output: Produces a logical formula in disjunctive normal form (OR of ANDs).

#### `verify_simplification(original::TwoLevelDNFFormula, simplified::TwoLevelDNFFormula)`

Verifies that the simplified formula is logically equivalent to the original.

### 3. Utility Functions

- `generate_binary_row`: Generates a binary row representation.
- `get_alphabet_to_rules`: Extracts the alphabet from the rules.
- `extract_rules_by_label`: Extracts rules for each label from the decision forest.
- `rules_of_the_forest`: Gets all rules from the forest.
- `concat_results`: Concatenates results and creates OR formulas of AND formulas.

### 4. Reporting Functions

- `genera_report_statistiche_ott`: Generates a statistical report for the optimized version.
- `genera_report_statistiche`: Generates a statistical report for the standard version.
- `stampa_dnf`: Prints OR of AND formulas.

## Workflow

1. Load and prepare the dataset.
2. Generate a decision forest using DecisionTree.jl.
3. Convert the forest to SOLE format.
4. Extract rules and atoms from the forest.
5. Generate the alphabet from the rules.
6. Generate truth combinations (either optimized or standard method).
7. Evaluate and simplify the resulting formulas.
8. Verify the simplification.
9. Generate a comprehensive statistical report.

## Optimizations

The script includes an optimized version (`truth_combinations_ott`) for generating truth combinations, which can be significantly faster for larger datasets or more complex forests.

## Parallel Processing

The script utilizes Julia's parallel processing capabilities to speed up computations, especially in the truth table generation and formula simplification steps.

## Output

The script generates a detailed statistical report, including:

- Information about rules and atoms
- Distribution of labels
- Simplified formulas
- Performance metrics

## Usage

To use the script, call the `experimenter` function with appropriate parameters. For example:

```julia
experimenter(10, "iris", -1, true)
```

This will run the analysis on the iris dataset, using 10 trees, no maximum depth limit, and the optimized truth table generation mode.

> [!NOTE]
> ### Note
>🚧 This script is designed for research and analysis purposes. For very large datasets or complex forests, computational resources and time requirements can be significant. 🚧
>1. **Research Prototype:** Lumen 🕯️ is primarily a research prototype aimed at exploring and validating complex algorithms and methodologies in symbolic learning and decision forest analysis.
>2. **Rapid Iteration:** The monolithic structure in the function allows for quick iterations and modifications during the research phase, facilitating agile experimentation and hypothesis testing.
>3. **Performance Analysis:** This structure enables comprehensive performance profiling and optimization of the entire workflow without the overhead of inter-module communications.
>4. **Conceptual Integrity:** By maintaining all components in a single block, we ensure a cohesive implementation of the underlying theoretical framework.

> [!TIP]
> ### Future Development
>While this monolithic approach is suitable for the current research phase, we acknowledge that it may pose challenges for long-term maintenance, scalability, and collaboration. Future versions of Lumen 🕯️ may be refactored into a more modular architecture to address these concerns.

> [!IMPORTANT]
> ### For Contributors and Users
>When working with or extending this code:
>- Exercise caution when making changes, as modifications may have wide-ranging effects.
>- Consider documenting any significant changes or additions thoroughly.
>- Be aware that future refactoring efforts may substantially alter the code structure.


