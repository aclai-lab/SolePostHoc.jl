# Documentation Project Lumen 🕯️

## Function Index:

```.jl
📦 Lumen
┣━━ NEW STRUCTURE CREATED
┃   ┣━━ TwoLevelDNFFormula
┃   ┣━━ generate_disjunct
┃   ┣━━ generate_custom_or_formula_from_mask
┃   ┗━━ evaluate (for TwoLevelDNFFormula, Atom, ScalarCondition, Connectives, SyntaxBranch)
┣━━ MAIN FUNCTIONS
┃   ┣━━ truth_combinations
┃   ┣━━ truth_combinations_opt (optimized version)
┃   ┣━━ process_combination
┃   ┣━━ get_alphabet_to_rules
┃   ┣━━ extract_rules_by_label
┃   ┣━━ rules_of_the_forest
┃   ┣━━ concat_results
┃   ┣━━ print_dnf
┃   ┣━━ print_disjunct
┃   ┗━━ minimize_dnf
┣━━ AUXILIARY FUNCTIONS
┃   ┣━━ print_rules_of_the_forest
┃   ┣━━ create_table
┃   ┣━━ print_estimated_time
┃   ┣━━ print_results_summary
┃   ┣━━ print_detailed_results
┃   ┣━━ print_progress_bar
┃   ┗━━ generate_binary_row
┣━━ TESTING AND VERIFICATION FUNCTIONS
┃   ┣━━ verify_simplification
┃   ┗━━ compare_truth_combinations
┣━━ REPORTING FUNCTIONS
┃   ┣━━ generate_statistics_report
┃   ┗━━ generate_statistics_report_opt (optimized version)
┗━━ MAIN FUNCTION
    ┗━━ experimenter
```
## Overview

Lumen🕯️  
**L**: Logic-driven
**U**: Unified
**M**: Minimal
**E**: Extractor of
**N**: Notions.

is a framework of Julia packet using the SOLE (SymbOlic LEarning) framework. It includes functionality for generating truth combinations, simplifying logical formulas, and producing detailed statistical reports.
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

### Main Function
```julia
lumen(
    dt,                            # input model
    :mitespresso;                  # minimization scheme
    vertical = 0.8,                # vertical accuracy (0.0-1.0)
    horizontal = 0.9,              # horizontal accuracy (0.0-1.0)
    silent = true                  # suppress output
)
```

### 1. Data Structures

- `MultivariateScalarAlphabet`: A type for handling multiple scalar alphabets.
- `TwoLevelDNFFormula`: A custom type representing a large OR of AND formulas.


#### Description of TwoLevelDNFFormula Struct: In-Depth Analysis
The `TwoLevelDNFFormula` struct is a specialized data structure designed to represent and manipulate large OR formulas composed of AND terms. This structure is particularly useful for handling complex logical expressions in decision tree analysis and rule-based systems.

#### Structure Definition and Visualization

##### Code Structure
```julia
struct TwoLevelDNFFormula <: Formula
    combinations::Vector{Int}  # Trit vector (-1, 0, 1)
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
end
```

##### Component Hierarchy
```ascii
TwoLevelDNFFormula
│
├── combinations
│   └── [TritVector, TritVector, TritVector, ...]
│       Values: -1 (not present / not essential), 0 (present / essential, >= ), 1 (present /essential, <)
│
├── num_atoms: Int
│
├── thresholds_by_feature
│   ├── 1 => [Float64, Float64, ...]
│   ├── 2 => [Float64, Float64, ...]
│   └── ...
│
└── atoms_by_feature
    ├── 1 => [(Float64, Bool), (Float64, Bool), ...]
    ├── 2 => [(Float64, Bool), (Float64, Bool), ...]
    └── ...
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

##### 1. Combinations Vector (Trit Vector)
Ternary representation of terms:
```ascii
Vector[TritVector[-1,0,1,1], TritVector[1,-1,0,1], TritVector[0,1,-1,1], ...]
                  │││└─ Atom 4: 1 (essential, <)
                  ││└── Atom 3: 1 (essential, <)
                  │└─── Atom 2: 0 (essential, >=)
                  └──── Atom 1: -1 (not essential)
```

Extended visualization:
```ascii
Combination Layout
┌─────┬─────┬─────┬─────┐
│ A1  │ A2  │ A3  │ A4  │
├─────┼─────┼─────┼─────┤
│ -1  │  0  │  1  │  1  │
└─────┴─────┴─────┴─────┘
```

[Rest of the structure remains the same as your original, continuing with Feature Thresholds, Atomic Propositions, etc...]

[All other sections of your README remain unchanged, including Core Functions, Utility Functions, Workflow, etc.]

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