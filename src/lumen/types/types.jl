using SoleData
import SoleLogics: atoms

# Definizione di costanti per la formattazione dell'output
const TITLE = "="^100
const COLORED_TITLE = "\033[1;32m"
const COLORED_INFO = "\033[1;34m"
const COLORED_ULTRA_OTT = "\033[1;35m"
const RESET = "\033[0m"


"""
A reentrant lock used to synchronize access to shared resources.
"""
const results_lock = ReentrantLock()


"""
# TwoLevelDNFFormula

## Overview
The `TwoLevelDNFFormula` struct is a specialized data structure designed to represent and manipulate large OR formulas composed of AND terms. This structure is particularly useful for handling complex logical expressions in decision tree analysis and rule-based systems.

## Structure Definition and Visualization

### Code Structure
```julia
struct TwoLevelDNFFormula <: Formula
    combinations::Vector{TritVector}
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
end
```

### Component Hierarchy
```ascii
TwoLevelDNFFormula
│
├── combinations
│   └── [TritVector, TritVector, TritVector, ...]
│
├── num_atoms: Int
│
├── thresholds_by_feature
│   ├── 1 => [Float64, Float64, ...]
│   ├── 2 => [Float64, Float64, ...]
│   └── ...
│
├── atoms_by_feature
    ├── 1 => [(Float64, Bool), (Float64, Bool), ...]
    ├── 2 => [(Float64, Bool), (Float64, Bool), ...]
    └── ...

```

### Formula Structure Visualization
```ascii
                       OR Formula
                           │
                 ┌─────────┴─────────┐
                AND                 AND
                 │                   │
        ┌────────┼────────┐     ┌────┴────┐
    Atom1     Atom2    Atom3   Atom1   Atom4
```

## Component Details

### 1. Combinations Vector
Trit representation of AND terms using TritVector:
```ascii
[TritVector("0101"), TritVector("1100"), TritVector("0011"), ...]
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

### 2. Feature Thresholds
```ascii
Thresholds Structure
┌─Feature 1─────┐  ┌─Feature 2──┐  ┌─Feature 3──────────┐
│ 2.5, 3.7, 4.2 │  │ 1.0, 2.0   │  │ 0.5, 1.5, 2.5, 3.5 │
└───────────────┘  └────────────┘  └────────────────────┘
```

### 3. Atomic Propositions
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

### 4. Prime Mask interpretation of combination Structure
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

### 5. Logical Flow
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
│   TritVector Combination   │
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

## Detailed Components

### Thresholds by Feature
Example structure:
```julia
{
    1 => [2.5, 3.7, 4.2],
    2 => [1.0, 2.0],
    3 => [0.5, 1.5, 2.5, 3.5]
}
```

### Atoms by Feature
Example structure:
```julia
{
    1 => [(2.5, true), (3.7, false), (4.2, true)],
    2 => [(1.0, true), (2.0, false)],
    3 => [(0.5, false), (1.5, true), (2.5, false)]
}
```

### Combination Structure
Example:
```julia
[
    [-1, 1, -1, 0],  # For combination 1
    [1, -1, 0, -1],  # For combination 2
    ...
]
```

## Logical Representation
The formula follows Disjunctive Normal Form (DNF):
```ascii
(A₁ ∧ A₂ ∧ ... ∧ Aₙ) ∨ (B₁ ∧ B₂ ∧ ... ∧ Bₙ) ∨ ... ∨ (Z₁ ∧ Z₂ ∧ ... ∧ Zₙ)

Where:
A₁, B₁, Z₁ = Atomic propositions
∧ = AND operator
∨ = OR operator
```

## Implementation Example
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

## Core Functionalities
1. **Formula Evaluation**
```ascii
Input → Threshold Check → TritVector Match → Result
  │           │               │                │
  └───────────┴───────────────┴────────────────┘
```

2. **Prime Analysis**
```ascii
Combinations → Minimization → Essential Bits → Optimization
     │            │            │              │
     └────────────┴────────────┴──────────────┘
```

## Advantages and Limitations

### Advantages
```ascii
┌─ Memory Efficiency (BitVector optimized)
├─ Fast Evaluation
├─ Flexible Structure
├─ Built-in Semplification
└─ Analysis Capabilities
```

### Limitations
```ascii
┌─ Fixed After Init
├─ Complex Management
└─ Careful Mapping Required
```
"""
struct TwoLevelDNFFormula <: Formula
    combinations::Vector{TritVector}
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
end

############################Constructs###########################################
# Costruttore che accetta atoms e combinations come TritVector
function TwoLevelDNFFormula(atoms::Vector, combinations::Vector{TritVector})
    # Validate inputs
    isempty(atoms) && throw(ArgumentError("Atoms vector cannot be empty"))

    # Process atoms
    my_atoms = unique(atoms)
    num_atoms = length(my_atoms)

    # Get number of features from the maximum feature index in atoms
    n_features = maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
    !isa(n_features, Integer) && error("Symbolic feature names not supported")

    # Create alphabet
    my_alphabet = process_alphabet(my_atoms, n_features)

    # Create thresholds dictionary from alphabet
    thresholds_by_feature = Dict(
        subalpha.featcondition[1].feature.i_variable => sort(subalpha.featcondition[2])
        for subalpha in my_alphabet.subalphabets
    )

    # Create atoms by feature dictionary
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
    for atom in my_atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(
            get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat),
            (threshold, true)
        )
    end

    # Sort atoms within each feature
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by=first)
    end

    return TwoLevelDNFFormula(combinations, num_atoms, thresholds_by_feature, atoms_by_feature)
end

# Costruttore che accetta atoms e combinations come BigInt
function TwoLevelDNFFormula(atoms::Vector, combinations::Vector{BigInt})
    num_atoms = length(atoms)
    # Converti i BigInt in TritVector usando la stessa logica di dict_to_tritvector
    trit_combinations = [bigint_to_tritvector(comb, num_atoms) for comb in combinations]

    # Usa il costruttore con Vector{TritVector}
    return TwoLevelDNFFormula(atoms, trit_combinations)
end
#################################################################################

############################ GETTER #############################################
"""
    eachcombination(f::TwoLevelDNFFormula) -> Vector{TritVector}

Returns the vector of all combinations (AND terms) in the DNF formula.
Each combination is represented as a TritVector where each position corresponds to an atom.

# Example
```julia
for combination in eachcombination(formula)
    println("Term complexity: \$(count(!iszero, combination))")
end
```
"""
eachcombination(f::TwoLevelDNFFormula) = f.combinations

"""
    eachthresholdsbyfeature(f::TwoLevelDNFFormula) -> Dict{Int,Vector{Float64}}

Returns a dictionary mapping feature indices to their sorted threshold values.
Each feature's thresholds are stored in ascending order.

# Example
```julia
thresholds = eachthresholdsbyfeature(formula)
for (feature, values) in thresholds
    println("Feature \$feature ranges from \$(minimum(values)) to \$(maximum(values))")
end
```
"""
eachthresholdsbyfeature(f::TwoLevelDNFFormula) = f.thresholds_by_feature

"""
    eachatomsbyfeature(f::TwoLevelDNFFormula) -> Dict{Int,Vector{Tuple{Float64,Bool}}}

Returns a dictionary mapping feature indices to their atomic propositions.
Each atom is represented as a tuple (threshold, is_geq) where:
- threshold: The numeric threshold value
- is_geq: true if the operator is ≥, false if the operator is <

# Example
```julia
atoms = eachatomsbyfeature(formula)
for (feature, atom_list) in atoms
    for (threshold, is_geq) in atom_list
        op = is_geq ? "≥" : "<"
        println("x\$feature \$op \$threshold")
    end
end
```
"""
eachatomsbyfeature(f::TwoLevelDNFFormula) = f.atoms_by_feature

"""
    nuberofatoms(f::TwoLevelDNFFormula) -> Int

Returns the total number of atomic propositions in the formula.

# Example
```julia
n = nuberofatoms(formula)
println("Formula contains \$n atomic propositions")
```
"""
nuberofatoms(f::TwoLevelDNFFormula) = f.num_atoms
#################################################################################

############################# SETTER ############################################ 
# TODO SOME TO RE - WATCH 
"""
    set_combinations(formula::TwoLevelDNFFormula, new_combinations::Vector{TritVector}) -> TwoLevelDNFFormula

Creates a new formula with modified combinations while preserving all other properties.
All combinations must be compatible with the existing atomic structure.

# Example
```julia
# Remove the first combination
new_formula = set_combinations(formula, eachcombination(formula)[2:end])
```

# Throws
- ArgumentError if new combinations are incompatible with existing atomic structure
"""
function set_combinations(formula::TwoLevelDNFFormula, new_combinations::Vector{TritVector})
    if any(length(c) != nuberofatoms(formula) for c in new_combinations)
        throw(ArgumentError("New combinations must have the same number of atoms as the original formula"))
    end

    return TwoLevelDNFFormula(
        new_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula)
    )
end

"""
    set_feature_thresholds(formula::TwoLevelDNFFormula, feature::Int, new_thresholds::Vector{Float64}) -> TwoLevelDNFFormula

Creates a new formula with modified thresholds for a specific feature.
Maintains consistency with atomic propositions.

# Example
```julia
# Update thresholds for feature 1
new_formula = set_feature_thresholds(formula, 1, [1.0, 2.0, 3.0])
```

# Throws
- KeyError if feature doesn't exist in the formula
- ArgumentError if new thresholds are incompatible with existing atoms
"""
function set_feature_thresholds(formula::TwoLevelDNFFormula, feature::Int, new_thresholds::Vector{Float64})
    if !haskey(eachatomsbyfeature(formula), feature)
        throw(KeyError("Feature $feature not found in formula"))
    end

    new_thresholds_by_feature = Dict(eachthresholdsbyfeature(formula))
    new_thresholds_by_feature[feature] = sort(new_thresholds)

    atoms = eachatomsbyfeature(formula)[feature]
    if !all(t -> t[1] ∈ new_thresholds, atoms)
        throw(ArgumentError("New thresholds must include all threshold values used in atomic propositions"))
    end

    return TwoLevelDNFFormula(
        eachcombination(formula),
        nuberofatoms(formula),
        new_thresholds_by_feature,
        eachatomsbyfeature(formula)
    )
end

"""
    set_feature_atoms(formula::TwoLevelDNFFormula, feature::Int, new_atoms::Vector{Tuple{Float64,Bool}}) -> TwoLevelDNFFormula

Creates a new formula with modified atomic propositions for a specific feature.
Updates combinations and thresholds accordingly.

# Example
```julia
# Update atoms for feature 1
new_atoms = [(1.0, true), (2.0, false)]
new_formula = set_feature_atoms(formula, 1, new_atoms)
```

# Throws
- KeyError if feature doesn't exist in the formula
"""
function set_feature_atoms(formula::TwoLevelDNFFormula, feature::Int, new_atoms::Vector{Tuple{Float64,Bool}})
    if !haskey(eachatomsbyfeature(formula), feature)
        throw(KeyError("Feature $feature not found in formula"))
    end

    new_atoms_by_feature = Dict(eachatomsbyfeature(formula))
    new_atoms_by_feature[feature] = sort(new_atoms, by=first)

    new_thresholds_by_feature = Dict(eachthresholdsbyfeature(formula))
    new_thresholds_by_feature[feature] = sort(unique([t for (t, _) in new_atoms]))

    old_atom_count = length(eachatomsbyfeature(formula)[feature])
    new_atom_count = nuberofatoms(formula) - old_atom_count + length(new_atoms)

    new_combinations = Vector{TritVector}()
    for comb in eachcombination(formula)
        new_comb = TritVector(zeros(Int8, new_atom_count))
        # TODO: Map old combination values to new structure
        push!(new_combinations, new_comb)
    end

    return TwoLevelDNFFormula(
        new_combinations,
        new_atom_count,
        new_thresholds_by_feature,
        new_atoms_by_feature
    )
end

"""
    add_feature(formula::TwoLevelDNFFormula, feature::Int, 
               atoms::Vector{Tuple{Float64,Bool}}) -> TwoLevelDNFFormula

Creates a new formula with an additional feature and its associated atomic propositions.

# Example
```julia
# Add new feature 3 with two atoms
new_atoms = [(1.5, true), (2.5, false)]
new_formula = add_feature(formula, 3, new_atoms)
```

# Throws
- ArgumentError if feature already exists
"""
function add_feature(formula::TwoLevelDNFFormula, feature::Int,
    atoms::Vector{Tuple{Float64,Bool}})
    if haskey(eachatomsbyfeature(formula), feature)
        throw(ArgumentError("Feature $feature already exists in formula"))
    end

    new_atoms_by_feature = Dict(eachatomsbyfeature(formula))
    new_atoms_by_feature[feature] = sort(atoms, by=first)

    new_thresholds_by_feature = Dict(eachthresholdsbyfeature(formula))
    new_thresholds_by_feature[feature] = sort(unique([t for (t, _) in atoms]))

    new_num_atoms = nuberofatoms(formula) + length(atoms)
    new_combinations = [
        TritVector(vcat(Vector{Int8}(c), zeros(Int8, length(atoms))))
        for c in eachcombination(formula)
    ]

    return TwoLevelDNFFormula(
        new_combinations,
        new_num_atoms,
        new_thresholds_by_feature,
        new_atoms_by_feature
    )
end

"""
    remove_feature(formula::TwoLevelDNFFormula, feature::Int) -> TwoLevelDNFFormula

Creates a new formula with the specified feature removed.
Updates all related structures accordingly.

# Example
```julia
# Remove feature 2
new_formula = remove_feature(formula, 2)
```

# Throws
- KeyError if feature doesn't exist
"""
function remove_feature(formula::TwoLevelDNFFormula, feature::Int)
    if !haskey(eachatomsbyfeature(formula), feature)
        throw(KeyError("Feature $feature not found in formula"))
    end

    new_atoms_by_feature = Dict(
        k => v for (k, v) in eachatomsbyfeature(formula)
        if k != feature
    )

    new_thresholds_by_feature = Dict(
        k => v for (k, v) in eachthresholdsbyfeature(formula)
        if k != feature
    )

    feature_atoms = eachatomsbyfeature(formula)[feature]
    remove_count = length(feature_atoms)
    new_num_atoms = nuberofatoms(formula) - remove_count

    new_combinations = Vector{TritVector}()
    for comb in eachcombination(formula)
        new_comb = TritVector(zeros(Int8, new_num_atoms))
        # TODO: Map remaining combination values correctly
        push!(new_combinations, new_comb)
    end

    return TwoLevelDNFFormula(
        new_combinations,
        new_num_atoms,
        new_thresholds_by_feature,
        new_atoms_by_feature
    )
end
#################################################################################

################# Utils for TwoLevelDNFFormula ##################################
# Utility functions for manipulating and analyzing TwoLevelDNFFormula

"""
    count_terms(formula::TwoLevelDNFFormula) -> Int

Returns the number of AND terms (combinations) in the DNF formula.

# Example
```julia
formula = TwoLevelDNFFormula(atoms, combinations)
num_terms = count_terms(formula)  # Returns number of AND terms
```
"""
function count_terms(formula::TwoLevelDNFFormula)
    return length(eachcombination(formula))
end

"""
    count_active_atoms(formula::TwoLevelDNFFormula) -> Dict{Int, Int}

Returns a dictionary mapping feature indices to the number of atoms using that feature.

# Example
```julia
formula = TwoLevelDNFFormula(atoms, combinations)
atom_counts = count_active_atoms(formula)  # Returns {1 => 3, 2 => 2, ...}
```
"""
function count_active_atoms(formula::TwoLevelDNFFormula)
    return Dict(k => length(v) for (k, v) in eachatomsbyfeature(formula))
end

"""
    get_term_complexity(formula::TwoLevelDNFFormula) -> Vector{Int}

Returns a vector containing the number of active atoms in each AND term.
Useful for analyzing the complexity of each combination.

# Example
```julia
complexities = get_term_complexity(formula)  # Returns [3, 2, 4, ...] for each term
```
"""
function get_term_complexity(formula::TwoLevelDNFFormula)
    return [count(!iszero, combination) for combination in eachcombination(formula)]
end

"""
    find_most_used_features(formula::TwoLevelDNFFormula, top_n::Int=3) -> Vector{Tuple{Int, Int}}

Returns the top N most frequently used features in the formula, sorted by usage count.
Returns tuples of (feature_index, count).

# Example
```julia
top_features = find_most_used_features(formula, 2)  # Returns [(1, 5), (3, 3)]
```
"""
function find_most_used_features(formula::TwoLevelDNFFormula, top_n::Int=3)
    feature_counts = count_active_atoms(formula)
    sorted_features = sort(collect(feature_counts), by=x -> x[2], rev=true)
    return sorted_features[1:min(top_n, length(sorted_features))]
end

"""
    get_feature_ranges(formula::TwoLevelDNFFormula) -> Dict{Int, Tuple{Float64, Float64}}

Returns a dictionary mapping feature indices to their (min, max) threshold values.

# Example
```julia
ranges = get_feature_ranges(formula)  # Returns {1 => (2.5, 4.2), 2 => (1.0, 2.0), ...}
```
"""
function get_feature_ranges(formula::TwoLevelDNFFormula)
    ranges = Dict{Int,Tuple{Float64,Float64}}()
    for (feature, thresholds) in eachthresholdsbyfeature(formula)
        ranges[feature] = (minimum(thresholds), maximum(thresholds))
    end
    return ranges
end


"""
    is_consistent(formula::TwoLevelDNFFormula) -> Bool

Checks if the formula's structure is internally consistent by verifying:
1. All referenced atoms exist in atoms_by_feature
2. Combination lengths match num_atoms
3. Thresholds are properly ordered
4. No duplicate thresholds exist

# Example
```julia
if !is_consistent(formula)
    @warn "Formula structure may be corrupted"
end
```
"""
function is_consistent(formula::TwoLevelDNFFormula)
    # Check combination lengths
    if any(length(c) != nuberofatoms(formula) for c in eachcombination(formula))
        return false
    end

    # Check threshold ordering and uniqueness
    for thresholds in values(eachthresholdsbyfeature(formula))
        if !issorted(thresholds) || length(unique(thresholds)) != length(thresholds)
            return false
        end
    end

    # Check atoms match thresholds
    for (feature, atoms) in eachatomsbyfeature(formula)
        if !haskey(eachthresholdsbyfeature(formula), feature)
            return false
        end
        thresholds = Set(t for (t, _) in atoms)
        if !issubset(thresholds, Set(eachthresholdsbyfeature(formula)[feature]))
            return false
        end
    end

    return true
end

"""
    get_formula_statistics(formula::TwoLevelDNFFormula) -> NamedTuple

Returns comprehensive statistics about the formula structure.

# Example
```julia
stats = get_formula_statistics(formula)
println("Average term complexity: \$(stats.avg_term_complexity)")
```
"""
function get_formula_statistics(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)

    # Calculate various statistics
    term_complexities = get_term_complexity(formula)
    return (
        num_terms=length(combinations),
        num_features=length(atoms_by_feat),
        total_atoms=nuberofatoms(formula),
        avg_term_complexity=mean(term_complexities),
        max_term_complexity=maximum(term_complexities),
        min_term_complexity=minimum(term_complexities),
        feature_usage=count_active_atoms(formula),
        complexity_distribution=countmap(term_complexities)
    )
end

"""
    visualize_formula_structure(formula::TwoLevelDNFFormula) -> String

Creates an ASCII visualization of the formula's structure.

# Example
```julia
println(visualize_formula_structure(formula))
```
"""
function visualize_formula_structure(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)

    # Create header
    result = "Formula Structure Visualization\n"
    result *= "═══════════════════════════\n\n"

    # Add formula overview
    result *= "┌── Formula Overview ──┐\n"
    result *= "│ Terms: $(length(combinations))".ljust(20) * "│\n"
    result *= "│ Features: $(length(atoms_by_feat))".ljust(20) * "│\n"
    result *= "│ Atoms: $(nuberofatoms(formula))".ljust(20) * "│\n"
    result *= "└────────────────────┘\n\n"

    # Visualize terms
    result *= "Term Structure:\n"
    for (i, combination) in enumerate(combinations)
        result *= "Term $i: "
        active_count = count(!iszero, combination)
        result *= "[$("█" ^ active_count)$(" " ^ (nuberofatoms(formula) - active_count))] ($active_count atoms)\n"
    end

    # Feature distribution
    result *= "\nFeature Distribution:\n"
    for (feature, atoms) in sort(collect(atoms_by_feat))
        result *= "F$feature: $("•" ^ length(atoms)) ($(length(atoms)) atoms)\n"
    end

    return result
end

"""
    analyze_feature_interactions(formula::TwoLevelDNFFormula) -> Dict

Analyzes how features interact with each other in the formula's terms.

# Example
```julia
interactions = analyze_feature_interactions(formula)
```
"""
function analyze_feature_interactions(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)
    features = collect(keys(atoms_by_feat))

    interactions = Dict{Tuple{Int,Int},Int}()

    for combination in combinations
        active_features = Set{Int}()
        for (feature, atoms) in atoms_by_feat
            if any(!iszero, combination[1:length(atoms)])
                push!(active_features, feature)
            end
        end

        # Count feature pairs
        for f1 in active_features
            for f2 in active_features
                if f1 < f2
                    pair = (f1, f2)
                    interactions[pair] = get(interactions, pair, 0) + 1
                end
            end
        end
    end

    return interactions
end

"""
    get_feature_importance(formula::TwoLevelDNFFormula) -> Vector{Tuple{Int, Float64}}

Calculates feature importance scores based on usage and interactions.
Returns sorted list of (feature_id, importance_score) tuples.

# Example
```julia
importance = get_feature_importance(formula)
```
"""
function get_feature_importance(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)

    # Initialize scores
    scores = Dict{Int,Float64}()

    # Base score from direct usage
    usage_counts = count_active_atoms(formula)
    for (feature, count) in usage_counts
        scores[feature] = count
    end

    # Add interaction scores
    interactions = analyze_feature_interactions(formula)
    for ((f1, f2), count) in interactions
        scores[f1] = get(scores, f1, 0.0) + count * 0.5
        scores[f2] = get(scores, f2, 0.0) + count * 0.5
    end

    # Normalize scores
    max_score = maximum(values(scores))
    scores = Dict(k => v / max_score for (k, v) in scores)

    return sort(collect(scores), by=x -> x[2], rev=true)
end

"""
    generate_latex_representation(formula::TwoLevelDNFFormula) -> String

Generates a LaTeX representation of the formula.

# Example
```julia
latex_str = generate_latex_representation(formula)
```
"""
function generate_latex_representation(formula::TwoLevelDNFFormula)
    thresholds_by_feat = eachthresholdsbyfeature(formula)
    atoms_by_feat = eachatomsbyfeature(formula)
    combinations = eachcombination(formula)

    terms = String[]

    for combination in combinations
        atoms = String[]
        for (feature, atom_list) in atoms_by_feat
            for (i, (threshold, is_geq)) in enumerate(atom_list)
                if !iszero(combination[i])
                    op = is_geq ? "\\geq" : "<"
                    push!(atoms, "x_{$(feature)} $op $(threshold)")
                end
            end
        end
        term = length(atoms) > 1 ? "\\left($(join(atoms, " \\wedge "))\\right)" : atoms[1]
        push!(terms, term)
    end

    return join(terms, " \\vee ")
end

"""
    plot_feature_distribution(formula::TwoLevelDNFFormula) -> String

Creates an ASCII plot showing the distribution of features and their thresholds.

# Example
```julia
println(plot_feature_distribution(formula))
```
"""
function plot_feature_distribution(formula::TwoLevelDNFFormula)
    thresholds_by_feat = eachthresholdsbyfeature(formula)

    result = "Feature Threshold Distribution\n"
    result *= "═══════════════════════════\n\n"

    max_thresholds = maximum(length(t) for t in values(thresholds_by_feat))
    width = 40

    for (feature, thresholds) in sort(collect(thresholds_by_feat))
        # Feature label
        result *= "F$(feature): "

        # Create scale
        min_val, max_val = minimum(thresholds), maximum(thresholds)
        range_val = max_val - min_val

        # Plot thresholds
        positions = Int[]
        for threshold in thresholds
            pos = round(Int, ((threshold - min_val) / range_val) * (width - 1)) + 1
            push!(positions, pos)
        end

        # Draw the line
        line = fill('─', width)
        for pos in positions
            line[pos] = '│'
        end

        result *= join(line)
        result *= " [$(round(min_val, digits=2)), $(round(max_val, digits=2))]\n"
    end

    return result
end

"""
    visualize_feature_coverage(formula::TwoLevelDNFFormula) -> String

Creates an ASCII visualization showing how features are covered across terms.
Uses a matrix-like representation where rows are terms and columns are features.

# Example
```julia
println(visualize_feature_coverage(formula))
```
"""
function visualize_feature_coverage(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)
    features = sort(collect(keys(atoms_by_feat)))

    result = "Feature Coverage Matrix\n"
    result *= "════════════════════\n\n"

    # Header
    result *= "Term│"
    for feature in features
        result *= " F$(lpad(string(feature), 2)) "
    end
    result *= "\n"

    # Separator
    result *= "────┼"
    result *= "────"^length(features)
    result *= "\n"

    # Matrix
    for (i, combination) in enumerate(combinations)
        result *= "$(lpad(string(i), 3)) │"

        atom_idx = 1
        for feature in features
            atoms = atoms_by_feat[feature]
            has_active = any(!iszero, combination[atom_idx:atom_idx+length(atoms)-1])
            result *= " $(has_active ? "■" : "□") "
            atom_idx += length(atoms)
        end
        result *= "\n"
    end

    return result
end

"""
    visualize_threshold_density(formula::TwoLevelDNFFormula) -> String

Creates an ASCII histogram showing the density of thresholds for each feature.
Uses variable-width bars to represent different concentrations of thresholds.

# Example
```julia
println(visualize_threshold_density(formula))
```
"""
function visualize_threshold_density(formula::TwoLevelDNFFormula)
    thresholds_by_feat = eachthresholdsbyfeature(formula)
    max_width = 40

    result = "Threshold Density Distribution\n"
    result *= "══════════════════════════\n\n"

    # Find overall range for normalization
    all_thresholds = vcat(values(thresholds_by_feat)...)
    global_min, global_max = minimum(all_thresholds), maximum(all_thresholds)
    range_size = global_max - global_min

    for (feature, thresholds) in sort(collect(thresholds_by_feat))
        result *= "F$(lpad(string(feature), 2)) │"

        # Create bins
        n_bins = 10
        bin_edges = range(global_min, global_max, length=n_bins + 1)
        bin_counts = zeros(Int, n_bins)

        # Count thresholds in each bin
        for threshold in thresholds
            bin = findlast(edge -> edge <= threshold, bin_edges)
            bin === nothing && continue
            bin > n_bins && continue
            bin_counts[bin] += 1
        end

        # Normalize and display
        max_count = maximum(bin_counts)
        for count in bin_counts
            bar_width = round(Int, (count / max_count) * max_width)
            result *= "█"^bar_width * " "
        end

        result *= " ($(length(thresholds)) thresholds)\n"
    end

    # Add scale
    result *= "\nScale: "
    result *= "$(round(global_min, digits=2)) "
    result *= "─"^(max_width - 20)
    result *= " $(round(global_max, digits=2))\n"

    return result
end

"""
    visualize_term_network(formula::TwoLevelDNFFormula) -> String

Creates an ASCII network visualization showing how terms are connected through shared features.
Terms that share more features are drawn closer together in the visualization.

# Example
```julia
println(visualize_term_network(formula))
```
"""
function visualize_term_network(formula::TwoLevelDNFFormula)
    combinations = eachcombination(formula)
    atoms_by_feat = eachatomsbyfeature(formula)

    result = "Term Relationship Network\n"
    result *= "═════════════════════\n\n"

    # Build adjacency matrix
    n_terms = length(combinations)
    adjacency = zeros(Int, n_terms, n_terms)

    for i in 1:n_terms
        for j in i+1:n_terms
            shared_features = 0
            atom_idx = 1

            for (feature, atoms) in atoms_by_feat
                range = atom_idx:atom_idx+length(atoms)-1
                if any(!iszero, combinations[i][range]) && any(!iszero, combinations[j][range])
                    shared_features += 1
                end
                atom_idx += length(atoms)
            end

            adjacency[i, j] = adjacency[j, i] = shared_features
        end
    end

    # Create simple ASCII visualization
    for i in 1:n_terms
        result *= "T$i "
        for j in 1:n_terms
            if i == j
                result *= "●"
            else
                strength = adjacency[i, j]
                result *= strength == 0 ? "  " :
                          strength == 1 ? "─ " :
                          strength == 2 ? "═ " : "▓ "
            end
        end
        result *= "\n"
    end

    # Add legend
    result *= "\nLegend:\n"
    result *= "● : Term\n"
    result *= "─ : 1 shared feature\n"
    result *= "═ : 2 shared features\n"
    result *= "▓ : 3+ shared features\n"

    return result
end

"""
    visualize_formula_complexity(formula::TwoLevelDNFFormula) -> String

Creates a detailed ASCII visualization of the formula's complexity metrics,
including term sizes, feature distributions, and interaction patterns.

# Example
```julia
println(visualize_formula_complexity(formula))
```
"""
function visualize_formula_complexity(formula::TwoLevelDNFFormula)
    stats = get_formula_statistics(formula)

    result = "Formula Complexity Analysis\n"
    result *= "════════════════════════\n\n"

    # Complexity metrics
    result *= "Complexity Metrics:\n"
    result *= "├── Terms: $(stats.num_terms)\n"
    result *= "├── Features: $(stats.num_features)\n"
    result *= "├── Total Atoms: $(stats.total_atoms)\n"
    result *= "├── Avg Term Size: $(round(stats.avg_term_complexity, digits=2))\n"
    result *= "└── Size Range: [$(stats.min_term_complexity), $(stats.max_term_complexity)]\n\n"

    # Term size distribution
    result *= "Term Size Distribution:\n"
    max_count = maximum(values(stats.complexity_distribution))
    for (size, count) in sort(collect(stats.complexity_distribution))
        bar = "█"^round(Int, (count / max_count) * 20)
        result *= "Size $size: $bar ($count terms)\n"
    end

    # Feature usage
    result *= "\nFeature Usage:\n"
    max_usage = maximum(values(stats.feature_usage))
    for (feature, usage) in sort(collect(stats.feature_usage))
        bar = "█"^round(Int, (usage / max_usage) * 20)
        result *= "F$feature: $bar ($usage atoms)\n"
    end

    return result
end
#################################################################################
################## MAN Two Level DNF Formula ###################################
#=
function manTwoLevelDNFFormula(use_ascii_art::Bool = false)
    # ANSI escape codes for blinking and colors
    BLINK = "\e[5m"
    BOLD = "\e[1m"
    RESET = "\e[0m"
    BLUE = "\e[34m"
    CYAN = "\e[36m"
    MAGENTA = "\e[35m"
    GREEN = "\e[32m"

    # Helper function for centered text
    center(text, width) = lpad(text, div(width + length(text), 2)) * " "^(width - div(width + length(text), 2))

    if use_ascii_art
        # Title banner - each line printed separately for proper formatting
        println(BLINK, BLUE, "╔════════════════════════════════════════════════════════════════════════════╗")
        println(BLUE, "║                                                                            ║")
        println(BLUE, "║   ████████╗██╗    ██╗ ██████╗ ██╗     ███████╗██╗   ██╗███████╗██╗         ║")
        println(BLUE, "║   ╚══██╔══╝██║    ██║██╔═══██╗██║     ██╔════╝██║   ██║██╔════╝██║         ║")
        println(BLUE, "║      ██║   ██║ █╗ ██║██║   ██║██║     █████╗  ██║   ██║█████╗  ██║         ║")
        println(BLUE, "║      ██║   ██║███╗██║██║   ██║██║     ██╔══╝  ╚██╗ ██╔╝██╔══╝  ██║         ║")
        println(BLUE, "║      ██║   ╚███╔███╔╝╚██████╔╝███████╗███████╗ ╚████╔╝ ███████╗███████╗    ║")
        println(BLUE, "║      ╚═╝    ╚══╝╚══╝  ╚═════╝ ╚══════╝╚══════╝  ╚═══╝  ╚══════╝╚══════╝    ║")
        println(BLUE, "║                                                                            ║")
        println(CYAN, "║                    DNF Formula Documentation Manual                        ║")
        println(BLUE, "║                                                                            ║")
        println(BLUE, "╚════════════════════════════════════════════════════════════════════════════╝", RESET)
    else
        println(BOLD, BLUE, repeat("=", 80))
        println(CYAN, center("TwoLevelDNFFormula Manual", 80))
        println(BLUE, repeat("=", 80), RESET)
    end

    println()
    # Structure Section - each line printed separately
    println(MAGENTA, "┌─────────────────── STRUCTURE ────────────────────┐")
    println(MAGENTA, "│                                                  │")
    println(MAGENTA, "│   TwoLevelDNFFormula                             │")
    println(MAGENTA, "│   ├── combinations: Vector{TritVector}           │")
    println(MAGENTA, "│   ├── num_atoms: Int                             │")
    println(MAGENTA, "│   ├── thresholds_by_feature: Dict{Int,Vector}    │")
    println(MAGENTA, "│   └── atoms_by_feature: Dict{Int,Vector}         │")
    println(MAGENTA, "│                                                  │")
    println(MAGENTA, "└──────────────────────────────────────────────────┘", RESET)
    
    println()
    # Analysis Functions - each line printed separately
    println(BLUE, "┌────────── ANALYSIS FUNCTIONS ───────────────┐")
    println(BLUE, "│                                             │")
    println(BLUE, "│   Feature Interactions:                     │")
    println(BLUE, "│   F1 ⟷ F2 ⟷ F3                              │")
    println(BLUE, "│   ↕    ↕    ↕                               │")
    println(BLUE, "│   F4 ⟷ F5 ⟷ F6                              │")
    println(BLUE, "│                                             │")
    println(BLUE, "│   Feature Importance:                       │")
    println(BLINK, "│   F1: ████████ 80%                          │")
    println(BLINK, "│   F2: ██████   60%                          │")
    println(BLINK, "│   F3: ████     40%                          │")
    println(BLUE, "│                                             │")
    println(BLUE, "└─────────────────────────────────────────────┘", RESET)

    println()
    # Evaluation Process - each line printed separately
    println(CYAN, "┌─────────── EVALUATION PROCESS ─────────────────┐")
    println(CYAN, "│                                                │")
    println(CYAN, "│   Input ➜ Evaluate Atoms ➜ Process ➜ Output    │")
    println(CYAN, "│     ↓           ↓            ↓        ↑        │")
    println(CYAN, "│   Parse       Check        Combine    │        │")
    println(CYAN, "│     ↓           ↓            ↓        │        │")
    println(CYAN, "│   Validate    Apply     Optimize ─────┘        │")
    println(CYAN, "│                                                │")
    println(CYAN, "└────────────────────────────────────────────────┘", RESET)

    println()
    # Visualization Functions - each line printed separately
    println(GREEN, "┌────────── VISUALIZATION FUNCTIONS ─────────┐")
    println(GREEN, "│                                            │")
    println(GREEN, "│   Structure View:                          │")
    println(BLINK, "│   ┌─────┐                                  │")
    println(BLINK, "│   │ F1  ├─────┐                            │")
    println(BLINK, "│   └─────┘     │                            │")
    println(BLINK, "│   ┌─────┐   ┌─┴───┐                        │")
    println(BLINK, "│   │ F2  ├───┤ AND ├───┐                    │")
    println(BLINK, "│   └─────┘   └─────┘   │                    │")
    println(BLINK, "│                     ┌─┴───┐                │")
    println(BLINK, "│                     │ OR  │                │")
    println(BLINK, "│                     └─────┘                │")
    println(GREEN, "│                                            │")
    println(GREEN, "└────────────────────────────────────────────┘", RESET)

    println()
    # Conversion Methods - each line printed separately
    println(MAGENTA, "┌──────────── CONVERSION METHODS ────────────┐")
    println(MAGENTA, "│                                            │")
    println(MAGENTA, "│    TwoLevelDNFFormula ⟺ SoleLogics.DNF     │")
    println(MAGENTA, "│                                            │")
    println(MAGENTA, "│    ┌──────────┐        ┌──────────┐        │")
    println(MAGENTA, "│    │ TwoLevel │ ─────⟶ │  Logic   │        │")
    println(MAGENTA, "│    │   DNF    │ ⟵───── │   DNF    │        │")
    println(MAGENTA, "│    └──────────┘        └──────────┘        │")
    println(MAGENTA, "│                                            │")
    println(MAGENTA, "└────────────────────────────────────────────┘", RESET)

    println()
    # Footer - each line printed separately
    println(BLUE, "╔════════════════════════════════════════════════════════════════════════════╗")
    println(BLUE, "║                                                                            ║")
    println(CYAN, "║                     End of Documentation Manual                            ║")
    println(BLUE, "║                                                                            ║")
    println(BLUE, "╚════════════════════════════════════════════════════════════════════════════╝", RESET)
end
=#
################################################################################
function SoleLogics.atoms(f::TwoLevelDNFFormula)
    Atom.(conditions(f))
end

nterms(f::TwoLevelDNFFormula) = length(eachcombination(f))
function natomsperterm(f::TwoLevelDNFFormula)
    sum = 0
    for combination in eachcombination(f)
        for i in 1:length(combination)
            if combination[i] != -1  # Count non-don't-care trits
                sum += 1
            end
        end
    end
    return sum / nterms(f)
end

function conditions(f::TwoLevelDNFFormula)
    return collect(
        Iterators.flatten(
            map(
                ((i_feature, atom_list),) -> map(
                    ((threshold, is_less_than),) -> ScalarCondition(
                        VariableValue(i_feature),
                        is_less_than ? (>) : (<=),
                        threshold,
                    ),
                    atom_list,
                ),
                sort(collect(pairs(eachatomsbyfeature(f)))),
            ),
        ),
    )
end

function Base.show(io::IO, f::TwoLevelDNFFormula)
    print(io, "TwoLevelDNFFormula($(nterms(f)) combinations)")
end

function Base.show(io::IO, ::MIME"text/plain", f::TwoLevelDNFFormula)
    stampa_dnf(io, f)
end

function generate_disjunct(
    combination::TritVector,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    feature_conditions = Dict{Int,Dict{Bool,Float64}}()
    current_atom_index = 1

    # Process each feature and its atoms
    for (feature, atoms) in atoms_by_feature
        for (threshold, _) in atoms
            # Make sure we don't exceed TritVector length
            if current_atom_index <= length(combination)
                trit_value = combination[current_atom_index]

                # Only process if not a don't-care value (-1)
                if trit_value != -1
                    # Initialize nested dict if needed
                    if !haskey(feature_conditions, feature)
                        feature_conditions[feature] = Dict{Bool,Float64}()
                    end

                    # trit_value: 1 => "< threshold" (is_greater_equal = false)
                    # trit_value: 0 => "≥ threshold" (is_greater_equal = true)
                    is_greater_equal = trit_value == 0

                    # Keep only the most restrictive threshold for each operation:
                    # For "<" (is_greater_equal = false) we want the smallest threshold
                    # For "≥" (is_greater_equal = true) we want the largest threshold
                    if !haskey(feature_conditions[feature], is_greater_equal) ||
                       (!is_greater_equal && threshold < feature_conditions[feature][is_greater_equal]) ||
                       (is_greater_equal && threshold > feature_conditions[feature][is_greater_equal])
                        feature_conditions[feature][is_greater_equal] = threshold
                    end
                end
            end
            current_atom_index += 1
        end
    end

    # Convert the conditions to Atoms
    atoms = Vector{Atom}()
    for (feature, conditions) in feature_conditions
        for (is_greater_equal, threshold) in conditions
            mc = ScalarMetaCondition(VariableValue(feature), is_greater_equal ? (≥) : (<))
            condition = ScalarCondition(mc, threshold)
            push!(atoms, Atom(condition))
        end
    end

    return isempty(atoms) ? ⊤ : LeftmostConjunctiveForm(atoms)
end



import Base: convert
using SoleLogics

function Base.convert(::Type{SoleLogics.DNF}, f::TwoLevelDNFFormula)
    conjuncts = [
        generate_disjunct(comb, nuberofatoms(f), eachthresholdsbyfeature(f), eachatomsbyfeature(f))
        for comb in eachcombination(f)
    ]
    filter!(!istop, conjuncts)
    return LeftmostDisjunctiveForm{LeftmostConjunctiveForm{Atom}}(conjuncts, true)
end

function Base.convert(::Type{TwoLevelDNFFormula}, f::SoleLogics.DNF)
    atoms = unique(SoleLogics.atoms(f))
    conds = SoleLogics.value.(atoms)
    num_atoms = length(atoms)
    # TODO:
    thresholds_by_feature = nothing
    atoms_by_feature = nothing
    combinations = nothing
    return TwoLevelDNFFormula(
        combinations,
        num_atoms,
        thresholds_by_feature,
        atoms_by_feature,
    )
end


"""
Evaluates the given `TwoLevelDNFFormula` by recursively evaluating its AND formulas and returning `true` if any of them evaluate to `true` for the given `assignment`.

This function assumes that the `TwoLevelDNFFormula` represents a disjunction of AND formulas, where each AND formula is represented by a `BitVector` in the `combinations` field.
"""
function evaluate(f::TwoLevelDNFFormula, assignment)
    for combination in eachcombination(f)
        disjunct = generate_disjunct(
            combination,
            nuberofatoms(f),
            eachthresholdsbyfeature(f),
            eachatomsbyfeature(f),
        )
        if evaluate(disjunct, assignment)
            return true
        end
    end
    return false
end

"""
Evaluates the given `Atom` by checking the condition associated with it and the values in the provided `assignment`.

If the condition is a `ScalarCondition`, it checks if the value for the corresponding feature is present in the `assignment`. If so, it applies the specified operator to the feature value and the threshold, and returns the result. If the feature is not present in the `assignment`, it raises an error.

If the condition is not a `ScalarCondition`, it raises an error indicating that the condition type is not supported.
"""
function evaluate(atom::Atom, assignment)
    condition = atom.value
    if condition isa ScalarCondition
        feature = condition.metacond.feature.i_variable
        threshold = condition.threshold
        operator = condition.metacond.test_operator

        if haskey(assignment, feature)
            return operator(assignment[feature], threshold)
        else
            error("L'assegnazione non contiene un valore per la feature $feature")
        end
    else
        error("Tipo di condizione non supportato: $(typeof(condition))")
    end
end

"""
Evaluates the given `ScalarCondition` by checking the value for the corresponding feature in the provided `assignment`. If the feature is present in the `assignment`, it applies the specified operator to the feature value and the threshold, and returns the result. If the feature is not present in the `assignment`, it raises an error.

This function is a helper function used within the `evaluate` method for `Atom` objects.
"""
function evaluate(condition::ScalarCondition, assignment)
    feature = condition.metacond.feature.i_variable
    threshold = condition.threshold
    operator = condition.metacond.test_operator

    if haskey(assignment, feature)
        return operator(assignment[feature], threshold)
    else
        error("L'assegnazione non contiene un valore per la feature $feature")
    end
end

"""
Evaluates the given `SoleLogics.NamedConnective{:∧}` (conjunction) by recursively evaluating each of its arguments and returning `true` if all of the arguments evaluate to `true`.

This function is used as part of the overall evaluation of logical expressions represented by the `SyntaxBranch` and `Atom` types.
"""
function evaluate(conjunction::SoleLogics.NamedConnective{:∧}, assignment)
    return all(evaluate(arg, assignment) for arg in conjunction.args)
end

"""
Evaluates the given `SoleLogics.NamedConnective{:∨}` (disjunction) by recursively evaluating each of its arguments and returning `true` if any of the arguments evaluate to `true`.

This function is used as part of the overall evaluation of logical expressions represented by the `SyntaxBranch` and `Atom` types.
"""
function evaluate(disjunction::SoleLogics.NamedConnective{:∨}, assignment)
    return any(evaluate(arg, assignment) for arg in disjunction.args)
end

"""
Evaluates the given `SoleLogics.NamedConnective{:¬}` (negation) by recursively evaluating its single argument and returning the opposite boolean value.

This function is used as part of the overall evaluation of logical expressions represented by the `SyntaxBranch` and `Atom` types.
"""
function evaluate(negation::SoleLogics.NamedConnective{:¬}, assignment)
    return !evaluate(negation.args[1], assignment)
end


"""
Evaluate the given `SyntaxBranch` by recursively evaluating its left and right children, and returning the result of the AND operation on the two child evaluations.

This function assumes that the `SyntaxBranch` represents a binary tree structure with an AND operation between the left and right child nodes.
"""
function evaluate(branch::SyntaxBranch, assignment)
    left_result = evaluate(branch.children[1], assignment)
    right_result = evaluate(branch.children[2], assignment)
    # Assumiamo che l'operatore sia AND (∧)
    return left_result && right_result
end


"""
Prints a human-readable representation of a `TwoLevelDNFFormula` to the specified output stream `io`.

The representation includes the number of combinations in the formula, and then prints each OR-AND combination using the `stampa_disjunct` function.
"""
function stampa_dnf(
    io::IO,
    formula::TwoLevelDNFFormula,
    max_combinations::Int=nterms(formula),
)
    println(io, "TwoLevelDNFFormula with $(nterms(formula)) combinations:")
    for (i, combination) in enumerate(eachcombination(formula)[1:min(max_combinations, end)])
        disjunct = generate_disjunct(
            combination,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
        )
        print(io, "  OR[$i]: ")
        stampa_disjunct(io, disjunct)
        println(io)

        if i == max_combinations && nterms(formula) > max_combinations
            println(
                io,
                "  ... (other $(nterms(formula) - max_combinations) combinations not shown)",
            )
        end
    end
end

function stampa_disjunct(io::IO, formula::Formula)
    if formula isa Atom
        print(io, formula)
    elseif formula isa SoleLogics.NamedConnective{:∧}
        print(io, "(")
        for (i, atom) in enumerate(formula.args)
            i > 1 && print(io, " ∧ ")
            print(io, atom)
        end
        print(io, ")")
    else
        print(io, formula)
    end
end