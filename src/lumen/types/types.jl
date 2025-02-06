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

function TwoLevelDNFFormula(combinations::AbstractVector, atoms::AbstractVector{<:Atom} = [])
    if isempty(f)
        # TODO @Marco
        num_atoms = length(atoms)
        thresholds_by_feature = Dict{Int,Vector{Float64}}()
        atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
        # TODO check that atoms are all < (or >, >=, <=?)
        conds = SoleLogics.value.(atoms)
        feats = unique(SoleLogics.feature.(atoms))
        # TODO check that there are no dual conditions
        for cond in conds
            push!(thresholds_by_feature[findfirst(x->x==feature(cond), feats)], TODO)
            push!(atoms_by_feature[findfirst(x->x==feature(cond), feats)], TODO)
        end
        TwoLevelDNFFormula(combinations, num_atoms, thresholds_by_feature, atoms_by_feature)
        # return TwoLevelDNFFormula(BitVector[], 0, Dict{Int,Vector{Float64}}(), Dict{Int,Vector{Tuple{Float64,Bool}}}(), Vector{Int}[])
    else
        error("TODO implement constructor.")
    end
end

eachcombination(f::TwoLevelDNFFormula) = f.combinations
#eachmaskedcombination(f::TwoLevelDNFFormula) = isempty(f.prime_mask) ? eachcombination(f) : f.prime_mask
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
                sort(collect(pairs(f.atoms_by_feature))),
            ),
        ),
    )
end

function SoleLogics.atoms(f::TwoLevelDNFFormula)
    Atom.(conditions(f))
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
        generate_disjunct(comb, f.num_atoms, f.thresholds_by_feature, f.atoms_by_feature)
        for comb in eachcombination(f)
    ]
    filter!(!istop, conjuncts)
    return LeftmostDisjunctiveForm{LeftmostConjunctiveForm{Atom}}(conjuncts, true)
end

function Base.convert(::Type{TwoLevelDNFFormula}, f::SoleLogics.Formula)
    if SoleLogics.istop(f)
        TwoLevelDNFFormula([])
    else
        f = SoleLogics.dnf(f)
        atoms = unique(SoleLogics.atoms(f))
        conds = SoleLogics.value.(atoms)
        combinations = begin
            if f isa SoleLogics.LeftmostDisjunctiveForm{<:Union{Atom,Literal,LeftmostConjunctiveForm{<:Union{Atom,Literal}}}}
                disjs = SoleLogics.disjuncts(f)
                # 
                function disjunct_to_combination!(combination, disj::Atom, conds)
                    combination[findall(==(disj), conds)] .= 1
                    if SoleLogics.hasdual(disj)
                        combination[findall(==(SoleLogics.dual(disj)), conds)] .= 0
                    end
                    combination
                end
                function disjunct_to_combination!(combination, disj::Literal, conds)
                    if SoleLogics.ispos(disj)
                        disjunct_to_combination!(combination, atom(disj), conds)
                    else
                        combination[findall(==(disj), conds)] .= 0
                        if SoleLogics.hasdual(disj)
                            combination[findall(==(SoleLogics.dual(disj)), conds)] .= 1
                        end
                        combination
                    end
                end
                function disjunct_to_combination(disj, conds)
                    combination = fill(-1, length(conds))
                    disjunct_to_combination!(combination, disj, conds)
                end
                disjunct_to_combination(disj, conds) = error("Cannot convert disjunct of type $(typeof(disj)) to combination.")
                function disjunct_to_combination(disj::LeftmostConjunctiveForm, conds)
                    combination = fill(-1, length(conds))
                    for conj in SoleLogics.conjuncts(disj)
                        disjunct_to_combination!(combination, conj, conds)
                    end
                    combination
                end
                combinations = [disjunct_to_combination(disj, conds) for disj in disjs]
            else
                error("Could not convert formula of type $(typeof(f)) to TwoLevelDNFFormula.")
            end
        end
        return TwoLevelDNFFormula(combinations, atoms, )
    end
end


"""
Evaluates the given `TwoLevelDNFFormula` by recursively evaluating its AND formulas and returning `true` if any of them evaluate to `true` for the given `assignment`.

This function assumes that the `TwoLevelDNFFormula` represents a disjunction of AND formulas, where each AND formula is represented by a `BitVector` in the `combinations` field.
"""
function evaluate(f::TwoLevelDNFFormula, assignment)
    for combination in eachcombination(f)
        disjunct = generate_disjunct(
            combination,
            f.num_atoms,
            f.thresholds_by_feature,
            f.atoms_by_feature,
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
    max_combinations::Int = nterms(formula),
)
    println(io, "TwoLevelDNFFormula with $(nterms(formula)) combinations:")
    for (i, combination) in enumerate(eachcombination(formula)[1:min(max_combinations, end)])
        disjunct = generate_disjunct(
            combination,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
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