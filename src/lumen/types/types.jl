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
    combinations::Vector{BitVector}
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
    prime_mask::Vector{Vector{Int}}
end
```

### Component Hierarchy
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
Binary representation of AND terms using BitVectors:
```ascii
[BitVector("0101"), BitVector("1100"), BitVector("0011"), ...]
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

### 4. Prime Mask Structure
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
│    BitVector Combination   │
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

### Prime Mask Structure
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
Input → Threshold Check → BitVector Match → Result
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
    combinations::Vector{BitVector}
    num_atoms::Int
    thresholds_by_feature::Dict{Int,Vector{Float64}}
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}
    prime_mask::Vector{Vector{Int}} 
end

function TwoLevelDNFFormula(f::AbstractVector)
    if isempty(f)
        # return TwoLevelDNFFormula(BitVector[], 0, Dict{Int,Vector{Float64}}(), Dict{Int,Vector{Tuple{Float64,Bool}}}(), Vector{Int}[])
    else
        error("TODO implement constructor.")
    end
end

eachcombination(f::TwoLevelDNFFormula) = f.combinations
eachmaskedcombination(f::TwoLevelDNFFormula) = isempty(f.prime_mask) ? eachcombination(f) : f.prime_mask
function SoleLogics.atoms(f::TwoLevelDNFFormula)
    Atom.(conditions(f))
end

nterms(f::TwoLevelDNFFormula) = length(eachcombination(f))
function natomsperterm(f::TwoLevelDNFFormula)
    # print(eachcombination(f))
    # sum(c->sum(!=(-1), c), eachcombination(f))/nterms(f) |> print
    # readline()
    sum(c->sum(!=(-1), c), eachmaskedcombination(f))/nterms(f)
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

"""
Provides a custom string representation for the `TwoLevelDNFFormula` struct, displaying the number of combinations it contains.
"""
Base.show(io::IO, f::TwoLevelDNFFormula) =
    print(io, "TwoLevelDNFFormula($(nterms(f)) combinations)")

"""
Prints a human-readable representation of a `TwoLevelDNFFormula` to the specified output stream `io`.

The representation includes the number of combinations in the formula, and then prints each OR-AND combination using the `stampa_dnf` function.
"""
function Base.show(io::IO, ::MIME"text/plain", f::TwoLevelDNFFormula)
    stampa_dnf(io, f)
end


"""
Generate an AND formula based on the provided binary `combination` representing a combination of atoms. 
Calculate conditions for each feature based on the thresholds and atom properties. Construct atoms with corresponding conditions and return the resulting AND formula.
"""
function generate_disjunct(
    combination::BitVector,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    comb, _ =
        process_combination(combination, num_atoms, thresholds_by_feature, atoms_by_feature)
    conditions_by_feature = Dict{Int,Dict{Bool,Float64}}()

    for (feat, values) in comb
        if haskey(atoms_by_feature, feat)
            for (threshold, _) in atoms_by_feature[feat]
                if values[1] < threshold
                    # Condizione "<"
                    if !haskey(conditions_by_feature, feat) ||
                       !haskey(conditions_by_feature[feat], false) ||
                       threshold < conditions_by_feature[feat][false]
                        conditions_by_feature[feat] =
                            get(conditions_by_feature, feat, Dict{Bool,Float64}())
                        conditions_by_feature[feat][false] = threshold
                    end
                else
                    # Condizione "≥"
                    if !haskey(conditions_by_feature, feat) ||
                       !haskey(conditions_by_feature[feat], true) ||
                       threshold > conditions_by_feature[feat][true]
                        conditions_by_feature[feat] =
                            get(conditions_by_feature, feat, Dict{Bool,Float64}())
                        conditions_by_feature[feat][true] = threshold
                    end
                end
            end
        end
    end

    atoms = Vector{Atom}()
    for (feat, conditions) in conditions_by_feature
        for (is_greater_equal, threshold) in conditions
            mc = ScalarMetaCondition(VariableValue(feat), is_greater_equal ? (≥) : (<))
            condition = ScalarCondition(mc, threshold)
            push!(atoms, Atom(condition))
        end
    end

    # return isempty(atoms) ? ⊤ : (length(atoms) == 1 ? first(atoms) : ∧(atoms...))
    return isempty(atoms) ? ⊤ : LeftmostConjunctiveForm(atoms)
end



import Base: convert
using SoleLogics

function Base.convert(::Type{SoleLogics.DNF}, f::TwoLevelDNFFormula)
    conjuncts = [
        generate_disjunct(comb, f.num_atoms, f.thresholds_by_feature, f.atoms_by_feature) for comb in eachcombination(f)
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
        if f isa SoleLogics.LeftmostDisjunctiveForm{<:Union{Atom,Literal,LeftmostConjunctive{<:Union{Atom,Literal}}}}
            disjs = SoleLogics.disjuncts(f)
            # 
            # TODO @Marco
            combinations, num_atoms, thresholds_by_feature, atoms_by_feature, prime_mask = begin
                num_atoms = length(atoms)
                function disjunct_to_combination!(combination, disj::Atom, conds)
                    combination[findall(==(disj), conds)] = 1
                    if SoleLogics.hasdual(disj)
                        combination[findall(==(SoleLogics.dual(disj)), conds)] = 0
                    end
                    combination
                end
                function disjunct_to_combination!(combination, disj::Literal, conds)
                    if SoleLogics.ispos(disj)
                        disjunct_to_combination!(combination, disj(disj), conds)
                    else
                        combination[findall(==(disj), conds)] = 0
                        if SoleLogics.hasdual(disj)
                            combination[findall(==(SoleLogics.dual(disj)), conds)] = 1
                        end
                        combination
                    end
                end
                function disjunct_to_combination(disj, conds)
                    combination = fill(-1, length(conds))
                    disjunct_to_combination!(combination, disj, conds)
                end
                disjunct_to_combination(disj, conds) = error("Cannot convert disjunct of type $(typeof(disj)) to combination.")
                function disjunct_to_combination(disj::LeftmostConjunctive, conds)
                    combination = fill(-1, length(conds))
                    for conj in conjuncts(disj)
                        disjunct_to_combination!(combination, conj, conds)
                    end
                end
                # TODO test!
                combinations = [disjunct_to_combination(disj, conds) for disj in disjs]
                thresholds_by_feature = nothing
                atoms_by_feature = nothing
                prime_mask = nothing
                # 
                combinations, num_atoms, thresholds_by_feature, atoms_by_feature, prime_mask
            end
        else
            error("Could not convert formula of type $(typeof(f)) to TwoLevelDNFFormula.")
        end
        return TwoLevelDNFFormula(
            combinations,
            num_atoms,
            thresholds_by_feature,
            atoms_by_feature,
            prime_mask,
        )
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
                "  ... (altre $(nterms(formula) - max_combinations) combinazioni non mostrate)",
            )
        end
    end
end


"""
Prints the AND formula represented by the given `Formula` to the specified output stream `io`.

If the `formula` is an `Atom`, it is printed as-is. If the `formula` is a conjunction (`SoleLogics.NamedConnective{:∧}`), it is printed with each conjunct separated by `∧`. For any other type of `Formula`, it is printed as-is.
"""
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