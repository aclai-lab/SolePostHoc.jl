
"""
    generate_binary_row(num::Int, length::Int)

Generates a binary representation of the given integer `num` as a vector of booleans with the specified `length`.

Args:
    num::Int: The integer to convert to a binary representation.
    length::Int: The desired length of the binary representation.

Returns:
    Vector{Bool}: A vector of booleans representing the binary digits of `num`, padded to the specified `length`.
"""
function generate_binary_row(num::Int, length::Int)
    return [(num >> j) & 1 == 1 for j = 0:(length-1)]
end


"""
    print_progress_bar(progress::BigFloat, width::Int=50)

Prints a progress bar to the console.

Args:
    progress::BigFloat: The current progress as a value between 0 and 1.
    width::Int: The width of the progress bar in characters (default is 50).

Returns:
    Nothing. Prints the progress bar to the console.
"""
function print_progress_bar(progress::BigFloat, width::Int = 50)
    progress_float = Float64(progress)
    filled = Int(round(progress_float * width))
    print(
        "\r[",
        "⚪"^filled,
        "⚫"^(width - filled),
        "] ",
        lpad(round(progress_float * 100, digits = 1), 5),
        "%",
    )
    flush(stdout)
end


"""
    dict_to_bitvector(d::Dict{Any, Vector{BigInt}}, np::Int)

Convert a dictionary of `Vector{BigInt}` values to a dictionary of `BitVector` values.

Args:
    d::Dict{Any, Vector{BigInt}}: The input dictionary to be converted.
    np::Int: The number of bits to pad each `BigInt` value to.

Returns:
    Dict{Any, Vector{BitVector}}: A new dictionary with the same keys as `d`, but with the values converted to `BitVector`.
"""
function dict_to_bitvector(d::Dict{Any,Vector{BigInt}}, np::Int)
    # Use dictionary comprehension and digits2 for direct conversion
    Dict(k => [digits(v, base = 2, pad = np) |> BitVector for v in vals] for (k, vals) in d)
end

function dict_to_tritvector(d::Dict{Any,Vector{BigInt}}, np::Int)
    # Similar to BitVector version but creates TritVectors
    Dict(k => [
        let bits = digits(v, base = 2, pad = np)
            trit = TritVector(np)
            for (i, bit) in enumerate(bits)
                trit[i] = bit
            end
            trit
        end for v in vals
    ] for (k, vals) in d)
end

# """
#     create_table(labels::Vector{String}, features::Matrix{Float64})

# Create a DataFrame from the given feature matrix and labels.

# Args:
#     labels::Vector{String}: A vector of labels for each row in the feature matrix.
#     features::Matrix{Float64}: A matrix of feature values.

# Returns:
#     DataFrame: A DataFrame with the feature values and labels.
# """
# function create_table(labels::Vector{String}, features::Matrix{Float64})
#     df = DataFrame(features, [:SepalLength, :SepalWidth, :PetalLength, :PetalWidth])
#     df.Label = labels
#     return df
# end


"""
    process_alphabet(atom_prop::Vector{<:Atom{<:ScalarCondition{Float64,<:VariableValue,<:ScalarMetaCondition{<:VariableValue,typeof(<)}}}}, n::Int)

Constructs a `MultivariateScalarAlphabet` from a vector of `Atom` objects with scalar conditions.

The function processes the `Atom` objects to extract the thresholds and test operators for each feature, and then creates a `UnivariateScalarAlphabet` for each feature. The resulting `UnivariateScalarAlphabet` objects are then combined into a `MultivariateScalarAlphabet`.

If a feature is not present in the `Atom` objects, a default threshold of 42.0 is used, and the test operator is set to the first `Atom`'s test operator.
"""
function process_alphabet(
    # atom_prop::Vector{
    #     <:Atom{
    #         <:ScalarCondition{
    #             <:AbstractFloat,
    #             <:VariableNamedValue,
    #             <:ScalarMetaCondition{<:VariableNamedValue,typeof(<)},
    #         },
    #     },
    # },
    atom_prop::Vector{<:Atom{<:ScalarCondition}},
    n::Int,
)
    thresholds_by_feature = Dict{Int64,Vector{Float64}}()
    test_operator_by_feature = Dict{Int64,Any}()

    for atom in atom_prop
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        test_op = atom.value.metacond.test_operator

        if !haskey(thresholds_by_feature, feat)
            thresholds_by_feature[feat] = Float64[]
            test_operator_by_feature[feat] = test_op
        end
        push!(thresholds_by_feature[feat], threshold)
    end

    for feat = 1:n
        if !haskey(thresholds_by_feature, feat)
            thresholds_by_feature[feat] = [42.0]
            test_operator_by_feature[feat] = atom_prop[1].value.metacond.test_operator
        end
    end

    univariate_alphabets = []

    for (feat, thresholds) in sort(collect(thresholds_by_feature))
        test_op = test_operator_by_feature[feat]
        mc = ScalarMetaCondition(VariableValue(feat), test_op)

        min_threshold = minimum(thresholds)
        new_min = prevfloat(min_threshold)
        diff = abs(min_threshold - new_min)
        if diff >= 0.1
            @warn "Grande differenza nella soglia rilevata"
        end
        push!(thresholds, new_min)
        sort!(thresholds)

        alphabet = UnivariateScalarAlphabet((mc, thresholds))
        push!(univariate_alphabets, alphabet)
    end
    return MultivariateScalarAlphabet{ScalarCondition}(univariate_alphabets)
end


"""
    extract_atoms_bfs_order(tree)

Extracts atoms (decisional conditions) from a decisuion tree, using
a Breadth-First Search (BFS) approach.

# Functioning
The function implements the BFS algorithm using a FIFO queue:
1. Initialize a queue with the root of the tree
2. For each extracted node from the queue:
   - If the node is a Branch (decision node), extract the atom from the condition (`antecedent`)
   - Add the node's children to the queue (posconsequent, negconsequent)
   - If it's a CostantModel leaf, ignore it (it has no atoms)
3. Continue until the queue is empty

# Args
- `tree`: The root of the decision tree (Branch{String})

# Returns
- A vector iof atoms in BFS order (level by level, left to right)

# Example
```julia
atoms_bfs = extract_atoms_bfs_order(tree)
# Returns: [V2<1.5, V1<1.5, V3<3.5, V1<2.5, V4<1.5]
```

# Notes
- Uses generic types (Any[]) to handle both Branch nodes and ConstantModel leaves
- The BFS orders ensures that the nodes from the same level appear consequtively
- Only Branch nodes contribute to the result (contain atoms)
"""
function extract_atoms_bfs_order(tree)
    bfs_atoms = []  # Use generic [], not typified
    queue = Any[tree]  # Uses Any[] to contain both Branch and ConstantModel

    while !isempty(queue)
        current = popfirst!(queue)

        if typeof(current).name.name == :Branch
            push!(bfs_atoms, current.antecedent)
            push!(queue, current.posconsequent)
            push!(queue, current.negconsequent)
        end
        # If it's a ConstantModel, it doesn't contain antecedent so we do not add bfs_atoms
    end

    return bfs_atoms
end


"""
    take_first_percentage(atoms, c)

Takes only the first c% atoms from a vectore, where c is a value between 0 and 1.

# Args
- `atoms`: An atoms vector
- `c`: Percentage as a decimal value (0.0 = 0%, 1.0 = 100%)

# Returns
- Subvector containing only the first c% atoms

# Example
```julia
atoms = [V2<1.5, V1<1.5, V3<3.5, V1<2.5]
take_first_percentage(atoms, 0.5)  # Returns [V2<1.5, V1<1.5]
take_first_percentage(atoms, 0.25) # Returns [V2<1.5]
take_first_percentage(atoms, 1.0)  # Returns the whole vector
```
"""
function take_first_percentage(atoms, c)
    n_total = length(atoms)
    n_to_take = Int(ceil(n_total * c))  # Rounds up

    return atoms[1:min(n_to_take, n_total)]
end

# """
# Extracts the rules from a `DecisionForest` and organizes them by label.

# For each tree in the forest, this function extracts the rules, groups them by the label (consequent) of the rule, and prints the antecedents for each label.

# This function is primarily used for debugging and inspecting the rules learned by the decision forest model.
# """
# function extract_rules_by_label(f::DecisionForest)
#     for (i, tree) in enumerate(f.trees)
#         println("___________ ALBERO #$i ___________")
#         rules_by_label = Dict{String,Vector{Any}}()
#         rules = listrules(tree)
#         for rule in rules
#             label = string(consequent(rule))
#             antecedent_rule = antecedent(rule)
#             push!(get!(Vector{Any}, rules_by_label, label), antecedent_rule)
#         end
#         for (label, antecedents) in rules_by_label
#             println("Etichetta: ", label)
#             for antecedent in antecedents
#                 println(antecedent)
#             end
#         end
#     end
# end
