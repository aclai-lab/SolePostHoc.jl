using Random

"""
Prints a summary of the estimated time for a computation.

Args:
    total_estimated_time::BigFloat: The total estimated time for the computation, in seconds.

This function prints a message indicating the estimated time for the computation, formatted as seconds, minutes, or hours as appropriate.
"""
function print_estimated_time(total_estimated_time::BigFloat)
    if total_estimated_time < 60
        println("Estimated time: approx $(round(total_estimated_time, digits=2)) seconda")
    elseif total_estimated_time < 3600
        println("Estimated time: approx $(round(total_estimated_time / 60, digits=2)) minut")
    else
        println("Estimated time: approx $(round(total_estimated_time / 3600, digits=2)) hours")
    end
end


"""
    print_results_summary(real_time_per_combination::BigFloat, contradictions::Int64, num_combinations::BigInt, label_count::Dict{Any,Int64})

Prints a summary of the results, including the total number of combinations, the number of valid and invalid combinations, and the distribution of labels.

Args:
    real_time_per_combination::BigFloat: The precision of the time per combination.
    contradictions::Int64: The number of logical contradictions found.
    num_combinations::BigInt: The total number of combinations.
    label_count::Dict{Any,Int64}: A dictionary mapping labels to their counts.
"""
function print_results_summary(
    real_time_per_combination::BigFloat,
    contradictions::Int64,
    num_combinations::BigInt,
    label_count::Dict{Any,Int64},
)
    valid_combinations = num_combinations - BigInt(contradictions)
    valid_percentage = (Float64(valid_combinations) / Float64(num_combinations)) * 100
    invalid_percentage = (Float64(contradictions) / Float64(num_combinations)) * 100

    println()
    println("Time precision: $real_time_per_combination")
    println("Combination generation completed.")
    println("Total combinations: $num_combinations")
    println("Total logical contradictions found: $contradictions out of $num_combinations")
    println("Total valid combinations: $valid_combinations")

    println("Percentage of valid combinations: $(round(valid_percentage, digits=2))%")
    println("Percentage of invalid combinations: $(round(invalid_percentage, digits=2))%")

    println("\nLabel distribution:")
    for (label, count) in sort(collect(label_count), by = x -> x[2], rev = true)
        percentage = (Float64(count) / Float64(valid_combinations)) * 100
        println("$label: $count ($(round(percentage, digits=2))%)")
    end
end


"""
    print_detailed_results(results::Dict{Any,Vector{BigInt}})

Prints detailed results of the computation, including the number of combinations for each result and the first 10 combinations (or all combinations if there are 10 or fewer).

Args:
    results::Dict{Any,Vector{BigInt}}: A dictionary mapping results to the corresponding combinations.
"""
function print_detailed_results(results::Dict{Any,Vector{BigInt}})
    println("\nDetailed results:")
    for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
        println("[$result] ($(length(combinations)) combinations):")
        if length(combinations) > 10
            println("  First 10 combinations: $(combinations[1:10])")
            println("  ... ($(length(combinations)-10) more combinations)")
        else
            println("  Combinations: $combinations")
        end
    end
end

"""
    my_syntaxstring(formula::TwoLevelDNFFormula, orizontal::Float64 = 1.0)

Given a `TwoLevelDNFFormula` which uses TritVectors (`1` for <, `0` for ≥, `-1` for don't-care),
this function returns an array of conjunction-strings. Each entry in the returned array is one
term of the DNF (the conjunction of conditions for that combination).
"""
function my_syntaxstring(
    formula::TwoLevelDNFFormula,
    orizontal::Float64 = 1.0
    )
    # Number of features (horizontally) we want to consider
    num_orizontal = floor(Int, eachthresholdsbyfeature(formula).count * orizontal)
    
    # We'll accumulate each distinct conjunction-string in a Set to avoid duplicates
    result = Set{String}()

    # For each TritVector (one per conjunction)
    for combination in eachcombination(formula)
        feature_conditions = Dict{Int,Dict{String,Float64}}()
        current_atom_index = 1
        
        # For each feature and its sorted list of (threshold, bool), in ascending threshold order
        for (feature, atoms) in eachatomsbyfeature(formula)
            if feature <= num_orizontal
                for (threshold, _) in atoms
                    # Make sure we do not exceed the length of the TritVector
                    if current_atom_index <= length(combination)
                        trit_value = combination[current_atom_index]
                        
                        # Only proceed if not -1 (i.e. not a don't-care)
                        if trit_value != -1
                            # Initialize nested dict if needed
                            if !haskey(feature_conditions, feature)
                                feature_conditions[feature] = Dict{String,Float64}()
                            end
                            
                            # 1 => "< threshold"
                            # 0 => "≥ threshold"
                            op = (trit_value == 1) ? "<" : "≥"
                            
                            # We keep only the "most restrictive" threshold for each op:
                            #   for "<"  we want the *smallest* threshold
                            #   for "≥" we want the *largest* threshold
                            if !haskey(feature_conditions[feature], op) ||
                               (op == "<" && threshold < feature_conditions[feature][op]) ||
                               (op == "≥" && threshold > feature_conditions[feature][op])
                                feature_conditions[feature][op] = threshold
                            end
                        end
                    end
                    current_atom_index += 1
                end
            end
        end

        # Build the string for this conjunction
        current_atoms = String[]
        for (feature, conditions) in feature_conditions
            for (op, value) in conditions
                push!(current_atoms, "V$(feature) $(op) $(value)")
            end
        end

        if !isempty(current_atoms)
            # Sort the atoms for a stable display, then join them with " ∧ "
            push!(result, join(sort(current_atoms), " ∧ "))
        end
    end

    # Convert the Set back to a Vector
    return collect(result)
end


"""
    convert_DNF_formula(formula::TwoLevelDNFFormula, outcome, orizontal::Float64 = 1.0)

Constructs a `Rule` from the DNF formula, parsing the resulting string into a `SoleLogics.Formula`.
The `outcome` is any label you want to attach to the rule (e.g., a class label).
"""
function convert_DNF_formula(
    formula::TwoLevelDNFFormula,
    outcome,
    orizontal::Float64 = 1.0,
)
    # Build the array of conjunction-strings
    formulas = my_syntaxstring(formula, orizontal)
    
    # Join the conjunctions with " ∨ " to get the full DNF
    result = join(formulas, " ∨ ")
    
    # Parse the DNF string into a SoleLogics formula
    φ = SoleLogics.parseformula(
        result;
        atom_parser = a -> Atom(
            parsecondition(
                SoleData.ScalarCondition, 
                a; 
                featuretype = SoleData.VariableValue,
                featvaltype = Real
            )
        )
    )

    #= 
        We create the Rule using the outcome passed as a parameter
        println("mask:", formula.prime_mask)
        println("combination", formula.combinations)
        print_dnf(stdout, formula)
        println("formulas:", dump(formula)) 
    =#
    #println("comb:", formula.combinations)
    #print_dnf(stdout, formula)
    return Rule(φ, outcome)
end


"""
    leftmost_disjunctive_form_to_string(ldf)

Given a `LeftmostDisjunctiveForm{LeftmostConjunctiveForm}`, returns a string in DNF form.
"""
function leftmost_disjunctive_form_to_string(ldf, c=1.0, vetImportance=[])
    # Helper: extract the list of atoms from either a ConjunctiveForm or a single Atom
    function get_atoms(cf)
        if cf isa LeftmostConjunctiveForm
            return cf.grandchildren
        elseif cf isa Atom
            return [cf]
        else
            error("Unknown type in leftmost_disjunctive_form_to_string: $(typeof(cf))")
        end
    end

    # If c is 1.0 (100%) or vetImportance is empty, use all variables
    include_all_vars = (c == 1.0 || isempty(vetImportance))
    
    # Calculate which variables to include based on importance
    included_vars = Set{Int}()
    if !include_all_vars
        # Determine how many variables to include based on c
        num_vars_to_include = max(1, ceil(Int, length(vetImportance) * c))
        
        # Get the most important variables
        sorted_indices = sortperm(1:length(vetImportance), by=i -> findfirst(==(i), vetImportance))
        for i in 1:num_vars_to_include
            push!(included_vars, sorted_indices[i])
        end
    end

    conjunctive_strings = String[]

    # Each child is one term in the top-level disjunction
    for child in ldf.grandchildren
        these_atoms = get_atoms(child)

        atom_strings = String[]
        for atom in these_atoms
            cond = atom.value
            i_var = cond.metacond.feature.i_variable
            
            # Skip this atom if the variable is not in the included set
            if !include_all_vars && !(i_var in included_vars)
                continue
            end
            
            op_fun = cond.metacond.test_operator
            thr = cond.threshold

            # Turn the operator function into a standard symbol
            op_str = op_fun === (<)  ? "<"  :
                     op_fun === (<=) ? "≤"  :
                     op_fun === (>)  ? ">"  :
                     op_fun === (>=) ? "≥"  :
                     string(op_fun)

            push!(atom_strings, "V$(i_var) $(op_str) $(thr)")
        end

        # Skip empty conjunctions (when all atoms were filtered out)
        if isempty(atom_strings)
            continue
        end
        
        # Join all atoms in this conjunction
        conj_str = "(" * join(atom_strings, " ∧ ") * ")"
        push!(conjunctive_strings, conj_str)
    end
    
    # Skip if all conjunctions were filtered out
    if isempty(conjunctive_strings)
        return "∅"  # Empty set symbol to indicate no conditions remain
    end

    # Finally, join all conjunctions with " ∨ "
    return join(conjunctive_strings, " ∨ ")
end
