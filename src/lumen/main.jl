module Lumen

# using Revise
using Logging
using Dates
using DataStructures
using SoleModels
using SoleModels: DecisionSet
using AbstractTrees
using SoleData
using SoleLogics
using SoleLogics: normalize_formula, dnf
using SoleData: MultivariateScalarAlphabet, UnivariateScalarAlphabet
using DecisionTree: load_data, build_forest, apply_forest
using ModalDecisionTrees
using SoleLogics
using BenchmarkTools, StatProfilerHTML
using Base: intersect
using Base.Threads: @threads
using Base.Threads: Atomic, atomic_add!
using Profile
using ConcurrentCollections
using ProgressMeter
using DataFrames

import DecisionTree as DT

export lumen

"""
    lumen(model, minimization_scheme = :mitespresso;
        start_time = time(),
        vertical = 1.0,
        horizontal = 1.0,
        ott_mode = false,
        controllo = false,
        minimization_kwargs = (;),
        filteralphabetcallback! = identity,
        solemodel = nothing,
        apply_function = SoleModels.apply,
        silent = false,
        return_info = true,
        vetImportance = [],
        kwargs...
    )

Logic-driven Unified Minimal Extractor of Notions (LUMEN): A function that extracts and minimizes
logical rules from a decision tree ensemble model into DNF (Disjunctive Normal Form) formulas.

# Arguments
- `model`: The decision tree ensemble model to analyze
- `minimization_scheme::Symbol`: DNF minimization algorithm to use (default: :mitespresso)

# Keyword arguments
- `start_time`: Start time for performance measurement
- `vertical::Real`: Vertical coverage parameter (α) for rule extraction (0.0 < vertical ≤ 1.0)
- `horizontal::Real`: Horizontal coverage parameter (β) for rule extraction (0.0 < horizontal ≤ 1.0)
- `ott_mode::Bool`: Flag to enable optimized processing mode
- `controllo::Bool`: Flag to enable validation mode for comparing results
- `minimization_kwargs::NamedTuple`: Additional arguments for the minimization algorithm
- `filteralphabetcallback!`: Callback function for alphabet processing (default: identity)
- `solemodel`: Pre-constructed SOLE model (if `nothing`, will be created from `model`)
- `apply_function`: Function used to apply the model (default: SoleModels.apply)
- `silent::Bool`: Flag to suppress output messages
- `return_info::Bool`: Flag to include additional information in the return value (default: true)
- `vetImportance`: Vector for feature importance values
- `kwargs...`: Additional keyword arguments

# Processing Steps
1. Rule Extraction:
   - Extracts logical rules from the decision trees in the model
   - Validates and processes input parameters

2. Alphabet Construction:
   - Extracts atoms from the rules
   - Constructs and processes the logical alphabet
   - Applies custom alphabet filtering if provided

3. Truth Table Generation:
   - Creates truth combinations based on the extracted atoms
   - Handles both standard and optimized processing modes
   - Generates L-class path forest for rule representation

4. DNF Minimization:
   - Applies the specified minimization scheme to simplify the DNF formulas
   - Verifies the correctness of the simplification
   - Reports the reduction in term count

5. Result Validation (when controllo = true):
   - Compares results between different processing modes
   - Validates optimization correctness

# Return
- When `return_info = false`: Returns a `DecisionSet` containing minimized rules
- When `return_info = true`: Returns a tuple with:
  - `DecisionSet` containing minimized rules
  - Named tuple with additional information:
    - `vectPrePostNumber`: Vector of tuples showing (pre-minimization, post-minimization) term counts
    - `unminimized_ds`: `DecisionSet` containing unminimized rules

# Notes
- Parameters `vertical` and `horizontal` must be in range (0.0, 1.0]
- Setting both parameters to 1.0 enforces strong rule extraction
- The function aligns with Algorithm 1 from the reference implementation
- Performance statistics are generated based on processing time and rule reduction

# Example
```julia
model = load_decision_tree_model()
start_time = time()
decision_set, info = lumen(model, :mitespresso)
```

See also
[`AbstractModel`](@ref),
[`DecisionList`](@ref),
[`listrules`](@ref),
[`rulemetrics`](@ref).
"""
function lumen(
    modelJ, # actualy truth_combinations usa model
    minimization_scheme::Symbol=:abc;
    vertical::Real=1.0,
    horizontal::Real=1.0,
    ott_mode::Bool=false,
    controllo::Bool=false,
    start_time=time(),
    minimization_kwargs::NamedTuple=(;),
    filteralphabetcallback=identity,
    solemodel=nothing,
    apply_function=SoleModels.apply,
    silent=false,
    return_info=true, # TODO must default to `false`.
    vetImportance=[],
    merge_negated_elements=false,
    testott=nothing,
    alphabetcontroll=nothing,
    kwargs...,
)
    if vertical <= 0.0 || vertical > 1.0 || horizontal <= 0.0 || horizontal > 1.0
        @warn "Invalid parameters, setting both to 1"
        vertical = horizontal = 1.0
    end
    model = isnothing(solemodel) ? SoleModels.solemodel(modelJ) : solemodel

    if isnothing(apply_function)
        if modelJ isa DT.Ensemble
            apply_function = DT.apply_forest
        else
            apply_function = SoleModels.apply
        end
    end

    is_minimization_scheme_espresso =
        minimization_scheme == :mitespresso ||
        minimization_scheme == :boom ||
        minimization_scheme == :abc # || minimization_scheme == :texasespresso

    # PART 2.a: Starter Ruleset Extraction
    silent || println(
        "\n\n$COLORED_TITLE$TITLE\n PART 2.a STARTER RULESET ESTRACTION \n$TITLE$RESET",
    )

    ruleset = @time begin
        if isensemble(model)
            rs = unique([
                listrules(tree; use_shortforms=true) for tree in SoleModels.models(model)
            ])
            # TODO maybe also sort?
            rs isa Vector{<:Vector{<:Any}} ? reduce(vcat, rs) : rs
        else
            listrules(model; use_shortforms=true)
        end
    end
    silent || println(ruleset)

    silent || println("\n\n$COLORED_TITLE$TITLE\n PART 2.b ATOM EXTRACTION \n$TITLE$RESET")

    num_all_atoms, my_atoms, my_alphabet = begin
        all_atoms = collect(atoms(SoleModels.alphabet(model, false)))
        num_all_atoms = length(all_atoms)
        my_atoms = unique(all_atoms)

        silent || println(my_atoms)

        silent || println(
            "\n\n$COLORED_TITLE$TITLE\n PART 2.c ALPHABET EXTRACTION pt2 \n$TITLE$RESET",
        )

        # Get number of features from the maximum feature index in atoms
        n_features =
            maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
        !isa(n_features, Integer) && error("Symbolic feature names not supported")

        my_alphabet = filteralphabetcallback(process_alphabet(my_atoms, n_features))

        all(x -> (x == (<)), SoleData.test_operator.(subalphabets(my_alphabet))) ||
        error("Only < operator supported")

        num_all_atoms, my_atoms, my_alphabet
    end

    if testott == nothing && alphabetcontroll == nothing
        silent ||
            println("\n\n$COLORED_TITLE$TITLE\n PART 3 TABLE GENERATION \n$TITLE$RESET")

        # PART 3: Table Generation

        results, label_count = @time "Lumen: computing combinations" begin
            if ott_mode

                silent || println("my_alphabet ", my_alphabet)
                silent || println("my_atoms", my_atoms)

                truth_combinations_ott(
                    modelJ,
                    my_alphabet,
                    my_atoms,
                    vertical,
                    vetImportance;
                    silent,
                    apply_function,
                ) # Combination generate with cartesian join
            else
                truth_combinations(
                    modelJ,
                    my_alphabet,
                    my_atoms,
                    vertical;
                    silent,
                    apply_function,
                ) # Combination generate with classic Binary clock table
            end
        end


        # Feature Processing TODO REMOVE
        #= OLD CODE were we haven't Constructor

            num_atoms = length(my_atoms)
            thresholds_by_feature = Dict(
                subalpha.featcondition[1].feature.i_variable => sort(subalpha.featcondition[2])
                for subalpha in my_alphabet.subalphabets
            )

            atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
            for atom in my_atoms
                feat = atom.value.metacond.feature.i_variable
                threshold = atom.value.threshold
                push!(
                    get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat),
                    (threshold, true),
                    )
            end
            for (_, atom_list) in atoms_by_feature
                sort!(atom_list, by=first)
            end

            # Results Processing
            combined_results =
                Lumen.concat_results(results, num_atoms, thresholds_by_feature, atoms_by_feature)
        =#


        # I take combinations with the same class L and disjoin them with logical OR
        combined_results = Lumen.concat_results(results, my_atoms)

        return_info && (unminimized_rules = Rule[])
        minimized_rules = Rule[]
        vectPrePostNumber = Vector{Tuple{Int,Int}}()

        # merge_negated_elements version:
        # To create a minimized version of the majority class, we can minimize and then negate the disjunction of all other classes.

        silent ||
            println("combined_results prima dell'ordinamento: ", collect(combined_results))

        sorted_results = sort(collect(combined_results), by=p -> length(p[2].combinations))
        max_key = sorted_results[end][1]  # Prendi la chiave dell'ultimo elemento


        silent || println("sorted_results: ", sorted_results)
        silent || println("max_key: ", max_key)

        # Process each result
        for (result, formula) in sorted_results
            silent || println("Performing minimization for: $result")
            
            if result == max_key && merge_negated_elements
                silent || println("Computing last rule as negation of disjunction of all others")

                try
                    #println("yet minimized rules: ",minimized_rules[0].antecedent)
                    # Extract antecedents from all previously minimized rules
                    antecedents =
                        [el_to_string(rule.antecedent) for rule in minimized_rules]

                    # Convert each antecedent to string directly
                    #println("string: ",antecedents)
                    antecedent_strings = [string("¬($ant)") for ant in antecedents]

                    #println(antecedent_strings)

                    # Create big disjunction string
                    big_dnf_str = join(antecedent_strings, " ∧ ")

                    silent || println("Debug: Negated string: $big_dnf_str")

                    # Parse negated formula
                    φ = SoleLogics.parseformula(
                        big_dnf_str;
                        atom_parser=a -> Atom(
                            parsecondition(
                                SoleData.ScalarCondition,
                                a;
                                featuretype=SoleData.VariableValue,
                                featvaltype=Real,
                            ),
                        ),
                    )

                    silent || println("Debug: Parsed formula: $φ")


                    # Apply dnf
                    #final_antecedent_dnf = φ                                               # EITHER THE ¬(DNF)
                    final_antecedent_dnf = normalize_formula(φ, :cnf)        # OR THE CNF
                    #println(dump(final_antecedent_dnf.grandchildren))

                    #final_antecedent_dnf = normalize_formula(final_antecedent_dnf, :dnf)  # TODO REMOVE
                    silent || println("Debug: DNF result: $final_antecedent_dnf")                               #    THIS FOR TEST
                    new_rule = Rule(final_antecedent_dnf, result)

                    # One-line per estrarre solo la formula da string(new_rule)
                    big_dnf_str = strip(split(replace(replace(string(new_rule), r"\e\[[0-9;]*m" => ""), r"^▣\s+" => ""), "↣")[1])
                    # Parse negated formula
                    φ = SoleLogics.parseformula(
                        big_dnf_str;
                        atom_parser=a -> Atom(
                            parsecondition(
                                SoleData.ScalarCondition,
                                a;
                                featuretype=SoleData.VariableValue,
                                featvaltype=Real,
                            ),
                        ),
                    )
                    new_rule = Rule(φ, result)
                    silent || println("new_rule:", new_rule)
                    push!(minimized_rules, new_rule)

                catch e
                    @error "Error creating negated rule: $e"
                    # FALLBACK

                    # TODO RESTORE CODE REMOVE COMMENT ..
                    #if return_info .
                    #    push!(unminimized_rules, convert_DNF_formula(formula, result, 1.0))
                    #end

                    formula_semplificata = Lumen.minimizza_dnf(
                        Val(minimization_scheme),
                        formula;
                        minimization_kwargs...,
                    )
                    try

                        if !is_minimization_scheme_espresso
                            silent || println("==========================")
                            silent ||
                                println("comb:", eachcombination(formula_semplificata))
                            silent || println("==========================")

                            ntermpresemp = nterms(formula)
                            ntermpostsemp = nterms(formula_semplificata)
                            push!(vectPrePostNumber, (ntermpresemp, ntermpostsemp))

                            silent || println("Term origin: ", nterms(formula))
                            silent ||
                                println("Term after semp: ", nterms(formula_semplificata))
                            silent || println("Atoms/term orig: ", natomsperterm(formula))

                            silent || println(
                                "Atom/term post semp: ",
                                #natomsperterm(formula_semplificata),
                            )
                        end
                        silent || println()
                        if is_minimization_scheme_espresso
                            formula_string = leftmost_disjunctive_form_to_string(
                                formula_semplificata,
                                horizontal,
                                vetImportance,
                            )
                            φ = SoleLogics.parseformula(
                                formula_string;
                                atom_parser=a -> Atom(
                                    parsecondition(
                                        SoleData.ScalarCondition,
                                        a;
                                        featuretype=SoleData.VariableValue,
                                        featvaltype=Real,
                                    ),
                                ),
                            )
                            new_rule = Rule(φ, result)
                        else
                            new_rule = convert_DNF_formula(
                                formula_semplificata,
                                result,
                                horizontal,
                            )
                        end
                        #new_rule = Rule(ant, result)
                        silent || println(new_rule)
                        push!(minimized_rules, new_rule)
                    catch e
                        @error "Simplification error: $e"
                    end
                    # END FALLBACK
                end

            else    # !merge_negated_elements or result is not max // finaly, continue in classic minimization

                # TODO RESTORE CODE REMOVE COMMENT ..
                #if return_info .
                #    push!(unminimized_rules, convert_DNF_formula(formula, result, 1.0))
                #end

                formula_semplificata = Lumen.minimizza_dnf(
                    Val(minimization_scheme),
                    formula;
                    minimization_kwargs...,
                )
                try

                    if !is_minimization_scheme_espresso
                        silent || println("==========================")
                        silent || println("comb:", eachcombination(formula_semplificata))
                        silent || println("==========================")

                        ntermpresemp = nterms(formula)
                        ntermpostsemp = nterms(formula_semplificata)
                        push!(vectPrePostNumber, (ntermpresemp, ntermpostsemp))

                        silent || println("Term origin: ", nterms(formula))
                        silent || println("Term after semp: ", nterms(formula_semplificata))
                        silent || println("Atoms/term orig: ", natomsperterm(formula))

                        silent || println(
                            "Atom/term post semp: ",
                            #natomsperterm(formula_semplificata),
                        )
                    end
                    silent || println()
                    if is_minimization_scheme_espresso
                        formula_string = leftmost_disjunctive_form_to_string(
                            formula_semplificata,
                            horizontal,
                            vetImportance,
                        )
                        φ = SoleLogics.parseformula(
                            formula_string;
                            atom_parser=a -> Atom(
                                parsecondition(
                                    SoleData.ScalarCondition,
                                    a;
                                    featuretype=SoleData.VariableValue,
                                    featvaltype=Real,
                                ),
                            ),
                        )
                        new_rule = Rule(φ, result)
                    else
                        new_rule =
                            convert_DNF_formula(formula_semplificata, result, horizontal)
                    end
                    #new_rule = Rule(ant, result)
                    silent || println(new_rule)
                    push!(minimized_rules, new_rule)
                catch e
                    @error "Simplification error: $e"
                end
            end
        end

        # Return results
        ds = DecisionSet(minimized_rules)
        if !return_info
            return ds
        end

        info = (;)
        info = merge(info, (; vectPrePostNumber=vectPrePostNumber))
        return_info &&
            (info = merge(info, (; unminimized_ds=DecisionSet(unminimized_rules))))

        return ds, info
    end

    if testott != nothing
        silent || println(
            "\n\n$COLORED_INFO$TITLE\n PART 2.d IS THE OPTIMIZATION VALID?\n$TITLE$RESET",
        )
        testOttt(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function, testott)
        return nothing, nothing
    end

    if alphabetcontroll != nothing
        silent || println("\n\n$COLORED_INFO$TITLE\n ANALIZE ONLY ALPHABET\n$TITLE$RESET")
        debug_combinations(
            modelJ,
            my_alphabet,
            my_atoms,
            vertical;
            silent,
            apply_function,
            alphabetcontroll,
        )
        return nothing, nothing
    end
end

##
# Types
##

# TritVector TODO Choose and implement one ?
include("types/trit-vector.jl")
include("types/balanced-trit-vector.jl")
include("types/balanced-ternary-vector.jl")

# ConstVariable e myOwnTypes
include("types/types.jl")

##
# Utils
##

include("reportAlphabet.jl")

# File report management
include("utils/report.jl")

# IO-Algorithmic-Utils - Project's algorithmic IO
include("utils/IO.jl")

# Minor-Algorithmic-Utils - Project's algorithmic core
include("utils/minor.jl")

# Algorithmic-Utils - Project's algorithmic core
include("utils/core.jl")

# Algorithmic-Optimization-Utils - Project's algorithmic core when running in ott_mode
include("utils/coreOttMode.jl")
include("utils/minimization.jl")
include("deprecate.jl")

end
