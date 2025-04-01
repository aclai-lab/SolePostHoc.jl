module Lumen

using Revise
using Logging
using Dates
using DataStructures
using SoleModels
using SoleModels: DecisionSet
using AbstractTrees
using SoleData
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
    minimization_scheme::Symbol=:mitespresso;
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
    kwargs...
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

    is_minimization_scheme_espresso = minimization_scheme == :mitespresso # TODO we use this only for BYPASS in the code

    # PART 2.a: Starter Ruleset Extraction
    silent || println(
        "\n\n$COLORED_TITLE$TITLE\n PART 2.a STARTER RULESET ESTRACTION \n$TITLE$RESET"
    )

    ruleset = @time begin
        if isensemble(model)
            rs = unique([listrules(tree; use_shortforms=true) for tree in SoleModels.models(model)])
            # TODO maybe also sort?
            rs isa Vector{<:Vector{<:Any}} ? reduce(vcat, rs) : rs
        else
            listrules(model; use_shortforms=true)
        end
    end
    silent || println(ruleset)

    silent || println(
        "\n\n$COLORED_TITLE$TITLE\n PART 2.b ATOM EXTRACTION \n$TITLE$RESET",
    )

    num_all_atoms, my_atoms, my_alphabet = begin
        all_atoms = collect(atoms(SoleModels.alphabet(model, false)))
        num_all_atoms = length(all_atoms)
        my_atoms = unique(all_atoms)

        silent || println(my_atoms)

        silent || println(
            "\n\n$COLORED_TITLE$TITLE\n PART 2.c ALPHABET EXTRACTION pt2 \n$TITLE$RESET",
        )

        # Get number of features from the maximum feature index in atoms
        n_features = maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
        !isa(n_features, Integer) && error("Symbolic feature names not supported")

        my_alphabet = filteralphabetcallback(process_alphabet(my_atoms, n_features))

        all(x -> (x == (<)), SoleData.test_operator.(subalphabets(my_alphabet))) ||
        error("Only < operator supported")

        num_all_atoms, my_atoms, my_alphabet
    end

    if !controllo
        silent || println(
            "\n\n$COLORED_TITLE$TITLE\n PART 3 TABLE GENERATION \n$TITLE$RESET"
        )

        # PART 3: Table Generation
        results, label_count = @time "Lumen: computing combinations" begin
            if ott_mode
                truth_combinations_ott(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
            else
                truth_combinations(modelJ, my_alphabet, my_atoms, vertical; silent, apply_function)
            end
        end

        # Feature Processing
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
        # Results Processing
        combined_results =
            Lumen.concat_results(results, my_atoms)

        return_info && (unminimized_rules = Rule[])
        minimized_rules = Rule[]
        vectPrePostNumber = Vector{Tuple{Int,Int}}()

        # Process each result
        for (result, formula) in combined_results
            silent || println("Risultato: $result")

            if return_info
                push!(unminimized_rules, convert_DNF_formula(formula, result, 1.0))
            end

            formula_semplificata_t = @timed Lumen.minimizza_dnf(
                Val(minimization_scheme),
                formula;
                minimization_kwargs...,
            )
            formula_semplificata = formula_semplificata_t.value



            try
                @info "Simplification completed in $(formula_semplificata_t.time) seconds"

                if !is_minimization_scheme_espresso
                    silent || println("==========================")
                    silent || println("comb:", eachcombination(formula_semplificata))
                    silent || println("==========================")

                    ntermpresemp = nterms(formula)
                    ntermpostsemp = nterms(formula_semplificata)
                    push!(vectPrePostNumber, (ntermpresemp, ntermpostsemp))

                    silent || println("Term origin: ", nterms(formula))
                    silent || println(
                        "Term after semp: ",
                        nterms(formula_semplificata),
                    )
                    silent || println("Atoms/term orig: ", natomsperterm(formula))

                    silent || println(
                        "Atom/term post semp: ",
                        #natomsperterm(formula_semplificata),
                    )
                end
                silent || println()
                if is_minimization_scheme_espresso
                    formula_string = leftmost_disjunctive_form_to_string(formula_semplificata, horizontal, vetImportance)
                    φ = SoleLogics.parseformula(
                        formula_string;
                        atom_parser=a -> Atom(
                            parsecondition(
                                SoleData.ScalarCondition,
                                a;
                                featuretype=SoleData.VariableValue,
                                featvaltype=Real
                            )
                        )
                    )
                    new_rule = Rule(φ, result)
                else
                    new_rule = convert_DNF_formula(
                        formula_semplificata_t.value,
                        result,
                        horizontal
                    )
                end
                #new_rule = Rule(ant, result)
                println(new_rule)
                push!(minimized_rules, new_rule)
            catch e
                @error "Simplification error: $e"
            end
        end

        # Return results
        ds = DecisionSet(minimized_rules)
        if !return_info
            return ds
        end

        info = (;)
        info = merge(info, (; vectPrePostNumber=vectPrePostNumber))
        return_info && (info = merge(info, (; unminimized_ds=DecisionSet(unminimized_rules))))

        return ds, info
    end

    if controllo
        silent || println(
            "\n\n$COLORED_INFO$TITLE\n PART 2.d IS THE OPTIMIZATION VALID?\n$TITLE$RESET",
        )
        are_results_equal =
         Lumen.compare_truth_combinations(modelJ, my_alphabet, my_atoms, vertical; apply_function, silent, kwargs...)
        if are_results_equal
            @info "\nOptimization valid: results are identical."
        else
            @warn "\nWARNING: Optimization might not be valid. Results differ."
        end
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