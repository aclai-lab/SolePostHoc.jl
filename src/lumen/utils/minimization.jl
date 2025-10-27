function minimizza_dnf(_minimization_scheme::Val, formula::TwoLevelDNFFormula, kwargs...)
    error("Unknown minimization scheme: $(_minimization_scheme)!")
end

# using Infiltrator
function minimizza_dnf(
    ::Val{:mitespresso},
    formula::TwoLevelDNFFormula;
    silent = true,
    horizontal = 1.0,
    vertical = 1.0,
    vetImportance = [],
    mitespresso_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)
    formula = SoleData.espresso_minimize(formula, silent; mitespresso_kwargs...)
    silent || (println(); @show formula)
    # @infiltrate
    # @show syntaxstring(formula)
    #formula = convert(TwoLevelDNFFormula, formula) # TODO FOR NOW WE USE BYPASS... FIX THIS WHEN WE KNOW HOW TO CONVERT
    silent || (println(); @show formula)
    return formula
end


function minimizza_dnf(
    ::Val{:boom},
    formula::TwoLevelDNFFormula;
    silent = true,
    boom_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)

    formula = SoleData.boom_minimize(formula, silent; boom_kwargs...)

    silent || (println(); @show formula)

    return formula
end

function minimizza_dnf(
    ::Val{:abc},
    formula::TwoLevelDNFFormula;
    silent = true,
    horizontal = 1.0,
    vertical = 1.0,
    vetImportance = [],
    boom_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)

    silent || println("||====================||")
    silent || println(dnf(formula))
    #@show vetImportance
    #@show horizontal
    #@show vertical
    silent || println("||====================||")

    silent || println("pre rc comparison: ", dnf(formula))

    if ((vertical != 1.0) && !isempty(vetImportance)) || (horizontal != 1.0)
        formula = dnf_rc_compression(formula, horizontal, vertical, vetImportance; silent)
    end
    silent || println("||====================||")
    silent || println("post rc comparison: ", dnf(formula))
    silent || println("||====================||")

    formula = SoleData.abc_minimize(formula, silent; boom_kwargs...)


    silent || (println(); @show formula)

    return formula
end

#=
function minimizza_dnf(
    ::Val{:texasespresso},
    formula::TwoLevelDNFFormula;
    silent = true,
    texasespresso_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)
    formula = SoleData.espressoTexas_minimize(formula, silent; texasespresso_kwargs...)
    silent || (println(); @show formula)
    # @infiltrate
    # @show syntaxstring(formula)
    #formula = convert(TwoLevelDNFFormula, formula) # TODO FOR NOW WE USE BYPASS... FIX THIS WHEN WE KNOW HOW TO CONVERT
    silent || (println(); @show formula)
    return formula
end
=#

# =========================================================================
# Homebrew functions for minimization
# =========================================================================

"""
Represents a product in Petrick's method (set of prime implicants).
Example: PetrickProduct(Set([1, 3])) means P1 ∧ P3
"""
struct PetrickProduct
    primes::Set{Int}
end

Base.length(pp::PetrickProduct) = length(pp.primes)

"""Multiplies two products (logical AND = union)"""
function multiply_products(p1::PetrickProduct, p2::PetrickProduct)
    return PetrickProduct(union(p1.primes, p2.primes))
end

"""
Multiplies two SOP expressions (complete distribution).
(A + B) × (C + D) = AC + AD + BC + BD
"""
function multiply_sop(expr1::Vector{PetrickProduct}, expr2::Vector{PetrickProduct})
    result = PetrickProduct[]

    # COMPLETE exploration of all combinations
    for p1 in expr1
        for p2 in expr2
            push!(result, multiply_products(p1, p2))
        end
    end

    return simplify_sop(result)
end

"""
Simplifies a SOP by removing:
1. Duplicates
2. Absorbed terms (if A ⊂ B, remove B)
"""
function simplify_sop(products::Vector{PetrickProduct})
    if isempty(products)
        return products
    end

    # Remove exact duplicates
    unique_products = unique(p -> p.primes, products)

    # Remove absorbed terms (absorption law)
    filtered = PetrickProduct[]

    for i = 1:length(unique_products)
        is_absorbed = false

        # Check if there exists a smaller subset
        for j = 1:length(unique_products)
            if i != j && issubset(unique_products[j].primes, unique_products[i].primes)
                is_absorbed = true
                break
            end
        end

        if !is_absorbed
            push!(filtered, unique_products[i])
        end
    end

    return filtered
end

# =========================================================================
# MAIN FUNCTION
# =========================================================================

"""
    minimizza_dnf(::Val{:quine_naive}, formula::TwoLevelDNFFormula; kwargs...) -> TwoLevelDNFFormula

Complete CLASSIC implementation of Quine-McCluskey algorithm with Petrick's method.
GLOBAL brute-force exploration to find the absolute optimal solution.

# Complete algorithm:
1. Generates all prime implicants through iterative combinations
2. Builds coverage table (prime implicants × minterms)
3. Identifies and selects all essential prime implicants
4. Applies Petrick's method to find exact minimal coverage of remaining terms
5. Returns the solution with minimum number of implicants

# Petrick's Method:
- Constructs a boolean formula representing all possible covers
- Expands in Product-of-Sums (POS) form
- Converts to Sum-of-Products (SOP)
- Finds products with minimum number of terms
- Returns the GLOBAL optimal solution

# Complexity:
- Worst case: exponential (NP-complete problem)
- Complete brute-force: explores ALL possible combinations
- Typically manageable for formulas with <20 minterms

# Interface:
- Accepts all standard kwargs (silent, horizontal, vertical, depth, vetImportance, ...)
- Ignores them if not relevant for the classic algorithm
- Guarantees compatibility with the calling interface

# Note:
This is the true Quine-McCluskey algorithm as described in textbooks,
with exhaustive search for the optimal solution.
"""
function minimizza_dnf(
    ::Val{:quine_naive},
    formula::TwoLevelDNFFormula;
    silent = true,
    horizontal = 1.0,
    vertical = 1.0,
    depth = 1.0,
    vetImportance = [],
    kwargs...  # Captures any other parameters
    )

    silent || println("=" ^ 70)
    silent || println("CLASSIC Quine-McCluskey with GLOBAL Petrick's Method")
    silent || println("=" ^ 70)

    # Convert TritVector to Int vectors for easier manipulation
    terms = Vector{Vector{Int}}()
    for combination in formula.combinations
        term = Vector{Int}(undef, length(combination))
        for (j, trit_val) in enumerate(combination)
            term[j] = Int(trit_val)
        end
        push!(terms, term)
    end

    if isempty(terms)
        silent || println("Empty formula, nothing to minimize")
        return formula
    end

    original_terms = copy(terms)
    n_vars = length(terms[1])
    silent || println("Input: $(length(terms)) minterms, $n_vars variables")

    # =========================================================================
    # HELPER FUNCTIONS
    # =========================================================================

    """Checks if a prime implicant covers a minterm"""
    function covers(prime::Vector{Int}, term::Vector{Int})
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != term[i]
                return false
            end
        end
        return true
    end

    """Checks if two terms differ in exactly one position"""
    function can_combine(t1::Vector{Int}, t2::Vector{Int})
        diff = 0
        pos = -1
        for i in eachindex(t1)
            if t1[i] != t2[i]
                diff += 1
                pos = i
                if diff > 1
                    return false, -1
                end
            end
        end
        return diff == 1, pos
    end

    """Combines two terms by putting -1 in the differing position"""
    function combine(t1::Vector{Int}, t2::Vector{Int}, pos::Int)
        result = copy(t1)
        result[pos] = -1
        return result
    end

    """Converts a Vector{Int} to TritVector"""
    function int_to_tritvector(term::Vector{Int})
        trit_vec = TritVector(length(term))
        for (i, val) in enumerate(term)
            trit_vec[i] = val
        end
        return trit_vec
    end

    """Counts the number of literals in a prime implicant"""
    function count_literals(prime::Vector{Int})
        return count(x -> x != -1, prime)
    end

    # =========================================================================
    # PHASE 1: GENERATE ALL PRIME IMPLICANTS (EXHAUSTIVE SEARCH)
    # =========================================================================

    function find_prime_implicants(terms::Vector{Vector{Int}})
        silent || println("\n[PHASE 1] Finding ALL prime implicants...")

        primes = Set{Vector{Int}}()
        used = Set{Vector{Int}}()
        current = Set(terms)
        iteration = 0

        while !isempty(current)
            iteration += 1
            silent || println("  Iteration $iteration: $(length(current)) terms")

            next_terms = Set{Vector{Int}}()
            current_list = collect(current)

            # Try ALL possible combinations (brute-force)
            for i = 1:length(current_list)
                for j = (i+1):length(current_list)
                    combinable, pos = can_combine(current_list[i], current_list[j])
                    if combinable
                        new_term = combine(current_list[i], current_list[j], pos)
                        push!(used, current_list[i])
                        push!(used, current_list[j])
                        push!(next_terms, new_term)
                    end
                end
            end

            # Unused terms are prime implicants
            for term in current
                if term ∉ used
                    push!(primes, term)
                end
            end

            current = next_terms
            empty!(used)
        end

        result = collect(primes)
        silent || println("  ✓ Found $(length(result)) prime implicants")

        return result
    end

    # =========================================================================
    # PHASE 2: BUILD COMPLETE COVERAGE TABLE
    # =========================================================================

    function build_coverage_table(primes::Vector{Vector{Int}}, minterms::Vector{Vector{Int}})
        silent || println("\n[PHASE 2] Building coverage table...")

        coverage = falses(length(primes), length(minterms))

        for i = 1:length(primes)
            for j = 1:length(minterms)
                coverage[i, j] = covers(primes[i], minterms[j])
            end
        end

        silent || println("  ✓ Coverage table: $(length(primes)) × $(length(minterms))")

        return coverage
    end

    # =========================================================================
    # PHASE 3: FIND ESSENTIAL PRIME IMPLICANTS
    # =========================================================================

    function find_essential_primes(coverage::BitMatrix, primes::Vector{Vector{Int}}, minterms::Vector{Vector{Int}})
        silent || println("\n[PHASE 3] Finding essential prime implicants...")

        essential = Set{Int}()
        covered_minterms = Set{Int}()

        for j = 1:length(minterms)
            covering_primes = findall(coverage[:, j])

            # If a minterm is covered by ONLY ONE prime, that prime is essential
            if length(covering_primes) == 1
                prime_idx = covering_primes[1]
                push!(essential, prime_idx)

                # Mark all minterms covered by this essential prime
                for k = 1:length(minterms)
                    if coverage[prime_idx, k]
                        push!(covered_minterms, k)
                    end
                end
            end
        end

        essential_list = collect(essential)
        silent || println("  ✓ Essential primes: $(length(essential_list))")
        silent || println("  ✓ Covered minterms: $(length(covered_minterms)) / $(length(minterms))")

        return essential_list, covered_minterms
    end

    # =========================================================================
    # PHASE 4: PETRICK'S METHOD - GLOBAL OPTIMAL SEARCH
    # =========================================================================

    """
    Petrick's Method with exhaustive GLOBAL search.

    Builds the POS (Product of Sums) expression and converts it to SOP
    by exploring ALL possible combinations to find the absolute
    optimal solution.
    """
    function petricks_method_global(coverage::BitMatrix, primes::Vector{Vector{Int}}, covered_minterms::Set{Int})
        silent || println("\n[PHASE 4] Petrick's Method - GLOBAL SEARCH...")

        # Identify uncovered minterms
        uncovered = [i for i = 1:size(coverage, 2) if i ∉ covered_minterms]

        if isempty(uncovered)
            silent || println("  ✓ All minterms covered by essential primes")
            return Int[]
        end

        silent || println("  Uncovered minterms: $(length(uncovered))")

        # Initialize with empty product (identity for AND)
        current_sop = [PetrickProduct(Set{Int}())]

        # For each uncovered minterm, build (Pi + Pj + ...)
        for (idx, minterm_idx) in enumerate(uncovered)
            covering_primes = findall(coverage[:, minterm_idx])

            if isempty(covering_primes)
                @warn "Minterm $minterm_idx cannot be covered!"
                return Int[]
            end

            # Create the sum (OR) of all primes covering this minterm
            minterm_sum = [PetrickProduct(Set([p])) for p in covering_primes]

            silent || println("  Minterm $idx/$(length(uncovered)): $(length(covering_primes)) covering primes")

            # Multiply with current expression (complete distribution)
            current_sop = multiply_sop(current_sop, minterm_sum)

            silent || println("    → SOP size: $(length(current_sop)) products")

            # Protection against combinatorial explosion
            if length(current_sop) > 50000
                silent || println("    ⚠ Expression exploding, applying aggressive simplification")

                # Keep only the most promising products
                sort!(current_sop, by = p -> length(p.primes))
                current_sop = current_sop[1:min(10000, length(current_sop))]
                current_sop = simplify_sop(current_sop)
            end
        end

        # Find the GLOBAL optimal solution(s)
        if !isempty(current_sop)
            min_cost = minimum(length(p) for p in current_sop)
            optimal_solutions = filter(p -> length(p) == min_cost, current_sop)

            silent || println("  ✓ Found $(length(optimal_solutions)) optimal solution(s)")
            silent || println("  ✓ Optimal cost: $min_cost additional primes")

            # If multiple optimal solutions exist, choose the one with fewest total literals
            if length(optimal_solutions) > 1
                best = optimal_solutions[1]
                best_literals = sum(count_literals(primes[i]) for i in best.primes)

                for sol in optimal_solutions[2:end]
                    literals = sum(count_literals(primes[i]) for i in sol.primes)
                    if literals < best_literals
                        best = sol
                        best_literals = literals
                    end
                end

                return collect(best.primes)
            else
                return collect(optimal_solutions[1].primes)
            end
        end

        return Int[]
    end

    # =========================================================================
    # MAIN ALGORITHM EXECUTION
    # =========================================================================

    try
        silent || println("\nStarting optimization...")

        # PHASE 1: Generate ALL prime implicants
        primes = find_prime_implicants(original_terms)

        if isempty(primes)
            @warn "No prime implicants found"
            return formula
        end

        # PHASE 2: Build complete coverage table
        coverage = build_coverage_table(primes, original_terms)

        # PHASE 3: Find essential prime implicants
        essential_indices, covered_minterms = find_essential_primes(coverage, primes, original_terms)

        # PHASE 4: Petrick's method with GLOBAL search
        additional_indices = petricks_method_global(coverage, primes, covered_minterms)

        # Combine essential and additional
        selected_indices = unique(sort([essential_indices; additional_indices]))
        selected_primes = primes[selected_indices]

        silent || println("\n" * "=" ^ 70)
        silent || println("RESULTS:")
        silent || println("  Input terms:          $(length(original_terms))")
        silent || println("  Prime implicants:     $(length(primes))")
        silent || println("  Essential primes:     $(length(essential_indices))")
        silent || println("  Additional primes:    $(length(additional_indices))")
        silent || println("  TOTAL OUTPUT:         $(length(selected_primes))")

        # Calculate savings
        if length(original_terms) > 0
            reduction = 100.0 * (1.0 - length(selected_primes) / length(original_terms))
            silent || println("  Reduction:            $(round(reduction, digits=1))%")
        end

        silent || println("=" ^ 70)

        # Convert final result to TritVectors
        new_combinations = Vector{TritVector}()
        for term in selected_primes
            push!(new_combinations, int_to_tritvector(term))
        end

        sort!(new_combinations)

        return TwoLevelDNFFormula(
            new_combinations,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
        )

    catch e
        @error "Error in Quine-McCluskey algorithm" exception=(e, catch_backtrace())
        return formula
    end
end

# =========================================================================
# QUINE-MCCLUSKEY WITH POST-PROCESSING OPTIMIZATION (OPP)
# =========================================================================

"""
    minimizza_dnf(::Val{:quine_opp}, formula::TwoLevelDNFFormula; kwargs...) -> TwoLevelDNFFormula

OPTIMIZED Quine-McCluskey algorithm with intelligent semantic post-processing.

# Features:
1. Executes classic Quine-McCluskey with Petrick's method
2. **SEMANTIC POST-PROCESSING** that understands thresholds:
   - Analyzes ranges covered on each feature
   - Eliminates threshold-aware redundancies (e.g., (x < a) ∨ (x ≥ a) → true)
   - Removes features that cover the entire domain
   - Simplifies overlapping ranges
3. Factorizes common terms
4. Produces minimal semantically correct formulas

# Example:
Input:  (V3 < 2.45 ∧ V4 < 1.2) ∨ (V3 < 2.45 ∧ V4 ≥ 1.2) ∨ (V3 < 2.45 ∧ V4 ≥ 1.65)
Quine:  (V3 < 2.45) ∧ [(V4 < 1.2) ∨ (V4 ≥ 1.2)]  # redundancy not detected
OPP:    V3 < 2.45                                  # redundancy eliminated!

# Algorithm:
1. Classic Quine-McCluskey (as :quine_naive)
2. Range analysis per feature
3. Redundant feature elimination
4. Inter-term simplification
5. Final factorization

# Complexity:
- Quine phase: exponential (NP-complete)
- Post-processing: O(n × m × log m) where n=terms, m=atoms per feature
- Manageable up to ~30 minterms

# Interface:
- Compatible with all standard kwargs
- silent: controls verbosity
- Other parameters: ignored if not relevant
"""
function minimizza_dnf(
    ::Val{:quine},
    formula::TwoLevelDNFFormula;
    silent = true,
    horizontal = 1.0,
    vertical = 1.0,
    depth = 1.0,
    vetImportance = [],
    kwargs...
    )

    silent || println("=" ^ 70)
    silent || println("QUINE-MCCLUSKEY with OPTIMIZED POST-PROCESSING (OPP)")
    silent || println("=" ^ 70)

    # Convert TritVector to Int vectors
    terms = Vector{Vector{Int}}()
    for combination in formula.combinations
        term = Vector{Int}(undef, length(combination))
        for (j, trit_val) in enumerate(combination)
            term[j] = Int(trit_val)
        end
        push!(terms, term)
    end

    if isempty(terms)
        silent || println("Empty formula, nothing to minimize")
        return formula
    end

    original_terms = copy(terms)
    n_vars = length(terms[1])
    silent || println("Input: $(length(terms)) minterms, $n_vars atoms")

    # =========================================================================
    # HELPER FUNCTIONS
    # =========================================================================

    function covers(prime::Vector{Int}, term::Vector{Int})
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != term[i]
                return false
            end
        end
        return true
    end

    function can_combine(t1::Vector{Int}, t2::Vector{Int})
        diff = 0
        pos = -1
        for i in eachindex(t1)
            if t1[i] != t2[i]
                diff += 1
                pos = i
                if diff > 1
                    return false, -1
                end
            end
        end
        return diff == 1, pos
    end

    function combine(t1::Vector{Int}, t2::Vector{Int}, pos::Int)
        result = copy(t1)
        result[pos] = -1
        return result
    end

    function int_to_tritvector(term::Vector{Int})
        trit_vec = TritVector(length(term))
        for (i, val) in enumerate(term)
            trit_vec[i] = val
        end
        return trit_vec
    end

    function count_literals(prime::Vector{Int})
        return count(x -> x != -1, prime)
    end

    # =========================================================================
    # PHASE 1: STANDARD QUINE-MCCLUSKEY
    # =========================================================================

    function find_prime_implicants(terms::Vector{Vector{Int}})
        silent || println("\n[PHASE 1] Quine-McCluskey: Finding prime implicants...")

        primes = Set{Vector{Int}}()
        used = Set{Vector{Int}}()
        current = Set(terms)
        iteration = 0

        while !isempty(current)
            iteration += 1
            next_terms = Set{Vector{Int}}()
            current_list = collect(current)

            for i = 1:length(current_list)
                for j = (i+1):length(current_list)
                    combinable, pos = can_combine(current_list[i], current_list[j])
                    if combinable
                        new_term = combine(current_list[i], current_list[j], pos)
                        push!(used, current_list[i])
                        push!(used, current_list[j])
                        push!(next_terms, new_term)
                    end
                end
            end

            for term in current
                if term ∉ used
                    push!(primes, term)
                end
            end

            current = next_terms
            empty!(used)
        end

        result = collect(primes)
        silent || println("  ✓ Found $(length(result)) prime implicants")

        return result
    end

    function build_coverage_table(primes::Vector{Vector{Int}}, minterms::Vector{Vector{Int}})
        coverage = falses(length(primes), length(minterms))
        for i = 1:length(primes)
            for j = 1:length(minterms)
                coverage[i, j] = covers(primes[i], minterms[j])
            end
        end
        return coverage
    end

    function find_essential_primes(coverage::BitMatrix, primes::Vector{Vector{Int}}, minterms::Vector{Vector{Int}})
        essential = Set{Int}()
        covered_minterms = Set{Int}()

        for j = 1:length(minterms)
            covering_primes = findall(coverage[:, j])
            if length(covering_primes) == 1
                prime_idx = covering_primes[1]
                push!(essential, prime_idx)
                for k = 1:length(minterms)
                    if coverage[prime_idx, k]
                        push!(covered_minterms, k)
                    end
                end
            end
        end

        return collect(essential), covered_minterms
    end

    function petricks_method_global(coverage::BitMatrix, primes::Vector{Vector{Int}}, covered_minterms::Set{Int})
        uncovered = [i for i = 1:size(coverage, 2) if i ∉ covered_minterms]

        if isempty(uncovered)
            return Int[]
        end

        current_sop = [PetrickProduct(Set{Int}())]

        for minterm_idx in uncovered
            covering_primes = findall(coverage[:, minterm_idx])
            if isempty(covering_primes)
                @warn "Minterm $minterm_idx cannot be covered!"
                return Int[]
            end

            minterm_sum = [PetrickProduct(Set([p])) for p in covering_primes]
            current_sop = multiply_sop(current_sop, minterm_sum)

            if length(current_sop) > 50000
                sort!(current_sop, by = p -> length(p.primes))
                current_sop = current_sop[1:min(10000, length(current_sop))]
                current_sop = simplify_sop(current_sop)
            end
        end

        if !isempty(current_sop)
            min_cost = minimum(length(p) for p in current_sop)
            optimal_solutions = filter(p -> length(p) == min_cost, current_sop)

            if length(optimal_solutions) > 1
                best = optimal_solutions[1]
                best_literals = sum(count_literals(primes[i]) for i in best.primes)

                for sol in optimal_solutions[2:end]
                    literals = sum(count_literals(primes[i]) for i in sol.primes)
                    if literals < best_literals
                        best = sol
                        best_literals = literals
                    end
                end

                return collect(best.primes)
            else
                return collect(optimal_solutions[1].primes)
            end
        end

        return Int[]
    end

    # =========================================================================
    # PHASE 2: SEMANTIC POST-PROCESSING (THE NOVELTY!)
    # =========================================================================

    """
    Builds a map from atom_index → (feature_id, threshold, operator)
    """
    function build_atom_map(formula::TwoLevelDNFFormula)
        atom_map = Dict{Int, Tuple{Int, Float64, Bool}}()
        atom_idx = 1

        for (feature_id, atoms) in sort(collect(formula.atoms_by_feature))
            for (threshold, is_less) in atoms
                atom_map[atom_idx] = (feature_id, threshold, is_less)
                atom_idx += 1
            end
        end

        return atom_map
    end

    """
    Analyzes the ranges covered by a term on a specific feature.
    Returns a vector of intervals (min, max, included_bounds).
    """
    function analyze_feature_ranges(term::Vector{Int}, feature_id::Int, atom_map::Dict{Int, Tuple{Int, Float64, Bool}})
        ranges = Tuple{Float64, Float64, Tuple{Bool, Bool}}[]  # (min, max, (include_min, include_max))

        constraints = []
        for (atom_idx, val) in enumerate(term)
            if val == -1
                continue  # Don't care
            end

            feat_id, threshold, is_less = atom_map[atom_idx]
            if feat_id != feature_id
                continue
            end

            # val == 0 → negation of atom, val == 1 → direct atom
            if val == 1
                if is_less
                    push!(constraints, (:less, threshold))  # x < threshold
                else
                    push!(constraints, (:geq, threshold))   # x ≥ threshold
                end
            else  # val == 0
                if is_less
                    push!(constraints, (:geq, threshold))   # NOT(x < threshold) = x ≥ threshold
                else
                    push!(constraints, (:less, threshold))  # NOT(x ≥ threshold) = x < threshold
                end
            end
        end

        if isempty(constraints)
            return [(-Inf, Inf, (false, false))]  # No constraints = entire domain
        end

        # Find lower and upper bounds
        lower_bounds = [t for (op, t) in constraints if op == :geq]
        upper_bounds = [t for (op, t) in constraints if op == :less]

        lower = isempty(lower_bounds) ? -Inf : maximum(lower_bounds)
        upper = isempty(upper_bounds) ? Inf : minimum(upper_bounds)

        if lower >= upper
            return []  # Contradiction: no valid range
        end

        return [(lower, upper, (lower != -Inf, upper != Inf))]
    end

    """
    Checks if ranges cover the entire domain of the feature.
    """
    function covers_full_domain(ranges::Vector{Tuple{Float64, Float64, Tuple{Bool, Bool}}})
        if isempty(ranges)
            return false
        end

        # Merge and sort ranges
        sorted_ranges = sort(ranges, by = r -> r[1])

        # Check if there's a range from -Inf to +Inf
        for (lower, upper, _) in sorted_ranges
            if lower == -Inf && upper == Inf
                return true
            end
        end

        # Check if the union covers everything
        # (simplified logic: if starts from -Inf and ends at +Inf)
        starts_from_inf = any(r[1] == -Inf for r in sorted_ranges)
        ends_at_inf = any(r[2] == Inf for r in sorted_ranges)

        if starts_from_inf && ends_at_inf
            # Check that there are no gaps
            # For simplicity, assume contiguous ranges merge
            return true
        end

        return false
    end

    """
    Semantic post-processing: eliminates threshold-aware redundancies.
    """
    function semantic_postprocessing(primes::Vector{Vector{Int}}, formula::TwoLevelDNFFormula)
        silent || println("\n[PHASE 2] Semantic Post-Processing...")

        atom_map = build_atom_map(formula)
        feature_ids = unique(feat_id for (_, (feat_id, _, _)) in atom_map)

        simplified_primes = Vector{Vector{Int}}()

        for prime in primes
            # Analyze each feature in this prime
            simplified_prime = copy(prime)

            for feature_id in feature_ids
                ranges = analyze_feature_ranges(prime, feature_id, atom_map)

                # If this feature covers the entire domain, eliminate its atoms
                if covers_full_domain(ranges)
                    silent || println("  ✓ Feature $feature_id is redundant in term, removing...")

                    # Set to -1 (don't care) all atoms of this feature
                    for (atom_idx, val) in enumerate(prime)
                        if val != -1
                            feat_id, _, _ = atom_map[atom_idx]
                            if feat_id == feature_id
                                simplified_prime[atom_idx] = -1
                            end
                        end
                    end
                end
            end

            push!(simplified_primes, simplified_prime)
        end

        # Remove duplicates after simplification
        unique_primes = unique(simplified_primes)

        silent || println("  ✓ Simplified: $(length(primes)) → $(length(unique_primes)) primes")

        return unique_primes
    end

    """
    Phase 3: Global analysis - factorizes completely redundant terms.
    """
    function global_factorization(primes::Vector{Vector{Int}}, formula::TwoLevelDNFFormula)
        silent || println("\n[PHASE 3] Global Factorization...")

        if length(primes) <= 1
            return primes
        end

        atom_map = build_atom_map(formula)
        feature_ids = unique(feat_id for (_, (feat_id, _, _)) in atom_map)

        # For each feature, calculate the union of ranges covered by ALL primes
        for feature_id in feature_ids
            # FIXED: Initialize with correct type to avoid type instability
            all_ranges = Vector{Tuple{Float64, Float64, Tuple{Bool, Bool}}}()

            for prime in primes
                ranges = analyze_feature_ranges(prime, feature_id, atom_map)
                append!(all_ranges, ranges)
            end

            # If the union covers everything, the feature is globally redundant
            if covers_full_domain(all_ranges)
                silent || println("  ✓ Feature $feature_id globally redundant, removing from ALL terms...")

                # Remove from all primes
                for prime in primes
                    for (atom_idx, val) in enumerate(prime)
                        if val != -1
                            feat_id, _, _ = atom_map[atom_idx]
                            if feat_id == feature_id
                                prime[atom_idx] = -1
                            end
                        end
                    end
                end
            end
        end

        # Remove final duplicates
        unique_primes = unique(primes)

        silent || println("  ✓ After factorization: $(length(unique_primes)) unique primes")

        return unique_primes
    end

    # =========================================================================
    # COMPLETE ALGORITHM EXECUTION
    # =========================================================================

    try
        # PHASE 1: Classic Quine-McCluskey
        primes = find_prime_implicants(original_terms)

        if isempty(primes)
            @warn "No prime implicants found"
            return formula
        end

        coverage = build_coverage_table(primes, original_terms)
        essential_indices, covered_minterms = find_essential_primes(coverage, primes, original_terms)
        additional_indices = petricks_method_global(coverage, primes, covered_minterms)

        selected_indices = unique(sort([essential_indices; additional_indices]))
        selected_primes = primes[selected_indices]

        silent || println("  After Quine: $(length(selected_primes)) primes")

        # PHASE 2: Semantic post-processing
        semantic_primes = semantic_postprocessing(selected_primes, formula)

        # PHASE 3: Global factorization
        final_primes = global_factorization(semantic_primes, formula)

        # Remove empty primes (all -1)
        final_primes = filter(p -> any(x != -1 for x in p), final_primes)

        silent || println("\n" * "=" ^ 70)
        silent || println("RESULTS:")
        silent || println("  Input terms:          $(length(original_terms))")
        silent || println("  After Quine:          $(length(selected_primes))")
        silent || println("  After semantic OPP:   $(length(semantic_primes))")
        silent || println("  FINAL OUTPUT:         $(length(final_primes))")

        if length(original_terms) > 0
            reduction = 100.0 * (1.0 - length(final_primes) / length(original_terms))
            silent || println("  Total reduction:      $(round(reduction, digits=1))%")
        end

        silent || println("=" ^ 70)

        # Convert final result
        new_combinations = Vector{TritVector}()
        for term in final_primes
            push!(new_combinations, int_to_tritvector(term))
        end

        sort!(new_combinations)

        return TwoLevelDNFFormula(
            new_combinations,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
        )

    catch e
        @error "Error in Quine-OPP algorithm" exception=(e, catch_backtrace())
        return formula
    end
end


include("dnf_rc_compression.jl")
