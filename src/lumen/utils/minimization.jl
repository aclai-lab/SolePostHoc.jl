function minimizza_dnf(_minimization_scheme::Val, formula::TwoLevelDNFFormula, kwargs...)
    error("Unknown minimization scheme: $(_minimization_scheme)!")
end

# using Infiltrator
function minimizza_dnf(
    ::Val{:mitespresso},
    formula::TwoLevelDNFFormula;
    silent = true,
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
    boom_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)

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


"""
    minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula, horizontal = 1.0)

Simplifies a custom OR formula using the Quine algorithm.

This function takes a `TwoLevelDNFFormula` object and applies the Quine algorithm to minimize the number of combinations in the formula. It returns a new `TwoLevelDNFFormula` object with the simplified combinations.
"""
function minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula)
    # Initialize the starting terms
    terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end

    if isempty(terms)
        return formula  # If there are no terms, return the empty formula directly
    end

    # Function to check if one cube covers another
    function covers(cube1::Vector{Int8}, cube2::Vector{Int8})
        for i in eachindex(cube1)
            if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                return false
            end
        end
        return true
    end

    # Function to combine two cubes into a generalized cube
    function combine_cubes(cube1::Vector{Int8}, cube2::Vector{Int8})
        result = Vector{Int8}(undef, length(cube1))
        for i in eachindex(cube1)
            result[i] = ifelse(cube1[i] == cube2[i], cube1[i], Int8(-1))
        end
        return result
    end

    # Minimization step: combining terms
    function minimize_step!(terms::Vector{Vector{Int8}})
        new_terms = Vector{Vector{Int8}}()
        used = fill(false, length(terms))

        # Iterate over all term pairs for combinations
        for i in 1:(length(terms) - 1)
            for j in (i + 1):length(terms)
                diff_count = 0
                for k in eachindex(terms[i])
                    if terms[i][k] != terms[j][k] && terms[i][k] != Int8(-1) && terms[j][k] != Int8(-1)
                        diff_count += 1
                        if diff_count > 1
                            break
                        end
                    end
                end
                # Combine terms only if they differ in a single position
                if diff_count == 1
                    push!(new_terms, combine_cubes(terms[i], terms[j]))
                    used[i] = true
                    used[j] = true
                end
            end
        end

        # Add unused terms
        for i in 1:length(terms)
            if !used[i]
                push!(new_terms, terms[i])
            end
        end

        # Remove duplicates
        unique!(new_terms)

        if length(new_terms) == length(terms)
            return new_terms  # No new combinations, terminate
        end

        return minimize_step!(new_terms)  # Continue iterating
    end

    # Find the minimal coverage of original terms using prime implicants
    function find_minimal_cover(terms::Vector{Vector{Int8}}, primes::Vector{Vector{Int8}})
        coverage = falses(length(primes), length(terms))

        for i in 1:length(primes)
            for j in 1:length(terms)
                coverage[i, j] = covers(primes[i], terms[j])
            end
        end

        # Find minimal cover using backtracking
        function find_cover_backtrack(covered::BitVector, selected_implicants::Vector{Int}, candidates::Vector{Int})
            if all(covered)
                return selected_implicants
            end

            best_solution = nothing

            for candidate in candidates
                new_covered = covered .| coverage[candidate, :]
                new_solution = find_cover_backtrack(new_covered, [selected_implicants; candidate], setdiff(candidates, [candidate]))

                if new_solution !== nothing && (best_solution === nothing || length(new_solution) < length(best_solution))
                    best_solution = new_solution
                end
            end

            return best_solution
        end

        candidates = collect(1:length(primes))
        covered = falses(length(terms))

        minimal_solution = find_cover_backtrack(covered, Int[], candidates)
        return primes[minimal_solution]
    end

    try
        # Step 1: Generate prime implicants
        prime_implicants = minimize_step!(terms)

        # Step 2: Find the minimal coverage of terms
        minimized_terms = find_minimal_cover(terms, prime_implicants)

        # Convert minimized terms into readable combinations
        new_combinations = BitVector[]
        seen = Set{BitVector}()  # Set to avoid duplicates

        for term in minimized_terms
            combo = falses(nuberofatoms(formula))
            for (i, val) in enumerate(term)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            if combo ∉ seen
                push!(seen, combo)
                push!(new_combinations, combo)
            end
        end

        sort!(new_combinations)
        return TwoLevelDNFFormula(
            new_combinations,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            collect(minimized_terms),
        )
    catch e
        @warn "Error during minimization: $e"
        return formula  # Return the original formula in case of error
    end
end


#=MINIMUMS THAT DON'T SEEM GLOBAL TO ME=#
function minimizza_dnf(::Val{:quine_naive}, formula::TwoLevelDNFFormula)
    # Convert terms into binary vectors
    terms = [Vector{Int}([x ? 1 : 0 for x in term]) for term in eachcombination(formula)]
    
    if isempty(terms)
        return formula
    end
    
    # 1. Find all possible prime implicants
    function find_prime_implicants(terms)
        function can_combine(t1, t2)
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
        
        function combine(t1, t2, pos)
            result = copy(t1)
            result[pos] = -1
            return result
        end
        
        primes = Set{Vector{Int}}()
        used = Set{Vector{Int}}()
        current = copy(terms)
        
        while !isempty(current)
            next_terms = Vector{Vector{Int}}()
            
            for i in 1:length(current)
                for j in (i+1):length(current)
                    combinable, pos = can_combine(current[i], current[j])
                    if combinable
                        push!(used, current[i])
                        push!(used, current[j])
                        push!(next_terms, combine(current[i], current[j], pos))
                    end
                end
            end
            
            # Add unused terms as prime implicants
            for term in current
                if term ∉ used
                    push!(primes, term)
                end
            end
            
            current = unique(next_terms)
            empty!(used)
        end
        
        return collect(primes)
    end
    
    # 2. Check if a prime implicant covers a term
    function covers(prime, term)
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != term[i]
                return false
            end
        end
        return true
    end
    
    # 3. Generate all possible combinations of k elements from an array
    function generate_combinations(arr, k)
        n = length(arr)
        if k > n
            return Vector{Vector{eltype(arr)}}()
        end
        
        result = Vector{Vector{eltype(arr)}}()
        
        # Recursive function to generate combinations
        function recursive_combine(start, current)
            if length(current) == k
                push!(result, copy(current))
                return
            end
            
            for i in start:n
                push!(current, arr[i])
                recursive_combine(i + 1, current)
                pop!(current)
            end
        end
        
        recursive_combine(1, eltype(arr)[])
        return result
    end
    
    # 4. Generate all possible combinations of prime implicants
    function all_combinations(primes)
        result = Vector{Vector{Vector{Int}}}()
        for k in 1:length(primes)
            append!(result, generate_combinations(primes, k))
        end
        return result
    end
    
    # 5. Check if a combination covers all terms
    function covers_all(combo, terms)
        for term in terms
            if !any(prime -> covers(prime, term), combo)
                return false
            end
        end
        return true
    end
    
    # Execute the full algorithm
    try
        # Find all prime implicants
        primes = find_prime_implicants(terms)
        
        # Find the minimum combination that covers all terms
        min_cover = primes  # Fallback to the complete solution
        min_size = length(primes)
        
        # Try all possible combinations of prime implicants
        for combo in all_combinations(primes)
            if covers_all(combo, terms) && length(combo) < min_size
                min_cover = combo
                min_size = length(combo)
            end
        end
        
        # Convert the result into a BitVector
        result_combinations = BitVector[]
        seen = Set{BitVector}()
        
        for term in min_cover
            combo = falses(nuberofatoms(formula))
            for (i, val) in enumerate(term)
                if val == 1
                    combo[i] = true
                end
            end
            
            if combo ∉ seen
                push!(seen, combo)
                push!(result_combinations, combo)
            end
        end
        
        sort!(result_combinations)
        return TwoLevelDNFFormula(
            result_combinations,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            collect(min_cover)
        )
    catch e
        @warn "Error during minimization: $e"
        return formula
    end
end


#=GLOBAL MINIMUMS BUT INCORRECT=#
function minimizza_dnf(::Val{:quine_oldstyle}, formula::TwoLevelDNFFormula)
    # Initialization of starting terms
    terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end
    
    if isempty(terms)
        return formula  # If there are no terms, return the empty formula directly
    end

    # Function to check if one cube covers another
    function covers(cube1::Vector{Int8}, cube2::Vector{Int8})
        for i in eachindex(cube1)
            if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                return false
            end
        end
        return true
    end

    # Function to combine two cubes into a new generalized cube
    function combine_cubes(cube1::Vector{Int8}, cube2::Vector{Int8})
        result = Vector{Int8}(undef, length(cube1))
        for i in eachindex(cube1)
            result[i] = ifelse(cube1[i] == cube2[i], cube1[i], Int8(-1))
        end
        return result
    end

    # Minimization phase: term combination
    function minimize_step!(terms::Vector{Vector{Int8}})
        if length(terms) <= 1
            return terms
        end
        
        new_terms = Vector{Vector{Int8}}()
        used = fill(false, length(terms))
        
        # Iterate over all term pairs for combinations
        for i in 1:(length(terms)-1)
            for j in (i+1):length(terms)
                diff_count = 0
                for k in eachindex(terms[i])
                    if terms[i][k] != terms[j][k] && terms[i][k] != Int8(-1) && terms[j][k] != Int8(-1)
                        diff_count += 1
                        if diff_count > 1
                            break
                        end
                    end
                end
                # Combine terms only if they differ in a single position
                if diff_count == 1
                    push!(new_terms, combine_cubes(terms[i], terms[j]))
                    used[i] = true
                    used[j] = true
                end
            end
        end
        
        # Add unused terms
        for i in 1:length(terms)
            if !used[i]
                push!(new_terms, terms[i])
            end
        end
        
        # Remove duplicates and sort
        unique!(new_terms)
        
        if length(new_terms) == length(terms)
            return new_terms  # No new combinations, terminate
        end
        
        return minimize_step!(new_terms)  # Continue iterating
    end

    # Find the minimal coverage of original terms using prime implicants
    function find_minimal_coverage(terms::Vector{Vector{Int8}}, primes::Vector{Vector{Int8}})
        if isempty(primes)
            return terms
        end
        
        # Construct coverage matrix
        coverage = falses(length(primes), length(terms))
        for i in 1:length(primes)
            for j in 1:length(terms)
                coverage[i, j] = covers(primes[i], terms[j])
            end
        end
        
        # Find essential implicants
        selected_primes = Int[]
        essential_coverage = falses(size(coverage, 2))
        
        for j in 1:size(coverage, 2)
            covering_primes = findall(coverage[:, j])
            if length(covering_primes) == 1
                push!(selected_primes, covering_primes[1])
                essential_coverage .|= coverage[covering_primes[1], :]
            end
        end
        
        # If needed, select additional implicants using a greedy approach
        if !all(essential_coverage)
            remaining_terms = .!essential_coverage
            selected = Set(selected_primes)
            
            while any(remaining_terms)
                best_prime = -1
                max_coverage = 0
                
                for (i, prime) in enumerate(primes)
                    if i ∉ selected
                        coverage_count = count(j -> remaining_terms[j] && coverage[i, j], 1:size(coverage, 2))
                        if coverage_count > max_coverage
                            max_coverage = coverage_count
                            best_prime = i
                        end
                    end
                end
                
                if best_prime == -1 || max_coverage == 0
                    break
                end
                
                push!(selected, best_prime)
                remaining_terms .&= .!coverage[best_prime, :]
            end
            
            return primes[collect(selected)]
        end
        
        return primes[selected_primes]
    end

    try
        # Phase 1: Generate prime implicants
        prime_implicants = minimize_step!(terms)
        # Phase 2: Find minimal coverage of terms
        minimized_terms = find_minimal_coverage(terms, prime_implicants)
        
        # Convert minimized terms into readable combinations
        new_combinations = BitVector[]
        seen = Set{BitVector}()  # Set to avoid duplicates

        for term in minimized_terms
            combo = falses(nuberofatoms(formula))
            for (i, val) in enumerate(term)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            if combo ∉ seen
                push!(seen, combo)
                push!(new_combinations, combo)
            end
        end
        
        sort!(new_combinations)
        return TwoLevelDNFFormula(
            new_combinations,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            collect(minimized_terms),
        )
    catch e
        @warn "Error during minimization: $e"
        return formula  # Return the original formula in case of error
    end
end


 #CLASSIC ESPRESSO STABLE (LOCAL MINIMAL)
 function minimizza_dnf(::Val{:espresso}, formula::TwoLevelDNFFormula; minimization_method_kwargs...)
    # Convert TritVectors to the internal representation used by Espresso
    # We'll map: 1 -> 1, 0 -> 0, -1 -> 0 (since we only care about positive terms)
    terms = Vector{Vector{Int}}()
    for tritvec in eachcombination(formula)
        term = Vector{Int}()
        for i in 1:length(tritvec)
            val = tritvec[i]
            push!(term, val == 1 ? 1 : 0)
        end
        push!(terms, term)
    end

    function covers(cube1, cube2)
        for (b1, b2) in zip(cube1, cube2)
            if b1 != -1 && b2 != -1 && b1 != b2
                return false
            end
        end
        return true
    end

    function can_combine(cube1, cube2)
        diff_count = 0
        diff_pos = -1

        for (i, (b1, b2)) in enumerate(zip(cube1, cube2))
            if b1 != b2
                diff_count += 1
                diff_pos = i
                if diff_count > 1
                    return false, -1
                end
            end
        end

        return diff_count == 1, diff_pos
    end

    function combine_cubes(cube1, cube2, pos)
        result = copy(cube1)
        result[pos] = -1
        return result
    end

    function find_combinations(cubes)
        if isempty(cubes)
            return Vector{Vector{Int}}()
        end

        result = Set{Vector{Int}}()
        combined = Set{Vector{Int}}()
        current_cubes = copy(cubes)

        while true
            found_new = false
            for i = 1:length(current_cubes)
                for j = (i+1):length(current_cubes)
                    can_combine, pos = can_combine(current_cubes[i], current_cubes[j])
                    if can_combine
                        new_cube = combine_cubes(current_cubes[i], current_cubes[j], pos)
                        if new_cube ∉ result
                            push!(result, new_cube)
                            push!(combined, current_cubes[i])
                            push!(combined, current_cubes[j])
                            found_new = true
                        end
                    end
                end
            end

            for cube in current_cubes
                if cube ∉ combined
                    push!(result, cube)
                end
            end

            if !found_new || length(result) == 0
                break
            end

            current_cubes = collect(result)
            empty!(result)
            empty!(combined)
        end

        if isempty(result)
            return current_cubes
        end

        return collect(result)
    end

    function find_essential_cubes(cubes, terms)
        if isempty(cubes) || isempty(terms)
            return terms
        end

        essential = Vector{Vector{Int}}()
        remaining_terms = Set(terms)

        while !isempty(remaining_terms)
            best_cube = nothing
            max_coverage = 0

            for cube in cubes
                coverage = count(term -> covers(cube, term), collect(remaining_terms))
                if coverage > max_coverage
                    max_coverage = coverage
                    best_cube = cube
                end
            end

            if best_cube === nothing || max_coverage == 0
                append!(essential, collect(remaining_terms))
                break
            end

            push!(essential, best_cube)
            filter!(term -> !covers(best_cube, term), remaining_terms)
        end

        return essential
    end

    function espresso_minimize(terms; minimization_method_kwargs...)
        if isempty(terms)
            return Vector{Vector{Int}}()
        end

        combined_terms = find_combinations(terms)
        if isempty(combined_terms)
            return terms
        end

        essential = find_essential_cubes(combined_terms, terms)
        if isempty(essential)
            return terms
        end

        return essential
    end

    # Perform minimization
    minimized_terms = espresso_minimize(terms; minimization_method_kwargs...)

    # Convert result to TritVector
    new_combinations = TritVector[]
    seen = Set{TritVector}()  # Set to track already seen terms

    for term in minimized_terms
        # Create a new TritVector for this term
        trit_combo = TritVector(nuberofatoms(formula))
        
        for (i, val) in enumerate(term)
            if val == 1
                trit_combo[i] = 1
            elseif val == -1
                trit_combo[i] = -1
            else
                trit_combo[i] = 0
            end
        end

        # Add the term only if it hasn't been seen before
        if trit_combo ∉ seen
            push!(seen, trit_combo)
            push!(new_combinations, trit_combo)
        end
    end

    sort!(new_combinations)  # Assuming you have defined sorting for TritVector
    return TwoLevelDNFFormula(
        new_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
    )
end


#=
function minimizza_dnf(::Val{:espresso}, formula::TwoLevelDNFFormula; minimization_method_kwargs...)
    terms = [Vector{Int}([x ? 1 : 0 for x in term]) for term in eachcombination(formula)]

    function covers(cube1, cube2)
        for (b1, b2) in zip(cube1, cube2)
            if b1 != -1 && b2 != -1 && b1 != b2
                return false
            end
        end
        return true
    end

    function can_combine(cube1, cube2)
        diff_count = 0
        diff_pos = -1

        for (i, (b1, b2)) in enumerate(zip(cube1, cube2))
            if b1 != b2
                diff_count += 1
                diff_pos = i
                if diff_count > 1
                    return false, -1
                end
            end
        end

        return diff_count == 1, diff_pos
    end

    function combine_cubes(cube1, cube2, pos)
        result = copy(cube1)
        result[pos] = -1
        return result
    end

    function find_combinations(cubes)
        if isempty(cubes)
            return Vector{Vector{Int}}()
        end

        result = Set{Vector{Int}}()
        combined = Set{Vector{Int}}()
        current_cubes = copy(cubes)

        while true
            found_new = false
            for i = 1:length(current_cubes)
                for j = (i+1):length(current_cubes)
                    can_combine, pos =
                        can_combine(current_cubes[i], current_cubes[j])
                    if can_combine
                        new_cube = combine_cubes(current_cubes[i], current_cubes[j], pos)
                        if new_cube ∉ result
                            push!(result, new_cube)
                            push!(combined, current_cubes[i])
                            push!(combined, current_cubes[j])
                            found_new = true
                        end
                    end
                end
            end

            # Add non-combined cubes to the result
            for cube in current_cubes
                if cube ∉ combined
                    push!(result, cube)
                end
            end

            if !found_new || length(result) == 0
                break
            end

            current_cubes = collect(result)
            empty!(result)
            empty!(combined)
        end

        # If no combinations were found, return the original cubes
        if isempty(result)
            return current_cubes
        end

        return collect(result)
    end

    function find_essential_cubes(cubes, terms)
        if isempty(cubes) || isempty(terms)
            return terms  # Return original terms if there are no cubes
        end

        essential = Vector{Vector{Int}}()
        remaining_terms = Set(terms)

        while !isempty(remaining_terms)
            best_cube = nothing
            max_coverage = 0

            for cube in cubes
                coverage = count(term -> covers(cube, term), collect(remaining_terms))
                if coverage > max_coverage
                    max_coverage = coverage
                    best_cube = cube
                end
            end

            if best_cube === nothing || max_coverage == 0
                # If no more coverage is found, add remaining terms
                append!(essential, collect(remaining_terms))
                break
            end

            push!(essential, best_cube)
            filter!(term -> !covers(best_cube, term), remaining_terms)
        end

        return essential
    end

    function espresso_minimize(terms; minimization_method_kwargs...)
        if isempty(terms)
            return Vector{Vector{Int}}()
        end

        combined_terms = find_combinations(terms)
        if isempty(combined_terms)
            return terms  # Return original terms if no combinations exist
        end

        essential = find_essential_cubes(combined_terms, terms)
        if isempty(essential)
            return terms  # Return original terms if no essential terms exist
        end

        return essential
    end

    # Execute minimization
    minimized_terms = espresso_minimize(terms; minimization_method_kwargs...)

    # Convert result to BitVector and remove duplicates using unique!
    new_combinations = BitVector[]
    seen = Set{BitVector}()  # Set to track already seen terms

    for term in minimized_terms
        combo = falses(formula.num_atoms)
        for (i, val) in enumerate(term)
            if val == 1
                combo[i] = true
            end
        end

        # Add the term only if it hasn't been seen before
        if combo ∉ seen
            push!(seen, combo)
            push!(new_combinations, combo)
        end
    end

    sort!(new_combinations)
    return TwoLevelDNFFormula(
        new_combinations,
        formula.num_atoms,
        formula.thresholds_by_feature,
        eachatomsbyfeature(formula),
        collect(minimized_terms),
    )
end
=#

#=quine resolution...=#
function minimizza_dnf(::Val{:quine_strict}, formula::TwoLevelDNFFormula)
    # Initialize terms
    terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end
    
    if isempty(terms)
        return formula
    end

    # Function to check if a point is covered by a term
    function is_point_covered(point::Vector{Int8}, term::Vector{Int8})
        for i in eachindex(point)
            if term[i] != Int8(-1) && term[i] != point[i]
                return false
            end
        end
        return true
    end

    # Function to check if two terms can be combined
    # without introducing new coverage points
    function can_combine_safely(cube1::Vector{Int8}, cube2::Vector{Int8}, original_terms::Vector{Vector{Int8}})
        # First check if they differ in only one position
        diff_count = 0
        diff_pos = -1
        
        for i in eachindex(cube1)
            if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                diff_count += 1
                diff_pos = i
                if diff_count > 1
                    return false, -1
                end
            end
        end
        
        if diff_count != 1
            return false, -1
        end
        
        # Create the combined term for testing
        combined = copy(cube1)
        combined[diff_pos] = Int8(-1)
        
        # For each original term
        for original in original_terms
            # If the combined term covers this point
            if is_point_covered(original, combined)
                # Verify that it was covered by at least one of the original terms
                if !is_point_covered(original, cube1) && !is_point_covered(original, cube2)
                    return false, -1
                end
            end
        end
        
        return true, diff_pos
    end

    function minimize_step!(terms::Vector{Vector{Int8}}, original_terms::Vector{Vector{Int8}})
        if length(terms) <= 1
            return terms
        end
        
        new_terms = Vector{Vector{Int8}}()
        used = fill(false, length(terms))
        
        for i in 1:length(terms)-1
            for j in i+1:length(terms)
                can_combine, pos = can_combine_safely(terms[i], terms[j], original_terms)
                
                if can_combine
                    new_term = copy(terms[i])
                    new_term[pos] = Int8(-1)
                    push!(new_terms, new_term)
                    used[i] = used[j] = true
                end
            end
        end
        
        # Add unused terms
        for i in 1:length(terms)
            if !used[i]
                push!(new_terms, terms[i])
            end
        end
        
        unique!(new_terms)
        
        if length(new_terms) == length(terms)
            return new_terms
        end
        
        return minimize_step!(new_terms, original_terms)
    end

    function find_minimum_coverage(terms::Vector{Vector{Int8}}, primes::Vector{Vector{Int8}})
        if isempty(primes)
            return terms
        end
        
        coverage = falses(length(primes), length(terms))
        for i in 1:length(primes)
            for j in 1:length(terms)
                coverage[i,j] = is_point_covered(terms[j], primes[i])
            end
        end
        
        # Find essential implicants
        selected_primes = Int[]
        for j in 1:size(coverage, 2)
            covering_primes = findall(coverage[:, j])
            if length(covering_primes) == 1
                push!(selected_primes, covering_primes[1])
            end
        end
        
        if !isempty(selected_primes)
            covered_terms = vec(any(coverage[selected_primes, :], dims=1))
            if all(covered_terms)
                return primes[unique(selected_primes)]
            end
        end
        
        # Greedy approach for remaining terms
        selected = Set(selected_primes)
        uncovered = Set(1:size(coverage, 2))
        
        for i in selected_primes
            filter!(t -> !coverage[i, t], uncovered)
        end
        
        while !isempty(uncovered)
            best_coverage = 0
            best_prime = 0
            
            for (i, prime) in enumerate(primes)
                if i ∉ selected
                    coverage_count = count(j -> j in uncovered && coverage[i, j], 1:size(coverage, 2))
                    if coverage_count > best_coverage
                        best_coverage = coverage_count
                        best_prime = i
                    end
                end
            end
            
            if best_coverage == 0
                break
            end
            
            push!(selected, best_prime)
            filter!(j -> !coverage[best_prime, j], uncovered)
        end
        
        return primes[collect(selected)]
    end

    try
        # Keep a copy of the original terms for checking
        original_terms = copy(terms)
        primes = minimize_step!(terms, original_terms)
        minimized = find_minimum_coverage(original_terms, primes)
        
        result = Vector{BitVector}(undef, length(minimized))
        for (i, term) in enumerate(minimized)
            combo = falses(numberofatoms(formula))
            for (j, val) in enumerate(term)
                if val == Int8(1)
                    combo[j] = true
                end
            end
            result[i] = combo
        end
        
        return TwoLevelDNFFormula(
            result,
            numberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            minimized
        )
        
    catch e
        @warn "Error during minimization: $e"
        @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
        return formula
    end
end

#resolution test#
function minimizza_dnf(::Val{:quine_petrick}, formula::TwoLevelDNFFormula)
    terms = eachcombination(formula)

    if isempty(terms)
        return formula
    end

    # Find prime implicants
    function find_prime_implicants(terms::Vector{BitVector})
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
            result[pos] = -1  # Add a don't care
            return result
        end

        primes = Set{Vector{Int}}()
        used = Set{Vector{Int}}()
        current = [Vector{Int}([x ? 1 : 0 for x in term]) for term in terms]

        while !isempty(current)
            next_terms = Vector{Vector{Int}}()

            for i in 1:length(current)
                for j in (i+1):length(current)
                    combinable, pos = can_combine(current[i], current[j])
                    if combinable
                        push!(used, current[i])
                        push!(used, current[j])
                        push!(next_terms, combine(current[i], current[j], pos))
                    end
                end
            end

            for term in current
                if term ∉ used
                    push!(primes, term)
                end
            end

            current = unique(next_terms)
            empty!(used)
        end

        return collect(primes)
    end

    # Check if a term covers another term
    function covers(prime::Vector{Int}, term::BitVector)
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != (term[i] ? 1 : 0)
                return false
            end
        end
        return true
    end

    # Petrick's method to minimize coverage
    function petrick_method(terms::Vector{BitVector}, primes::Vector{Vector{Int}})
        equation = Vector{Vector{Int}}()
        for term in terms
            row = Vector{Int}()
            for (idx, prime) in enumerate(primes)
                if covers(prime, term)
                    push!(row, idx)
                end
            end
            push!(equation, row)
        end

        function combine_equations(equation::Vector{Vector{Int}})
            if isempty(equation)
                return [[]]
            end

            first_row = first(equation)
            rest_combinations = combine_equations(equation[2:end])
            return [vcat([x], combo) for x in first_row for combo in rest_combinations]
        end

        solutions = combine_equations(equation)

        min_solution = solutions[1]
        for sol in solutions
            if length(unique(sol)) < length(unique(min_solution))
                min_solution = sol
            end
        end

        return [primes[idx] for idx in unique(min_solution)]
    end

    primes = find_prime_implicants(terms)
    minimized_cover = petrick_method(terms, primes)

    # Create a new list of minimized combinations
    minimized_combinations = BitVector[]
    for prime in minimized_cover
        bitvec = falses(numberofatoms(formula))
        for (i, val) in enumerate(prime)
            if val == 1
                bitvec[i] = true
            end
        end
        push!(minimized_combinations, bitvec)
    end

    return TwoLevelDNFFormula(
        minimized_combinations,
        numberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        minimized_cover
    )
end

# quine fuso with espresso [DEFINITIVE !?!?!?!1?!?!?!1?!1?!]
function minimizza_dnf(::Val{:pina}, formula::TwoLevelDNFFormula)
    let spinner_state = Ref(1)  # Keep track of spinner state between calls
        function show_progress(desc::String, current::Int, total::Int)
            spinner = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
            percentage = round(Int, current * 100 / total)
            
            # Update spinner state
            spinner_state[] = spinner_state[] % length(spinner) + 1
            
            # Use spinner_state instead of current for spinner animation
            print("\r$desc: $(spinner[spinner_state[]]) $percentage%")
            flush(stdout)  # Ensure output is displayed immediately
            
            current == total && println()
        end

        function modified_espresso(terms)
            function can_combine(cube1, cube2)
                diff_count = 0
                diff_pos = -1
                
                for (i, (b1, b2)) in enumerate(zip(cube1, cube2))
                    if b1 != -1 && b2 != -1 && b1 != b2
                        diff_count += 1
                        diff_pos = i
                        if diff_count > 1
                            return false, -1
                        end
                    end
                end
                
                return diff_count == 1, diff_pos
            end

            function combine_cubes(cube1, cube2, pos)
                result = copy(cube1)
                result[pos] = -1
                return result
            end

            function find_combinations_for_quine(terms)
                result = Vector{Vector{Int8}}()
                combined = Set{Tuple{Vector{Int8}, Vector{Int8}}}()
                iteration = 1
                
                while true
                    found_new = false
                    total_combinations = div(length(terms) * (length(terms) - 1), 2)
                    progress_counter = 0
                    
                    println("\nEspresso - Iteration $iteration")
                    
                    for i in 1:length(terms)
                        for j in (i+1):length(terms)
                            progress_counter += 1
                            if progress_counter % 1000 == 0
                                show_progress("Combinations", progress_counter, total_combinations)
                            end
                            
                            if (terms[i], terms[j]) ∉ combined
                                can_combine, pos = can_combine(terms[i], terms[j])
                                if can_combine
                                    new_cube = combine_cubes(terms[i], terms[j], pos)
                                    push!(result, new_cube)
                                    push!(combined, (terms[i], terms[j]))
                                    found_new = true
                                end
                            end
                        end
                    end
                    
                    if !found_new
                        break
                    end
                    
                # Add new terms to terms for the next round
                    append!(terms, result)
                    empty!(result)
                    iteration += 1
                end
                
            # Filter redundant terms
                return unique(terms)
            end

            return find_combinations_for_quine(terms)
        end

    # Step 2: Modified version of Quine that works on the results of modified_espresso
        function modified_quine(terms)
            function covers(cube1::Vector{Int8}, cube2::Vector{Int8})
                for i in eachindex(cube1)
                    if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                        return false
                    end
                end
                return true
            end

            function find_essential_coverage(terms, primes)
                essential = Vector{Vector{Int8}}()
                covered = Set{Vector{Int8}}()
                
                println("\nQuine - Essential implicants phase")
                for (i, term) in enumerate(terms)
                    if i % 1000 == 0
                        show_progress("Term analysis", i, length(terms))
                    end
                    covering_primes = filter(p -> covers(p, term), primes)
                    if length(covering_primes) == 1
                        push!(essential, first(covering_primes))
                        push!(covered, term)
                    end
                end
                
            # Then use a greedy approach for the remaining
                remaining_terms = setdiff(Set(terms), covered)
                progress_counter = 0
                
                println("\nQuine - Greedy phase")
                while !isempty(remaining_terms)
                    progress_counter += 1
                    if progress_counter % 1000 == 0
                        show_progress("Optimization", progress_counter, length(terms))
                    end
                    
                    best_prime = nothing
                    max_coverage = 0
                    
                    for prime in primes
                        coverage = count(term -> term ∈ remaining_terms && covers(prime, term), terms)
                        if coverage > max_coverage
                            max_coverage = coverage
                            best_prime = prime
                        end
                    end
                    
                    if best_prime === nothing || max_coverage == 0
                        break
                    end
                    
                    push!(essential, best_prime)
                    filter!(term -> !covers(best_prime, term), remaining_terms)
                end
                
                return essential
            end

            primes = terms
            return find_essential_coverage(terms, primes)
        end

        println("\nInitialization")
        initial_terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
        for (i, term) in enumerate(eachcombination(formula))
            if i % 1000 == 0
                show_progress("Converting terms", i, length(eachcombination(formula)))
            end
            for (j, x) in enumerate(term)
                initial_terms[i][j] = x ? Int8(1) : Int8(0)
            end
        end
        
        if isempty(initial_terms)
            return formula
        end

        try
            println("\nStarting DNF minimization")
            
            espresso_terms = modified_espresso(initial_terms)
            final_terms = modified_quine(espresso_terms)
            
            println("\nFinalization")
            result_combinations = Vector{BitVector}()
            for (i, term) in enumerate(final_terms)
                if i % 1000 == 0
                    show_progress("Converting results", i, length(final_terms))
                end
                combo = falses(nuberofatoms(formula))
                for (i, val) in enumerate(term)
                    if val == Int8(1)
                        combo[i] = true
                    end
                end
                push!(result_combinations, combo)
            end
            
            println("\nMinimization completed!")
            
            return TwoLevelDNFFormula(
                result_combinations,
                nuberofatoms(formula),
                eachthresholdsbyfeature(formula),
                eachatomsbyfeature(formula),
                final_terms
            )
            
        catch e
            @warn "Error during minimization: $e"
            @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
            return formula
        end
    end
end

function minimizza_dnf(::Val{:pina_interrupted}, formula::TwoLevelDNFFormula)
    let spinner_state = Ref(1)  # Spinner state
        function show_progress(desc::String, current::Int, total::Int)
            spinner = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
            percentage = round(Int, current * 100 / total)

            # Update spinner state
            spinner_state[] = spinner_state[] % length(spinner) + 1

            print("\r$desc: $(spinner[spinner_state[]]) $percentage%")
            flush(stdout)  # Show output immediately

            current == total && println()
        end

        function can_combine(cube1, cube2)
            diff_count = 0
            diff_pos = -1

            for (i, (b1, b2)) in enumerate(zip(cube1, cube2))
                if b1 != -1 && b2 != -1 && b1 != b2
                    diff_count += 1
                    diff_pos = i
                    if diff_count > 1
                        return false, -1
                    end
                end
            end

            return diff_count == 1, diff_pos
        end

        function combine_cubes(cube1, cube2, pos)
            result = copy(cube1)
            result[pos] = -1
            return result
        end

        function find_combinations_for_quine(terms)
            result = Vector{Vector{Int8}}()
            iteration = 1

            while true
                found_new = false
                combined = Set{Tuple{Vector{Int8}, Vector{Int8}}}()

                println("\nEspresso - Iteration $iteration")

                for i in 1:length(terms)
                    for j in (i + 1):length(terms)
                        if (terms[i], terms[j]) ∉ combined
                            can_comb, pos = can_combine(terms[i], terms[j])
                            if can_comb
                                new_cube = combine_cubes(terms[i], terms[j], pos)
                                push!(result, new_cube)
                                push!(combined, (terms[i], terms[j]))
                                found_new = true
                            end
                        end
                    end
                end

                if !found_new
                    break
                end

                append!(terms, result)
                terms = unique(sort(terms))
                empty!(result)
                iteration += 1
            end

            return terms
        end

        function covers(cube1::Vector{Int8}, cube2::Vector{Int8})
            for i in eachindex(cube1)
                if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                    return false
                end
            end
            return true
        end

        function find_essential_coverage(terms, primes)
            essential = Vector{Vector{Int8}}()
            covered = Set{Vector{Int8}}()

            println("\nQuine - Essential implicants phase")
            for (i, term) in enumerate(terms)
                if i % 1000 == 0
                    show_progress("Analyzing terms", i, length(terms))
                end
                covering_primes = filter(p -> covers(p, term), primes)
                if length(covering_primes) == 1
                    push!(essential, first(covering_primes))
                    push!(covered, term)
                end
            end

            remaining_terms = setdiff(Set(terms), covered)
            println("\nQuine - Greedy phase")

            while !isempty(remaining_terms)
                best_prime = nothing
                max_coverage = 0

                for prime in primes
                    coverage = count(term -> term ∈ remaining_terms && covers(prime, term), terms)
                    if coverage > max_coverage
                        max_coverage = coverage
                        best_prime = prime
                    end
                end

                if best_prime === nothing || max_coverage == 0
                    break
                end

                push!(essential, best_prime)
                filter!(term -> !covers(best_prime, term), remaining_terms)
            end

            return essential
        end

        println("\nInitialization")
        initial_terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
        for (i, term) in enumerate(eachcombination(formula))
            if i % 1000 == 0
                show_progress("Converting terms", i, length(eachcombination(formula)))
            end
            for (j, x) in enumerate(term)
                initial_terms[i][j] = x ? Int8(1) : Int8(0)
            end
        end

        if isempty(initial_terms)
            return formula
        end

        try
            println("\nStarting DNF minimization")

            espresso_terms = find_combinations_for_quine(initial_terms)
            final_terms = find_essential_coverage(initial_terms, espresso_terms)

            println("\nFinalization")
            result_combinations = Vector{BitVector}()
            for (i, term) in enumerate(final_terms)
                if i % 1000 == 0
                    show_progress("Converting results", i, length(final_terms))
                end
                combo = falses(nuberofatoms(formula))
                for (i, val) in enumerate(term)
                    if val == Int8(1)
                        combo[i] = true
                    end
                end
                push!(result_combinations, combo)
            end

            println("\nMinimization completed!")

            return TwoLevelDNFFormula(
                result_combinations,
                nuberofatoms(formula),
                eachthresholdsbyfeature(formula),
                eachatomsbyfeature(formula),
                final_terms
            )

        catch e
            @warn "Error during minimization: $e"
            @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
            return formula
        end
    end
end

# WORKING with same results as espresso
function minimizza_dnf(::Val{:espresso_quine}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Original terms in 0/1 format to guarantee the same coverage
    ########################################################################
    original_terms_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimization with Espresso
    ########################################################################
    espresso_result = minimizza_dnf(Val(:espresso), formula)

    ########################################################################
    # 3) Convert Espresso results into cubes with -1 (if prime_mask not empty)
    ########################################################################
    espresso_terms = if isempty(espresso_result.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_result)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_result.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combination for prime implicants
    ########################################################################
    function quine_combination(implicants::Vector{Vector{Int8}})
        if length(implicants) <= 1
            return implicants
        end

        function can_combine(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # If both are not -1 and different, count difference
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combine(c1, c2)
            result = copy(c1)
            for i in eachindex(c1)
                if c1[i] == c2[i]
                    result[i] = c1[i]
                else
                    result[i] = Int8(-1)
                end
            end
            return result
        end

        new_ones = Vector{Vector{Int8}}()
        used = fill(false, length(implicants))

        # Try to combine each pair
        for i in 1:length(implicants)-1
            for j in i+1:length(implicants)
                if can_combine(implicants[i], implicants[j])
                    push!(new_ones, combine(implicants[i], implicants[j]))
                    used[i] = true
                    used[j] = true
                end
            end
        end

        # Add uncombined ones
        for i in 1:length(implicants)
            if !used[i]
                push!(new_ones, implicants[i])
            end
        end

        unique!(new_ones)

        # If no more reduction, we have the prime implicants
        if length(new_ones) == length(implicants)
            return new_ones
        else
            return quine_combination(new_ones)
        end
    end

    ########################################################################
    # 5) Quine: minimum coverage selection (preserves coverage on original_terms_bit)
    ########################################################################
    function quine_minimum_coverage(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        # Returns true if cube1 covers cube2
        function covers(cube1, cube2)
            for i in eachindex(cube1)
                if cube1[i] != Int8(-1) && cube1[i] != cube2[i]
                    return false
                end
            end
            return true
        end

        coverage = falses(length(prime_impls), length(orig_terms))
        for i in 1:length(prime_impls)
            for j in 1:length(orig_terms)
                coverage[i, j] = covers(prime_impls[i], orig_terms[j])
            end
        end

        # Find essential prime implicants
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        covered = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    covered[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), covered)
        while !isempty(uncovered)
            best_prime = 0
            best_score = 0
            for i in 1:length(prime_impls)
                if i ∉ selected
                    c = count(j -> coverage[i, j], uncovered)
                    if c > best_score
                        best_score = c
                        best_prime = i
                    end
                end
            end
            if best_prime == 0
                break
            end
            push!(selected, best_prime)
            for j in uncovered
                if coverage[best_prime, j]
                    covered[j] = true
                end
            end
            uncovered = findall(!isequal(true), covered)
        end

        return prime_impls[collect(selected)]
    end

    # Get prime implicants and final coverage
    prime_implicants = quine_combination(espresso_terms)
    final_coverage = quine_minimum_coverage(original_terms_bit, prime_implicants)

    ########################################################################
    # 6) Final step: "pushed" unification on any redundant features
    ########################################################################
    function final_unify_cubes!(implicants::Vector{Vector{Int8}})
        changed = true
        while changed
            changed = false
            newset = Set{Vector{Int8}}()
            used = fill(false, length(implicants))

            for i in 1:length(implicants)-1
                if used[i]
                    continue
                end
                for j in i+1:length(implicants)
                    if used[j]
                        continue
                    end
                    new_cube = unify_entire_feature(implicants[i], implicants[j], formula)
                    if new_cube !== nothing
                        push!(newset, new_cube)
                        used[i] = true
                        used[j] = true
                        changed = true
                        break
                    end
                end
                if !used[i]
                    push!(newset, implicants[i])
                end
            end
            if !used[end]
                push!(newset, implicants[end])
            end

            implicants = collect(newset)
            unique!(implicants)
        end
        return implicants
    end

    function unify_entire_feature(c1::Vector{Int8}, c2::Vector{Int8}, formula::TwoLevelDNFFormula)
        diffpos = Int[]
        for k in eachindex(c1)
            if c1[k] != c2[k]
                push!(diffpos, k)
            end
        end
        if isempty(diffpos)
            return nothing
        end

        feats = Set{Int}()
        for pos in diffpos
            (f, _thr, _op) = find_feature_op_threshold(pos, formula)
            push!(feats, f)
            if length(feats) > 1
                return nothing
            end
        end

        f = first(feats)
        posf = positions_of_feature(f, formula)

        all_ones = trues(length(posf))
        for (k, p) in enumerate(posf)
            if !(c1[p] == Int8(1) || c2[p] == Int8(1))
                all_ones[k] = false
            end
        end

        if all(all_ones)
            newcube = copy(c1)
            for p in posf
                newcube[p] = Int8(-1)
            end
            return newcube
        end
        return nothing
    end

    function positions_of_feature(feat::Int, formula::TwoLevelDNFFormula)
        out = Int[]
        offset = 0
        for (f, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for _ in atoms
                offset += 1
                if f == feat
                    push!(out, offset)
                end
            end
        end
        return out
    end

    function find_feature_op_threshold(idx::Int, formula::TwoLevelDNFFormula)
        offset = 0
        for (feat, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: index $idx out of range")
    end

    # Final unification, if you want it
    final_coverage_unified = final_unify_cubes!(final_coverage)

    ########################################################################
    # 7) Building final formula
    ########################################################################
    new_combinations = Vector{BitVector}()
    for cube in final_coverage_unified
        combo = falses(nuberofatoms(formula))
        for (i, val) in enumerate(cube)
            if val == Int8(1)
                combo[i] = true
            end
        end
        push!(new_combinations, combo)
    end
    unique!(new_combinations)

    return TwoLevelDNFFormula(
        new_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        final_coverage_unified,
    )
end

# AI GENERATED (TODO TEST it)
function minimizza_dnf(::Val{:espresso_quine_2}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Original terms in 0/1 format to guarantee the same coverage
    ########################################################################
    original_terms_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimization with Espresso
    ########################################################################
    espresso_result = minimizza_dnf(Val(:espresso), formula)

    ########################################################################
    # 3) Convert Espresso results into cubes with -1 (if prime_mask not empty)
    ########################################################################
    espresso_terms = if isempty(espresso_result.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_result)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_result.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combination for prime implicants
    ########################################################################
    function quine_combination(implicants::Vector{Vector{Int8}})
        if length(implicants) <= 1
            return implicants
        end

        function can_combine(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # If both are not -1 and different, count difference
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combine(c1, c2)
            result = copy(c1)
            for i in eachindex(c1)
                if c1[i] == c2[i]
                    result[i] = c1[i]
                else
                    result[i] = Int8(-1)
                end
            end
            return result
        end

        new_ones = Vector{Vector{Int8}}()
        used = fill(false, length(implicants))

        # Try to combine each pair
        for i in 1:length(implicants)-1
            for j in i+1:length(implicants)
                if can_combine(implicants[i], implicants[j])
                    push!(new_ones, combine(implicants[i], implicants[j]))
                    used[i] = true
                    used[j] = true
                end
            end
        end

        # Add uncombined ones
        for i in 1:length(implicants)
            if !used[i]
                push!(new_ones, implicants[i])
            end
        end

        unique!(new_ones)

        # If no more reduction, we have the prime implicants
        if length(new_ones) == length(implicants)
            return new_ones
        else
            return quine_combination(new_ones)
        end
    end

    ########################################################################
    # 5) Quine: minimum coverage selection (preserves coverage on original_terms_bit)
    ########################################################################
    function quine_minimum_coverage(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        # Returns true if cube1 covers cube2
        function covers(cube1, cube2)
            for i in eachindex(cube1)
                if cube1[i] != Int8(-1) && cube1[i] != cube2[i]
                    return false
                end
            end
            return true
        end

        coverage = falses(length(prime_impls), length(orig_terms))
        for i in 1:length(prime_impls)
            for j in 1:length(orig_terms)
                coverage[i, j] = covers(prime_impls[i], orig_terms[j])
            end
        end

        # Find essential prime implicants
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        covered = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    covered[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), covered)
        while !isempty(uncovered)
            best_prime = 0
            best_score = 0
            for i in 1:length(prime_impls)
                if i ∉ selected
                    c = count(j -> coverage[i, j], uncovered)
                    if c > best_score
                        best_score = c
                        best_prime = i
                    end
                end
            end
            if best_prime == 0
                break
            end
            push!(selected, best_prime)
            for j in uncovered
                if coverage[best_prime, j]
                    covered[j] = true
                end
            end
            uncovered = findall(!isequal(true), covered)
        end

        return prime_impls[collect(selected)]
    end

    ########################################################################
    # 6) Final step: "pushed" unification on any redundant features
    ########################################################################
    # Here we want to detect situations where, FOR EXAMPLE,
    #   (V3 < 5.0) and (V3 ≥ 4.95) cover all of V3 => unify into "don't care on V3".
    #
    # We do this in a "final_unify_cubes!" function that repeatedly
    # searches for pairs of unifiable cubes and replaces them.
    ########################################################################

    function final_unify_cubes!(implicants::Vector{Vector{Int8}})
        changed = true
        while changed
            changed = false
            newset = Set{Vector{Int8}}()
            used = fill(false, length(implicants))

            # Try to unify pairs
            for i in 1:length(implicants)-1
                if used[i]
                    continue
                end
                for j in i+1:length(implicants)
                    if used[j]
                        continue
                    end
                    new_cube = unify_entire_feature(implicants[i], implicants[j], formula)
                    if new_cube !== nothing
                        # unified
                        push!(newset, new_cube)
                        used[i] = true
                        used[j] = true
                        changed = true
                        break
                    end
                end
                if !used[i]
                    push!(newset, implicants[i])
                end
            end

            # last possible unused one
            if !used[end]
                push!(newset, implicants[end])
            end

            implicants = collect(newset)
            unique!(implicants)
        end
        return implicants
    end

    # unify_entire_feature(c1, c2, formula):
    #   - If c1 and c2 are identical on ALL features except one
    #     and on THAT feature they differ in a way that covers the entire domain,
    #     then returns a cube that has -1 on that feature. Otherwise nothing.
    function unify_entire_feature(
        c1::Vector{Int8},
        c2::Vector{Int8},
        formula::TwoLevelDNFFormula
    )
        # 1) Understand how many features the two cubes differ on.
        #    If they differ on more than 1 feature, => nothing.
        #    If they differ on 1 single feature, check if "union" of atoms = domain.
        #
        # To do this, identify positions where c1 != c2
        # and then check if they ALL correspond to a single feature.
        diffpos = Int[]
        for k in eachindex(c1)
            if c1[k] != c2[k]
                push!(diffpos, k)
            end
        end
        if isempty(diffpos)
            # c1 == c2 => no further union
            return nothing
        end

        # 2) Discover if diffpos all belong to the SAME feature
        feats = Set{Int}()
        for pos in diffpos
            (f, _thr, _op) = find_feature_op_threshold(pos, formula)
            push!(feats, f)
            if length(feats) > 1
                return nothing  # they differ on more than 1 feature => don't unify
            end
        end

        # Now feats = {f} = a single feature.
        # 3) Check if the sub-atoms on which c1 and c2 diverge cover the entire
        #    domain of that feature f. That is:
        #    - Take ALL atoms of feature f = setpos = [positions...].
        #    - c1 and c2, on those atoms, must "complement" each other and unite all bits = 0/1
        #      in order to make the entire range. If yes, => put -1 on ALL bits of f.
        f = first(feats)  # the feature they diverge on
        posf = positions_of_feature(f, formula)  # all bit indices that correspond to f

        # Build the "bitmask" union on those posf: (c1[i] or c2[i]) => if it's 1 on all => complementary
        # But an or isn't enough. We want to say: do c1 and c2, united, set ALL atoms of f?
        # Or, if "≥4.95" and "<5.0" are only 2 atoms out of 2 total => then yes, they cover everything.
        #
        # To simplify, use the logic: "If the size of posf = number of atoms on f",
        # and on that set (c1[i], c2[i]) take values such that the OR = 1 on all posf,
        # then that set of atoms = feature f in its entirety. For the example (≥4.95, <5.0) it's 2 atoms total.
        all_ones = trues(length(posf))
        for (k, p) in enumerate(posf)
            # if c1[p] or c2[p] == 1 => union => 1
            if !(c1[p] == Int8(1) || c2[p] == Int8(1))
                all_ones[k] = false
            end
        end

        # If the OR is 1 on ALL bits of posf => we're covering the entire domain of f
        if all(all_ones)
            # => unify => cube with -1 in ALL posf
            newcube = copy(c1)
            # Copy also bits in c2 where c1 had -1 (or vice versa),
            # since we'll put -1 on the ENTIRE feature f anyway.
            for p in posf
                newcube[p] = Int8(-1)
            end
            return newcube
        end

        return nothing
    end

    # Find all bit indices (1-based) that correspond to a given feature.
    function positions_of_feature(feat::Int, formula::TwoLevelDNFFormula)
        out = Int[]
        offset = 0
        # Traverse formula.atoms_by_feature in order, as in find_feature_op_threshold
        for (f, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for _ in atoms
                offset += 1
                if f == feat
                    push!(out, offset)
                end
            end
        end
        return out
    end

    # Returns (feature, threshold, op) for a bit index `idx` (1-based).
    function find_feature_op_threshold(idx::Int, formula::TwoLevelDNFFormula)
        offset = 0
        for (feat, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: index $idx out of range")
    end

    ########################################################################
    # Apply final unification
    ########################################################################
    final_coverage_unified = final_unify_cubes!(final_coverage)

    ########################################################################
    # 7) Convert final cubes to BitVector and build final formula
    ########################################################################
    new_combinations = Vector{BitVector}()
    for cube in final_coverage_unified
        combo = falses(nuberofatoms(formula))
        for (i, val) in enumerate(cube)
            if val == Int8(1)
                combo[i] = true
            end
        end
        push!(new_combinations, combo)
    end
    unique!(new_combinations)

    # Final return
    return TwoLevelDNFFormula(
        new_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        final_coverage_unified,
    )
end

# AI GENERATED (TODO TEST it) (RAM HUNGRY)
function minimizza_dnf(::Val{:espresso_quine_pp}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Original terms in 0/1 format to guarantee the same coverage
    ########################################################################
    original_terms_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimization with Espresso (initial step)
    ########################################################################
    espresso_result = minimizza_dnf(Val(:espresso), formula)
    println(length(eachcombination(espresso_result.combination)))

    ########################################################################
    # 3) Convert Espresso results into cubes with -1 (if prime_mask is not empty)
    ########################################################################
    espresso_terms = if isempty(espresso_result.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_result)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_result.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combination for prime implicants
    ########################################################################
    function quine_combination(implicants::Vector{Vector{Int8}})

        # ----------------------------------------------------
        # AVOID combining if it produces a cube all -1
        # ----------------------------------------------------
        function can_combine(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # If both are not -1 and different, count difference
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combine(c1, c2)
            result = copy(c1)
            for i in eachindex(c1)
                if c1[i] == c2[i]
                    result[i] = c1[i]
                else
                    result[i] = Int8(-1)
                end
            end
            return result
        end

        function is_all_dontcare(cube::Vector{Int8})
            return all(x -> x == Int8(-1), cube)
        end

        new_ones = Vector{Vector{Int8}}()
        used = fill(false, length(implicants))

        # Try to combine each pair
        for i in 1:length(implicants)-1
            for j in i+1:length(implicants)
                if can_combine(implicants[i], implicants[j])
                    new_cube = combine(implicants[i], implicants[j])

                    # BLOCK cube all -1
                    if is_all_dontcare(new_cube)
                        # if you don't want to admit the tautology,
                        # skip insertion
                        continue
                    end

                    push!(new_ones, new_cube)
                    used[i] = true
                    used[j] = true
                end
            end
        end

        # Add the uncombined ones
        for i in 1:length(implicants)
            if !used[i]
                push!(new_ones, implicants[i])
            end
        end

        unique!(new_ones)

        # If it doesn't reduce more, we have the prime implicants
        if length(new_ones) == length(implicants)
            return new_ones
        else
            return quine_combination(new_ones)
        end
    end

    ########################################################################
    # 5) Quine: minimum coverage selection (preserves coverage on original_terms_bit)
    ########################################################################
    function quine_minimum_coverage(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        function covers(cube1, cube2)
            # Returns true if cube1 covers cube2
            for i in eachindex(cube1)
                if cube1[i] != Int8(-1) && cube1[i] != cube2[i]
                    return false
                end
            end
            return true
        end

        coverage = falses(length(prime_impls), length(orig_terms))
        for i in 1:length(prime_impls)
            for j in 1:length(orig_terms)
                coverage[i, j] = covers(prime_impls[i], orig_terms[j])
            end
        end

        # Find essential prime implicants
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        covered = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    covered[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), covered)
        while !isempty(uncovered)
            best_prime = 0
            best_score = 0
            for i in 1:length(prime_impls)
                if i ∉ selected
                    c = count(j -> coverage[i, j], uncovered)
                    if c > best_score
                        best_score = c
                        best_prime = i
                    end
                end
            end
            if best_prime == 0
                break
            end
            push!(selected, best_prime)
            for j in uncovered
                if coverage[best_prime, j]
                    covered[j] = true
                end
            end
            uncovered = findall(!isequal(true), covered)
        end

        return prime_impls[collect(selected)]
    end

    ########################################################################
    # 6) Final step: "pushed" unification on any redundant features
    ########################################################################
    function final_unify_cubes!(implicants::Vector{Vector{Int8}})
        changed = true
        while changed
            changed = false
            newset = Set{Vector{Int8}}()
            used = fill(false, length(implicants))

            for i in 1:length(implicants)-1
                if used[i]
                    continue
                end
                for j in i+1:length(implicants)
                    if used[j]
                        continue
                    end
                    new_cube = unify_entire_feature(implicants[i], implicants[j], formula)
                    if new_cube !== nothing
                        # BLOCK tautology: if new_cube is completely -1, discard it
                        if all(x -> x == Int8(-1), new_cube)
                            continue
                        end
                        push!(newset, new_cube)
                        used[i] = true
                        used[j] = true
                        changed = true
                        break
                    end
                end
                if !used[i]
                    push!(newset, implicants[i])
                end
            end
            if !used[end]
                push!(newset, implicants[end])
            end

            implicants = collect(newset)
            unique!(implicants)
        end
        return implicants
    end

    function unify_entire_feature(
        c1::Vector{Int8}, 
        c2::Vector{Int8}, 
        formula::TwoLevelDNFFormula
    )
        diffpos = Int[]
        for k in eachindex(c1)
            if c1[k] != c2[k]
                push!(diffpos, k)
            end
        end
        if isempty(diffpos)
            return nothing
        end

        feats = Set{Int}()
        for pos in diffpos
            (f, _thr, _op) = find_feature_op_threshold(pos, formula)
            push!(feats, f)
            if length(feats) > 1
                return nothing
            end
        end

        # At this point c1 and c2 differ only on a single feature
        f = first(feats)
        posf = positions_of_feature(f, formula)

        all_ones = trues(length(posf))
        for (k, p) in enumerate(posf)
            if !(c1[p] == Int8(1) || c2[p] == Int8(1))
                all_ones[k] = false
            end
        end

        if all(all_ones)
            newcube = copy(c1)
            for p in posf
                newcube[p] = Int8(-1)
            end
            return newcube
        end
        return nothing
    end

    function positions_of_feature(feat::Int, formula::TwoLevelDNFFormula)
        out = Int[]
        offset = 0
        for (f, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for _ in atoms
                offset += 1
                if f == feat
                    push!(out, offset)
                end
            end
        end
        return out
    end

    function find_feature_op_threshold(idx::Int, formula::TwoLevelDNFFormula)
        offset = 0
        for (feat, atoms) in sort(collect(eachatomsbyfeature(formula)))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: index $idx out of range")
    end

    ############ MAIN FLOW BEGINNING ############

    # 1) Get the basic prime implicants
    prime_implicants = quine_combination(espresso_terms)

    # 2) Find minimum coverage
    final_coverage = quine_minimum_coverage(original_terms_bit, prime_implicants)

    # 3) Final unification (pushed)
    final_coverage_unified = final_unify_cubes!(final_coverage)

    # 4) POSSIBLE check for "residual tautology"
    # If there is exactly 1 cube and it's all -1, discard it:
    if length(final_coverage_unified) == 1 && all(x -> x == Int8(-1), final_coverage_unified[1])
        @warn "minimizza_dnf: found 1 completely -1 cube (tautology). Forcibly removing it!"
        final_coverage_unified = Vector{Vector{Int8}}()  # e.g. empty formula ("false")
        # Or you could go back to prime implicants, etc. Depends on your logic.
    end

    ########################################################################
    # 7) Final formula construction
    ########################################################################
    new_combinations = Vector{BitVector}()
    for cube in final_coverage_unified
        combo = falses(nuberofatoms(formula))
        for (i, val) in enumerate(cube)
            if val == Int8(1)
                combo[i] = true
            end
        end
        push!(new_combinations, combo)
    end
    unique!(new_combinations)

    return TwoLevelDNFFormula(
        new_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        final_coverage_unified,
    )
end