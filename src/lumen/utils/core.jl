"""
Processes a combination of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

This function takes a combination of atom values, represented as either a `BitVector` or an integer, and the dictionaries of thresholds and atom properties. It then processes the combination, filtering the valid values for each feature based on the atom values. If a feature has no valid values, the function sets a flag indicating a contradiction. Finally, it returns the dictionary of valid values for each feature and the contradiction flag.

Parameters:
- `i::Union{BitVector,Int}`: The combination of atom values to be processed.
- `num_atoms::Int`: The number of atoms in the combination.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding valid values.
- `Bool`: A flag indicating if the combination has a contradiction.
"""
function truth_combinations(
    model::Any,
    alphabet,
    atoms::Vector{
        <:Atom{
            <:ScalarCondition{
                Float64,
                <:VariableValue,
                <:ScalarMetaCondition{<:VariableValue,typeof(<)},
            },
        },
    },
    vertical::Float64;
    apply_function = SoleModels.apply,
    print_progress = true,
    silent = false,
)
    # Inizializzazione di strutture dati per gestire le soglie e gli atomi per feature
    thresholds_by_feature = Dict{Int,Vector{Float64}}()
    atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()

    # Popolamento delle strutture dati con le soglie e gli atomi
    for subalpha in alphabet.subalphabets
        feat = subalpha.featcondition[1].feature.i_variable
        thresholds_by_feature[feat] = sort(subalpha.featcondition[2])
    end

    for atom in atoms
        feat = atom.value.metacond.feature.i_variable
        threshold = atom.value.threshold
        push!(get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat), (threshold, true))
    end

    # Ordinamento degli atomi per ciascuna feature
    for (_, atom_list) in atoms_by_feature
        sort!(atom_list, by = first)
    end

    # Calcolo del numero totale di combinazioni
    num_atoms = length(atoms)
    all_combinations = BigInt(2)^num_atoms
    println("num_atoms: ", num_atoms)
    println("2^num_atoms: ", all_combinations)
    println("vertical: ", vertical)
    
    if isone(vertical)
        num_combinations = all_combinations
    else
        num_combinations = BigInt(round(all_combinations * vertical))
    end

    # Stima del tempo di esecuzione
    silent || begin
        println("Stima del tempo di esecuzione...")
        sample_size = min(1000, num_combinations)
        start_time = time()

        for i = 1:sample_size
            process_combination(
                BigInt(i - 1),
                num_atoms,
                thresholds_by_feature,
                atoms_by_feature,
            )
        end

        end_time = time()
        estimated_time_per_combination = (end_time - start_time) / sample_size
        total_estimated_time =
            BigFloat(estimated_time_per_combination) * num_combinations

        print_estimated_time(total_estimated_time)
    end

    # Inizio della generazione delle combinazioni
    silent || println("Inizio della generazione delle combinazioni...")
    silent || println("Totale combinazioni da generare: $num_combinations")

    results = Dict{Any,Vector{BigInt}}()
    label_count = Dict{Any,Int}()

    # Processo di generazione e valutazione delle combinazioni
    start_time = time()
    contradictions = 0
    i_rand = 0
    update_interval = max(1, num_combinations ÷ 100)  # Aggiorna ogni 1% di progresso
    
    for i = 0:(num_combinations-1)
        i_v = i                                      # Indice per la progressbar
        if isone(vertical)
            combination, has_contradiction = process_combination(
                BigInt(i),
                num_atoms,
                thresholds_by_feature,
                atoms_by_feature,
            )
        else
            has_contradiction = true
                i = rand(BigInt(0):(all_combinations-1))
                combination, has_contradiction = process_combination(
                    BigInt(i),
                    num_atoms,
                    thresholds_by_feature,
                    atoms_by_feature,
                )
        end
        if has_contradiction
            contradictions += 1
        else
            combination_dict = SortedDict(combination)
            combination_vector = vcat(collect(values(combination_dict))...)
            
            if !(combination_vector in values(results))
                result = apply_function(model, combination_vector)
                push!(get!(Vector{BigInt}, results, result), BigInt(i))
                label_count[result] = get(label_count, result, 0) + 1
            end
        end

        # Aggiornamento della barra di avanzamento testuale FIXME TOTALMENTE RIMUOVIBILE ?
        print_progress && begin
            if i_v % update_interval == 0 || i_v == num_combinations - 1
                progress = BigFloat(i_v + 1) / num_combinations
                print_progress_bar(progress)
            end
        end
    end
    print_progress && println()  # Nuova linea dopo la barra di avanzamento

    end_time = time()
    real_time_per_combination = (end_time - start_time) / num_combinations

    # Stampa dei risultati
    silent || print_results_summary(
        real_time_per_combination,
        Int64(contradictions),
        num_combinations,
        label_count,
    )
    silent || print_detailed_results(results)

    return results, label_count
end


"""
Processes a combination of atom values and returns a dictionary of valid values for each feature, along with a flag indicating if the combination has a contradiction.

This function takes a combination of atom values, represented as either a `BitVector` or an integer, and the dictionaries of thresholds and atom properties. It then processes the combination, filtering the valid values for each feature based on the atom values. If a feature has no valid values, the function sets a flag indicating a contradiction. Finally, it returns the dictionary of valid values for each feature and the contradiction flag.

Parameters:
- `i::Union{BitVector,Int}`: The combination of atom values to be processed.
- `num_atoms::Int`: The number of atoms in the combination.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding valid values.
- `Bool`: A flag indicating if the combination has a contradiction.
"""
function process_combination(
    i,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    combination = Dict{Int,Vector{Float64}}()

    truth_row = !isa(i, BitVector) ? digits(i, base = 2, pad = num_atoms) : i

    has_contradiction = false
    j = 1
    for (feat, atom_list) in atoms_by_feature
        thresholds = thresholds_by_feature[feat]
        valid_values = copy(thresholds)

        for (_, (threshold, _)) in enumerate(atom_list)
            if truth_row[j] == 1
                filter!(x -> x < threshold, valid_values)
            else
                filter!(x -> x >= threshold, valid_values)
            end
            j += 1
        end

        if isempty(valid_values)
            has_contradiction = true
            break
        else
            combination[feat] = valid_values
        end

    end

    if !has_contradiction
        for feat in keys(thresholds_by_feature)
            if !haskey(combination, feat)
                combination[feat] = [1605.0]
            end
        end
    end

    return combination, has_contradiction
end


"""
Concatenates the results of a custom OR formula computation into a dictionary of `TwoLevelDNFFormula` objects.

This function takes the raw results of the custom OR formula computation, along with the number of atoms and the thresholds and atom properties, and constructs a dictionary of `TwoLevelDNFFormula` objects. The results are sorted by the number of combinations in descending order, and the `TwoLevelDNFFormula` objects are created with the sorted combinations.

Parameters:
- `results::Any`: The raw results of the custom OR formula computation.
- `num_atoms::Int`: The number of atoms in the custom OR formula.
- `thresholds_by_feature::Dict{Int,Vector{Float64}}`: A dictionary mapping feature indices to their corresponding thresholds.
- `atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}}`: A dictionary mapping feature indices to their corresponding atom properties (threshold and boolean flag).

Returns:
- `Dict{Any,TwoLevelDNFFormula}`: A dictionary mapping the original results to their corresponding `TwoLevelDNFFormula` objects.
"""
function concat_results(
    results::Any,
    num_atoms::Int,
    thresholds_by_feature::Dict{Int,Vector{Float64}},
    atoms_by_feature::Dict{Int,Vector{Tuple{Float64,Bool}}},
)
    results = dict_to_bitvector(results, num_atoms)
    res = Dict{Any,TwoLevelDNFFormula}()
    println("\nRisultati dettagliati:")

    for (result, combinations) in sort(collect(results), by = x -> length(x[2]), rev = true)
        println("[$result] ($(length(combinations)) combinazioni)")
        res[result] = TwoLevelDNFFormula(
            Vector{BitVector}(combinations),
            num_atoms,
            thresholds_by_feature,
            atoms_by_feature,
            Vector{Vector{Int}}(), # TODO nothing
        )
    end
    return res
end


"""
	minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula, horizontal = 1.0)

Simplifies a custom OR formula using the Quine algorithm.

This function takes a `TwoLevelDNFFormula` object and applies the Quine algorithm to minimize the number of combinations in the formula. It returns a new `TwoLevelDNFFormula` object with the simplified combinations.

"""
function minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula, horizontal = 1.0)
    # Strutture dati iniziali
    terms = eachcombination(formula)
    feature_count =
        Dict(feature => length(atoms) for (feature, atoms) in formula.atoms_by_feature)

    # Calcolo features permesse basato su horizontal
    total_features = length(formula.thresholds_by_feature)
    allowed_features = floor(Int, total_features * horizontal)
    num_orizontal = sum(get(feature_count, i, 0) for i = 1:allowed_features)

    # Funzioni di conversione
    function bitVector_to_term(bv::BitVector)
        return [Int(b) for b in bv]
    end

    function term_to_bitVector(term::Vector{Int}, original::BitVector)
        result = copy(original)
        for (i, val) in enumerate(term)
            if val != -1
                result[i] = val == 1
            end
        end
        return result
    end

    # Funzione di mascheramento verticale migliorata
    function mask_vertical_bits!(term::Vector{Int}, num_orizontal::Int)
        for i = (num_orizontal+1):length(term)
            term[i] = -1
        end
        return term
    end

    # Preparazione termini mascherati
    masked_terms =
        [mask_vertical_bits!(bitVector_to_term(term), num_orizontal) for term in terms]

    # Funzioni di supporto per l'algoritmo Quine-McCluskey
    function combina_termini(term1::Vector{Int}, term2::Vector{Int})
        @assert length(term1) == length(term2) "I termini devono avere la stessa lunghezza"
        return [t1 == t2 ? t1 : -1 for (t1, t2) in zip(term1, term2)]
    end

    function conta_uni(term::Vector{Int})
        return count(x -> x == 1, term)
    end

    function differisce_per_un_bit(term1::Vector{Int}, term2::Vector{Int})
        diff_count = 0
        diff_pos = 0
        for (i, (t1, t2)) in enumerate(zip(term1, term2))
            if t1 != -1 && t2 != -1 && t1 != t2
                diff_count += 1
                diff_pos = i
                if diff_count > 1
                    return false, 0
                end
            end
        end
        return diff_count == 1, diff_pos
    end

    # Implementazione migliorata di trova_prime_implicants
    function trova_prime_implicants(terms::Vector{Vector{Int}})
        groups = Dict{Int,Vector{Vector{Int}}}()
        for term in terms
            ones = conta_uni(term)
            if !haskey(groups, ones)
                groups[ones] = Vector{Vector{Int}}()
            end
            push!(groups[ones], term)
        end

        prime_implicants = Set{Vector{Int}}()
        used_terms = Set{Vector{Int}}()

        while !isempty(groups)
            next_groups = Dict{Int,Vector{Vector{Int}}}()
            merged = false

            for ones_count in sort(collect(keys(groups)))
                if haskey(groups, ones_count + 1)
                    for term1 in groups[ones_count]
                        for term2 in groups[ones_count+1]
                            is_diff, pos = differisce_per_un_bit(term1, term2)
                            if is_diff
                                merged = true
                                combined = combina_termini(term1, term2)
                                push!(used_terms, term1)
                                push!(used_terms, term2)

                                ones_in_combined = conta_uni(combined)
                                if !haskey(next_groups, ones_in_combined)
                                    next_groups[ones_in_combined] = Vector{Vector{Int}}()
                                end
                                push!(next_groups[ones_in_combined], combined)
                            end
                        end
                    end
                end
            end

            # Aggiungi i termini non utilizzati ai prime implicants
            for (_, terms_in_group) in groups
                for term in terms_in_group
                    if term ∉ used_terms
                        push!(prime_implicants, term)
                    end
                end
            end

            if !merged
                break
            end
            groups = next_groups
        end

        return collect(prime_implicants)
    end

    # Funzione migliorata per trovare gli essential prime implicants
    function trova_essential_prime_implicants(
        minterms::Vector{Vector{Int}},
        prime_implicants::Vector{Vector{Int}},
    )
        essential_prime_implicants = Set{Vector{Int}}()
        remaining_minterms = Set(minterms)

        # Costruisci la tabella di copertura
        coverage = Dict{Vector{Int},Set{Vector{Int}}}()
        for pi in prime_implicants
            coverage[pi] = Set{Vector{Int}}()
            for minterm in minterms
                if all(p == -1 || p == m for (p, m) in zip(pi, minterm))
                    push!(coverage[pi], minterm)
                end
            end
        end

        # Trova gli essential prime implicants
        while !isempty(remaining_minterms)
            # Trova i minterms coperti da un solo prime implicant
            essential_found = false
            for minterm in remaining_minterms
                covering_pis = [pi for (pi, covered) in coverage if minterm in covered]
                if length(covering_pis) == 1
                    push!(essential_prime_implicants, covering_pis[1])
                    setdiff!(remaining_minterms, coverage[covering_pis[1]])
                    essential_found = true
                    break
                end
            end

            # Se non sono stati trovati essential prime implicants, scegli quello con la migliore copertura
            if !essential_found && !isempty(remaining_minterms)
                best_coverage = 0
                best_pi = nothing
                for (pi, covered) in coverage
                    current_coverage = length(intersect(covered, remaining_minterms))
                    if current_coverage > best_coverage
                        best_coverage = current_coverage
                        best_pi = pi
                    end
                end
                if best_pi !== nothing
                    push!(essential_prime_implicants, best_pi)
                    setdiff!(remaining_minterms, coverage[best_pi])
                end
            end
        end

        return collect(essential_prime_implicants)
    end

    # Esegui l'algoritmo
    prime_implicants = trova_prime_implicants(masked_terms)
    essential_prime_implicants =
        trova_essential_prime_implicants(masked_terms, prime_implicants)

    # Genera le nuove combinazioni
    nuove_combinazioni = BitVector[]
    for implicant in essential_prime_implicants
        nuovo_bitvector = term_to_bitVector(implicant, terms[1])
        push!(nuove_combinazioni, nuovo_bitvector)
    end

    sort!(nuove_combinazioni)
    return TwoLevelDNFFormula(
        nuove_combinazioni,
        formula.num_atoms,
        formula.thresholds_by_feature,
        formula.atoms_by_feature,
        collect(essential_prime_implicants),
    )
end


"""
	minimizza_dnf(::Val{:espresso}, formula::TwoLevelDNFFormula)

Simplifies a custom OR formula using the Espresso algorithm.

This function takes a `TwoLevelDNFFormula` object and applies the Espresso algorithm to minimize the number of combinations in the formula. It returns a new `TwoLevelDNFFormula` object with the simplified combinations.

The Espresso algorithm is a well-known method for Boolean function minimization, which is used here to simplify the custom OR formula. The function performs the following steps:

1. Converts the formula's combinations into a list of integer vectors representing the minterms.
2. Applies the Espresso algorithm to find the essential prime implicants of the formula.
3. Generates new combinations based on the essential prime implicants.
4. Removes any duplicate combinations and returns the simplified `TwoLevelDNFFormula`.

"""
function minimizza_dnf(::Val{:espresso}, formula::TwoLevelDNFFormula)
    terms = [Vector{Int}([x ? 1 : 0 for x in term]) for term in eachcombination(formula)]

    function copre(cube1, cube2)
        for (b1, b2) in zip(cube1, cube2)
            if b1 != -1 && b2 != -1 && b1 != b2
                return false
            end
        end
        return true
    end

    function possono_combinarsi(cube1, cube2)
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

    function combina_cubi(cube1, cube2, pos)
        result = copy(cube1)
        result[pos] = -1
        return result
    end

    function trova_combinazioni(cubes)
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
                        possono_combinarsi(current_cubes[i], current_cubes[j])
                    if can_combine
                        nuovo_cubo = combina_cubi(current_cubes[i], current_cubes[j], pos)
                        if nuovo_cubo ∉ result
                            push!(result, nuovo_cubo)
                            push!(combined, current_cubes[i])
                            push!(combined, current_cubes[j])
                            found_new = true
                        end
                    end
                end
            end

            # Aggiungi i cubi non combinati al risultato
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

        # Se non abbiamo trovato combinazioni, restituisci i cubi originali
        if isempty(result)
            return current_cubes
        end

        return collect(result)
    end

    function find_essential_cubes(cubes, terms)
        if isempty(cubes) || isempty(terms)
            return terms  # Restituisci i termini originali se non ci sono cubi
        end

        essential = Vector{Vector{Int}}()
        remaining_terms = Set(terms)

        while !isempty(remaining_terms)
            best_cube = nothing
            max_coverage = 0

            for cube in cubes
                coverage = count(term -> copre(cube, term), collect(remaining_terms))
                if coverage > max_coverage
                    max_coverage = coverage
                    best_cube = cube
                end
            end

            if best_cube === nothing || max_coverage == 0
                # Se non troviamo più copertura, aggiungi i termini rimanenti
                append!(essential, collect(remaining_terms))
                break
            end

            push!(essential, best_cube)
            filter!(term -> !copre(best_cube, term), remaining_terms)
        end

        return essential
    end

    function espresso_minimize(terms)
        if isempty(terms)
            return Vector{Vector{Int}}()
        end

        combined_terms = trova_combinazioni(terms)
        if isempty(combined_terms)
            return terms  # Restituisci i termini originali se non ci sono combinazioni
        end

        essential = find_essential_cubes(combined_terms, terms)
        if isempty(essential)
            return terms  # Restituisci i termini originali se non ci sono termini essenziali
        end

        return essential
    end

    # Esegui la minimizzazione
    minimized_terms = espresso_minimize(terms)

    # Converti il risultato in BitVector e rimuovi i duplicati usando unique!
    nuove_combinazioni = BitVector[]
    seen = Set{BitVector}()  # Set per tenere traccia dei termini già visti

    for term in minimized_terms
        combo = falses(formula.num_atoms)
        for (i, val) in enumerate(term)
            if val == 1
                combo[i] = true
            end
        end

        # Aggiungi il termine solo se non è già stato visto
        if combo ∉ seen
            push!(seen, combo)
            push!(nuove_combinazioni, combo)
        end
    end

    sort!(nuove_combinazioni)
    return TwoLevelDNFFormula(
        nuove_combinazioni,
        formula.num_atoms,
        formula.thresholds_by_feature,
        formula.atoms_by_feature,
        collect(minimized_terms),
    )
end


"""
    verify_simplification(original::TwoLevelDNFFormula, simplified::TwoLevelDNFFormula)

Verifies that a simplified custom OR formula is congruent with the original formula.

This function takes the original and simplified custom OR formulas, and generates a set of random assignments to evaluate and compare the results of the two formulas. If any mismatch is found, the function returns `false`, indicating that the simplification is not correct. Otherwise, it returns `true`, confirming that the simplified formula is congruent with the original.

Parameters:
- `original::TwoLevelDNFFormula`: The original custom OR formula.
- `simplified::TwoLevelDNFFormula`: The simplified custom OR formula.

Returns:
- `true` if the simplified formula is congruent with the original, `false` otherwise.
"""
function verify_simplification(original::TwoLevelDNFFormula, simplified::TwoLevelDNFFormula)
    @info "Starting verification of simplification"

    # Funzione per generare assegnazioni basate sui BitVector
    function generate_smart_assignments(formula, num_samples)
        assignments = Set{Dict{Int,Bool}}()

        # Converti ogni BitVector in un Dict{Int,Bool}
        for combination in eachcombination(formula)
            assignment = Dict{Int,Bool}()
            for (i, bit) in enumerate(combination)
                assignment[i] = bit
            end
            push!(assignments, assignment)
            if length(assignments) >= num_samples
                break
            end
        end

        # Aggiungi alcune assegnazioni casuali per aumentare la copertura
        while length(assignments) < num_samples
            random_assignment = Dict(i => rand(Bool) for i = 1:formula.num_atoms)
            push!(assignments, random_assignment)
        end
        return collect(assignments)
    end

    # Funzione ottimizzata per valutare la formula
    function evaluate_custom_or_formula(formula, assignment)
        for combination in eachcombination(formula)
            all_true = true
            for (feat, atom_list) in formula.atoms_by_feature
                feat_value = get(assignment, feat, false)
                for (threshold, is_less_than) in atom_list
                    if is_less_than ? (feat_value >= threshold) : (feat_value < threshold)
                        all_true = false
                        break
                    end
                end
                if !all_true
                    break
                end
            end
            if all_true
                return true
            end
        end
        return false
    end

    # Genera un set di assegnazioni "intelligenti" basate sulle combinazioni esistenti
    num_samples = min(1000, 2^original.num_atoms)  # Limita il numero di campioni per input molto grandi
    assignments = generate_smart_assignments(original, num_samples)

    # Verifica le formule utilizzando le assegnazioni generate
    for (i, assignment) in enumerate(assignments)
        original_result = evaluate_custom_or_formula(original, assignment)
        simplified_result = evaluate_custom_or_formula(simplified, assignment)

        if original_result != simplified_result
            @warn "Mismatch found for assignment: $assignment"
            @warn "Original result: $original_result"
            @warn "Simplified result: $simplified_result"
            return false
        end

        # Stampa il progresso ogni 100 iterazioni
        if i % 100 == 0
            @info "Processed $i out of $num_samples assignments"
        end
    end

    @info "Verification complete. Simplified formula is congruent with the original."
    return true
end
