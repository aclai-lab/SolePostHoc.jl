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
# CLASSIC QUINE STABLE (RAM HUNGRY)
function minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula)
    # Inizializzazione dei termini
    terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end
    
    if isempty(terms)
        return formula
    end

    # Funzione per controllare se un cubo ne copre un altro
    function copre(cube1::Vector{Int8}, cube2::Vector{Int8})
        for i in eachindex(cube1)
            if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                return false
            end
        end
        return true
    end

    # Funzione per combinare due cubi
    function combina_cubi!(result::Vector{Int8}, cube1::Vector{Int8}, cube2::Vector{Int8}, pos::Int)
        for i in eachindex(cube1)
            result[i] = (i == pos) ? Int8(-1) : cube1[i]
        end
    end

    function minimizza_step!(termini::Vector{Vector{Int8}})
        if length(termini) <= 1
            return termini
        end
        
        nuovi_termini = Vector{Vector{Int8}}()
        usati = fill(false, length(termini))
        temp_cube = Vector{Int8}(undef, length(termini[1]))
        
        # Prima fase: cerca coppie di termini che differiscono per una sola posizione
        for i in 1:length(termini)-1
            for j in i+1:length(termini)
                diff_count = 0
                diff_pos = -1
                
                for k in eachindex(termini[i])
                    if termini[i][k] != Int8(-1) && termini[j][k] != Int8(-1)
                        if termini[i][k] != termini[j][k]
                            diff_count += 1
                            diff_pos = k
                            if diff_count > 1
                                break
                            end
                        end
                    end
                end
                
                if diff_count == 1
                    combina_cubi!(temp_cube, termini[i], termini[j], diff_pos)
                    push!(nuovi_termini, copy(temp_cube))
                    usati[i] = usati[j] = true
                end
            end
        end
        
        # Aggiungi i termini non usati
        for i in 1:length(termini)
            if !usati[i]
                push!(nuovi_termini, termini[i])
            end
        end
        
        unique!(nuovi_termini)
        
        if length(nuovi_termini) == length(termini)
            return nuovi_termini
        end
        
        return minimizza_step!(nuovi_termini)
    end

    function trova_copertura_minima(termini::Vector{Vector{Int8}}, primi::Vector{Vector{Int8}})
        if isempty(primi)
            return termini
        end
        
        # Costruisci matrice di copertura
        coverage = falses(length(primi), length(termini))
        for i in 1:length(primi)
            for j in 1:length(termini)
                coverage[i,j] = copre(primi[i], termini[j])
            end
        end
        
        # Trova implicanti essenziali
        selected_primes = Int[]
        for j in 1:size(coverage, 2)
            covering_primes = findall(coverage[:, j])
            if length(covering_primes) == 1
                push!(selected_primes, covering_primes[1])
            end
        end
        
        # Se gli implicanti essenziali coprono tutto, siamo finiti
        if !isempty(selected_primes)
            covered_terms = vec(any(coverage[selected_primes, :], dims=1))
            if all(covered_terms)
                return primi[unique(selected_primes)]
            end
        end
        
        # Altrimenti, usa un approccio greedy per coprire i termini rimanenti
        selected = Set(selected_primes)
        uncovered = Set(1:size(coverage, 2))
        
        # Rimuovi i termini già coperti
        for i in selected_primes
            filter!(t -> !coverage[i, t], uncovered)
        end
        
        while !isempty(uncovered)
            best_coverage = 0
            best_prime = 0
            
            for (i, prime) in enumerate(primi)
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
        
        return primi[collect(selected)]
    end

    try
        primi = minimizza_step!(terms)
        minimized = trova_copertura_minima(terms, primi)
        
        # Converti il risultato in BitVectors
        risultato = Vector{BitVector}(undef, length(minimized))
        for (i, term) in enumerate(minimized)
            combo = falses(formula.num_atoms)
            for (j, val) in enumerate(term)
                if val == Int8(1)
                    combo[j] = true
                end
            end
            risultato[i] = combo
        end
        
        return TwoLevelDNFFormula(
            risultato,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
            minimized
        )
        
    catch e
        @warn "Errore durante la minimizzazione: $e"
        @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
        return formula
    end
end

#CLASSIC ESPRESSO STABLE (LOCAL MINIMAL)
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

function minimizza_dnf(::Val{:pina}, formula::TwoLevelDNFFormula)
    # Step 1: Versione modificata di Espresso
    function modified_espresso(terms)
        function possono_combinarsi(cube1, cube2)
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

        function combina_cubi(cube1, cube2, pos)
            result = copy(cube1)
            result[pos] = -1
            return result
        end

        function trova_combinazioni_per_quine(terms)
            result = Vector{Vector{Int8}}()
            combined = Set{Tuple{Vector{Int8}, Vector{Int8}}}()
            
            while true
                found_new = false
                
                for i in 1:length(terms)
                    for j in (i+1):length(terms)
                        if (terms[i], terms[j]) ∉ combined
                            can_combine, pos = possono_combinarsi(terms[i], terms[j])
                            if can_combine
                                new_cube = combina_cubi(terms[i], terms[j], pos)
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
                
                # Aggiungi i nuovi termini a terms per il prossimo round
                append!(terms, result)
                empty!(result)
            end
            
            # Filtra i termini ridondanti
            return unique(terms)
        end

        return trova_combinazioni_per_quine(terms)
    end

    # Step 2: Versione modificata di Quine che lavora sui risultati di modified_espresso
    function modified_quine(terms)
        function copre(cube1::Vector{Int8}, cube2::Vector{Int8})
            for i in eachindex(cube1)
                if cube1[i] != Int8(-1) && cube2[i] != Int8(-1) && cube1[i] != cube2[i]
                    return false
                end
            end
            return true
        end

        function trova_copertura_essenziale(terms, primes)
            essential = Vector{Vector{Int8}}()
            covered = Set{Vector{Int8}}()
            
            # Prima trova gli implicanti essenziali
            for term in terms
                covering_primes = filter(p -> copre(p, term), primes)
                if length(covering_primes) == 1
                    push!(essential, first(covering_primes))
                    push!(covered, term)
                end
            end
            
            # Poi usa un approccio greedy per i rimanenti
            remaining_terms = setdiff(Set(terms), covered)
            while !isempty(remaining_terms)
                best_prime = nothing
                max_coverage = 0
                
                for prime in primes
                    coverage = count(term -> term ∈ remaining_terms && copre(prime, term), terms)
                    if coverage > max_coverage
                        max_coverage = coverage
                        best_prime = prime
                    end
                end
                
                if best_prime === nothing || max_coverage == 0
                    break
                end
                
                push!(essential, best_prime)
                filter!(term -> !copre(best_prime, term), remaining_terms)
            end
            
            return essential
        end

        primes = terms
        return trova_copertura_essenziale(terms, primes)
    end

    # Converti la formula iniziale in termini Int8
    initial_terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            initial_terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end
    
    if isempty(initial_terms)
        return formula
    end

    try
        # Applica prima Espresso modificato
        espresso_terms = modified_espresso(initial_terms)
        
        # Poi applica Quine modificato
        final_terms = modified_quine(espresso_terms)
        
        # Converti il risultato finale in BitVectors
        result_combinations = Vector{BitVector}()
        for term in final_terms
            combo = falses(formula.num_atoms)
            for (i, val) in enumerate(term)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            push!(result_combinations, combo)
        end
        
        return TwoLevelDNFFormula(
            result_combinations,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
            final_terms
        )
        
    catch e
        @warn "Errore durante la minimizzazione: $e"
        @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
        return formula
    end
end

# FUNZIONANTE? stessi risultati di espresso
function minimizza_dnf(::Val{:espresso_quine}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Termini originali in formato 0/1 per garantire la stessa copertura
    ########################################################################
    termini_originali_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimizzazione con Espresso
    ########################################################################
    espresso_risultato = minimizza_dnf(Val(:espresso), formula)

    ########################################################################
    # 3) Converti i risultati di Espresso in cubi con -1 (se prime_mask non vuota)
    ########################################################################
    termini_espresso = if isempty(espresso_risultato.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_risultato)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_risultato.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combinazione per prime implicants
    ########################################################################
    function quine_combinazione(implicanti::Vector{Vector{Int8}})
        if length(implicanti) <= 1
            return implicanti
        end

        function possono_combinarsi(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # Se entrambi non sono -1 e diversi, contiamo differenza
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combina(c1, c2)
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

        nuovi = Vector{Vector{Int8}}()
        usati = fill(false, length(implicanti))

        # Tenta di combinare ogni paio
        for i in 1:length(implicanti)-1
            for j in i+1:length(implicanti)
                if possono_combinarsi(implicanti[i], implicanti[j])
                    push!(nuovi, combina(implicanti[i], implicanti[j]))
                    usati[i] = true
                    usati[j] = true
                end
            end
        end

        # Aggiungi i non combinati
        for i in 1:length(implicanti)
            if !usati[i]
                push!(nuovi, implicanti[i])
            end
        end

        unique!(nuovi)

        # Se non si riduce più, abbiamo i prime implicants
        if length(nuovi) == length(implicanti)
            return nuovi
        else
            return quine_combinazione(nuovi)
        end
    end

    ########################################################################
    # 5) Quine: selezione copertura minima (conserve la copertura su termini_originali_bit)
    ########################################################################
    function quine_copertura_minima(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        # Ritorna true se cube1 copre cube2
        function copre(cube1, cube2)
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
                coverage[i, j] = copre(prime_impls[i], orig_terms[j])
            end
        end

        # Trova prime implicants essenziali
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        coperti = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    coperti[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), coperti)
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
                    coperti[j] = true
                end
            end
            uncovered = findall(!isequal(true), coperti)
        end

        return prime_impls[collect(selected)]
    end

    # Otteniamo prime implicants e copertura finale
    prime_implicants = quine_combinazione(termini_espresso)
    copertura_finale = quine_copertura_minima(termini_originali_bit, prime_implicants)

    ########################################################################
    # 6) Passaggio finale: unificazione “spinta” su eventuali feature ridondanti
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
        for (f, atoms) in sort(collect(formula.atoms_by_feature))
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
        for (feat, atoms) in sort(collect(formula.atoms_by_feature))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: indice $idx fuori range")
    end

    # Unificazione finale, se la vuoi
        copertura_finale_unified = final_unify_cubes!(copertura_finale)

        ########################################################################
        # 7) Costruzione formula finale
        ########################################################################
        nuove_combinazioni = Vector{BitVector}()
        for cube in copertura_finale_unified
            combo = falses(formula.num_atoms)
            for (i, val) in enumerate(cube)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            push!(nuove_combinazioni, combo)
        end
        unique!(nuove_combinazioni)

        return TwoLevelDNFFormula(
            nuove_combinazioni,
            formula.num_atoms,
            formula.thresholds_by_feature,
            formula.atoms_by_feature,
            copertura_finale_unified,
        )
end

function minimizza_dnf(::Val{:espresso_quine_2}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Termini originali in formato 0/1 per garantire la stessa copertura
    ########################################################################
    termini_originali_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimizzazione con Espresso
    ########################################################################
    espresso_risultato = minimizza_dnf(Val(:espresso), formula)

    ########################################################################
    # 3) Converti i risultati di Espresso in cubi con -1 (se prime_mask non vuota)
    ########################################################################
    termini_espresso = if isempty(espresso_risultato.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_risultato)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_risultato.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combinazione per prime implicants
    ########################################################################
    function quine_combinazione(implicanti::Vector{Vector{Int8}})
        if length(implicanti) <= 1
            return implicanti
        end

        function possono_combinarsi(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # Se entrambi non sono -1 e diversi, contiamo differenza
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combina(c1, c2)
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

        nuovi = Vector{Vector{Int8}}()
        usati = fill(false, length(implicanti))

        # Tenta di combinare ogni paio
        for i in 1:length(implicanti)-1
            for j in i+1:length(implicanti)
                if possono_combinarsi(implicanti[i], implicanti[j])
                    push!(nuovi, combina(implicanti[i], implicanti[j]))
                    usati[i] = true
                    usati[j] = true
                end
            end
        end

        # Aggiungi i non combinati
        for i in 1:length(implicanti)
            if !usati[i]
                push!(nuovi, implicanti[i])
            end
        end

        unique!(nuovi)

        # Se non si riduce più, abbiamo i prime implicants
        if length(nuovi) == length(implicanti)
            return nuovi
        else
            return quine_combinazione(nuovi)
        end
    end

    ########################################################################
    # 5) Quine: selezione copertura minima (conserve la copertura su termini_originali_bit)
    ########################################################################
    function quine_copertura_minima(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        # Ritorna true se cube1 copre cube2
        function copre(cube1, cube2)
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
                coverage[i, j] = copre(prime_impls[i], orig_terms[j])
            end
        end

        # Trova prime implicants essenziali
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        coperti = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    coperti[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), coperti)
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
                    coperti[j] = true
                end
            end
            uncovered = findall(!isequal(true), coperti)
        end

        return prime_impls[collect(selected)]
    end

    # Otteniamo prime implicants e copertura finale
    prime_implicants = quine_combinazione(termini_espresso)
    copertura_finale = quine_copertura_minima(termini_originali_bit, prime_implicants)

    ########################################################################
    # 6) Passaggio finale: unificazione “spinta” su eventuali feature ridondanti
    ########################################################################
    # Qui vogliamo rilevare situazioni in cui, per ESEMPIO,
    #   (V3 < 5.0) e (V3 ≥ 4.95)  coprono tutto V3 => unisci in "don't care su V3".
    #
    # Lo facciamo in una funzione "final_unify_cubes!" che ripete più volte
    # la ricerca di coppie di cubi unificabili e le sostituisce.
    ########################################################################

    function final_unify_cubes!(implicants::Vector{Vector{Int8}})
        changed = true
        while changed
            changed = false
            newset = Set{Vector{Int8}}()
            used = fill(false, length(implicants))

            # Proviamo a unire le coppie
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
                        # unificati
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

            # ultimo eventuale non usato
            if !used[end]
                push!(newset, implicants[end])
            end

            implicants = collect(newset)
            unique!(implicants)
        end
        return implicants
    end

    # unify_entire_feature(c1, c2, formula):
    #   - Se c1 e c2 sono identici su TUTTE le feature tranne una
    #     e su QUELLA feature differiscono in modo da coprire l'intero dominio,
    #     allora restituisce un cubo che ha -1 su quella feature. Altrimenti nothing.
    function unify_entire_feature(
        c1::Vector{Int8},
        c2::Vector{Int8},
        formula::TwoLevelDNFFormula
    )
        # 1) Capire su quante feature differiscono i due cubi.
        #    Se differiscono su più di 1 feature, => nothing.
        #    Se differiscono su 1 sola feature, verifichiamo se "unione" di atomi = dominio.
        #
        # Per farlo, individuiamo le posizioni in cui c1 != c2
        # e poi controlliamo se corrispondono TUTTE a un'unica feature.
        diffpos = Int[]
        for k in eachindex(c1)
            if c1[k] != c2[k]
                push!(diffpos, k)
            end
        end
        if isempty(diffpos)
            # c1 == c2 => nessuna unione ulteriore
            return nothing
        end

        # 2) Scopriamo se diffpos tutte appartengono alla STESSA feature
        feats = Set{Int}()
        for pos in diffpos
            (f, _thr, _op) = find_feature_op_threshold(pos, formula)
            push!(feats, f)
            if length(feats) > 1
                return nothing  # differiscono su più di 1 feature => non unifichiamo
            end
        end

        # Ora feats = {f} = una sola feature.
        # 3) Verifichiamo se i sub-atomi su cui c1 e c2 divergono coprono tutto
        #    il dominio di quella feature f. Vale a dire:
        #    - Prendiamo TUTTI gli atomi di feature f = setpos = [positions...].
        #    - c1 e c2, su quei atomi, devono “complementarsi” e unire tutti i bit = 0/1
        #      in modo da fare l'intero range. Se sì, => mettiamo -1 su TUTTI i bit di f.
        f = first(feats)  # la feature su cui divergono
        posf = positions_of_feature(f, formula)  # tutti gli indici di bit che corrispondono a f

        # Costruiamo l'unione “bitmask” su quei posf: (c1[i] or c2[i]) => se è 1 su tutti => complementare
        # Ma in realtà non basta un or. Vogliamo dire: c1 e c2, uniti, settano TUTTI gli atomi di f?
        # Oppure, se “≥4.95” e “<5.0” sono solo 2 atomi su 2 totali => allora sì, coprono tutto.
        #
        # Per semplificare, usiamo la logica: “Se la dimensione di posf = numero di atomi su f”,
        # e su quell'insieme (c1[i], c2[i]) assumono valori tali da far sì che l'OR = 1 su tutti i posf,
        # allora quell'insieme di atomi = feature f al completo. Per l'esempio (≥4.95, <5.0) è 2 atomi totali.
        all_ones = trues(length(posf))
        for (k, p) in enumerate(posf)
            # se c1[p] o c2[p] == 1 => unione => 1
            if !(c1[p] == Int8(1) || c2[p] == Int8(1))
                all_ones[k] = false
            end
        end

        # Se l'OR è 1 su TUTTI i bit di posf => stiamo coprendo l'intero dominio di f
        if all(all_ones)
            # => unify => cubo con -1 in TUTTI i posf
            newcube = copy(c1)
            # Copiamo anche i bit in c2 dove c1 aveva -1 (o viceversa),
            # tanto poi mettiamo -1 su TUTTA la feature f.
            for p in posf
                newcube[p] = Int8(-1)
            end
            return newcube
        end

        return nothing
    end

    # Trova tutti gli indici di bit (1-based) che corrispondono a una feature data.
    function positions_of_feature(feat::Int, formula::TwoLevelDNFFormula)
        out = Int[]
        offset = 0
        # Scorriamo formula.atoms_by_feature in ordine, come in find_feature_op_threshold
        for (f, atoms) in sort(collect(formula.atoms_by_feature))
            for _ in atoms
                offset += 1
                if f == feat
                    push!(out, offset)
                end
            end
        end
        return out
    end

    # Restituisce (feature, threshold, op) per un indice di bit `idx` (1-based).
    function find_feature_op_threshold(idx::Int, formula::TwoLevelDNFFormula)
        offset = 0
        for (feat, atoms) in sort(collect(formula.atoms_by_feature))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: indice $idx fuori range")
    end

    ########################################################################
    # Applichiamo la unificazione finale
    ########################################################################
    copertura_finale_unified = final_unify_cubes!(copertura_finale)

    ########################################################################
    # 7) Convertiamo i cubi finali in BitVector e costruiamo la formula finale
    ########################################################################
    nuove_combinazioni = Vector{BitVector}()
    for cube in copertura_finale_unified
        combo = falses(formula.num_atoms)
        for (i, val) in enumerate(cube)
            if val == Int8(1)
                combo[i] = true
            end
        end
        push!(nuove_combinazioni, combo)
    end
    unique!(nuove_combinazioni)

    # Ritorno finale
    return TwoLevelDNFFormula(
        nuove_combinazioni,
        formula.num_atoms,
        formula.thresholds_by_feature,
        formula.atoms_by_feature,
        copertura_finale_unified,
    )
end

# TODO VALUTARE SE RIMANE (RAM HUNGRY)
function minimizza_dnf(::Val{:espresso_quine_pp}, formula::TwoLevelDNFFormula)
    ########################################################################
    # 1) Termini originali in formato 0/1 per garantire la stessa copertura
    ########################################################################
    termini_originali_bit = [
        Vector{Int8}([b ? 1 : 0 for b in combo])
        for combo in eachcombination(formula)
    ]

    ########################################################################
    # 2) Minimizzazione con Espresso (passo iniziale)
    ########################################################################
    espresso_risultato = minimizza_dnf(Val(:espresso), formula)
    println(length(eachcombination(espresso_risultato.combination)))

    ########################################################################
    # 3) Converti i risultati di Espresso in cubi con -1 (se prime_mask non vuota)
    ########################################################################
    termini_espresso = if isempty(espresso_risultato.prime_mask)
        [
            Vector{Int8}([bit ? 1 : 0 for bit in combo])
            for combo in eachcombination(espresso_risultato)
        ]
    else
        [Vector{Int8}(mask) for mask in espresso_risultato.prime_mask]
    end

    ########################################################################
    # 4) Quine–McCluskey: combinazione per prime implicants
    ########################################################################
    function quine_combinazione(implicanti::Vector{Vector{Int8}})

        # ----------------------------------------------------
        # EVITIAMO di combinare se produce un cubo tutto -1
        # ----------------------------------------------------
        function possono_combinarsi(c1, c2)
            diff_count = 0
            for i in eachindex(c1)
                # Se entrambi non sono -1 e diversi, contiamo differenza
                if c1[i] != Int8(-1) && c2[i] != Int8(-1) && c1[i] != c2[i]
                    diff_count += 1
                    if diff_count > 1
                        return false
                    end
                end
            end
            return diff_count == 1
        end

        function combina(c1, c2)
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

        nuovi = Vector{Vector{Int8}}()
        usati = fill(false, length(implicanti))

        # Tenta di combinare ogni paio
        for i in 1:length(implicanti)-1
            for j in i+1:length(implicanti)
                if possono_combinarsi(implicanti[i], implicanti[j])
                    new_cube = combina(implicanti[i], implicanti[j])

                    # BLOCCO del cubo tutto -1
                    if is_all_dontcare(new_cube)
                        # se non vuoi ammettere la tautologia,
                        # saltiamo l’inserimento
                        continue
                    end

                    push!(nuovi, new_cube)
                    usati[i] = true
                    usati[j] = true
                end
            end
        end

        # Aggiungi i non combinati
        for i in 1:length(implicanti)
            if !usati[i]
                push!(nuovi, implicanti[i])
            end
        end

        unique!(nuovi)

        # Se non si riduce più, abbiamo i prime implicants
        if length(nuovi) == length(implicanti)
            return nuovi
        else
            return quine_combinazione(nuovi)
        end
    end

    ########################################################################
    # 5) Quine: selezione copertura minima (conserva la copertura su termini_originali_bit)
    ########################################################################
    function quine_copertura_minima(
        orig_terms::Vector{Vector{Int8}},
        prime_impls::Vector{Vector{Int8}}
    )
        if isempty(prime_impls)
            return Vector{Vector{Int8}}()
        end

        function copre(cube1, cube2)
            # Ritorna true se cube1 copre cube2
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
                coverage[i, j] = copre(prime_impls[i], orig_terms[j])
            end
        end

        # Trova prime implicants essenziali
        selected = Set{Int}()
        for j in 1:size(coverage, 2)
            covering = findall(coverage[:, j])
            if length(covering) == 1
                push!(selected, covering[1])
            end
        end

        coperti = falses(size(coverage, 2))
        for i in selected
            for j in 1:size(coverage, 2)
                if coverage[i, j]
                    coperti[j] = true
                end
            end
        end

        uncovered = findall(!isequal(true), coperti)
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
                    coperti[j] = true
                end
            end
            uncovered = findall(!isequal(true), coperti)
        end

        return prime_impls[collect(selected)]
    end

    ########################################################################
    # 6) Passaggio finale: unificazione “spinta” su eventuali feature ridondanti
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
                        # BLOCCO tautologia: se new_cube è completamente -1, lo scartiamo
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

        # A questo punto c1 e c2 differiscono solo su una singola feature
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
        for (f, atoms) in sort(collect(formula.atoms_by_feature))
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
        for (feat, atoms) in sort(collect(formula.atoms_by_feature))
            for (thr, op) in atoms
                offset += 1
                if offset == idx
                    return (feat, thr, op)
                end
            end
        end
        error("find_feature_op_threshold: indice $idx fuori range")
    end

    ############ INIZIO FLUSSO PRINCIPALE ############

    # 1) Otteniamo i prime implicants di base
    prime_implicants = quine_combinazione(termini_espresso)

    # 2) Troviamo la copertura minima
    copertura_finale = quine_copertura_minima(termini_originali_bit, prime_implicants)

    # 3) Unificazione finale (spinta)
    copertura_finale_unified = final_unify_cubes!(copertura_finale)

    # 4) EVENTUALE check di “tautologia residua”
    # Se c'è esattamente 1 cubo ed è tutto -1, lo scartiamo:
    if length(copertura_finale_unified) == 1 && all(x -> x == Int8(-1), copertura_finale_unified[1])
        @warn "minimizza_dnf: trovato 1 cubo completamente -1 (tautologia). Lo rimuovo forzatamente!"
        copertura_finale_unified = Vector{Vector{Int8}}()  # ad es. formula vuota (“false”)
        # Oppure potresti tornare indietro ai prime implicants, ecc. Dipende dalla tua logica.
    end

    ########################################################################
    # 7) Costruzione formula finale
    ########################################################################
    nuove_combinazioni = Vector{BitVector}()
    for cube in copertura_finale_unified
        combo = falses(formula.num_atoms)
        for (i, val) in enumerate(cube)
            if val == Int8(1)
                combo[i] = true
            end
        end
        push!(nuove_combinazioni, combo)
    end
    unique!(nuove_combinazioni)

    return TwoLevelDNFFormula(
        nuove_combinazioni,
        formula.num_atoms,
        formula.thresholds_by_feature,
        formula.atoms_by_feature,
        copertura_finale_unified,
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
