function minimizza_dnf(_minimization_scheme::Val, formula::TwoLevelDNFFormula, kwargs...)
    error("Unknown minimization scheme: $(_minimization_scheme)!")
end

function minimizza_dnf(
    ::Val{:mitespresso},
    formula::TwoLevelDNFFormula;
    silent = false,
    mitespresso_kwargs...,
)
    formula = convert(SoleLogics.DNF, formula)
    silent || (println(); @show formula)
    formula = SoleData.espresso_minimize(formula, silent; mitespresso_kwargs...)
    silent || (println(); @show formula)
    #silent || println("===========================")
    #silent || println("dump:" * string(dump(formula)))
    #silent || println("===========================")
    #formula = convert(TwoLevelDNFFormula, formula)   TODO
    #silent || (println(); @show formula)             TODO IN PROGRESS BY GIO
    return formula
end


"""
	minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula, horizontal = 1.0)

Simplifies a custom OR formula using the Quine algorithm.

This function takes a `TwoLevelDNFFormula` object and applies the Quine algorithm to minimize the number of combinations in the formula. It returns a new `TwoLevelDNFFormula` object with the simplified combinations.

"""
#==#
function minimizza_dnf(::Val{:quine}, formula::TwoLevelDNFFormula)
    # Inizializzazione dei termini di partenza
    terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
    for (i, term) in enumerate(eachcombination(formula))
        for (j, x) in enumerate(term)
            terms[i][j] = x ? Int8(1) : Int8(0)
        end
    end

    if isempty(terms)
        return formula  # Se non ci sono termini, restituisci direttamente la formula vuota
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

    # Funzione per combinare due cubi in un nuovo cubo generalizzato
    function combina_cubi(cube1::Vector{Int8}, cube2::Vector{Int8})
        result = Vector{Int8}(undef, length(cube1))
        for i in eachindex(cube1)
            result[i] = ifelse(cube1[i] == cube2[i], cube1[i], Int8(-1))
        end
        return result
    end

    # Fase di minimizzazione: combinazione dei termini
    function minimizza_step!(termini::Vector{Vector{Int8}})
        nuovi_termini = Vector{Vector{Int8}}()
        usati = fill(false, length(termini))

        # Itera su tutte le coppie di termini per combinazioni
        for i in 1:(length(termini) - 1)
            for j in (i + 1):length(termini)
                diff_count = 0
                for k in eachindex(termini[i])
                    if termini[i][k] != termini[j][k] && termini[i][k] != Int8(-1) && termini[j][k] != Int8(-1)
                        diff_count += 1
                        if diff_count > 1
                            break
                        end
                    end
                end
                # Combina termini solo se differiscono in una singola posizione
                if diff_count == 1
                    push!(nuovi_termini, combina_cubi(termini[i], termini[j]))
                    usati[i] = true
                    usati[j] = true
                end
            end
        end

        # Aggiungi termini non usati
        for i in 1:length(termini)
            if !usati[i]
                push!(nuovi_termini, termini[i])
            end
        end

        # Rimuovi duplicati
        unique!(nuovi_termini)

        if length(nuovi_termini) == length(termini)
            return nuovi_termini  # Nessuna nuova combinazione, termina
        end

        return minimizza_step!(nuovi_termini)  # Continua a iterare
    end

    # Trova la copertura minima dei termini originali usando i primi implicanti
    function trova_copertura_minima(termini::Vector{Vector{Int8}}, primi::Vector{Vector{Int8}})
        coverage = falses(length(primi), length(termini))

        for i in 1:length(primi)
            for j in 1:length(termini)
                coverage[i, j] = copre(primi[i], termini[j])
            end
        end

        # Trova la copertura minima tramite backtracking
        function trova_copertura_backtrack(coperti::BitVector, implicanti_selezionati::Vector{Int}, candidati::Vector{Int})
            if all(coperti)
                return implicanti_selezionati
            end

            best_soluzione = nothing

            for candidato in candidati
                nuovi_coperti = coperti .| coverage[candidato, :]
                nuova_soluzione = trova_copertura_backtrack(nuovi_coperti, [implicanti_selezionati; candidato], setdiff(candidati, [candidato]))

                if nuova_soluzione !== nothing && (best_soluzione === nothing || length(nuova_soluzione) < length(best_soluzione))
                    best_soluzione = nuova_soluzione
                end
            end

            return best_soluzione
        end

        candidati = collect(1:length(primi))
        coperti = falses(length(termini))

        soluzione_minima = trova_copertura_backtrack(coperti, Int[], candidati)
        return primi[soluzione_minima]
    end

    try
        # Fase 1: Genera i primi implicanti
        primi_implicanti = minimizza_step!(terms)

        # Fase 2: Trova la copertura minima dei termini
        minimized_terms = trova_copertura_minima(terms, primi_implicanti)

        # Conversione dei termini minimizzati in combinazioni leggibili
        nuove_combinazioni = BitVector[]
        seen = Set{BitVector}()  # Set per evitare duplicati

        for term in minimized_terms
            combo = falses(nuberofatoms(formula))
            for (i, val) in enumerate(term)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            if combo ∉ seen
                push!(seen, combo)
                push!(nuove_combinazioni, combo)
            end
        end

        sort!(nuove_combinazioni)
        return TwoLevelDNFFormula(
            nuove_combinazioni,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            collect(minimized_terms),
        )
    catch e
        @warn "Errore durante la minimizzazione: $e"
        return formula  # Restituisci la formula originale in caso di errore
    end
end


#=MINIMI CHE NON MI SEMBRANO GLOBALI=#
function minimizza_dnf(::Val{:quine_naive}, formula::TwoLevelDNFFormula)
    # Converti i termini in vettori binari
    terms = [Vector{Int}([x ? 1 : 0 for x in term]) for term in eachcombination(formula)]
    
    if isempty(terms)
        return formula
    end
    
    # 1. Trova tutti i primi implicanti possibili
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
            
            # Aggiungi i termini non utilizzati come primi implicanti
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
    
    # 2. Verifica se un primo implicante copre un termine
    function covers(prime, term)
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != term[i]
                return false
            end
        end
        return true
    end
    
    # 3. Genera tutte le possibili combinazioni di k elementi da un array
    function generate_combinations(arr, k)
        n = length(arr)
        if k > n
            return Vector{Vector{eltype(arr)}}()
        end
        
        result = Vector{Vector{eltype(arr)}}()
        
        # Funzione ricorsiva per generare le combinazioni
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
    
    # 4. Genera tutte le possibili combinazioni di primi implicanti
    function all_combinations(primes)
        result = Vector{Vector{Vector{Int}}}()
        for k in 1:length(primes)
            append!(result, generate_combinations(primes, k))
        end
        return result
    end
    
    # 5. Verifica se una combinazione copre tutti i termini
    function covers_all(combo, terms)
        for term in terms
            if !any(prime -> covers(prime, term), combo)
                return false
            end
        end
        return true
    end
    
    # Esegui l'algoritmo completo
    try
        # Trova tutti i primi implicanti
        primes = find_prime_implicants(terms)
        
        # Trova la combinazione minima che copre tutti i termini
        min_cover = primes  # Fallback alla soluzione completa
        min_size = length(primes)
        
        # Prova tutte le possibili combinazioni di primi implicanti
        for combo in all_combinations(primes)
            if covers_all(combo, terms) && length(combo) < min_size
                min_cover = combo
                min_size = length(combo)
            end
        end
        
        # Converti il risultato in BitVector
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
        @warn "Errore durante la minimizzazione: $e"
        return formula
    end
end

#=MINIMI GLOBALI MA ERRATI=#
function minimizza_dnf(::Val{:quine_oldstyle}, formula::TwoLevelDNFFormula)
    # Inizializzazione dei termini di partenza
        terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
        for (i, term) in enumerate(eachcombination(formula))
            for (j, x) in enumerate(term)
                terms[i][j] = x ? Int8(1) : Int8(0)
            end
    end
    
    if isempty(terms)
        return formula  # Se non ci sono termini, restituisci direttamente la formula vuota
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

    # Funzione per combinare due cubi in un nuovo cubo generalizzato
    function combina_cubi(cube1::Vector{Int8}, cube2::Vector{Int8})
            result = Vector{Int8}(undef, length(cube1))
            for i in eachindex(cube1)
                result[i] = ifelse(cube1[i] == cube2[i], cube1[i], Int8(-1))
            end
            return result
    end

    # Fase di minimizzazione: combinazione dei termini
    function minimizza_step!(termini::Vector{Vector{Int8}})
        if length(termini) <= 1
            return termini
        end
        
            nuovi_termini = Vector{Vector{Int8}}()
            usati = fill(false, length(termini))
            
        # Itera su tutte le coppie di termini per combinazioni
            for i in 1:(length(termini)-1)
                for j in (i+1):length(termini)
                    diff_count = 0
                    for k in eachindex(termini[i])
                        if termini[i][k] != termini[j][k] && termini[i][k] != Int8(-1) && termini[j][k] != Int8(-1)
                            diff_count += 1
                            if diff_count > 1
                                break
                            end
                        end
                    end
                # Combina termini solo se differiscono in una singola posizione
                    if diff_count == 1
                    push!(nuovi_termini, combina_cubi(termini[i], termini[j]))
                            usati[i] = true
                            usati[j] = true
                    end
                end
            end
            
            # Aggiungi termini non usati
            for i in 1:length(termini)
                if !usati[i]
                    push!(nuovi_termini, termini[i])
                end
            end
            
        # Rimuovi duplicati e ordina
            unique!(nuovi_termini)
            
            if length(nuovi_termini) == length(termini)
            return nuovi_termini  # Nessuna nuova combinazione, termina
            end
            
        return minimizza_step!(nuovi_termini)  # Continua a iterare
    end

    # Trova la copertura minima dei termini originali usando i primi implicanti
    function trova_copertura_minima(termini::Vector{Vector{Int8}}, primi::Vector{Vector{Int8}})
        if isempty(primi)
            return termini
        end
        
        # Costruzione della matrice di copertura
            coverage = falses(length(primi), length(termini))
            for i in 1:length(primi)
                for j in 1:length(termini)
                    coverage[i, j] = copre(primi[i], termini[j])
                end
            end
            
        # Trova gli implicanti essenziali
            selected_primes = Int[]
            essential_coverage = falses(size(coverage, 2))
            
            for j in 1:size(coverage, 2)
                covering_primes = findall(coverage[:, j])
                if length(covering_primes) == 1
                    push!(selected_primes, covering_primes[1])
                    essential_coverage .|= coverage[covering_primes[1], :]
                end
            end
            
        # Se necessario, seleziona altri implicanti con approccio greedy
            if !all(essential_coverage)
                remaining_terms = .!essential_coverage
                selected = Set(selected_primes)
                
                while any(remaining_terms)
                    best_prime = -1
                    max_coverage = 0
                    
                    for (i, prime) in enumerate(primi)
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
                
                return primi[collect(selected)]
            end
            
            return primi[selected_primes]
    end

    try
        # Fase 1: Genera i primi implicanti
        primi_implicanti = minimizza_step!(terms)
        # Fase 2: Trova la copertura minima dei termini
        minimized_terms = trova_copertura_minima(terms, primi_implicanti)
        
        # Conversione dei termini minimizzati in combinazioni leggibili
        nuove_combinazioni = BitVector[]
        seen = Set{BitVector}()  # Set per evitare duplicati

        for term in minimized_terms
            combo = falses(nuberofatoms(formula))
            for (i, val) in enumerate(term)
                if val == Int8(1)
                    combo[i] = true
                end
            end
            if combo ∉ seen
                push!(seen, combo)
                push!(nuove_combinazioni, combo)
            end
        end
        
        sort!(nuove_combinazioni)
        return TwoLevelDNFFormula(
            nuove_combinazioni,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            collect(minimized_terms),
        )
    catch e
        @warn "Errore durante la minimizzazione: $e"
        return formula  # Restituisci la formula originale in caso di errore
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
                    can_combine, pos = possono_combinarsi(current_cubes[i], current_cubes[j])
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
                coverage = count(term -> copre(cube, term), collect(remaining_terms))
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
            filter!(term -> !copre(best_cube, term), remaining_terms)
        end

        return essential
    end

    function espresso_minimize(terms; minimization_method_kwargs...)
        if isempty(terms)
            return Vector{Vector{Int}}()
        end

        combined_terms = trova_combinazioni(terms)
        if isempty(combined_terms)
            return terms
        end

        essential = find_essential_cubes(combined_terms, terms)
        if isempty(essential)
            return terms
        end

        return essential
    end

    # Esegui la minimizzazione
    minimized_terms = espresso_minimize(terms; minimization_method_kwargs...)

    # Converti il risultato in TritVector
    nuove_combinazioni = TritVector[]
    seen = Set{TritVector}()  # Set per tenere traccia dei termini già visti

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

        # Aggiungi il termine solo se non è già stato visto
        if trit_combo ∉ seen
            push!(seen, trit_combo)
            push!(nuove_combinazioni, trit_combo)
        end
    end

    sort!(nuove_combinazioni)  # Assuming you have defined sorting for TritVector
    return TwoLevelDNFFormula(
        nuove_combinazioni,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
    )
end

#=
function minimizza_dnf(::Val{:espresso}, formula::TwoLevelDNFFormula; minimization_method_kwargs...)
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

    function espresso_minimize(terms; minimization_method_kwargs...)
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
    minimized_terms = espresso_minimize(terms; minimization_method_kwargs...)

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
        eachatomsbyfeature(formula),
        collect(minimized_terms),
    )
end
=#
#=prova risoluzione...=#
function minimizza_dnf(::Val{:quine_strict}, formula::TwoLevelDNFFormula)
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

    # Funzione per verificare se un punto è coperto da un termine
    function punto_coperto(punto::Vector{Int8}, termine::Vector{Int8})
        for i in eachindex(punto)
            if termine[i] != Int8(-1) && termine[i] != punto[i]
                return false
            end
        end
        return true
    end

    # Funzione per verificare se due termini possono essere combinati
    # senza introdurre nuovi punti di copertura
    function possono_combinarsi_safe(cube1::Vector{Int8}, cube2::Vector{Int8}, termini_originali::Vector{Vector{Int8}})
        # Prima verifica se differiscono per una sola posizione
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
        
        # Crea il termine combinato per il test
        combined = copy(cube1)
        combined[diff_pos] = Int8(-1)
        
        # Per ogni termine originale
        for original in termini_originali
            # Se il termine combinato copre questo punto
            if punto_coperto(original, combined)
                # Verifica che fosse coperto da almeno uno dei termini originali
                if !punto_coperto(original, cube1) && !punto_coperto(original, cube2)
                    return false, -1
                end
            end
        end
        
        return true, diff_pos
    end

    function minimizza_step!(termini::Vector{Vector{Int8}}, termini_originali::Vector{Vector{Int8}})
        if length(termini) <= 1
            return termini
        end
        
        nuovi_termini = Vector{Vector{Int8}}()
        usati = fill(false, length(termini))
        
        for i in 1:length(termini)-1
            for j in i+1:length(termini)
                can_combine, pos = possono_combinarsi_safe(termini[i], termini[j], termini_originali)
                
                if can_combine
                    nuovo_termine = copy(termini[i])
                    nuovo_termine[pos] = Int8(-1)
                    push!(nuovi_termini, nuovo_termine)
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
        
        return minimizza_step!(nuovi_termini, termini_originali)
    end

    function trova_copertura_minima(termini::Vector{Vector{Int8}}, primi::Vector{Vector{Int8}})
        if isempty(primi)
            return termini
        end
        
        coverage = falses(length(primi), length(termini))
        for i in 1:length(primi)
            for j in 1:length(termini)
                coverage[i,j] = punto_coperto(termini[j], primi[i])
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
        
        if !isempty(selected_primes)
            covered_terms = vec(any(coverage[selected_primes, :], dims=1))
            if all(covered_terms)
                return primi[unique(selected_primes)]
            end
        end
        
        # Approccio greedy per i termini rimanenti
        selected = Set(selected_primes)
        uncovered = Set(1:size(coverage, 2))
        
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
        # Manteniamo una copia dei termini originali per il controllo
        termini_originali = copy(terms)
        primi = minimizza_step!(terms, termini_originali)
        minimized = trova_copertura_minima(termini_originali, primi)
        
        risultato = Vector{BitVector}(undef, length(minimized))
        for (i, term) in enumerate(minimized)
            combo = falses(nuberofatoms(formula))
            for (j, val) in enumerate(term)
                if val == Int8(1)
                    combo[j] = true
                end
            end
            risultato[i] = combo
        end
        
        return TwoLevelDNFFormula(
            risultato,
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            minimized
        )
        
    catch e
        @warn "Errore durante la minimizzazione: $e"
        @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
        return formula
    end
end

#prova risoluzione#
function minimizza_dnf(::Val{:quine_petrick}, formula::TwoLevelDNFFormula)
    terms = eachcombination(formula)

    if isempty(terms)
        return formula
    end

    # Trova i primi implicanti
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
            result[pos] = -1  # Aggiungi un don't care
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

    # Verifica se un termine copre un altro termine
    function covers(prime::Vector{Int}, term::BitVector)
        for i in eachindex(prime)
            if prime[i] != -1 && prime[i] != (term[i] ? 1 : 0)
                return false
            end
        end
        return true
    end

    # Metodo di Petrick per minimizzare la copertura
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

    # Crea una nuova lista di combinazioni minimizzate
    minimized_combinations = BitVector[]
    for prime in minimized_cover
        bitvec = falses(nuberofatoms(formula))
        for (i, val) in enumerate(prime)
            if val == 1
                bitvec[i] = true
            end
        end
        push!(minimized_combinations, bitvec)
    end

    return TwoLevelDNFFormula(
        minimized_combinations,
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        minimized_cover
    )
end


# quine fuso ad espresso [DEFINITIVO !?!?!?1?!?!?!1?!1?!]
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
                iteration = 1
                
                while true
                    found_new = false
                    total_combinations = div(length(terms) * (length(terms) - 1), 2)
                    progress_counter = 0
                    
                    println("\nEspresso - Iterazione $iteration")
                    
                    for i in 1:length(terms)
                        for j in (i+1):length(terms)
                            progress_counter += 1
                            if progress_counter % 1000 == 0
                                show_progress("Combinazioni", progress_counter, total_combinations)
                            end
                            
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
                    iteration += 1
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
                
                println("\nQuine - Fase implicanti essenziali")
                for (i, term) in enumerate(terms)
                    if i % 1000 == 0
                        show_progress("Analisi termini", i, length(terms))
                    end
                    covering_primes = filter(p -> copre(p, term), primes)
                    if length(covering_primes) == 1
                        push!(essential, first(covering_primes))
                        push!(covered, term)
                    end
                end
                
            # Poi usa un approccio greedy per i rimanenti
                remaining_terms = setdiff(Set(terms), covered)
                progress_counter = 0
                
                println("\nQuine - Fase greedy")
                while !isempty(remaining_terms)
                    progress_counter += 1
                    if progress_counter % 1000 == 0
                        show_progress("Ottimizzazione", progress_counter, length(terms))
                    end
                    
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

        println("\nInizializzazione")
        initial_terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
        for (i, term) in enumerate(eachcombination(formula))
            if i % 1000 == 0
                show_progress("Conversione termini", i, length(eachcombination(formula)))
            end
            for (j, x) in enumerate(term)
                initial_terms[i][j] = x ? Int8(1) : Int8(0)
            end
        end
        
        if isempty(initial_terms)
            return formula
        end

        try
            println("\nAvvio minimizzazione DNF")
            
            espresso_terms = modified_espresso(initial_terms)
            final_terms = modified_quine(espresso_terms)
            
            println("\nFinalizzazione")
            result_combinations = Vector{BitVector}()
            for (i, term) in enumerate(final_terms)
                if i % 1000 == 0
                    show_progress("Conversione risultati", i, length(final_terms))
                end
                combo = falses(nuberofatoms(formula))
                for (i, val) in enumerate(term)
                    if val == Int8(1)
                        combo[i] = true
                    end
                end
                push!(result_combinations, combo)
            end
            
            println("\nMinimizzazione completata!")
            
            return TwoLevelDNFFormula(
                result_combinations,
                nuberofatoms(formula),
                eachthresholdsbyfeature(formula),
                eachatomsbyfeature(formula),
                final_terms
            )
            
        catch e
            @warn "Errore durante la minimizzazione: $e"
            @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
            return formula
        end
    end
end

function minimizza_dnf(::Val{:pina_interrupted}, formula::TwoLevelDNFFormula)
    let spinner_state = Ref(1)  # Stato del "spinner"
        function show_progress(desc::String, current::Int, total::Int)
            spinner = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
            percentage = round(Int, current * 100 / total)

            # Aggiorna lo stato dello spinner
            spinner_state[] = spinner_state[] % length(spinner) + 1

            print("\r$desc: $(spinner[spinner_state[]]) $percentage%")
            flush(stdout)  # Mostra subito l'output

            current == total && println()
        end

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
            iteration = 1

            while true
                found_new = false
                combined = Set{Tuple{Vector{Int8}, Vector{Int8}}}()

                println("\nEspresso - Iterazione $iteration")

                for i in 1:length(terms)
                    for j in (i + 1):length(terms)
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

                append!(terms, result)
                terms = unique(sort(terms))
                empty!(result)
                iteration += 1
            end

            return terms
        end

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

            println("\nQuine - Fase implicanti essenziali")
            for (i, term) in enumerate(terms)
                if i % 1000 == 0
                    show_progress("Analisi termini", i, length(terms))
                end
                covering_primes = filter(p -> copre(p, term), primes)
                if length(covering_primes) == 1
                    push!(essential, first(covering_primes))
                    push!(covered, term)
                end
            end

            remaining_terms = setdiff(Set(terms), covered)
            println("\nQuine - Fase greedy")

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

        println("\nInizializzazione")
        initial_terms = [Vector{Int8}(undef, length(term)) for term in eachcombination(formula)]
        for (i, term) in enumerate(eachcombination(formula))
            if i % 1000 == 0
                show_progress("Conversione termini", i, length(eachcombination(formula)))
            end
            for (j, x) in enumerate(term)
                initial_terms[i][j] = x ? Int8(1) : Int8(0)
            end
        end

        if isempty(initial_terms)
            return formula
        end

        try
            println("\nAvvio minimizzazione DNF")

            espresso_terms = trova_combinazioni_per_quine(initial_terms)
            final_terms = trova_copertura_essenziale(initial_terms, espresso_terms)

            println("\nFinalizzazione")
            result_combinations = Vector{BitVector}()
            for (i, term) in enumerate(final_terms)
                if i % 1000 == 0
                    show_progress("Conversione risultati", i, length(final_terms))
                end
                combo = falses(nuberofatoms(formula))
                for (i, val) in enumerate(term)
                    if val == Int8(1)
                        combo[i] = true
                    end
                end
                push!(result_combinations, combo)
            end

            println("\nMinimizzazione completata!")

            return TwoLevelDNFFormula(
                result_combinations,
                nuberofatoms(formula),
                eachthresholdsbyfeature(formula),
                eachatomsbyfeature(formula),
                final_terms
            )

        catch e
            @warn "Errore durante la minimizzazione: $e"
            @warn "Stack trace: " * sprint(showerror, e, catch_backtrace())
            return formula
        end
    end
end

# FUNZIONANTI stessi risultati di espresso
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
        error("find_feature_op_threshold: indice $idx fuori range")
    end

    # Unificazione finale, se la vuoi
        copertura_finale_unified = final_unify_cubes!(copertura_finale)

        ########################################################################
        # 7) Costruzione formula finale
        ########################################################################
        nuove_combinazioni = Vector{BitVector}()
        for cube in copertura_finale_unified
            combo = falses(nuberofatoms(formula))
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
            nuberofatoms(formula),
            eachthresholdsbyfeature(formula),
            eachatomsbyfeature(formula),
            copertura_finale_unified,
        )
end

# IA GENERATED (TODO TEST it)
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

    # Restituisce (feature, threshold, op) per un indice di bit `idx` (1-based).
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
        combo = falses(nuberofatoms(formula))
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
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        copertura_finale_unified,
    )
end

# # IA GENERATED (TODO TEST it) (RAM HUNGRY)
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
        combo = falses(nuberofatoms(formula))
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
        nuberofatoms(formula),
        eachthresholdsbyfeature(formula),
        eachatomsbyfeature(formula),
        copertura_finale_unified,
    )
end