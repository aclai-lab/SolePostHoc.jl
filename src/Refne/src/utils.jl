using Random


function generate_univers_of_combinations(min_values::Vector, max_values::Vector, n_points::Int=100)
    # Creation of linspace for each min/max pair
    ranges = [LinRange(min, max, n_points) for (min, max) in zip(min_values, max_values)]

    # Generate the Cartesian product of the combinations
    all_combinations = collect(product(ranges...))

    # Convert the result into a matrix of Float64
    # Each row corresponds to a combination
    mat = hcat([collect(comb) for comb in all_combinations]...)'

    return mat
end

######################################################################################################
function generate_univers_of_combinations_ott(min_values::Vector, max_values::Vector, n_points::Int=100)
    n_dims = length(min_values)
    total_combinations = n_points^n_dims
    
    # Threshold per decidere strategia
    reasonable_limit = 50000  # Limite ragionevole per refne
    
    if total_combinations <= reasonable_limit
        # CASO NORMALE: Generazione diretta ottimizzata
        println("Generazione diretta: $total_combinations combinazioni")
        
        result = Matrix{Float64}(undef, total_combinations, n_dims)
        
        # Pre-calcola i valori per ogni dimensione
        dim_values = [collect(range(min_val, max_val, length=n_points)) 
                     for (min_val, max_val) in zip(min_values, max_values)]
        
        # Generazione ottimizzata
        Threads.@threads for i in 1:total_combinations
            linear_idx = i - 1
            for dim in 1:n_dims
                value_idx = linear_idx % n_points + 1
                linear_idx ÷= n_points
                result[i, dim] = dim_values[dim][value_idx]
            end
        end
        
        return result
        
    else
        # CASO ESTREMO: Campionamento RIDOTTO e INTELLIGENTE per refne
        println("Troppe combinazioni ($total_combinations)!")
        
        # Per refne, un numero fisso e ragionevole di campioni
        target_samples = min(10000, reasonable_limit)  # Max 10k campioni per velocità
        
        println("Generando $target_samples campioni strategici per symbolic learning...")
        
        result = Matrix{Float64}(undef, target_samples, n_dims)
        
        # Pre-calcola i valori possibili per ogni dimensione
        dim_values = [collect(range(min_val, max_val, length=n_points)) 
                     for (min_val, max_val) in zip(min_values, max_values)]
        
        # STRATEGIA OTTIMIZZATA PER SYMBOLIC LEARNING
        sample_idx = 1
        
        # 1. CAMPIONI AGLI ESTREMI (molto importanti per decision trees)
        n_extremes = target_samples ÷ 4
        for i in 1:n_extremes
            for dim in 1:n_dims
                # Alterna tra min e max per ogni dimensione
                if rand() < 0.5
                    result[sample_idx, dim] = dim_values[dim][1]      # min
                else
                    result[sample_idx, dim] = dim_values[dim][end]    # max
                end
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 2. CAMPIONI "PURI" (una dimensione alla volta agli estremi)
        n_pure = target_samples ÷ 6
        for i in 1:n_pure
            # Scegli una dimensione casuale da mettere all'estremo
            extreme_dim = rand(1:n_dims)
            extreme_value = rand([1, n_points])  # min o max
            
            for dim in 1:n_dims
                if dim == extreme_dim
                    result[sample_idx, dim] = dim_values[dim][extreme_value]
                else
                    # Altre dimensioni: valore casuale dal grid
                    result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
                end
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 3. CAMPIONI CENTRALI (per catturare comportamenti intermedi)
        n_center = target_samples ÷ 6
        center_indices = [max(1, (n_points + 1) ÷ 2) for _ in 1:n_dims]
        
        for i in 1:n_center
            for dim in 1:n_dims
                # Intorno al centro con piccola variazione
                center_idx = center_indices[dim]
                variation = rand(-1:1)
                value_idx = clamp(center_idx + variation, 1, n_points)
                result[sample_idx, dim] = dim_values[dim][value_idx]
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 4. CAMPIONI CASUALI (per robustezza)
        while sample_idx <= target_samples
            for dim in 1:n_dims
                value_idx = rand(1:n_points)
                result[sample_idx, dim] = dim_values[dim][value_idx]
            end
            sample_idx += 1
        end
        
        println("Generati $target_samples campioni strategici")
        return result
    end
end

# Funzione specializzata per casi estremi ad alta dimensionalità
function generate_symbolic_grid_samples(min_values::Vector, max_values::Vector, n_points::Int, target_samples::Int=5000)
    """
    Genera campioni ottimizzati per symbolic learning su grid discreto
    Privilegia combinazioni che aiutano i decision trees a trovare pattern
    """
    n_dims = length(min_values)
    result = Matrix{Float64}(undef, target_samples, n_dims)
    
    # Valori possibili per ogni dimensione
    dim_values = [collect(range(min_val, max_val, length=n_points)) 
                 for (min_val, max_val) in zip(min_values, max_values)]
    
    sample_idx = 1
    
    # Strategia 1: Tutti gli estremi globali (2^n_dims combinazioni, ma sample solo alcune)
    n_global_extremes = min(target_samples ÷ 3, 2^min(n_dims, 10))  # Max 1024 per evitare esplosione
    
    for i in 1:n_global_extremes
        binary_pattern = digits(i-1, base=2, pad=n_dims)
        for dim in 1:n_dims
            extreme_idx = binary_pattern[dim] == 0 ? 1 : n_points
            result[sample_idx, dim] = dim_values[dim][extreme_idx]
        end
        sample_idx += 1
        if sample_idx > target_samples break end
    end
    
    # Strategia 2: Focus su singole dimensioni
    n_single_focus = min(target_samples ÷ 3, n_dims * n_points)
    
    for focus_dim in 1:n_dims
        for value_idx in 1:n_points
            if sample_idx > target_samples break end
            
            # Dimensione focus: valore specifico
            result[sample_idx, focus_dim] = dim_values[focus_dim][value_idx]
            
            # Altre dimensioni: valori casuali
            for dim in 1:n_dims
                if dim != focus_dim
                    result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
                end
            end
            
            sample_idx += 1
            if sample_idx > target_samples break end
        end
    end
    
    # Strategia 3: Riempimento casuale
    while sample_idx <= target_samples
        for dim in 1:n_dims
            result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
        end
        sample_idx += 1
    end
    
    return result
end



