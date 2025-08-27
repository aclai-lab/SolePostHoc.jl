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

function generate_univers_of_combinations_ott(min_values::Vector, max_values::Vector, n_points::Int=100)
    n_dims = length(min_values)
    total_combinations = n_points^n_dims
    
    # Threshold to decide strategy
    reasonable_limit = 50000  # Reasonable limit for refne
    
    if total_combinations <= reasonable_limit
        # NORMAL CASE: Optimized direct generation
        println("Direct generation: $total_combinations combinations")
        
        result = Matrix{Float64}(undef, total_combinations, n_dims)
        
        # Pre-calculate values for each dimension
        dim_values = [collect(range(min_val, max_val, length=n_points)) 
                     for (min_val, max_val) in zip(min_values, max_values)]
        
        # Optimized generation
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
        # EXTREME CASE: REDUCED and INTELLIGENT sampling for refne
        println("Too many combinations ($total_combinations)!")
        
        # For refne, a fixed and reasonable number of samples
        target_samples = min(10000, reasonable_limit)  # Max 10k samples for speed
        
        println("Generating $target_samples strategic samples for symbolic learning...")
        
        result = Matrix{Float64}(undef, target_samples, n_dims)
        
        # Pre-calculate possible values for each dimension
        dim_values = [collect(range(min_val, max_val, length=n_points)) 
                     for (min_val, max_val) in zip(min_values, max_values)]
        
        # OPTIMIZED STRATEGY FOR SYMBOLIC LEARNING
        sample_idx = 1
        
        # 1. EXTREME SAMPLES (very important for decision trees)
        n_extremes = target_samples ÷ 4
        for i in 1:n_extremes
            for dim in 1:n_dims
                # Alternate between min and max for each dimension
                if rand() < 0.5
                    result[sample_idx, dim] = dim_values[dim][1]      # min
                else
                    result[sample_idx, dim] = dim_values[dim][end]    # max
                end
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 2. "PURE" SAMPLES (one dimension at a time at extremes)
        n_pure = target_samples ÷ 6
        for i in 1:n_pure
            # Choose a random dimension to put at extreme
            extreme_dim = rand(1:n_dims)
            extreme_value = rand([1, n_points])  # min or max
            
            for dim in 1:n_dims
                if dim == extreme_dim
                    result[sample_idx, dim] = dim_values[dim][extreme_value]
                else
                    # Other dimensions: random value from grid
                    result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
                end
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 3. CENTRAL SAMPLES (to capture intermediate behaviors)
        n_center = target_samples ÷ 6
        center_indices = [max(1, (n_points + 1) ÷ 2) for _ in 1:n_dims]
        
        for i in 1:n_center
            for dim in 1:n_dims
                # Around center with small variation
                center_idx = center_indices[dim]
                variation = rand(-1:1)
                value_idx = clamp(center_idx + variation, 1, n_points)
                result[sample_idx, dim] = dim_values[dim][value_idx]
            end
            sample_idx += 1
            if sample_idx > target_samples break end
        end
        
        # 4. RANDOM SAMPLES (for robustness)
        while sample_idx <= target_samples
            for dim in 1:n_dims
                value_idx = rand(1:n_points)
                result[sample_idx, dim] = dim_values[dim][value_idx]
            end
            sample_idx += 1
        end
        
        println("Generated $target_samples strategic samples")
        return result
    end
end

function generate_symbolic_grid_samples(min_values::Vector, max_values::Vector, n_points::Int, target_samples::Int=5000)
    """
    Generates optimized samples for symbolic learning on discrete grid
    Prioritizes combinations that help decision trees find patterns
    """
    n_dims = length(min_values)
    result = Matrix{Float64}(undef, target_samples, n_dims)
    
    # Possible values for each dimension
    dim_values = [collect(range(min_val, max_val, length=n_points)) 
                 for (min_val, max_val) in zip(min_values, max_values)]
    
    sample_idx = 1
    
    # Strategy 1: All global extremes (2^n_dims combinations, but sample only some)
    n_global_extremes = min(target_samples ÷ 3, 2^min(n_dims, 10))  # Max 1024 to avoid explosion
    
    for i in 1:n_global_extremes
        binary_pattern = digits(i-1, base=2, pad=n_dims)
        for dim in 1:n_dims
            extreme_idx = binary_pattern[dim] == 0 ? 1 : n_points
            result[sample_idx, dim] = dim_values[dim][extreme_idx]
        end
        sample_idx += 1
        if sample_idx > target_samples break end
    end
    
    # Strategy 2: Focus on single dimensions
    n_single_focus = min(target_samples ÷ 3, n_dims * n_points)
    
    for focus_dim in 1:n_dims
        for value_idx in 1:n_points
            if sample_idx > target_samples break end
            
            # Focus dimension: specific value
            result[sample_idx, focus_dim] = dim_values[focus_dim][value_idx]
            
            # Other dimensions: random values
            for dim in 1:n_dims
                if dim != focus_dim
                    result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
                end
            end
            
            sample_idx += 1
            if sample_idx > target_samples break end
        end
    end
    
    # Strategy 3: Random filling
    while sample_idx <= target_samples
        for dim in 1:n_dims
            result[sample_idx, dim] = dim_values[dim][rand(1:n_points)]
        end
        sample_idx += 1
    end
    
    return result
end

