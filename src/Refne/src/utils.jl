function generate_univers_of_combinations(min_values::Vector{Any}, max_values::Vector{Any}, n_points::Int=100)
    # Creation of linspace for each min/max pair
    ranges = [LinRange(min, max, n_points) for (min, max) in zip(min_values, max_values)]

    # Generate the Cartesian product of the combinations
    all_combinations = collect(product(ranges...))

    # Convert the result into a matrix of Float64
    # Each row corresponds to a combination
    mat = hcat([collect(comb) for comb in all_combinations]...)'

    return mat
end