function generate_univers_of_combinations(min_values::Vector{Any}, max_values::Vector{Any}, n_points::Int=100)
    # Creazione dei linspace per ogni coppia min/max
    ranges = [LinRange(min, max, n_points) for (min, max) in zip(min_values, max_values)]

    # Genera il prodotto cartesiano delle combinazioni
    all_combinations = collect(product(ranges...))

    # Converte il risultato in una matrice di Float64
    # Ogni riga corrisponde a una combinazione
    mat = hcat([collect(comb) for comb in all_combinations]...)'

    return mat
end