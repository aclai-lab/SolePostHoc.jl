module BATrees

# Includi i file con le funzioni secondarie
include("run.jl")
include("trad.jl")

# Importa i moduli
using .RunModule
using .TradModule

export batrees

"""
    main_function()

Funzione 'principale' che fa uso di run_function e trad_function
definite negli altri file.
"""
function batrees(f; dataset_name="iris", num_trees=5, max_depth=3)
    WRAP_batrees(f; dataset_name="iris", num_trees=5, max_depth=3)
    trad()
end

end 