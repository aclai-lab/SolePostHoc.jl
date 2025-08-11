module BATrees

# Includi i file con le funzioni secondarie
include("core.jl")
include("apiBaSole.jl")

# Importa i moduli
using .RunModule
using .TradModule

export batrees


"""
    batrees(f; dataset_name="iris", num_trees=10, max_depth=10, dsOutput=true)

Builds and trains a set of binary decision trees OR using the specified function `f`.

# Arguments
- `f`: An SoleForest.
- `dataset_name::String`: The name of the dataset to be used. Default is "iris".
- `num_trees::Int`: The number of trees to be built. Default is 10.
- `max_depth::Int`: The maximum depth of each tree. Default is 10.
- `dsOutput::Bool`: A flag indicating whether to return the dsStruct output. Default is true. if false, returns the result single tree.

# Returns
- If `dsOutput` is true, returns the result is in DecisionSet ds.
- If `dsOutput` is false, returns the result is SoleTree t`.

# Example
"""
function batrees(f=nothing; dataset_name="iris", num_trees=10, max_depth=10, dsOutput=true, mode_obj = 0)
    if (isnothing(f))
        if (dsOutput)
            WRAP_batrees(nothing, max_depth, dataset_name=dataset_name, num_trees=num_trees, mod = mode_obj)
            t = BAinSoleTree()
            return t
        end
    else
        if (dsOutput)
            class_map = WRAP_batrees(f,max_depth, mod = mode_obj)

            println("FINITO WRAPPPP======>",class_map)
            ds = BAinDS(class_map)
            t = BAinSoleTree()
            println("=========================")
            println("ds: ", ds)
            println("t: ", t)
            println("=========================")
            return ds
        else
            WRAP_batrees(f , max_depth, mod = mode_obj)
            t = BAinSoleTree()
            return t
        end
    end
end

end