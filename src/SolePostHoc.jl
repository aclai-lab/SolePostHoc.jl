module SolePostHoc

using Reexport
using SoleData
using SoleLogics
using SoleModels

# wrapper for SoleXplorer
using DataFrames
import SoleModels: listrules
listrules(m, X::AbstractDataFrame, y, args...; kwargs...) = listrules(m, args...; kwargs...)

# Function mse (from Metrics.jl)
mse = (y_pred, t_true) -> (sum((y_true .- y_pred) .^ 2) / length(y_true))

include("rule-extraction.jl")

#@reexport using .RuleExtraction

end
