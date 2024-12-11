module SolePostHoc

using SoleData
using SoleLogics
using SoleModels

# Function mse of Metrics.jl
mse = (y_pred,t_true) -> (sum((y_true .- y_pred).^2) / length(y_true))

include("rule-extraction.jl")

using .RuleExtraction

end
