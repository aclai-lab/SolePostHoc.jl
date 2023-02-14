module SolePostHoc

using Metrics: mse
using SoleBase
using SoleData
using SoleLogics
using SoleModels
# TODO using SoleFeatures: findcorrelation

export extract_rules

include("rule-extraction.jl")
# Write your package code here!

end
