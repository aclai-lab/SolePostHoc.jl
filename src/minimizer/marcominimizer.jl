module SoleMarcoMinimizer

using Revise
using Logging
using Dates
using DataStructures
using SoleModels
using AbstractTrees
using SoleData
using SoleData: UnivariateScalarAlphabet

using ModalDecisionTrees
using SoleLogics
using BenchmarkTools, StatProfilerHTML
using Base: intersect
using Base.Threads: @threads
using Base.Threads: Atomic, atomic_add!
using Profile
using ConcurrentCollections
using ProgressMeter

##
# Types
##

# TritVector TODO Choose and implement one 
include("types/trit-vector.jl")
include("types/balanced-trit-vector.jl")
include("types/balanced-ternary-vector.jl")

# ConstVariable e myOwnTypes
include("types/types.jl")

##
# Utils
##

# Gestione dei report su file
include("utils/report.jl")

# IO-Algoritmic-Utils - IO algoritmico del progetto
include("utils/IO.jl")

# Minor-Algoritmic-Utils - core algoritmico del progetto
include("utils/minor.jl")

# Algoritmic-Utils - core algoritmico del progetto
include("utils/core.jl")

# Algoritmic-Optimization-Utils - core algoritmico del progetto se avviato in ott_mode
include("utils/coreOttMode.jl")

include("utils/minimization.jl")

include("deprecate.jl")

end