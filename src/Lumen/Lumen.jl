module Lumen

using SoleLogics
const SL = SoleLogics
using SoleModels
const SM = SoleModels
using SoleData
const SD = SoleData

using CategoricalArrays
using DataFrames
using IterTools

using StatsBase: countmap

using ABC_jll

# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #
"""
    AbstractConfig

Abstract base type for all LUMEN configuration structs.

Concrete subtypes encapsulate the parameters needed to control a specific
algorithm variant. Using a common supertype allows generic code to accept
any configuration object without being tied to a particular implementation.

See also: [`LumenShannonConfig`](@ref)
"""
abstract type AbstractConfig end

abstract type AbstractMinimization end
struct Abc <: AbstractMinimization end
struct MitEspresso <: AbstractMinimization end
export Abc, MitEspresso

# ---------------------------------------------------------------------------- #
#                                   include                                    #
# ---------------------------------------------------------------------------- #
include("config.jl")
export LumenShannonConfig
include("ensemble.jl")
include("thresholds.jl")
include("atoms.jl")
include("extract.jl")

# include("atoms.jl")
# export ThresholdSpace, _extract_atoms_bfs_order, extract_atoms_bfs_order
# include("apply.jl")
# include("sequential_minimizer.jl")
# include("super_minimizer.jl")
include("lumen_shannon.jl")
export lumen_shannon

end