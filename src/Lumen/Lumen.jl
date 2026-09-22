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

abstract type AbstractMinimizeSetup end
struct Fast <: AbstractMinimizeSetup end
struct Balanced <: AbstractMinimizeSetup end
struct Sop <: AbstractMinimizeSetup end
export Fast, Balanced, Sop

abstract type AbstractEncoding end
struct Univariate <: AbstractEncoding end
struct Multivariate <: AbstractEncoding end
export Univariate, Multivariate

# ---------------------------------------------------------------------------- #
#                                   include                                    #
# ---------------------------------------------------------------------------- #
export LumenShannonConfig, minimization_scheme
include("config.jl")
export LumenEnsemble
include("ensemble.jl")
export ThresholdSpace
include("thresholds.jl")
export LumenDNF, PlaLayout, Minimizer
export write_pla, read_pla, run_minimization, refine!, compact!
include("minimize.jl")
export gather_atoms
include("extract.jl")
# export lumen_shannon, assign, decision_set, leafclasses
export lumen_shannon, assign, decision_set
include("lumen_shannon.jl")

end