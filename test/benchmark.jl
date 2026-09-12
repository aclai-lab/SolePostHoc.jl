using Test

using SolePostHoc
const SP = SolePostHoc

using SoleXplorer
const SX = SoleXplorer

using SoleModels
const SM = SoleModels

using BenchmarkTools
using RDatasets
using DataFrames
using Random

using InteractiveUtils
using JET

using CategoricalArrays

# ---------------------------------------------------------------------------- #
#                                 iris test                                    #
# ---------------------------------------------------------------------------- #
iris_df = dataset("datasets", "iris")
y = String.(iris_df.Species)
Xdf = select(iris_df, Not(:Species))

dsc = setup_dataset(Xdf, y; model=SX.RandomForestClassifier(n_trees=20, max_depth=5), rng=42)
modelc = solexplorer(dsc)
solem = modelc.sole[1]

# config = SP.Lumen.LumenConfig(; minimization_scheme=:mitespresso)
# SP.Lumen.super_lumen(config, solem, M=100000, N=10);

config = SP.Lumen.LumenConfig(; minimization_scheme=Abc, M=5000)
lumen = SP.Lumen.lumen_shannon(config, solem)
@btime SP.Lumen.lumen_shannon(config, solem);
# 3.915 s (15192644 allocations: 801.34 MiB)
# 2.043 s (15223660 allocations: 755.99 MiB)
# 2.087 s (14735030 allocations: 745.08 MiB)
# 2.059 s (14733296 allocations: 745.04 MiB)
# 2.054 s (14772318 allocations: 744.97 MiB)
# 2.064 s (14759051 allocations: 744.59 MiB)
# 2.057 s (14759046 allocations: 744.59 MiB)
# 3.507 s (12960860 allocations: 604.03 MiB) PILE!
# 3.450 s (8954959 allocations: 471.73 MiB) PILE!
# 2.158 s (7168008 allocations: 407.08 MiB) PILE!

# ensemble
# 36.778 μs (433 allocations: 30.41 KiB)
# atoms
# 631.257 ns (5 allocations: 7.30 KiB)
# apply
# 11.769 μs (12 allocations: 8.05 KiB)


@code_warntype SP.Lumen.lumen_shannon(config, solem)
result = @report_opt SP.Lumen.lumen_shannon(config, solem)
open("jet_report.txt", "w") do io
    show(io, MIME"text/plain"(), result)
end

ctx, per_class_terms, config = SP.Lumen.lumen_shannon(config, solem);

@code_warntype SP.Lumen._finalize_decision_set(ctx, per_class_terms, config)
result = @report_opt SP.Lumen._extract_atoms_bfs_order(model)
open("jet_report.txt", "w") do io
    show(io, MIME"text/plain"(), result)
end
# ---------------------------------------------------------------------------- #
#                            apply on matrix data                              #
# ---------------------------------------------------------------------------- #
using SoleLogics
const SL = SoleLogics
using SoleData
const SD = SoleData
using StatsBase
X = Matrix(Xdf)
m = solem
d = X
@btime apply(m, X)
@code_warntype apply(m, X)
# 728.493 μs (11063 allocations: 774.88 KiB)


function apply(
    m::DecisionEnsemble{U,Branch{S}},
    d::Matrix{T}
)::Vector{S} where {U,S<:CategoricalValue,T<:AbstractFloat}
    subm = models(m)
    preds = Vector{S}(undef, size(d, 1))
    pool = CategoricalArrays.pool(first(info(m).supporting_predictions))::CategoricalPool{String,UInt32}
    counts = zeros(Int, length(levels(pool)))

    @inbounds for i in axes(d, 1)
        row = @view d[i, :]
        fill!(counts, 0)
        for sm in subm
            counts[levelcode(apply(sm, row))] += 1
        end
        preds[i] = pool[argmax(counts)]
    end

    return preds
end

function apply(
    m::Branch{S},
    d::SubArray{T,1}
)::S where {S<:CategoricalValue,T<:AbstractFloat}
    return checkantecedent(m, d) ?
        apply(posconsequent(m), d) :
        apply(negconsequent(m), d)
end

function apply(
    m::ConstantModel{S},
    ::SubArray{T,1}
)::S where {S<:CategoricalValue,T<:AbstractFloat}
    outcome(m)
end

function checkantecedent(
    m::Branch{S},
    d::SubArray{T,1}
)::Bool where {S,T<:AbstractFloat}
    check(antecedent(m), d)
end

function check(
    atom::Atom{S}, d::SubArray{T,1}
)::Bool where {S<:SM.ScalarCondition,T<:AbstractFloat}
        cond = SL.value(atom)

        cond_threshold = SD.threshold(cond)
        cond_operator = SD.test_operator(cond)
        cond_feature = SD.feature(cond)

        col = SD.i_variable(cond_feature)
        
        return cond_operator(d[col], cond_threshold)
end

struct LumenEnsemble{T<:AbstractFloat,R<:Integer}
    feat::Vector{Int32}      # >0: split variable; -1: leaf
    thr::Vector{T}
    left::Vector{Int32}      # taken when antecedent holds
    right::Vector{Int32}
    leaf::Vector{R}          # level code, valid where feat < 0
    roots::Vector{Int32}
    pool::CategoricalPool{String,R}
end

function flatten(m::DecisionEnsemble{U,Branch{S}}) where {U,S<:CategoricalValue}
    T = Float64
    feat, thr = Int32[], T[]
    left, right, leaf = Int32[], Int32[], UInt32[]
    roots = Int32[]

    function push!!(node)::Int32
        push!(feat, 0); push!(thr, zero(T))
        push!(left, 0); push!(right, 0); push!(leaf, 0)
        id = Int32(length(feat))

        if node isa ConstantModel
            feat[id] = -1
            leaf[id] = UInt32(levelcode(outcome(node)))
        else
            cond = SL.value(antecedent(node))
            feat[id] = Int32(SD.i_variable(SD.feature(cond)))
            thr[id]  = T(SD.threshold(cond))
            left[id]  = push!!(posconsequent(node))
            right[id] = push!!(negconsequent(node))
        end
        return id
    end

    for sm in models(m)
        push!(roots, push!!(sm))
    end
    pool = CategoricalArrays.pool(first(info(m).supporting_predictions))::CategoricalPool{String,UInt32}
    return LumenEnsemble{T,UInt32}(feat, thr, left, right, leaf, roots, pool)
end

function apply(
    f::LumenEnsemble{T,R},
    d::Matrix{T}
) where {T<:AbstractFloat,R}
    n = size(d, 1)
    pool = f.pool
    preds  = Vector{CategoricalValue{String,R}}(undef, n)
    counts = zeros(Int, length(levels(pool)))

    @inbounds for i in 1:n
        fill!(counts, 0)
        for r in f.roots
            node = r
            while f.feat[node] > 0
                node = d[i, f.feat[node]] < f.thr[node] ? f.left[node] : f.right[node]
            end
            counts[f.leaf[node]] += 1
        end
        preds[i] = pool[argmax(counts)]
    end
    return preds
end

@btime begin
    r=flatten(solem)
    apply(r, X)
end