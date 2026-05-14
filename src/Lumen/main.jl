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


include("config.jl")

export lumen, LumenConfig, LumenResult

const Operators = Union{typeof(<),typeof(>),typeof(≤),typeof(≥)}
const Float = Union{Float32,Float64}

# ---------------------------------------------------------------------------- #
#                   initial minimization algorithms setup                      #
# ---------------------------------------------------------------------------- #
function setup_espresso()
    espressobinary = try
        joinpath(SD.load(SD.MITESPRESSOLoader()), "espresso")
    catch e
        error("Failed to setup espresso binary: $e")
    end
    isfile(espressobinary) ||
        error("espresso binary not found at $espressobinary")
    return espressobinary
end

function setup_boom() end

function ensure_abc_binary(; force_rebuild=false)
    abc_binary = joinpath(@__DIR__, "abc")
    if isfile(abc_binary) && !force_rebuild
        @info "ABC binary already exists at: $abc_binary (skipping download/compilation)"
        return abc_binary
    end
    @info "Setting up ABC binary..."
    abc_temp_dir = mktempdir(; prefix="abc_build_")
    abc_url = "https://github.com/berkeley-abc/abc/archive/refs/heads/master.tar.gz"
    try
        @info "Downloading ABC source code..."
        tarfile = joinpath(abc_temp_dir, "abc-master.tar.gz")
        run(`curl -L -o $tarfile $abc_url`)
        extract_dir = joinpath(abc_temp_dir, "extract")
        mkdir(extract_dir)
        @info "Extracting ABC source code..."
        if success(`which tar`)
            run(`tar -xzf $tarfile -C $extract_dir`)
        else
            error("System tar command not found. Please install tar utilities.")
        end
        abc_source_dir = joinpath(extract_dir, "abc-master")
        if !isdir(abc_source_dir)
            error("Failed to extract ABC source code - directory not found")
        end
        @info "Compiling ABC... This may take a few minutes."
        old_dir = pwd()
        cd(abc_source_dir)
        try
            if !success(`which make`)
                error("make command not found. Please install build tools (make, gcc, etc.)")
            end
            run(`make ABC_USE_NO_READLINE=1`)
            compiled_binary = joinpath(abc_source_dir, "abc")
            if isfile(compiled_binary)
                cp(compiled_binary, abc_binary; force=true)
                chmod(abc_binary, 0o755)
                @info "ABC compiled successfully at: $abc_binary"
            else
                error("ABC compilation completed but binary not found")
            end
        finally
            cd(old_dir)
        end
        return abc_binary
    catch e
        @error "Failed to download/compile ABC: $e. Consider downloading ABC manually from https://github.com/berkeley-abc/abc"
        rethrow(e)
    finally
        try
            rm(abc_temp_dir; recursive=true, force=true)
        catch cleanup_error
            @warn "Failed to cleanup temporary directory: $cleanup_error"
        end
    end
end

function setup_abc()
    abcbinary = try
        ensure_abc_binary()
    catch e
        error("Failed to setup ABC binary: $e")
    end
    isfile(abcbinary) || error("ABC binary not found at $abcbinary")
    try
        run(`$abcbinary -h`; wait=false)
    catch e
        error("ABC binary are not working properly: $e")
    end
    return abcbinary
end

function setup_quine() end

# ---------------------------------------------------------------------------- #
#                                 print utils                                  #
# ---------------------------------------------------------------------------- #
function _featurename(f::SD.VariableValue)
    return if isnothing(f.i_name)
        f.i_variable isa Feature ?
        "$(f.i_variable)" : "V$(f.i_variable)"
    else
        "$(f.i_name)"
    end
end

# ---------------------------------------------------------------------------- #
#                                 LumenResult                                  #
# ---------------------------------------------------------------------------- #
struct LumenResult
    decision_set::DecisionSet
    info::NamedTuple

    LumenResult(ds, info) = new(ds, info)
    LumenResult(ds) = new(ds, (;))
end

Base.length(lr::LumenResult) = length(lr.decision_set)

# ---------------------------------------------------------------------------- #
#                      extra methods for SoleLogics Atom                       #
# ---------------------------------------------------------------------------- #
@inline get_operator(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.test_operator

@inline get_feature(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.feature

@inline get_threshold(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.threshold

@inline get_i_variable(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.feature.i_variable

# ---------------------------------------------------------------------------- #
#                            operator family utils                             #
# ---------------------------------------------------------------------------- #
const _supported_operators = ((<), (≥), (>), (≤))

@inline _is_lt_family(op) = op === (<) || op === (≤)
@inline _is_gt_family(op) = op === (>) || op === (≥)

# ---------------------------------------------------------------------------- #
#                         sampling threshold                                   #
# ---------------------------------------------------------------------------- #
# Convert any atom's (operator, threshold) into a sampling point in the
# universal :lt / descending convention used by _truths_by_thresholds:
#
#   < t  →  t               boundary at t; sample t (t < t = false, region above)
#   ≤ t  →  t               same boundary
#   ≥ t  →  nextfloat(t)    first value satisfying ≥ t; now  s < s = false  above
#   > t  →  nextfloat(t)    same shift
#
# After sorting all sampling points descending, the Gray-code truth table and
# generate_disjunct (which always emits < and ≥) work correctly because:
#   - an atom "x < t"  fires when the sample is < t, i.e. below boundary t
#   - an atom "x ≥ t"  fires when the sample is ≥ t = nextfloat(t) - ε,
#     which in the descending table corresponds to the same slot as "x < nextfloat(t)"
@inline function _sampling_threshold(op, thr::T) where {T<:AbstractFloat}
    return (op === (≥) || op === (>)) ? nextfloat(thr) : thr
end

# ---------------------------------------------------------------------------- #
#                                 depth utils                                  #
# ---------------------------------------------------------------------------- #
function _extract_atoms_bfs_order(tree::SM.AbstractModel)
    bfs_atoms = SL.Atom{SD.AbstractCondition}[]
    queue = SM.AbstractModel[tree]
    while !isempty(queue)
        current = popfirst!(queue)
        if current isa SM.Branch
            push!(bfs_atoms, antecedent(current))
            push!(queue, SM.posconsequent(current))
            push!(queue, SM.negconsequent(current))
        end
    end
    return bfs_atoms
end

function _take_first_percentage(
    atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
    depth::Float64
)
    n_total = length(atoms)
    n_to_take = Int(ceil(n_total * depth))
    return @view atoms[1:min(n_to_take, n_total)]
end

# ---------------------------------------------------------------------------- #
#                              thresholds utils                                #
# ---------------------------------------------------------------------------- #
@inline _atoms_for_feature(
    atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
    feat::Symbol
) = filter(a -> Symbol(SM.featurename(get_feature(a))) == feat, atoms)

function _truths_by_thresholds(thresholds::Vector{<:Float})
    ntruths = length(thresholds)
    truths = Vector{BitVector}(undef, ntruths + 1)
    @inbounds for i = 1:ntruths+1
        truths[i] = BitVector(undef, ntruths)
        val = 2^(i - 1) - 1
        for j = 1:ntruths
            truths[i][j] = !((val >> (j - 1)) & 1 == 1)
        end
    end
    return truths
end

@inline _truths_by_thresholds(
    thresholds::Vector{T}
) where {T<:Vector{<:Float}} = _truths_by_thresholds.(thresholds)

function _truths_by_thresholds(value::Float, thresholds::Vector{<:Float})
    isnan(value) && return BitVector()
    idx = findfirst(==(value), thresholds)
    return isnothing(idx) ?
           falses(length(thresholds)) :
           _truths_by_thresholds(thresholds)[idx]
end

@inline _truths_by_thresholds(
    values::Tuple{Vararg{<:Float}},
    thresholds::Vector{T}
) where {T<:Vector{<:Float}} = _truths_by_thresholds.(values, thresholds)

# All thresholds are now always :lt / descending sampling points.
# The boundary point covers the region x < t_min (below the smallest threshold).
function _thrs_with_boundary(thresholds::Vector{T}) where {T<:Float}
    isempty(thresholds) && return T[NaN]
    nthrs = length(thresholds)
    result = Vector{T}(undef, nthrs + 1)
    result[1:nthrs] .= thresholds
    result[end] = prevfloat(last(thresholds))
    return result
end

@inline _thrs_with_boundary(
    thresholds::Vector{T}
) where {T<:Vector{<:Float}} = _thrs_with_boundary.(thresholds)

# keep old 2-arg signature for callers that still pass op_families
@inline _thrs_with_boundary(
    thresholds::Vector{T},
    ::Vector{Symbol}
) where {T<:Vector{<:Float}} = _thrs_with_boundary(thresholds)

# ---------------------------------------------------------------------------- #
#                              generate disjuncts                              #
# ---------------------------------------------------------------------------- #
function push_disjunct!(
    disjuncts::Vector{SL.Atom},
    i::Int,
    featurename::Symbol,
    operator,
    threshold::Real
)
    feature = SD.VariableValue(i, featurename)
    mc = SD.ScalarMetaCondition(feature, operator)
    condition = SD.ScalarCondition(mc, threshold)
    push!(disjuncts, SL.Atom(condition))
end

# All sampling thresholds are descending; emit < and ≥ only.
function generate_disjunct(
    truths::Vector{BitVector},
    thresholds::Vector{T},
    features::Vector{Symbol},
    op_families::Vector{Symbol}
) where {T<:Vector{<:Float}}
    disjuncts = Vector{SL.Atom}()

    @inbounds for i in eachindex(thresholds)
        isempty(thresholds[i]) && continue
        idx0 = findall(x -> !x, truths[i])
        idx1 = findall(identity, truths[i])

        # upper bound: x < thresholds[max(idx0)]
        isempty(idx0) ||
            push_disjunct!(disjuncts, i, features[i], <, thresholds[i][maximum(idx0)])
        # lower bound: ¬(x < thresholds[min(idx1)]) → emitted as ≥
        # but we want only <, so we emit x < thresholds[min(idx1)-1] negated...
        # actually just keep ≥ for lower bound — ABC handles it fine
        isempty(idx1) ||
            push_disjunct!(disjuncts, i, features[i], ≥, thresholds[i][minimum(idx1)])
    end

    return disjuncts
end

# ---------------------------------------------------------------------------- #
#              Lazy column view over Iterators.product(thrs...)                #
# ---------------------------------------------------------------------------- #
struct _ProductColumn{T,V<:AbstractVector{<:AbstractVector{T}}} <: AbstractVector{T}
    levels::V
    lens::Vector{Int}
    strides::Vector{Int}
    j::Int
    nrows::Int
end

Base.IndexStyle(::Type{<:_ProductColumn}) = IndexLinear()
Base.size(c::_ProductColumn) = (c.nrows,)
Base.axes(c::_ProductColumn) = (Base.OneTo(c.nrows),)
Base.eltype(::Type{_ProductColumn{T,V}}) where {T,V} = T

@inline function Base.getindex(c::_ProductColumn, i::Int)
    @boundscheck checkbounds(c, i)
    idx = ((i - 1) ÷ c.strides[c.j]) % c.lens[c.j] + 1
    @inbounds return c.levels[c.j][idx]
end

function _product_columntable(
    thrs_with_p::Vector{<:AbstractVector{T}},
    featurenames::Vector{Symbol},
) where {T}
    n = length(thrs_with_p)
    lens = length.(thrs_with_p)
    strides = ones(Int, n)
    @inbounds for j in 2:n
        strides[j] = strides[j-1] * lens[j-1]
    end
    nrows = prod(lens)
    names = Tuple(featurenames)
    cols = ntuple(j -> _ProductColumn{T,typeof(thrs_with_p)}(
            thrs_with_p, lens, strides, j, nrows
        ), n)
    return NamedTuple{names}(cols)
end

# ---------------------------------------------------------------------------- #
#                          extract rules data struct                           #
# ---------------------------------------------------------------------------- #
struct ExtractRulesData{
    P,
    C<:Base.Iterators.ProductIterator,
    T<:Vector{<:Float},
    F<:SM.Label,
    L<:SM.Label
}
    predictions::P
    combinations::C
    thresholds::Vector{T}
    featurenames::Vector{F}
    classnames::AbstractVector{L}
    op_families::Vector{Symbol}   # always :lt for every feature

    ExtractRulesData(
        predictions::P,
        combinations::C,
        thresholds::Vector{T},
        featurenames::Vector{F},
        classnames::AbstractVector{L},
        op_families::Vector{Symbol}
    ) where {
        P,
        C<:Base.Iterators.ProductIterator,
        T<:Vector{<:Float},
        F<:SM.Label,
        L<:SM.Label
    } = new{P,C,T,F,L}(
        predictions, combinations, thresholds, featurenames, classnames, op_families
    )

    function ExtractRulesData(extractor::LumenConfig, model::SM.AbstractModel)
        depth = get_depth(extractor)

        # STEP 2 — extract and deduplicate atoms
        raw_atoms = if depth < 1.0
            mapreduce(
                vcat, SM.models(model); init=SL.Atom{SD.AbstractCondition}[]
            ) do t
                all_atoms_bfs = _extract_atoms_bfs_order(t)
                _take_first_percentage(all_atoms_bfs, depth)
            end
        else
            SL.atoms(SM.alphabet(model, false))
        end

        atoms = unique(raw_atoms)

        # STEP 3 — validate operators
        let unsupported = unique(
                op for op in get_operator.(atoms)
                if op ∉ _supported_operators
            )
            isempty(unsupported) || throw(ArgumentError(
                "Only '<', '≥', '>', '≤' operators are currently supported. " *
                "Found unsupported operators: $(unsupported).",
            ))
        end

        # STEP 4 — feature names present in atoms
        features = Symbol.(SM.featurename.(unique!(get_feature.(atoms))))

        # STEP 5 — canonical feature list and class labels from model
        featurenames = SM.info(model, :featurenames)
        classnames = unique!(SM.info(model, :supporting_labels))

        type = get_float_type(extractor)

        # STEP 6 — build sampling thresholds, always :lt / descending
        #
        # Each atom (op, t) is mapped to a sampling point s:
        #   op ∈ {<, ≤}  →  s = t
        #   op ∈ {≥, >}  →  s = nextfloat(t)
        #
        # All sampling points for the same feature are deduplicated and sorted
        # descending.  The existing _truths_by_thresholds Gray-code and
        # generate_disjunct (emitting < and ≥) then work correctly for every
        # original operator, because the sampling points encode the correct
        # region boundaries in the :lt sense.
        thresholds = Vector{Vector{type}}(undef, length(featurenames))
        op_families = fill(:lt, length(featurenames))

        @inbounds for i in eachindex(featurenames)
            idx = findfirst(f -> f == featurenames[i], features)
            if isnothing(idx)
                thresholds[i] = type[]
            else
                feat_atoms = _atoms_for_feature(atoms, features[idx])
                sampling_thrs = Set{type}()
                for a in feat_atoms
                    t = type(get_threshold(a))
                    op = get_operator(a)
                    if op === (<) || op === (≤)
                        push!(sampling_thrs, t)
                    else  # ≥ or >
                        push!(sampling_thrs, prevfloat(t))
                        push!(sampling_thrs, t)
                        push!(sampling_thrs, nextfloat(t))
                    end
                end
                thresholds[i] = sort!(collect(sampling_thrs); rev=true)
            end
        end

        # STEP 7 — augment with boundary point (prevfloat of minimum sampling thr)
        thrs_with_p = _thrs_with_boundary(thresholds)

        # STEP 8 — Cartesian product of all sampling vectors
        combinations = Iterators.product(thrs_with_p...)

        # STEP 9 — apply model to all combinations
        tbl = _product_columntable(thrs_with_p, Symbol.(featurenames))
        d = PropositionalLogiset(tbl)

        predictions = get_apply_function(extractor)(
            model,
            d;
            use_multithreads=get_use_multithreads(extractor),
            suppress_parity_warning=true
        )

        #@show predictions


        return ExtractRulesData(
            predictions, combinations, thresholds, featurenames, classnames, op_families
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                   methods                                    #
# ---------------------------------------------------------------------------- #
function get_thresholds(
    e::ExtractRulesData;
    prev_float::Bool=false,
    float_type::Type=Float64
)
    thresholds = [float_type.(t) for t in e.thresholds]
    return prev_float ? _thrs_with_boundary(thresholds) : thresholds
end

@inline get_featurenames(e::ExtractRulesData) = e.featurenames
@inline get_classnames(e::ExtractRulesData) = e.classnames
@inline get_op_families(e::ExtractRulesData) = e.op_families

@inline function get_truths(e::ExtractRulesData, i::Int)
    _truths_by_thresholds(IterTools.nth(e.combinations, i), e.thresholds)
end

function truths_by_groups(e::ExtractRulesData, i::Int)
    idxs = findall(==(e.classnames[i]), e.predictions)
    truths = [get_truths(e, i) for i in idxs]
    return isempty(truths) ? Vector{BitVector}[] : truths
end

# ---------------------------------------------------------------------------- #
#                                  get atoms                                   #
# ---------------------------------------------------------------------------- #
function get_atoms(
    e::ExtractRulesData;
    grouped::Bool=false,
    float_type::Type=Float64
)
    thresholds = get_thresholds(e; prev_float=true, float_type)
    features = get_featurenames(e)
    op_families = get_op_families(e)

    return if grouped
        truths = vcat(get_truths(e)...)
        get_atoms(truths, thresholds, features, op_families)
    else
        [get_atoms(truths, thresholds, features, op_families)
         for truths in get_truths(e)]
    end
end

function get_atoms(e::ExtractRulesData, c::SM.Label; float_type::Type=Float64)
    i = findfirst(==(c), get_classnames(e))
    isnothing(i) ? nothing : get_atoms(e, i; float_type)
end

function get_atoms(e::ExtractRulesData, i::Int; float_type::Type=Float64)
    truths = truths_by_groups(e, i)
    thresholds = get_thresholds(e; prev_float=false, float_type)
    featurenames = get_featurenames(e)
    op_families = get_op_families(e)
    get_atoms(truths, thresholds, featurenames, op_families)
end

function get_atoms(
    truths::Vector{Vector{BitVector}},
    thresholds::Vector{T},
    featurenames::Vector{Symbol},
    op_families::Vector{Symbol}
) where {T<:Vector{<:Float}}
    conjuncts = Vector{Vector{SL.Atom}}(undef, length(truths))
    Threads.@threads for i in eachindex(truths)
        conjuncts[i] =
            generate_disjunct(truths[i], thresholds, featurenames, op_families)
    end
    return conjuncts
end

# ---------------------------------------------------------------------------- #
#                                get conjuncts                                 #
# ---------------------------------------------------------------------------- #
function get_conjuncts(e::ExtractRulesData, i::Int)
    atoms = get_atoms(e, i)
    [get_conjuncts(atom) for atom in atoms]
end

function get_conjuncts(e::ExtractRulesData, c::SM.Label)
    i = findfirst(==(c), get_classnames(e))
    isnothing(i) ? nothing : get_conjuncts(e, i)
end

@inline get_conjuncts(a::Vector{Vector{SL.Atom}}) = get_conjuncts.(a)
@inline get_conjuncts(a::Vector{SL.Atom}) = isempty(a) ?
                                            ⊤ : SL.LeftmostConjunctiveForm{SL.Literal}(SL.Literal.(a))

# ---------------------------------------------------------------------------- #
#                                get formulas                                  #
# ---------------------------------------------------------------------------- #
@inline get_formula(e::ExtractRulesData, i::Int) =
    get_formula(get_conjuncts(e, i))

function get_formula(e::ExtractRulesData, c::SM.Label)
    i = findfirst(==(c), get_classnames(e))
    isnothing(i) ? nothing : get_formula(e, i)
end

@inline get_formula(grouped_conj::Vector{SL.LeftmostConjunctiveForm{SL.Atom}}) =
    SL.LeftmostDisjunctiveForm{SL.LeftmostConjunctiveForm{SL.Literal}}(
        grouped_conj, true
    )

# ---------------------------------------------------------------------------- #
#                           dnf minimization refine                            #
# ---------------------------------------------------------------------------- #
function _refine_dnf(
    terms::Vector{<:Union{SL.LeftmostConjunctiveForm{SL.Atom},SyntaxStructure}}
)
    length(terms) ≤ 1 && return terms
    all_bounds = map(term -> SD.extract_term_bounds(term; silent=true), terms)
    keep_mask = map(enumerate(all_bounds)) do (i, bounds_i)
        !any(j -> i ≠ j && SD.strictly_dominates(
                all_bounds[j], bounds_i), eachindex(all_bounds))
    end
    kept_terms = terms[keep_mask]
    return isempty(kept_terms) ? terms : kept_terms
end

# ---------------------------------------------------------------------------- #
#                              minimization core                               #
# ---------------------------------------------------------------------------- #
function run_minimization(
    ::Val{:abc},
    extractor::LumenConfig,
    atoms::Vector{Vector{SL.Atom}}
)
    minimized_formula =
        SD.abc_minimize(
            atoms,
            get_binary(extractor);
            fast=1,
            depth=get_depth(extractor),
            float_type=get_float_type(extractor)
        )
    return _refine_dnf(minimized_formula)
end

function run_minimization(
    ::Val{:mitespresso},
    extractor::LumenConfig,
    atoms::Vector{Vector{SL.Atom}}
)
    minimized_formula =
        SD.espresso_minimize(
            atoms,
            get_binary(extractor);
            depth=get_depth(extractor),
            float_type=get_float_type(extractor)
        )
    return _refine_dnf(minimized_formula)
end

# ---------------------------------------------------------------------------- #
#                                    lumen                                     #
# ---------------------------------------------------------------------------- #
function lumen(
    config::LumenConfig,
    model::SM.AbstractModel
)
    float_type = get_float_type(config)

    extractrulesdata = ExtractRulesData(config, model)

    classes = get_classnames(extractrulesdata)
    nclasses = length(classes)

    formulas =
        Vector{Vector{Union{
            SL.LeftmostConjunctiveForm{SL.Atom{float_type}},
            SyntaxStructure
        }}}(undef, nclasses)

    Threads.@threads for i in 1:nclasses
        atoms = get_atoms(extractrulesdata, i; float_type)
        formulas[i] = isempty(atoms) ?
                      SL.Atom{SD.AbstractCondition}[] :
                      run_minimization(
            Val(get_minimization_scheme(config)), config, atoms
        )
    end

    valid_mask = .!isempty.(formulas)
    formulas = formulas[valid_mask]
    classes = classes[valid_mask]

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end

function lumen(
    config::LumenConfig,
    model::Vector{SM.AbstractModel}
)
    ds = map(model) do m
        lumen(config, m)
    end
    return LumenResult(ds)
end

function lumen(
    model::SM.AbstractModel,
    args...;
    kwargs...
)
    lumen(LumenConfig(; kwargs...), model)
end

function lumen(
    model::Vector{SM.AbstractModel},
    args...;
    kwargs...
)
    ds = map(model) do m
        lumen(m, args...; kwargs...)
    end
    return LumenResult(ds)
end

end