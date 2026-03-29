using Test
using SoleXplorer

using MLJ
using DataFrames
using Random
# using Downloads

Xc, yc = @load_iris
Xc = DataFrame(Xc)

seed=11

resampling = Holdout(fraction_train=0.8, shuffle=true)
# resampling = pCV(nfolds=2, fraction_train=0.6, shuffle=true)

# model = DecisionTreeClassifier()
model = RandomForestClassifier(n_trees=5)
# model = RandomForestClassifier(n_trees=20)
# model = RandomForestClassifier(n_trees=100)

modelc = solexplorer(Xc, yc; model, resampling, seed)

# ---------------------------------------------------------------------------- #
# using SolePostHoc
using SoleModels
using SoleData

# using SoleData: PLA

# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #
abstract type AbstractRuleExtractorInfo end

const AbstractCondition       = SoleData.AbstractCondition
const ScalarCondition         = SoleData.ScalarCondition
const ScalarMetaCondition     = SoleData.ScalarMetaCondition
const VariableValue           = SoleData.VariableValue

const Atom                    = SoleLogics.Atom
const LeftmostConjunctiveForm = SoleLogics.LeftmostConjunctiveForm
const LeftmostDisjunctiveForm = SoleLogics.LeftmostDisjunctiveForm
const Literal                 = SoleLogics.Literal

const AbstractModel           = SoleModels.AbstractModel
const Label                   = SoleModels.Label

# ---------------------------------------------------------------------------- #
#                   initial minimization algorithms setup                      #
# ---------------------------------------------------------------------------- #
function setup_espresso() end
function setup_boom() end

function setup_abc() 
    # auto setup ABC binary if not specified
    abcbinary = try
        joinpath(SoleData.load(SoleData.ABCLoader()), "abc")
    catch e
        error("Failed to setup ABC binary: $e")
    end

    # verify that binary exists and is executable
    isfile(abcbinary) || error("ABC binary not found at $abcbinary")

    #abcbinary = ensure_abc_binary(; force_rebuild = force_rebuild_abc) emergency version

    # test that ABC binary is working
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
_featurename(f::SoleData.VariableValue) = isnothing(f.i_name) ? "V$(f.i_variable)" : "[$(f.i_name)]"

# ---------------------------------------------------------------------------- #
#                                 Lumen struct                                 #
# ---------------------------------------------------------------------------- #
struct LuminoRuleExtractor <: AbstractRuleExtractorInfo
    minimization_scheme :: Symbol
    binary              :: String
    depth               :: Float64
    vertical            :: Float64
    horizontal          :: Float64
    minimization_kwargs :: NamedTuple 
    filt_alphabet       :: Base.Callable
    apply_function      :: Base.Callable
    importance          :: Vector
    check_opt           :: Bool
    check_alphabet      :: Bool
    rng                 :: AbstractRNG

    function LuminoRuleExtractor(;
        minimization_scheme :: Symbol        = :abc,
        depth               :: Float64       = 1.0,
        vertical            :: Float64       = 1.0,
        horizontal          :: Float64       = 1.0,
        minimization_kwargs :: NamedTuple    = (;),
        filt_alphabet       :: Base.Callable = identity,
        apply_function      :: Base.Callable = identity,
        importance          :: Vector        = Float64[],
        check_opt           :: Bool          = false,
        check_alphabet      :: Bool          = false,
        rng                 :: AbstractRNG   = Random.TaskLocalRNG()
    )
        # validate coverage parameters - must be positive and ≤ 1.0
        # these parameters control the proportion of instances that must be covered by rules
        if vertical ≤ 0.0 || vertical > 1.0 || horizontal ≤ 0.0 || horizontal > 1.0 || depth ≤ 0.0 || depth > 1.0
            throw(
                ArgumentError(
                    "vertical, depth and horizontal parameters must be in range (0.0, 1.0]. " *
                    "Got vertical=$(vertical), depth=$(depth), horizontal=$(horizontal). " *
                    "These parameters control rule coverage and must be meaningful proportions.",
                ),
            )
        end

        # validate minimization scheme
        valid_schemes = Dict(
            :mitespresso  => setup_espresso(),
            :boom         => setup_boom(),
            :abc          => setup_abc(),
            :abc_balanced => setup_abc(),
            :abc_thorough => setup_abc(),
            :quine        => setup_quine(),
            :quine_naive  => setup_quine()
        )
        
        if minimization_scheme ∉ keys(valid_schemes)
            throw(
                ArgumentError(
                    "minimization_scheme must be one of: $(keys(valid_schemes) |> collect). " *
                    "Got: $(minimization_scheme). "
                )
            )
        end

        binary = valid_schemes[minimization_scheme]

        new(
            minimization_scheme,
            binary,
            depth,
            vertical,
            horizontal,
            minimization_kwargs,
            filt_alphabet,
            apply_function,
            importance,
            check_opt,
            check_alphabet,
            rng
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
@inline get_minimization_scheme(r::LuminoRuleExtractor) = r.minimization_scheme
@inline get_binary(r::LuminoRuleExtractor)              = r.binary
@inline get_depth(r::LuminoRuleExtractor)               = r.depth
@inline get_vertical(r::LuminoRuleExtractor)            = r.vertical
@inline get_horizontal(r::LuminoRuleExtractor)          = r.horizontal
@inline get_minimization_kwargs(r::LuminoRuleExtractor) = r.minimization_kwargs
@inline get_filt_alphabet(r::LuminoRuleExtractor)       = r.filt_alphabet
@inline get_apply_function(r::LuminoRuleExtractor)      = r.apply_function
@inline get_importance(r::LuminoRuleExtractor)          = r.importance
@inline get_check_opt(r::LuminoRuleExtractor)           = r.check_opt
@inline get_check_alphabet(r::LuminoRuleExtractor)      = r.check_alphabet
@inline get_rng(r::LuminoRuleExtractor)                 = r.rng

# ---------------------------------------------------------------------------- #
#                      extra methods for SoleLogics Atom                       #
# ---------------------------------------------------------------------------- #
@inline get_operator(atom::Atom{<:AbstractCondition})   = atom.value.metacond.test_operator
@inline get_feature(atom::Atom{<:AbstractCondition})    = atom.value.metacond.feature
@inline get_threshold(atom::Atom{<:AbstractCondition})  = atom.value.threshold
@inline get_i_variable(atom::Atom{<:AbstractCondition}) = atom.value.metacond.feature.i_variable

# ---------------------------------------------------------------------------- #
#                                 depth utils                                  #
# ---------------------------------------------------------------------------- #
function _extract_atoms_bfs_order(tree::AbstractModel)
    bfs_atoms = Atom{AbstractCondition}[]
    queue     = AbstractModel[tree]

    while !isempty(queue)
        current = popfirst!(queue)

        if current isa SoleModels.Branch
            push!(bfs_atoms, antecedent(current))
            push!(queue, SoleModels.posconsequent(current))
            push!(queue, SoleModels.negconsequent(current))
        end
    end

    return bfs_atoms
end

function _take_first_percentage(
    atoms :: Vector{<:Atom{<:ScalarCondition}},
    depth :: Float64
)
    n_total   = length(atoms)
    n_to_take = Int64(ceil(n_total * depth))

    return atoms[1:min(n_to_take, n_total)]
end

# ---------------------------------------------------------------------------- #
#                              thresholds utils                                #
# ---------------------------------------------------------------------------- #
@inline _atoms_for_feature(
    atoms :: Vector{<:Atom{<:ScalarCondition}},
    feat  :: Symbol
) = filter(a -> SoleModels.featurename(get_feature(a)) == feat, atoms)

function _truths_by_thresholds(thresholds::Vector{Float64})
    ntruths = length(thresholds)
    truths  = Vector{BitVector}(undef, ntruths+1)

    @inbounds for i = 1:ntruths+1
        truths[i] = BitVector(undef, ntruths)
        val = 2^(i-1) - 1
        for j = 1:ntruths
            truths[i][j] = !((val >> (j-1)) & 1 == 1)
        end
    end

    return truths
end

@inline  _truths_by_thresholds(thresholds::Vector{Vector{Float64}}) = _truths_by_thresholds.(thresholds)

function _truths_by_thresholds(value::Float64, thresholds::Vector{Float64})
    isnan(value) && return BitVector()
    
    idx = findfirst(==(value), thresholds)
    return isnothing(idx) ? falses(length(thresholds)) : _truths_by_thresholds(thresholds)[idx]
end
@inline  _truths_by_thresholds(values::Tuple{Vararg{Float64}}, thresholds::Vector{Vector{Float64}}) =
    _truths_by_thresholds.(values, thresholds)

function _thrs_with_prevfloat(thresholds::Vector{Float64})
    isempty(thresholds) && return [NaN]
    
    nthrs  = length(thresholds)
    result = Vector{Float64}(undef, nthrs + 1)
    @inbounds for i = 1:nthrs
        result[i] = thresholds[i]
    end
    result[end] = prevfloat(last(thresholds))
    return result
end
@inline  _thrs_with_prevfloat(thresholds::Vector{Vector{Float64}}) = _thrs_with_prevfloat.(thresholds)

# ---------------------------------------------------------------------------- #
#                              generate disjunts                               #
# ---------------------------------------------------------------------------- #
function push_disjunct!(
    disjuncts   :: Vector{Atom},
    i           :: Int64,
    featurename :: Symbol,
    operator    :: Union{typeof(<), typeof(≥)},
    threshold   :: Real
)
    feature   = VariableValue(i, featurename)
    mc        = ScalarMetaCondition(feature, operator)
    condition = ScalarCondition(mc, threshold)

    push!(disjuncts, Atom(condition))
end

function generate_disjunct(
    truths     :: Vector{BitVector},
    thresholds :: Vector{Vector{Float64}},
    features   :: Vector{Symbol}
)
    disjuncts = Vector{Atom}()

    @inbounds for i in eachindex(thresholds)
        idx0 = findall(x -> !x,  truths[i])
        idx1 = findall(identity, truths[i])

        isempty(idx0) ||
            push_disjunct!(disjuncts, i, features[i], <, thresholds[i][maximum(idx0)])
        isempty(idx1) ||
            push_disjunct!(disjuncts, i, features[i], ≥, thresholds[i][minimum(idx1)])
    end

    return disjuncts
end

# ---------------------------------------------------------------------------- #
#                          extract rules data struct                           #
# ---------------------------------------------------------------------------- #
struct ExtractRulesData
    grp_truths :: Vector{Vector{Vector{BitVector}}}
    thresholds :: Vector{Vector{Float64}}
    features   :: Vector{<:Label}
    classnames :: Vector{<:Label}

    ExtractRulesData(
        grp_truths :: Vector{Vector{Vector{BitVector}}},
        thresholds :: Vector{Vector{Float64}},
        features   :: Vector{<:Label},
        classnames :: Vector{<:Label}
    ) = new(grp_truths, thresholds, features, classnames)

    function ExtractRulesData(extractor::LuminoRuleExtractor, model::AbstractModel)
        depth = get_depth(extractor)

        atoms = unique!(if depth < 1.0
            mapreduce(vcat, SoleModels.models(model); init=Atom{AbstractCondition}[]) do t
                all_atoms_bfs = _extract_atoms_bfs_order(t)
                _take_first_percentage(all_atoms_bfs, depth)
            end
        else
            SoleLogics.atoms(SoleModels.alphabet(model, false))
        end)

        # validate supported operators
        isempty(unique(op for op in get_operator.(atoms) if op != (<))) || throw(
            ArgumentError(
                "Only '<' operator is currently supported. " *
                "Found unsupported operators: $(unsupported). " *
                "This limitation may be addressed in future versions. " *
                "Consider preprocessing your model to use only '<' conditions.",
            ),
        )

        features     = SoleModels.featurename.(unique!(get_feature.(atoms)))
        featurenames = SoleModels.info(model, :featurenames)
        classnames   = unique!(SoleModels.info(model, :supporting_labels))

        thresholds   = Vector{Vector{Float64}}(undef, length(featurenames))

        @inbounds for i in eachindex(featurenames)
            idx = findfirst(f -> f == featurenames[i], features)
            thresholds[i] = isnothing(idx) ? 
                Float64[] :
                sort!(get_threshold.(_atoms_for_feature(atoms, features[idx])), rev=true)
        end

        thrs_with_p  = _thrs_with_prevfloat(thresholds)

        combinations = collect(Iterators.product(thrs_with_p...))
        predictions  = get_apply_function(extractor)(
            model,
            scalarlogiset(DataFrame(combinations, featurenames), allow_propositional=true);
            suppress_parity_warning=true
        )

        truths = Vector{Vector{BitVector}}(undef, length(combinations))
        Threads.@threads for i in eachindex(combinations)
            truths[i] = _truths_by_thresholds(combinations[i], thresholds)
        end

        grp_truths = Vector{Vector{Vector{BitVector}}}(undef, length(classnames))
        Threads.@threads for i in eachindex(classnames)
            indices = findall(==(classnames[i]), predictions)
            grp_truths[i] = truths[indices]
        end

        return ExtractRulesData(grp_truths, thresholds, featurenames, classnames)
    end
end

# ---------------------------------------------------------------------------- #
#                                   methods                                    #
# ---------------------------------------------------------------------------- #
@inline  get_grouped_truths(e::ExtractRulesData) = e.grp_truths

function get_thresholds(e::ExtractRulesData; prev_float::Bool=false)
    thresholds = e.thresholds
    return prev_float ? _thrs_with_prevfloat(thresholds) : thresholds
end

@inline  get_features(e::ExtractRulesData)   = e.features
@inline  get_classnames(e::ExtractRulesData) = e.classnames

function get_grouped_truths(e::ExtractRulesData, c::Label)
    i = findfirst(g -> get_classnames(g) == c, e.grp_truths)
    isnothing(i) ? nothing : get_grouped_truths(e, i)
end

@inline  get_grouped_truths(e::ExtractRulesData, i::Int64) = e.grp_truths[i]

@inline  get_truths(e::ExtractRulesData) = [get_truths(e, i) for i in eachindex(get_classnames(e))]
@inline  get_truths(e::ExtractRulesData, i::Int64) = get_grouped_truths(e, i)

@inline  get_truth(e::ExtractRulesData, i::Int64, j::Int64) = get_grouped_truths(e, i)[j]
@inline  get_truth(e::ExtractRulesData, i::Int64) = get_truths(e)[i]

# ---------------------------------------------------------------------------- #
#                                  get atoms                                   #
# ---------------------------------------------------------------------------- #
function get_atoms(e::ExtractRulesData; grouped::Bool=false)
    thresholds = get_thresholds(e, prev_float=true)
    features   = get_features(e)

    return if grouped
        truths = vcat(get_truths(e)...)
        get_atoms(truths, thresholds, features)
    else
        [get_atoms(truths, thresholds, features) for truths in get_truths(e)]
    end
end

function get_atoms(e::ExtractRulesData, c::Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_atoms(e, i)
end

function get_atoms(e::ExtractRulesData, i::Int64)
    truths     = get_truths(e, i)
    thresholds = get_thresholds(e, prev_float=false)
    features   = get_features(e)

    get_atoms(truths, thresholds, features)
end

function get_atoms(
    truths     :: Vector{Vector{BitVector}},
    thresholds :: Vector{Vector{Float64}},
    features   :: Vector{Symbol}
)
    conjuncts = Vector{Vector{Atom}}(undef, length(truths))

    Threads.@threads for i in eachindex(truths)
        conjuncts[i] = generate_disjunct(truths[i], thresholds, features)
    end

    return conjuncts
end

# ---------------------------------------------------------------------------- #
#                                get conjuncts                                 #
# ---------------------------------------------------------------------------- #
function get_conjuncts(e::ExtractRulesData, i::Int64)
    atoms = get_atoms(e, i)
    [get_conjuncts(atom) for atom in atoms]
end

function get_conjuncts(e::ExtractRulesData, c::Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_conjuncts(e, i)
end

@inline  get_conjuncts(a::Vector{Vector{Atom}}) = get_conjuncts.(a)
@inline  get_conjuncts(a::Vector{Atom}) = isempty(a) ? ⊤ : LeftmostConjunctiveForm{Literal}(Literal.(a))

# ---------------------------------------------------------------------------- #
#                                get formulas                                  #
# ---------------------------------------------------------------------------- #
@inline  get_formula(e::ExtractRulesData, i::Int64) = get_formula(get_conjuncts(e, i))

function get_formula(e::ExtractRulesData, c::Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_formula(e, i)
end

@inline  get_formula(grouped_conj::Vector{LeftmostConjunctiveForm{Atom}}) =
    LeftmostDisjunctiveForm{LeftmostConjunctiveForm{Literal}}(grouped_conj, true)

# ---------------------------------------------------------------------------- #
#                           dnf minimization refine                            #
# ---------------------------------------------------------------------------- #
function _refine_dnf(terms::Vector{<:Union{LeftmostConjunctiveForm{Atom}, SyntaxStructure}})  
    length(terms) ≤ 1 && return terms
    
    all_bounds = map(term -> SoleData.extract_term_bounds(term; silent=true), terms)
    
    # find terms not strictly dominated by any other term
    keep_mask = map(enumerate(all_bounds)) do (i, bounds_i)
        !any(j -> i ≠ j && SoleData.strictly_dominates(all_bounds[j], bounds_i), eachindex(all_bounds))
    end
    
    kept_terms = terms[keep_mask]
    
    # safety check: never return empty formula
    return isempty(kept_terms) ? terms : kept_terms
end

# ---------------------------------------------------------------------------- #
#                              minimization core                               #
# ---------------------------------------------------------------------------- #
function run_minimization(extractor::LuminoRuleExtractor, atoms::Vector{Vector{Atom}})
    minimized_formula = abc_minimize(extractor, atoms; fast = 1, depth=get_depth(extractor))

    if get_minimization_scheme(extractor) in
        [:mitespresso, :boom, :abc, :abc_balanced, :abc_thorough]
        minimized_formula = _refine_dnf(minimized_formula)
    end

    minimized_formula
end

# ---------------------------------------------------------------------------- #
#                                    lumen                                     #
# ---------------------------------------------------------------------------- #
function lumen(
    extractor :: LuminoRuleExtractor,
    model     :: AbstractModel
)
    # handle special test modes -> TODO
    # if !isnothing(config.testott) || !isnothing(config.alphabetcontroll)
    #     return handle_test_modes(model, config)
    # end

    # Determine apply function -> TODO
    # apply_function = determine_apply_function(model, config.apply_function)
    # config = @set config.apply_function = apply_function

    # extract conjuncts
    extractrulesdata = ExtractRulesData(extractor, model)
    nclasses         = length(get_classnames(extractrulesdata))

    # formulas         = Vector{SoleLogics.LeftmostLinearForm}(undef, nclasses)
    formulas         = Vector{Vector{Union{LeftmostConjunctiveForm{Atom}, SyntaxStructure}}}(undef, nclasses)

    Threads.@threads for i in 1:nclasses
        atoms        = get_atoms(extractrulesdata, i)
        formulas[i]  = run_minimization(extractor, atoms)
    end

    return LeftmostDisjunctiveForm.(formulas)
end

function lumen(
    extractor :: LuminoRuleExtractor,
    model     :: Vector{AbstractModel}
)
    map(enumerate(model)) do (_, m)
        lumen(extractor, m)
    end
end

# ---------------------------------------------------------------------------- #
#                                  abc utils                                   #
# ---------------------------------------------------------------------------- #
function clean_abc_output(raw_pla::String)
    lines = split(raw_pla, '\n')
    pla_lines = filter(l -> !isempty(strip(l)) &&
        (startswith(l, '.') || occursin(r"^[01\-]+ ", l)), lines)
    return join(pla_lines, '\n')
end

_patchnothing(v, d) = isnothing(v) ? d : v

# ---------------------------------------------------------------------------- #
#                                abc minimize                                  #
# ---------------------------------------------------------------------------- #
function abc_minimize(
    :: LuminoRuleExtractor,
    atoms      :: Vector{Vector{Atom}};
    fast::Int64 = 1,
    allow_scalar_range_conditions::Bool = false,
    depth::Float64=1.0
)
    abcbinary = get_binary(extractor)

    # convert formula to pla string format
    pla_string, fnames = PLA.formula_to_pla(
        atoms;
        allow_scalar_range_conditions,
        removewhitespaces=true,
        pretty_op=false
    )

    # Create temporary files for input/output
    mktempdir() do tmp
        inputfile  = joinpath(tmp, "in.pla")
        outputfile = joinpath(tmp, "out.pla")
        
        write(inputfile, pla_string)

        abc_commands = if fast == 1
            "read $inputfile; strash; collapse; write $outputfile"
        elseif fast == 0
            "read $inputfile; strash; balance; rewrite; refactor; balance; rewrite -z; collapse; sop; fx; strash; balance; collapse; write $outputfile"
        else
            "read $inputfile; sop; strash; dc2; collapse; strash; dc2; collapse; sop; write $outputfile"
        end

        run(`$abcbinary -c $abc_commands`)

        # isfile(outputfile) || return conjuncts

        minimized_pla_raw = read(outputfile, String)
        minimized_pla_raw = replace(minimized_pla_raw, ">=" => "≥")
        isempty(strip(minimized_pla_raw)) && return atoms

        minimized_pla = clean_abc_output(minimized_pla_raw)
        conditionstype = allow_scalar_range_conditions ? SoleData.RangeScalarCondition : ScalarCondition

        return PLA.pla_to_formula(minimized_pla, fnames; conditionstype)
    end
end

# ---------------------------------------------------------------------------- #
# test con depth < 1.0 e > 1.0 per extract_atoms
seed = 11
scalar_range_conds = false
encoding=:univariate

extractor = LuminoRuleExtractor(
    minimization_scheme = :abc,
    depth               = 1.0,
    vertical            = 1.0,
    horizontal          = 1.0,
    minimization_kwargs = (;),
    filt_alphabet       = identity,
    apply_function      = SoleModels.apply,
    importance          = Float64[],
    check_opt           = false,
    check_alphabet      = false,
    rng                 = Random.Xoshiro(seed)
)

# @btime begin
    formulas = lumen(extractor, modelc.sole);
# end;

test = SoleXplorer.extractrules(LumenRuleExtractor(), (;), modelc.ds, modelc.sole)

# Lumen 20 trees
# 28.552 s (334369836 allocations: 26.52 GiB)
# 11.468 s (90994927  allocations: 22.66 GiB)
# 10.867 s (89848762  allocations: 22.61 GiB)
# 7.185  s (39312210  allocations: 13.78 GiB)
# 5.259  s (15941330  allocations: 6.22 GiB)
# 5.104  s (15872751  allocations: 6.21 GiB)
# 4.922  s (15846072  allocations: 6.21 GiB)
# 3.339  s (15032124  allocations: 6.20 GiB)
# 2.935  s (14954190  allocations: 6.19 GiB)

# Reference result
# [LeftmostDisjunctiveForm with 6 grandchildren:
#         ([V3] ≥ 4.85) ∧ ([V4] ≥ 1.75)
#         ([V3] ≥ 4.95) ∧ ([V4] < 1.55)
#         ([V3] ≥ 4.85) ∧ ([V3] < 4.95) ∧ ([V4] ≥ 1.65)
#         ([V3] ≥ 2.45) ∧ ([V3] < 4.85) ∧ ([V4] ≥ 1.65) ∧ ([V4] < 1.75)
#         [V3] ≥ 5.449999999999999
#         ([V3] ≥ 2.45) ∧ ([V4] ≥ 1.75) ∧ ([V2] < 3.1)
# , DNF with 5 disjuncts and literals of type Atom:
#         ([V3] ≥ 4.85) ∧ ([V3] < 4.95) ∧ ([V4] < 1.65)
#         ([V3] ≥ 2.45) ∧ ([V3] < 4.85) ∧ ([V4] < 1.65)
#         ([V3] ≥ 4.85) ∧ ([V3] < 5.449999999999999) ∧ ([V4] ≥ 1.55) ∧ ([V4] < 1.65)
#         ([V3] ≥ 4.95) ∧ ([V3] < 5.449999999999999) ∧ ([V4] ≥ 1.55) ∧ ([V4] < 1.75)
#         ([V3] ≥ 2.45) ∧ ([V3] < 4.85) ∧ ([V4] ≥ 1.75) ∧ ([V2] ≥ 3.1)
# , LeftmostDisjunctiveForm with 1 Atom{typeof(<)}}} grandchildren:
#         [V3] < 2.45
# ]

#  LumenResult(▣
# ├[1/3] (V3 < 2.45)  ↣  CategoricalArrays.CategoricalValue{String, UInt32}["setosa"] : NamedTuple()
# ├[2/3] 
    # ((V3 ≥ 4.85) ∧ (V4 ≥ 1.75)) ∨ 
    # ((V3 ≥ 4.95) ∧ (V4 < 1.55)) ∨ 
    # ((V3 ≥ 4.85) ∧ (V3 < 4.95) ∧ (V4 ≥ 1.65)) ∨ 
    # ((V3 ≥ 2.45) ∧ (V3 < 4.85) ∧ (V4 ≥ 1.65) ∧ (V4 < 1.75)) ∨ 
    # (V3 ≥ 5.449999999999999) ∨ 
    # ((V2 < 3.1) ∧ (V3 ≥ 2.45) ∧ (V4 ≥ 1.75))  ↣  CategoricalArrays.CategoricalValue{String, UInt32}["virginica"] : NamedTuple()
# └[3/3] 
    # ((V3 ≥ 4.85) ∧ (V3 < 4.95) ∧ (V4 < 1.65)) ∨ 
    # ((V3 ≥ 2.45) ∧ (V3 < 4.85) ∧ (V4 < 1.65)) ∨ 
    # ((V3 ≥ 4.85) ∧ (V3 < 5.449999999999999) ∧ (V4 ≥ 1.55) ∧ (V4 < 1.65)) ∨ 
    # ((V3 ≥ 4.95) ∧ (V3 < 5.449999999999999) ∧ (V4 ≥ 1.55) ∧ (V4 < 1.75)) ∨ 
    # ((V2 ≥ 3.1) ∧ (V3 ≥ 2.45) ∧ (V3 < 4.85) ∧ (V4 ≥ 1.75))  ↣  CategoricalArrays.CategoricalValue{String, UInt32}["versicolor"] : NamedTuple()