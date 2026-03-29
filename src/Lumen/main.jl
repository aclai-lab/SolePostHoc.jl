module Lumen

using SoleLogics
const SL = SoleLogics
using SoleModels
const SM = SoleModels
using SoleData
const SD = SoleData

export lumen, LumenConfig, LumenResult

# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #
abstract type AbstractConfig end

# ---------------------------------------------------------------------------- #
#                   initial minimization algorithms setup                      #
# ---------------------------------------------------------------------------- #
function setup_espresso() end
function setup_boom() end

function setup_abc() 
    # auto setup ABC binary if not specified
    abcbinary = try
        joinpath(SD.load(SD.ABCLoader()), "abc")
    catch e
        error("Failed to setup ABC binary: $e")
    end

    # verify that binary exists and is executable
    isfile(abcbinary) || error("ABC binary not found at $abcbinary")

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
_featurename(f::SD.VariableValue) = isnothing(f.i_name) ?
    "V$(f.i_variable)" :
    "[$(f.i_name)]"

# ---------------------------------------------------------------------------- #
#                                 Lumen struct                                 #
# ---------------------------------------------------------------------------- #
struct LumenConfig <: AbstractConfig
    minimization_scheme::Symbol
    binary::String
    depth::Float64
    vertical::Float64
    horizontal::Float64
    minimization_kwargs::NamedTuple
    filt_alphabet::Base.Callable
    apply_function::Base.Callable
    importance::Vector
    check_opt::Bool
    check_alphabet::Bool

    function LumenConfig(;
        minimization_scheme::Symbol=:abc,
        depth::Float64=1.0,
        vertical::Float64=1.0,
        horizontal::Float64=1.0,
        minimization_kwargs::NamedTuple=(;),
        filt_alphabet::Base.Callable=identity,
        apply_function::Base.Callable=identity,
        importance::Vector=Float64[],
        check_opt::Bool=false,
        check_alphabet::Bool=false,
    )
        # validate coverage parameters - must be positive and ≤ 1.0
        # these parameters control the proportion of instances that must be covered by rules
        if vertical ≤ 0.0 || vertical > 1.0 ||
            horizontal ≤ 0.0 || horizontal > 1.0 ||
            depth ≤ 0.0 || depth > 1.0
            throw(ArgumentError(
                "vertical, depth and horizontal parameters must be in range (0.0, 1.0]. " *
                "Got vertical=$(vertical), depth=$(depth), horizontal=$(horizontal). " *
                "These parameters control rule coverage and must be meaningful proportions.",
            ),)
        end

        # validate minimization scheme
        valid_schemes = Dict(
            :mitespresso => setup_espresso(),
            :boom => setup_boom(),
            :abc => setup_abc(),
            :abc_balanced => setup_abc(),
            :abc_thorough => setup_abc(),
            :quine => setup_quine(),
            :quine_naive => setup_quine()
        )
        
        if minimization_scheme ∉ keys(valid_schemes)
            throw(ArgumentError(
                "minimization_scheme must be one of: $(keys(valid_schemes) |> collect). " *
                "Got: $(minimization_scheme). "
            ))
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
            check_alphabet
        )
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
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
@inline get_minimization_scheme(r::LumenConfig) = r.minimization_scheme
@inline get_binary(r::LumenConfig)              = r.binary
@inline get_depth(r::LumenConfig)               = r.depth
@inline get_vertical(r::LumenConfig)            = r.vertical
@inline get_horizontal(r::LumenConfig)          = r.horizontal
@inline get_minimization_kwargs(r::LumenConfig) = r.minimization_kwargs
@inline get_filt_alphabet(r::LumenConfig)       = r.filt_alphabet
@inline get_apply_function(r::LumenConfig)      = r.apply_function
@inline get_importance(r::LumenConfig)          = r.importance
@inline get_check_opt(r::LumenConfig)           = r.check_opt
@inline get_check_alphabet(r::LumenConfig)      = r.check_alphabet
@inline get_rng(r::LumenConfig)                 = r.rng

# ---------------------------------------------------------------------------- #
#                      extra methods for SoleLogics SL.Atom                       #
# ---------------------------------------------------------------------------- #
@inline get_operator(atom::SL.Atom{<:SD.AbstractCondition})   = atom.value.metacond.test_operator
@inline get_feature(atom::SL.Atom{<:SD.AbstractCondition})    = atom.value.metacond.feature
@inline get_threshold(atom::SL.Atom{<:SD.AbstractCondition})  = atom.value.threshold
@inline get_i_variable(atom::SL.Atom{<:SD.AbstractCondition}) = atom.value.metacond.feature.i_variable

# ---------------------------------------------------------------------------- #
#                                 depth utils                                  #
# ---------------------------------------------------------------------------- #
function _extract_atoms_bfs_order(tree::SM.AbstractModel)
    bfs_atoms = SL.Atom{SD.AbstractCondition}[]
    queue     = SM.AbstractModel[tree]

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
    atoms :: Vector{<:SL.Atom{<:SD.ScalarCondition}},
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
    atoms :: Vector{<:SL.Atom{<:SD.ScalarCondition}},
    feat  :: Symbol
) = filter(a -> SM.featurename(get_feature(a)) == feat, atoms)

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
    disjuncts   :: Vector{SL.Atom},
    i           :: Int64,
    featurename :: Symbol,
    operator    :: Union{typeof(<), typeof(≥)},
    threshold   :: Real
)
    feature   = SD.VariableValue(i, featurename)
    mc        = SD.ScalarMetaCondition(feature, operator)
    condition = SD.ScalarCondition(mc, threshold)

    push!(disjuncts, SL.Atom(condition))
end

function generate_disjunct(
    truths     :: Vector{BitVector},
    thresholds :: Vector{Vector{Float64}},
    features   :: Vector{Symbol}
)
    disjuncts = Vector{SL.Atom}()

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
    features   :: Vector{<:SM.Label}
    classnames :: Vector{<:SM.Label}

    ExtractRulesData(
        grp_truths :: Vector{Vector{Vector{BitVector}}},
        thresholds :: Vector{Vector{Float64}},
        features   :: Vector{<:SM.Label},
        classnames :: Vector{<:SM.Label}
    ) = new(grp_truths, thresholds, features, classnames)

    function ExtractRulesData(extractor::LumenConfig, model::SM.AbstractModel)
        depth = get_depth(extractor)

        atoms = unique!(if depth < 1.0
            mapreduce(vcat, SM.models(model); init=SL.Atom{SD.AbstractCondition}[]) do t
                all_atoms_bfs = _extract_atoms_bfs_order(t)
                _take_first_percentage(all_atoms_bfs, depth)
            end
        else
            SL.atoms(SM.alphabet(model, false))
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

        features     = SM.featurename.(unique!(get_feature.(atoms)))
        featurenames = SM.info(model, :featurenames)
        classnames   = unique!(SM.info(model, :supporting_labels))

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

function get_grouped_truths(e::ExtractRulesData, c::SM.Label)
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

function get_atoms(e::ExtractRulesData, c::SM.Label)
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
    conjuncts = Vector{Vector{SL.Atom}}(undef, length(truths))

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

function get_conjuncts(e::ExtractRulesData, c::SM.Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_conjuncts(e, i)
end

@inline  get_conjuncts(a::Vector{Vector{SL.Atom}}) = get_conjuncts.(a)
@inline  get_conjuncts(a::Vector{SL.Atom}) = isempty(a) ? ⊤ : SL.LeftmostConjunctiveForm{SL.Literal}(SL.Literal.(a))

# ---------------------------------------------------------------------------- #
#                                get formulas                                  #
# ---------------------------------------------------------------------------- #
@inline  get_formula(e::ExtractRulesData, i::Int64) = get_formula(get_conjuncts(e, i))

function get_formula(e::ExtractRulesData, c::SM.Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_formula(e, i)
end

@inline  get_formula(grouped_conj::Vector{SL.LeftmostConjunctiveForm{SL.Atom}}) =
    SL.LeftmostDisjunctiveForm{SL.LeftmostConjunctiveForm{SL.Literal}}(grouped_conj, true)

# ---------------------------------------------------------------------------- #
#                           dnf minimization refine                            #
# ---------------------------------------------------------------------------- #
function _refine_dnf(terms::Vector{<:Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure}})  
    length(terms) ≤ 1 && return terms
    
    all_bounds = map(term -> SD.extract_term_bounds(term; silent=true), terms)
    
    # find terms not strictly dominated by any other term
    keep_mask = map(enumerate(all_bounds)) do (i, bounds_i)
        !any(j -> i ≠ j && SD.strictly_dominates(all_bounds[j], bounds_i), eachindex(all_bounds))
    end
    
    kept_terms = terms[keep_mask]
    
    # safety check: never return empty formula
    return isempty(kept_terms) ? terms : kept_terms
end

# ---------------------------------------------------------------------------- #
#                              minimization core                               #
# ---------------------------------------------------------------------------- #
function run_minimization(extractor::LumenConfig, atoms::Vector{Vector{SL.Atom}})
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
    extractor :: LumenConfig,
    model     :: SM.AbstractModel
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

    # formulas         = Vector{SL.LeftmostLinearForm}(undef, nclasses)
    formulas         = Vector{Vector{Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure}}}(undef, nclasses)

    Threads.@threads for i in 1:nclasses
        atoms        = get_atoms(extractrulesdata, i)
        formulas[i]  = run_minimization(extractor, atoms)
    end

    return SL.LeftmostDisjunctiveForm.(formulas)
end

function lumen(
    extractor :: LumenConfig,
    model     :: Vector{SM.AbstractModel}
)
    ds = map(enumerate(model)) do (_, m)
        lumen(extractor, m)
    end

    return LumenResult(ds)
end

end