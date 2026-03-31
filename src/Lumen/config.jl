# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #
abstract type AbstractConfig end

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
        apply_function::Base.Callable=SM.apply,
        importance::Vector=Float64[],
        check_opt::Bool=false,
        check_alphabet::Bool=false,
    )
        # validate coverage parameters - must be positive and ≤ 1.0
        # these parameters control the proportion of instances
        # that must be covered by rules
        if vertical ≤ 0.0 || vertical > 1.0 ||
            horizontal ≤ 0.0 || horizontal > 1.0 ||
            depth ≤ 0.0 || depth > 1.0
            throw(ArgumentError(
                "vertical, depth and horizontal parameters must be in range " *
                "(0.0, 1.0]. Got vertical=$(vertical), depth=$(depth), " *
                "horizontal=$(horizontal). These parameters control " *
                "rule coverage and must be meaningful proportions.",
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
                "minimization_scheme must be one of: " *
                "$(keys(valid_schemes) |> collect). " *
                "Got: $(minimization_scheme)."
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
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
@inline get_minimization_scheme(r::LumenConfig) = r.minimization_scheme
@inline get_binary(r::LumenConfig) = r.binary
@inline get_depth(r::LumenConfig) = r.depth
@inline get_vertical(r::LumenConfig) = r.vertical
@inline get_horizontal(r::LumenConfig) = r.horizontal
@inline get_minimization_kwargs(r::LumenConfig) = r.minimization_kwargs
@inline get_filt_alphabet(r::LumenConfig) = r.filt_alphabet
@inline get_apply_function(r::LumenConfig) = r.apply_function
@inline get_importance(r::LumenConfig) = r.importance
@inline get_check_opt(r::LumenConfig) = r.check_opt
@inline get_check_alphabet(r::LumenConfig) = r.check_alphabet