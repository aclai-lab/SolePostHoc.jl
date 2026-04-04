# ---------------------------------------------------------------------------- #
#                                    types                                     #
# ---------------------------------------------------------------------------- #

"""
    AbstractConfig

Abstract base type for all LUMEN configuration structs.

Concrete subtypes encapsulate the parameters needed to control a specific
algorithm variant. Using a common supertype allows generic code to accept
any configuration object without being tied to a particular implementation.

See also: [`LumenConfig`](@ref)
"""
abstract type AbstractConfig end

# ---------------------------------------------------------------------------- #
#                                 Lumen struct                                 #
# ---------------------------------------------------------------------------- #

"""
    LumenConfig <: AbstractConfig

Configuration object for the LUMEN rule-extraction algorithm.

Bundles every tunable parameter into a single, validated, immutable struct.
All fields are set through the keyword constructor, which performs range
validation and resolves the correct minimizer binary before storing anything.

# Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `minimization_scheme` | `Symbol` | `:abc` | DNF minimization algorithm to use. |
| `binary` | `String` | *(auto)* | Absolute path to the minimizer executable, resolved automatically from `minimization_scheme`. |
| `depth` | `Float64` | `1.0` | Fraction of each tree's BFS-ordered atoms to include ∈ (0, 1]. `1.0` uses the full alphabet. |
| `vertical` | `Float64` | `1.0` | Instance-coverage parameter α ∈ (0, 1]. |
| `horizontal` | `Float64` | `1.0` | Feature-coverage parameter β ∈ (0, 1]. |
| `minimization_kwargs` | `NamedTuple` | `(;)` | Extra keyword arguments forwarded verbatim to the chosen minimizer. |
| `filt_alphabet` | `Base.Callable` | `identity` | Optional callback applied to the logical alphabet before rule extraction. |
| `apply_function` | `Base.Callable` | `SM.apply` | Function used to evaluate the model on generated input combinations. |
| `importance` | `Vector` | `Float64[]` | Feature-importance weights; influences rule construction when non-empty. |
| `check_opt` | `Bool` | `false` | When `true`, validates the OTT optimisation against the standard algorithm. |
| `check_alphabet` | `Bool` | `false` | When `true`, runs alphabet-analysis diagnostics instead of full extraction. |

# Supported minimization schemes

| Scheme | Backend | Notes |
|--------|---------|-------|
| `:mitespresso` | MIT Espresso | Balanced speed / quality. |
| `:boom` | BOOM | Aggressive minimisation. |
| `:abc` | Berkeley ABC | Fast, moderate compression. |
| `:abc_balanced` | Berkeley ABC | Balanced ABC variant. |
| `:abc_thorough` | Berkeley ABC | Thorough ABC variant. |
| `:quine` | Quine–McCluskey | Exact minimisation. |
| `:quine_naive` | Quine–McCluskey | Naïve variant, educational use. |

# Validation

The constructor throws `ArgumentError` when:
- Any of `vertical`, `depth`, or `horizontal` is outside (0.0, 1.0].
- `minimization_scheme` is not one of the supported symbols listed above.

# Examples

```julia
# Default configuration
cfg = LumenConfig()

# Custom scheme and coverage parameters
cfg = LumenConfig(
    minimization_scheme = :mitespresso,
    depth               = 0.7,
    vertical            = 0.9,
    horizontal          = 0.8,
)

# Pass extra kwargs to the minimizer and use a custom alphabet filter
cfg = LumenConfig(
    minimization_scheme  = :abc,
    minimization_kwargs  = (timeout = 30,),
    filt_alphabet        = alph -> my_filter(alph),
)
```

See also: [`lumen`](@ref), [`LumenResult`](@ref), [`AbstractConfig`](@ref)
"""
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
    float_type::Type

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
        float_type::Type=Float64
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
            check_alphabet,
            float_type
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #

"""
    get_minimization_scheme(r::LumenConfig) -> Symbol

Return the DNF minimization algorithm identifier stored in `r`.
"""
@inline get_minimization_scheme(r::LumenConfig) = r.minimization_scheme

"""
    get_binary(r::LumenConfig) -> String

Return the absolute path to the minimizer executable stored in `r`.
"""
@inline get_binary(r::LumenConfig) = r.binary

"""
    get_depth(r::LumenConfig) -> Float64

Return the depth coverage parameter δ ∈ (0, 1] stored in `r`.
"""
@inline get_depth(r::LumenConfig) = r.depth

"""
    get_vertical(r::LumenConfig) -> Float64

Return the instance-coverage parameter α ∈ (0, 1] stored in `r`.
"""
@inline get_vertical(r::LumenConfig) = r.vertical

"""
    get_horizontal(r::LumenConfig) -> Float64

Return the feature-coverage parameter β ∈ (0, 1] stored in `r`.
"""
@inline get_horizontal(r::LumenConfig) = r.horizontal

"""
    get_minimization_kwargs(r::LumenConfig) -> NamedTuple

Return the extra keyword arguments forwarded to the minimizer stored in `r`.
"""
@inline get_minimization_kwargs(r::LumenConfig) = r.minimization_kwargs

"""
    get_filt_alphabet(r::LumenConfig) -> Base.Callable

Return the alphabet-filter callback stored in `r`.
"""
@inline get_filt_alphabet(r::LumenConfig) = r.filt_alphabet

"""
    get_apply_function(r::LumenConfig) -> Base.Callable

Return the model-application function stored in `r`.
"""
@inline get_apply_function(r::LumenConfig) = r.apply_function

"""
    get_importance(r::LumenConfig) -> Vector

Return the feature-importance weight vector stored in `r`.
"""
@inline get_importance(r::LumenConfig) = r.importance

"""
    get_check_opt(r::LumenConfig) -> Bool

Return `true` if OTT-optimisation validation is enabled in `r`.
"""
@inline get_check_opt(r::LumenConfig) = r.check_opt

"""
    get_check_alphabet(r::LumenConfig) -> Bool

Return `true` if alphabet-analysis diagnostics are enabled in `r`.
"""
@inline get_check_alphabet(r::LumenConfig) = r.check_alphabet

"""
    get_float_type(r::LumenConfig) -> Type

Return the floating-point type stored in `r`.
"""
@inline get_float_type(r::LumenConfig) = r.float_type
