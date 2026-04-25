module Lumen

using SoleLogics
const SL = SoleLogics
using SoleModels
const SM = SoleModels
using SoleData
const SD = SoleData

using CategoricalArrays
using DataFrames

include("config.jl")

export lumen, LumenConfig, LumenResult

const Operators = Union{typeof(<),typeof(>),typeof(≤),typeof(≥)}
const Float = Union{Float32,Float64}

# ---------------------------------------------------------------------------- #
#                   initial minimization algorithms setup                      #
# ---------------------------------------------------------------------------- #
"""
    setup_espresso() -> String

Automatically locate and validate the Espresso logic minimizer binary.

Attempts to load the MIT Espresso binary via `SoleData.MITESPRESSOLoader`.
If the binary cannot be found or loaded, an informative error is raised.

# Returns
- `String`: Absolute path to the verified Espresso executable.

# Throws
- `ErrorException`: If the loader fails or the binary is not found
  at the expected path.

# Notes
This function is called internally by the LUMEN pipeline when `:mitespresso` is
selected as the minimization scheme.

See also: [`setup_abc`](@ref), [`lumen`](@ref)
"""
function setup_espresso()
    # auto setup espresso binary if not specified
    espressobinary = try
        joinpath(SD.load(SD.MITESPRESSOLoader()), "espresso")
    catch e
        error("Failed to setup espresso binary: $e")
    end

    # verify that binary exists and is executable
    isfile(espressobinary) ||
        error("espresso binary not found at $espressobinary")

    return espressobinary
end

"""
    setup_boom() -> Nothing

Placeholder for the BOOM minimizer setup routine.

This function is reserved for future integration of the BOOM logic minimization
tool. Currently a no-op pending evaluation of the minimizer.

# Notes
- Not yet implemented.
- TODO: evaluate and implement this minimizer.
"""
function setup_boom() end # TODO: evaluate this minimizer

"""
    setup_abc() -> String

Automatically locate, validate, and smoke-test the ABC logic synthesis binary.

Attempts to load the ABC binary via `SoleData.ABCLoader`. After locating the
binary it performs a basic health-check by running `abc -h` to ensure the
executable is functional.

# Returns
- `String`: Absolute path to the verified ABC executable.

# Throws
- `ErrorException`: If the loader fails, the binary is missing, 
  or the health-check invocation raises an exception.

# Notes
This function is called internally when `:abc` (or its variants) is selected as
the minimization scheme.

See also: [`setup_espresso`](@ref), [`lumen`](@ref)
"""
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

"""
    setup_quine() -> Nothing

Placeholder for the Quine–McCluskey minimizer setup routine.

Reserved for future integration of the Quine–McCluskey algorithm. Currently a
no-op pending evaluation.

# Notes
- Not yet implemented.
- TODO: evaluate and implement this minimizer.
"""
function setup_quine() end # TODO: evaluate this minimizer

# ---------------------------------------------------------------------------- #
#                                 print utils                                  #
# ---------------------------------------------------------------------------- #
"""
    _featurename(f::SD.VariableValue) -> String

Return a human-readable name for a `VariableValue` feature.

If the feature carries an explicit name (`i_name`), it is returned wrapped in
square brackets. Otherwise the fallback `"V<index>"` string is produced using
the feature's integer index (`i_variable`).

# Arguments
- `f::SD.VariableValue`: The feature descriptor to format.

# Returns
- `String`: Either `"[<name>]"` or `"V<index>"`.
"""
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
"""
    LumenResult

Lightweight container for the output produced by [`lumen`](@ref).

# Fields
- `decision_set::DecisionSet`: The minimized rule set extracted from the model.
- `info::NamedTuple`: Auxiliary metadata. Empty `(;)` when not requested.

# Constructors

```julia
LumenResult(decision_set, info) # Full construction with metadata.
LumenResult(decision_set)       # Convenience constructor; info defaults to (;).
```

# Examples
```julia
result = lumen(model)
rules  = result.decision_set
meta   = result.info            # NamedTuple – may be empty
```

See also: [`lumen`](@ref), [`LumenConfig`](@ref)
"""
struct LumenResult
    decision_set::DecisionSet
    info::NamedTuple

    LumenResult(ds, info) = new(ds, info)
    LumenResult(ds) = new(ds, (;))
end

"""
    Base.length(lr::LumenResult) -> Int

Return the number of rules contained in the result's `decision_set`.
"""
Base.length(lr::LumenResult) = length(lr.decision_set)

# ---------------------------------------------------------------------------- #
#                      extra methods for SoleLogics Atom                       #
# ---------------------------------------------------------------------------- #
"""
    get_operator(atom::SL.Atom{<:SD.AbstractCondition}) -> Function

Return the test operator (e.g. `<`, `≥`) embedded in a scalar condition atom.
"""
@inline get_operator(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.test_operator

"""
    get_feature(atom::SL.Atom{<:SD.AbstractCondition}) -> SD.AbstractFeature

Return the feature object embedded in a scalar condition atom.
"""
@inline get_feature(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.feature

"""
    get_threshold(atom::SL.Atom{<:SD.AbstractCondition}) -> Real

Return the threshold value stored in a scalar condition atom.
"""
@inline get_threshold(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.threshold

"""
    get_i_variable(atom::SL.Atom{<:SD.AbstractCondition}) -> Int

Return the integer variable index of the feature inside a scalar condition atom.
"""
@inline get_i_variable(atom::SL.Atom{<:SD.AbstractCondition}) =
    atom.value.metacond.feature.i_variable

# ---------------------------------------------------------------------------- #
#                            operator family utils                             #
# ---------------------------------------------------------------------------- #
"""
    _supported_operators

The set of scalar comparison operators that LUMEN currently supports.

Supported:
- `<`  and its negation `≥`  (strictly-less family)
- `>`  and its negation `≤`  (strictly-greater family)

All four operators are accepted when extracting atoms from a model; after
normalization via [`_normalize_atom`](@ref) only `<` and `≥` will survive in
practice. The full set is kept here as a safety-net to catch genuinely
unsupported operators (e.g. `==`, `!=`) that normalization does not handle.
"""
const _supported_operators = ((<), (≥), (>), (≤))

"""
    _is_lt_family(op) -> Bool

Return `true` when `op` belongs to the strictly-less family (`<` or `≤`).

The `<`/`≥` family encodes conditions as `value < threshold` and sorts
thresholds in **descending** order (largest first) so that the Gray-code
bit pattern in [`_truths_by_thresholds`](@ref) is consistent.

See also: [`_is_gt_family`](@ref)
"""
@inline _is_lt_family(op) = op === (<) || op === (≤)

"""
    _is_gt_family(op) -> Bool

Return `true` when `op` belongs to the strictly-greater family (`>` or `≥`).

The `>`/`≤` family encodes conditions as `value > threshold` and sorts
thresholds in **ascending** order (smallest first) so that the Gray-code
bit pattern in [`_truths_by_thresholds`](@ref) is consistent with the
reversed ordering direction.

See also: [`_is_lt_family`](@ref)
"""
@inline _is_gt_family(op) = op === (>) || op === (≥)

"""
    _feature_op_family(atoms, feat) -> Symbol

Determine the operator family used by the atoms belonging to `feat`.

Inspects all atoms whose feature name matches `feat` and returns:
- `:lt` if every operator in that group belongs to the `<`/`≤` family.
- `:gt` if every operator in that group belongs to the `>`/`≥` family.

# Throws
- `ArgumentError`: If the feature's atoms mix both families, which would make
  the threshold encoding ambiguous.
"""
function _feature_op_family(
    atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
    feat::Symbol
)
    feat_atoms = _atoms_for_feature(atoms, feat)
    ops = unique(get_operator.(feat_atoms))

    has_lt = any(_is_lt_family, ops)
    has_gt = any(_is_gt_family, ops)

    (has_lt && has_gt) && throw(ArgumentError(
        "Feature '$feat' mixes '<'/'≤' and '>'/'≥' operators. " *
        "Each feature must use operators from a single comparison family."
    ))

    return has_lt ? :lt : :gt
end

# ---------------------------------------------------------------------------- #
#                            atom normalization                                #
# ---------------------------------------------------------------------------- #
"""
    _normalize_atom(atom::SL.Atom{<:SD.ScalarCondition})
        -> SL.Atom{<:SD.ScalarCondition}

Rewrite atoms that use `>` or `≤` into the canonical `<`/`≥` family so that
every feature ends up using a single comparison direction.

The rewrites are lossless under IEEE 754 floating-point arithmetic:
- `value > t`  →  `value ≥ nextfloat(t)`
- `value ≤ t`  →  `value < nextfloat(t)`
- `value < t`  →  unchanged
- `value ≥ t`  →  unchanged

This resolves the mixed-family error that arises when a `DecisionList` contains
rules whose antecedents use both operator families on the same feature (e.g.
`sepal_length ≤ 5.6` in one rule and `sepal_length ≥ 5.9` in another).

# Arguments
- `atom::SL.Atom{<:SD.ScalarCondition}`: The atom to normalize.

# Returns
- `SL.Atom{<:SD.ScalarCondition}`: An equivalent atom whose operator belongs to
  the canonical `<`/`≥` family.

# Notes
`nextfloat(t)` is used rather than an arbitrary epsilon because `x > t` is
satisfied by exactly the IEEE 754 doubles that are also `≥ nextfloat(t)` —
no information is lost and no ad-hoc constant is introduced.

See also: [`_feature_op_family`](@ref), [`ExtractRulesData`](@ref)
"""
function _normalize_atom(
    atom::SL.Atom{<:SD.ScalarCondition}
)::SL.Atom{SD.ScalarCondition}
    op = get_operator(atom)

    # fast path: already in the canonical family
    op in ((<), (≥)) && return atom

    thr  = get_threshold(atom)
    feat = get_feature(atom)

    new_op, new_thr = if op === (>)
        (≥), nextfloat(thr)   # > t  →  ≥ nextfloat(t)
    else                       # op === (≤)
        (<), nextfloat(thr)   # ≤ t  →  <  nextfloat(t)
    end

    i_name = isnothing(feat.i_name) ? feat.i_variable : feat.i_name

    new_mc = SD.ScalarMetaCondition(
        SD.VariableValue(feat.i_variable, i_name), new_op
    )

    SL.Atom(SD.ScalarCondition(new_mc, new_thr))
end

# ---------------------------------------------------------------------------- #
#                                 depth utils                                  #
# ---------------------------------------------------------------------------- #
"""
    _extract_atoms_bfs_order(tree::SM.AbstractModel)
        -> Vector{SL.Atom{SD.AbstractCondition}}

Traverse a decision-tree model in breadth-first order and return the antecedent
atoms encountered at each `Branch` node.

The traversal visits left (positive) and right (negative) sub-trees in BFS order.
Only `SM.Branch` nodes contribute atoms; leaf nodes are silently skipped.

# Arguments
- `tree::SM.AbstractModel`: Root of the decision tree (or sub-tree) to traverse.

# Returns
- `Vector{SL.Atom{SD.AbstractCondition}}`: Atoms in BFS visitation order.
"""
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

"""
    _take_first_percentage(
        atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
        depth::Float64
    ) -> Vector{<:SL.Atom{<:SD.ScalarCondition}}

Return the first `ceil(length(atoms) × depth)` elements of `atoms`.

Used to implement partial-depth extraction: only atoms from the upper portion
of a decision tree (as visited in BFS order) are retained.

# Arguments
- `atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}}`: Ordered list of atoms.
- `depth::Float64`: Fraction ∈ (0, 1] of atoms to keep.

# Returns
- A sub-vector containing at most `ceil(n × depth)` elements.
"""
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
"""
    _atoms_for_feature(
        atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
        feat::Symbol
    ) -> Vector{<:SL.Atom{<:SD.ScalarCondition}}

Filter `atoms` to only those whose feature name matches `feat`.

# Arguments
- `atoms`: Collection of scalar-condition atoms.
- `feat::Symbol`: Target feature name (as returned by `SM.featurename`).

# Returns
- Sub-vector of atoms whose feature matches `feat`.
"""
@inline _atoms_for_feature(
    atoms::Vector{<:SL.Atom{<:SD.ScalarCondition}},
    feat::Symbol
) = filter(a -> SM.featurename(get_feature(a)) == feat, atoms)

"""
    _truths_by_thresholds(thresholds::Vector{Float64}) -> Vector{BitVector}

Build a lookup table mapping each "threshold region" to the truth-value
assignment it implies for the `length(thresholds)` binary conditions.

For `n` thresholds there are `n + 1` possible ordinal regions. The `i`-th entry
encodes, for each condition `j`, whether `value < thresholds[j]` is true in that
region using a Gray-code–inspired bit pattern.

When thresholds are sorted **descending** (the `<`/`≥` family default) the bit
pattern is consistent with the `<` semantics.  When thresholds are sorted
**ascending** (the `>`/`≤` family) the same bit pattern is interpreted as
`value > thresholds[j]` by the callers that supply ascending-sorted vectors.

# Arguments
- `thresholds::Vector{Float64}`: Sorted threshold values.

# Returns
- `Vector{BitVector}` of length `n + 1`, one entry per ordinal region.

---

    function _truths_by_thresholds(thresholds::Vector{<:Float})
        -> Vector{Vector{BitVector}}

Broadcast version: applies `_truths_by_thresholds` element-wise to a vector of
threshold vectors (one per feature).

---

    _truths_by_thresholds(
        thresholds::Vector{T}
    ) where {T<:Vector{<:Float}} -> BitVector

Return the truth-value assignment for a single concrete `value`
against `thresholds`.

Returns an empty `BitVector` when `value` is `NaN`, and `falses(n)` when `value`
is not found in `thresholds`.

---

    _truths_by_thresholds(value::Float, thresholds::Vector{<:Float})
        -> Vector{BitVector}

Broadcast version: element-wise application for a tuple of values paired with a
vector of per-feature threshold vectors.
"""
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

"""
    _thrs_with_boundary(
        thresholds::Vector{T},
        family::Symbol
    ) where {T<:Float} -> Vector{Float64}

Append the appropriate boundary point to `thresholds` depending on the operator
family, ensuring that all `n + 1` ordinal regions induced by `n` thresholds
are sampled.

- `:lt` family (`<`/`≤`): thresholds sorted **descending** → appends
  `prevfloat(last(thresholds))` 
    to cover the region **below the smallest threshold**
  (i.e. values smaller than every condition).

- `:gt` family (`>`/`≥`): thresholds sorted **ascending** → appends
  `nextfloat(last(thresholds))` 
    to cover the region **above the largest threshold**
  (i.e. values larger than every condition).

Returns `[NaN]` for an empty input vector.

# Examples

```julia
# :lt  — thresholds [4.8, 4.7, 1.9] (descending)
# regions: x < 1.9 | 1.9 ≤ x < 4.7 | 4.7 ≤ x < 4.8 | x ≥ 4.8
# boundary needed: prevfloat(1.9)  ← covers  x < 1.9
_thrs_with_boundary([4.8, 4.7, 1.9], :lt)
# → [4.8, 4.7, 1.9, prevfloat(1.9)]

# :gt  — thresholds [1.9, 4.7, 4.8] (ascending)
# regions: x ≤ 1.9 | 1.9 < x ≤ 4.7 | 4.7 < x ≤ 4.8 | x > 4.8
# boundary needed: nextfloat(4.8)  ← covers  x > 4.8
_thrs_with_boundary([1.9, 4.7, 4.8], :gt)
# → [1.9, 4.7, 4.8, nextfloat(4.8)]
```

---

    _thrs_with_boundary(
        thresholds::Vector{T},
        op_families::Vector{Symbol}
    ) where {T<:Vector{<:Float}} -> Vector{Vector{T}}

Element-wise version: applies `_thrs_with_boundary`
to each per-feature threshold
vector using the corresponding operator family.
"""
function _thrs_with_boundary(
    thresholds::Vector{T},
    family::Symbol
) where {T<:Float}
    isempty(thresholds) && return [NaN]

    nthrs  = length(thresholds)
    result = Vector{T}(undef, nthrs + 1)
    result[1:nthrs] .= thresholds

    # :lt (descending) → boundary point is BELOW the minimum threshold
    #                    prevfloat(last) because last is the smallest value
    # :gt (ascending)  → boundary point is ABOVE the maximum threshold
    #                    nextfloat(last) because last is the largest value
    result[end] = family === :lt ?
        prevfloat(last(thresholds)) :
        nextfloat(last(thresholds))

    return result
end

@inline _thrs_with_boundary(
    thresholds::Vector{T},
    op_families::Vector{Symbol}
) where {T<:Vector{<:Float}} = _thrs_with_boundary.(thresholds, op_families)

# ---------------------------------------------------------------------------- #
#                              generate disjunts                               #
# ---------------------------------------------------------------------------- #
"""
    push_disjunct!(
        disjuncts::Vector{SL.Atom},
        i::Int,
        featurename::Symbol,
        operator,
        threshold::Real
    ) -> Nothing

Construct a new `SL.Atom` encoding the scalar condition
`feature[i] <operator> threshold` and append it to `disjuncts` in-place.

# Arguments
- `disjuncts`: Accumulator vector to append the new atom to.
- `i::Int`: Integer index identifying the feature (used for `SD.VariableValue`).
- `featurename::Symbol`: Human-readable feature name.
- `operator`: Any of `<`, `≥`, `>`, `≤`.
- `threshold::Real`: The numeric comparison threshold.
"""
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

"""
    generate_disjunct(
        truths::Vector{BitVector},
        thresholds::Vector{T},
        features::Vector{<:Feature},
        op_families::Vector{Symbol}
    ) where {T<:Vector{<:Float}} -> Vector{SL.Atom}

Derive the tightest bounding atoms for a single truth-value assignment.

For each feature `i` the function inspects `truths[i]` and the feature's
operator family (`op_families[i]`) to determine the bounding conditions:

- `:lt` family (`<`/`≥`): thresholds are sorted **descending**.
  - `idx0` (false bits) → emit `value < max_threshold_where_false`
  - `idx1` (true bits)  → emit `value ≥ min_threshold_where_true`

- `:gt` family (`>`/`≤`): thresholds are sorted **ascending**.
  - `idx0` (false bits) → emit `value ≤ min_threshold_where_false`
  - `idx1` (true bits)  → emit `value > max_threshold_where_true`

# Arguments
- `truths::Vector{BitVector}`: Per-feature truth assignments over the
  sorted thresholds.
- `thresholds::Vector{Vector{Float64}}`: Per-feature sorted threshold vectors.
- `features::Vector{Symbol}`: Feature names aligned with `thresholds`.
- `op_families::Vector{Symbol}`: Per-feature operator family (`:lt` or `:gt`).

# Returns
- `Vector{SL.Atom}`: Atoms encoding the implied scalar conditions.
"""
function generate_disjunct(
    truths::Vector{BitVector},
    thresholds::Vector{T},
    features::Vector{Symbol},
    op_families::Vector{Symbol}
) where {T<:Vector{<:Float}}
    disjuncts = Vector{SL.Atom}()

    @inbounds for i in eachindex(thresholds)
        idx0 = findall(x -> !x,  truths[i])
        idx1 = findall(identity, truths[i])

        if op_families[i] === :lt
            # thresholds sorted descending
            #   → same logic as the original pipeline:
            # false bits (idx0) → upper bound via 
            # true  bits (idx1) → lower bound via ≥
            isempty(idx0) ||
                push_disjunct!(
                    disjuncts, i, features[i], <, thresholds[i][maximum(idx0)]
                )
            isempty(idx1) ||
                push_disjunct!(
                    disjuncts, i, features[i], ≥, thresholds[i][minimum(idx1)]
                )
        else  # :gt family — thresholds sorted ascending
            # false bits (idx0) → upper bound via ≤
            #   (value is NOT > any of these)
            # true  bits (idx1) → lower bound via >
            #   (value IS > all of these)
            isempty(idx0) ||
                push_disjunct!(
                    disjuncts, i, features[i], ≤, thresholds[i][minimum(idx0)]
                )
            isempty(idx1) ||
                push_disjunct!(
                    disjuncts, i, features[i], >, thresholds[i][maximum(idx1)]
                )
        end
    end

    return disjuncts
end

# ---------------------------------------------------------------------------- #
#                          extract rules data struct                           #
# ---------------------------------------------------------------------------- #
"""
    ExtractRulesData

Intermediate data structure that aggregates all information needed to build and
minimize per-class DNF formulas.

# Fields
- `grp_truths::Vector{Vector{Vector{BitVector}}}`: For each class, a list of
  truth-value assignments (one per input combination that the model assigns to
  that class). Each assignment is a `Vector{BitVector}` – one `BitVector` per
  feature.
- `thresholds::Vector{Vector{Float64}}`: Per-feature sorted threshold vectors
  derived from the model's alphabet.
- `features::Vector{<:SM.Label}`: Ordered feature names aligned with
  `thresholds`.
- `classnames::Vector{<:SM.Label}`: Unique class labels in the model.
- `op_families::Vector{Symbol}`: Per-feature operator family (`:lt` or `:gt`),
  used by [`generate_disjunct`](@ref) to emit the correct comparison operators.

# Constructors

```julia
# Low-level constructor: supply all fields directly.
ExtractRulesData(grp_truths, thresholds, features, classnames, op_families)

# High-level constructor: derive everything from a LumenConfig and a model.
ExtractRulesData(extractor::LumenConfig, model::SM.AbstractModel)
```

The high-level constructor:
1. Extracts atoms from the model (respecting the `depth` parameter).
2. Normalizes all atoms to the canonical `<`/`≥` family via
   [`_normalize_atom`](@ref).
3. Builds per-feature sorted threshold vectors.
4. Enumerates all threshold-induced input combinations.
5. Applies the model to those combinations to obtain class labels.
6. Groups truth assignments by predicted class.

See also: [`lumen`](@ref), [`LumenConfig`](@ref), [`get_atoms`](@ref)
"""
struct ExtractRulesData{T<:Vector{<:Float},F<:SM.Label,L<:SM.Label}
    # grp_truths::Vector{Vector{Vector{BitVector}}}
    predictions::Vector{L}
    combinations::Array{<:NTuple}
    thresholds::Vector{T}
    featurenames::Vector{F}
    classnames::AbstractVector{L}
    op_families::Vector{Symbol}

    ExtractRulesData(
        # grp_truths::Vector{Vector{Vector{BitVector}}},
        predictions::Vector{L},
        combinations::Array{<:NTuple},
        thresholds::Vector{T},
        featurenames::Vector{F},
        classnames::AbstractVector{L},
        op_families::Vector{Symbol}
    ) where {T<:Vector{<:Float},F<:SM.Label,L<:SM.Label} =
        new{T,F,L}(
            predictions, combinations, thresholds, featurenames, classnames, op_families
        )

    function ExtractRulesData(extractor::LumenConfig, model::SM.AbstractModel)
        # -------------------------------------------------------------------- #
        # STEP 1 — Read the depth parameter from the configuration.
        # `depth ∈ (0, 1]`: if < 1.0, only atoms from the upper levels of the
        # tree are used (partial extraction); if == 1.0, all atoms are used.
        # -------------------------------------------------------------------- #
        depth = get_depth(extractor)

        # -------------------------------------------------------------------- #
        # STEP 2 — Extract the atoms (scalar conditions) from the model,
        # normalize them to the canonical </>= family, and deduplicate.
        #
        # Two strategies depending on `depth`:
        #
        #   depth < 1.0 → partial extraction by depth
        #       - Iterates over every tree in the model (SM.models).
        #       - For each tree, visits nodes in BFS order
        #         (_extract_atoms_bfs_order),
        #         yielding atoms ordered from root to leaves.
        #       - Retains only the first `depth`% of BFS atoms for that tree
        #         (_take_first_percentage),
        #         simulating a cut at a relative depth.
        #       - Concatenates all atoms collected across trees
        #         (mapreduce + vcat).
        #
        #   depth == 1.0 → full extraction
        #       - Directly retrieves the alphabet of the entire model
        #         (SM.alphabet(model, false)) and extracts all its atoms.
        #
        # In both cases `_normalize_atom` is broadcast over the raw atom list to
        # rewrite any `>` or `≤` operator into the canonical `<`/`≥` family
        # (see [`_normalize_atom`](@ref)), and `unique!` removes duplicates
        # in-place. Normalization ensures that `_feature_op_family` never
        # encounters a mixed-family feature, which would otherwise arise when a
        # DecisionList mixes operator families across its rules.
        # -------------------------------------------------------------------- #
        atoms = unique!(_normalize_atom.(if depth < 1.0
            mapreduce(
                vcat, SM.models(model); init=SL.Atom{SD.AbstractCondition}[]
            ) do t
                all_atoms_bfs = _extract_atoms_bfs_order(t)
                _take_first_percentage(all_atoms_bfs, depth)
            end
        else
            SL.atoms(SM.alphabet(model, false))
        end))

        # -------------------------------------------------------------------- #
        # STEP 3 — Validate that every operator present in the extracted atoms
        # belongs to the supported set: `<`, `≥`, `>`, `≤`.
        #
        # After normalization only `<` and `≥` should remain; this check acts
        # as a safety-net for genuinely unsupported operators (e.g. `==`, `!=`)
        # that `_normalize_atom` does not handle.
        # -------------------------------------------------------------------- #
        let unsupported = unique(
                op for op in get_operator.(atoms)
                if op ∉ _supported_operators
            )
            isempty(unsupported) || throw(ArgumentError(
                "Only '<', '≥', '>', '≤' operators are currently supported. " *
                "Found unsupported operators: $(unsupported). " *
                "This limitation may be addressed in future versions. " *
                "Consider preprocessing your model to use only " *
                "supported conditions.",
            ))
        end

        # -------------------------------------------------------------------- #
        # STEP 4 — Derive the names of the features present in the
        #          extracted atoms.
        #
        # - get_feature.(atoms)  → feature object for each atom.
        # - unique!              → removes duplicate feature objects in-place.
        # - SM.featurename.      → converts each feature object to its 
        #                          symbolic name.
        #
        # `features` therefore contains the names of only the features actually
        # referenced by the atoms (a subset of the model's full feature set).
        # -------------------------------------------------------------------- #
        features = SM.featurename.(unique!(get_feature.(atoms))) 
        # TODO: if we dont have featurename ?

        # -------------------------------------------------------------------- #
        # STEP 5 — Retrieve the canonical feature name list and class labels
        # from the model.
        #
        # - featurenames: the model's canonical feature ordering (may include
        #   features absent from the extracted atoms, e.g. when depth < 1.0).
        # - classnames: unique class labels present in the model's leaves.
        # -------------------------------------------------------------------- #
        featurenames = SM.info(model, :featurenames)
        classnames = unique!(SM.info(model, :supporting_labels))

        # -------------------------------------------------------------------- #
        # STEP 6 — Build, for each feature in `featurenames`, its sorted
        # threshold vector and determine its operator family.
        #
        # For each feature i:
        #   - Look up its name in `features` (the extracted atoms).
        #   - If not found (idx == nothing): the feature does not appear in the
        #     extracted atoms → assign an empty threshold vector [] and default
        #     family :lt (unused, since no thresholds means no conditions).
        #   - If found:
        #       * Determine the operator family via _feature_op_family.
        #         Because atoms have already been normalized by _normalize_atom,
        #         every feature is guaranteed to use a single family here.
        #       * Filter atoms belonging to that feature (_atoms_for_feature),
        #         extract their threshold values (get_threshold.), and sort:
        #           - descending for the :lt family (consistent with `value < t`
        #             encoding used by _truths_by_thresholds).
        #           - ascending  for the :gt family (consistent with `value > t`
        #             encoding, where larger thresholds correspond to
        #             later bits).
        # -------------------------------------------------------------------- #
        type = get_float_type(extractor)
        thresholds = Vector{Vector{type}}(undef, length(featurenames))
        op_families = Vector{Symbol}(undef, length(featurenames))

        @inbounds for i in eachindex(featurenames)
            idx = findfirst(f -> f == featurenames[i], features)
            if isnothing(idx)
                thresholds[i]  = type[]
                op_families[i] = :lt # default (irrelevant: no thresholds)
            else
                family = _feature_op_family(atoms, features[idx])
                op_families[i] = family
                thresholds[i] = sort!(
                    get_threshold.(_atoms_for_feature(atoms, features[idx]));
                    rev=(family === :lt) # descending for :lt, ascending for :gt
                )
            end
        end

        # -------------------------------------------------------------------- #
        # STEP 7 — Augment each threshold vector with the correct
        #          boundary point.
        #
        # For each feature, one extra sampling point is appended to cover the
        # ordinal region that lies beyond the extreme threshold:
        #
        #   :lt family (descending, e.g. [4.8, 4.7, 1.9]):
        #     → appends prevfloat(last) = prevfloat(1.9)
        #     → covers the region x < 1.9  (below the smallest threshold)
        #
        #   :gt family (ascending, e.g. [1.9, 4.7, 4.8]):
        #     → appends nextfloat(last) = nextfloat(4.8)
        #     → covers the region x > 4.8  (above the largest threshold)
        #
        # This ensures that all n+1 ordinal regions induced by n thresholds
        # are represented in the Cartesian product generated in STEP 8.
        # Without this extra point, the class occupying the extreme region
        # would never appear in `predictions` and no rules would be extracted
        # for it (e.g. virginica in a ≤-only iris tree).
        # -------------------------------------------------------------------- #
        thrs_with_p = _thrs_with_boundary(thresholds, op_families)

        # -------------------------------------------------------------------- #
        # STEP 8 — Generate the Cartesian product of all augmented threshold
        # vectors.
        #
        # Iterators.product(thrs_with_p...) produces every possible combination
        # of threshold values (one per feature), systematically covering all
        # input-space regions induced by the model's thresholds.
        # collect() materialises the lazy iterator into an array of tuples.
        # --
        # THIS IS THE CORE OF OUR Algorithm NOTICE, IF WE OPTIMIZE HERE WE HAVE
        # HUGE BOOST !!!
        # -------------------------------------------------------------------- #
        combinations = collect(Iterators.product(thrs_with_p...))

        # -------------------------------------------------------------------- #
        # STEP 9 — Apply the model to all generated combinations.
        #
        # - The combinations are packed into a DataFrame with the canonical
        #   feature names and converted into a scalar logiset (scalarlogiset),
        #   which is the format expected by the model.
        # - get_apply_function(extractor) returns the configured application
        #   function (e.g. SoleModels.apply or DT.apply_forest).
        # - `predictions` is a vector of class labels, one per combination.
        # -------------------------------------------------------------------- #
        predictions = get_apply_function(extractor)(
            model,
            PropositionalLogiset(DataFrame(combinations, featurenames));
            suppress_parity_warning=true
        )

        # -------------------------------------------------------------------- #
        # STEP 10 — Construct and return the instance with all computed data.
        #
        # Note: `featurenames` (canonical model ordering) is used instead of
        # `features` (atom-extraction ordering) to guarantee alignment with
        # `thresholds` and `op_families`,
        # which were both built over `featurenames`.
        # -------------------------------------------------------------------- #
        return new{Vector{<:type},eltype(featurenames),eltype(classnames)}(
            predictions, combinations, thresholds, featurenames, classnames, op_families
        )
    end
end

# ---------------------------------------------------------------------------- #
#                                   methods                                    #
# ---------------------------------------------------------------------------- #
# """
#     get_grouped_truths(e::ExtractRulesData) -> Vector{Vector{Vector{BitVector}}}

# Return the full per-class grouped truth assignments stored in `e`.

# ---

#     get_grouped_truths(e::ExtractRulesData, i::Int) -> Vector{Vector{BitVector}}

# Return the truth assignments for the `i`-th class.

# ---

#     get_grouped_truths(e::ExtractRulesData, c::SM.Label)
#         -> Union{Vector{Vector{BitVector}}, Nothing}

# Return the truth assignments for the class whose label equals `c`, or `nothing`
# if `c` is not found.
# """
# @inline get_grouped_truths(e::ExtractRulesData) = e.grp_truths

# function truths_by_thresholds(t::TableTruths, i::Int)
#     _truths_by_thresholds(t.combinations[i], t.thresholds)
# end

# function truths_by_groups(t::TableTruths, i::Int)
#     indices = findall(==(t.classnames[i]), t.predictions)
#     [truths_by_thresholds(t, i) for i in indicies]
# end

"""
    get_thresholds(
        e::ExtractRulesData;
        prev_float::Bool=false
        float_type::Type=Float64
    )
        -> Vector{Vector{float_type}}

Return the per-feature threshold vectors stored in `e`.

When `prev_float=true` the boundary-augmented form produced by
[`_thrs_with_boundary`](@ref) is returned, using each feature's operator family
to determine the correct boundary point:
- `:lt` family → `prevfloat` of the last (smallest) threshold.
- `:gt` family → `nextfloat` of the last (largest) threshold.
"""
function get_thresholds(
    e::ExtractRulesData;
    prev_float::Bool=false,
    float_type::Type=Float64
)
    thresholds = [float_type.(t) for t in e.thresholds]
    op_families = e.op_families
    return prev_float ?
        _thrs_with_boundary(thresholds, op_families) : thresholds
end

"""
    get_featurenames(e::ExtractRulesData) -> Vector{<:SM.Label}

Return the ordered feature-name vector stored in `e`.
"""
@inline get_featurenames(e::ExtractRulesData) = e.featurenames

"""
    get_classnames(e::ExtractRulesData) -> Vector{<:SM.Label}

Return the unique class-label vector stored in `e`.
"""
@inline get_classnames(e::ExtractRulesData) = e.classnames

"""
    get_op_families(e::ExtractRulesData) -> Vector{Symbol}

Return the per-feature operator family vector stored in `e`.
Each entry is either `:lt` (for `<`/`≥` models) or `:gt` (for `>`/`≤` models).
"""
@inline get_op_families(e::ExtractRulesData) = e.op_families

# function get_grouped_truths(e::ExtractRulesData, c::SM.Label)
#     i = findfirst(g -> get_classnames(g) == c, e.grp_truths)
#     isnothing(i) ? nothing : get_grouped_truths(e, i)
# end

# @inline get_grouped_truths(e::ExtractRulesData, i::Int) = e.grp_truths[i]

"""
    get_truths(e::ExtractRulesData) -> Vector{Vector{Vector{BitVector}}}

Return all per-class truth-assignment lists.

---

    get_truths(e::ExtractRulesData, i::Int) -> Vector{Vector{BitVector}}

Return the truth-assignment list for class `i`.
"""
# @inline get_truths(e::ExtractRulesData) =
#     [get_truths(e, i) for i in eachindex(get_classnames(e))]
@inline function get_truths(e::ExtractRulesData, i::Int)
    _truths_by_thresholds(e.combinations[i], e.thresholds)
end

function truths_by_groups(e::ExtractRulesData, i::Int)
    idxs = findall(==(e.classnames[i]), e.predictions)
    [get_truths(e, i) for i in idxs]
end

# """
#     get_truth(e::ExtractRulesData, i::Int, j::Int) -> Vector{BitVector}

# Return the `j`-th truth assignment for class `i`.

# ---

#     get_truth(e::ExtractRulesData, i::Int) -> Vector{Vector{BitVector}}

# Return all truth assignments for class `i` (alias for `get_truths(e, i)`).
# """
# @inline get_truth(e::ExtractRulesData, i::Int, j::Int) =
#     get_grouped_truths(e, i)[j]
# @inline get_truth(e::ExtractRulesData, i::Int) = get_truths(e)[i]

# ---------------------------------------------------------------------------- #
#                                  get atoms                                   #
# ---------------------------------------------------------------------------- #
"""
    get_atoms(e::ExtractRulesData; grouped::Bool=false)
        -> Union{Vector{Vector{Vector{SL.Atom}}}, Vector{Vector{SL.Atom}}}

Derive bounding atoms for all classes stored in `e`.

When `grouped=false` (default) returns one `Vector{Vector{SL.Atom}}` per class.
When `grouped=true` all class truth assignments are pooled before
deriving atoms,
returning a single `Vector{Vector{SL.Atom}}`.

---

    get_atoms(e::ExtractRulesData, c::SM.Label)
        -> Union{Vector{Vector{SL.Atom}}, Nothing}

Derive atoms for the class labelled `c`, or `nothing` if not found.

---

    get_atoms(e::ExtractRulesData, i::Int) -> Vector{Vector{SL.Atom}}

Derive atoms for the `i`-th class.

---

    get_atoms(
        truths::Vector{Vector{BitVector}},
        thresholds::Vector{T},
        features::Vector{<:Feature},
        op_families::Vector{Symbol}
    ) where {T<:Vector{<:Float}} -> Vector{Vector{SL.Atom}}

Low-level worker: given a flat list of truth assignments, thresholds, feature
names, and per-feature operator families, parallelise the call to
[`generate_disjunct`](@ref) over all assignments and return
the resulting atom vectors.
"""
function get_atoms(
    e::ExtractRulesData;
    grouped::Bool=false,
    float_type::Type=Float64
)
    thresholds  = get_thresholds(e; prev_float=true, float_type)
    features    = get_features(e)
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
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_atoms(e, i)
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
"""
    get_conjuncts(e::ExtractRulesData, i::Int)
        -> Vector{SL.LeftmostConjunctiveForm{SL.Literal}}

Build a conjunctive form for each input combination assigned to class `i`.

Delegates atom retrieval to [`get_atoms`](@ref) and wraps each atom vector using
the scalar `get_conjuncts(::Vector{SL.Atom})` overload.

---

    get_conjuncts(e::ExtractRulesData, c::SM.Label) -> Union{..., Nothing}

Build conjunctive forms for class `c`, or `nothing` if the class is not found.

---

    get_conjuncts(a::Vector{Vector{SL.Atom}})
        -> Vector{SL.LeftmostConjunctiveForm{SL.Literal}}

Map `get_conjuncts` over every atom vector in `a`.

---

    get_conjuncts(a::Vector{SL.Atom})
        -> Union{⊤, SL.LeftmostConjunctiveForm{SL.Literal}}

Wrap a single atom vector in a `LeftmostConjunctiveForm`.
Returns `⊤` (tautology) for an empty vector.
"""
function get_conjuncts(e::ExtractRulesData, i::Int)
    atoms = get_atoms(e, i)
    [get_conjuncts(atom) for atom in atoms]
end

function get_conjuncts(e::ExtractRulesData, c::SM.Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_conjuncts(e, i)
end

@inline get_conjuncts(a::Vector{Vector{SL.Atom}}) = get_conjuncts.(a)
@inline get_conjuncts(a::Vector{SL.Atom}) = isempty(a) ?
    ⊤ : SL.LeftmostConjunctiveForm{SL.Literal}(SL.Literal.(a))

# ---------------------------------------------------------------------------- #
#                                get formulas                                  #
# ---------------------------------------------------------------------------- #
"""
    get_formula(e::ExtractRulesData, i::Int) -> SL.LeftmostDisjunctiveForm

Build the full DNF formula for class `i` by wrapping its conjunctive forms in a
disjunction.

---

    get_formula(e::ExtractRulesData, c::SM.Label)
        -> Union{SL.LeftmostDisjunctiveForm, Nothing}

Build the DNF formula for class `c`, or `nothing` if the class is not found.

---

    get_formula(
        grouped_conj::Vector{SL.LeftmostConjunctiveForm{SL.Atom}}
    ) -> SL.LeftmostDisjunctiveForm{SL.LeftmostConjunctiveForm{SL.Literal}}

Wrap a vector of conjunctive forms into a disjunctive (DNF) formula.
"""
@inline get_formula(e::ExtractRulesData, i::Int) =
    get_formula(get_conjuncts(e, i))

function get_formula(e::ExtractRulesData, c::SM.Label)
    i = findfirst(g -> get_classname(g) == c, get_grouped_truths(e))
    isnothing(i) ? nothing : get_formula(e, i)
end

@inline get_formula(grouped_conj::Vector{SL.LeftmostConjunctiveForm{SL.Atom}}) =
    SL.LeftmostDisjunctiveForm{SL.LeftmostConjunctiveForm{SL.Literal}}(
        grouped_conj, true
    )

# ---------------------------------------------------------------------------- #
#                           dnf minimization refine                            #
# ---------------------------------------------------------------------------- #
"""
    _refine_dnf(terms::Vector{
        <:Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure
    }}) -> Vector{...}

Remove DNF terms that are strictly dominated by another term
in the same formula.

A term `t_i` is strictly dominated by `t_j` (i ≠ j) when the hyper-rectangle
described by `t_j`'s bounds is entirely contained within that of `t_i`, making
`t_i` logically redundant.

The function:
1. Extracts bounds for every term via `SD.extract_term_bounds`.
2. Marks and removes dominated terms.
3. Returns the original vector unchanged as a safety fallback if all terms would
   be removed.

# Arguments
- `terms`: Non-empty vector of conjunctive terms forming a DNF formula.

# Returns
- Pruned vector of terms; never empty
  (returns `terms` if pruning would empty it).

# Notes
Requires at least two terms to perform any pruning; single-term inputs are
returned immediately.
"""
function _refine_dnf(
    terms::Vector{<:Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure}}
)
    length(terms) ≤ 1 && return terms

    all_bounds = map(term -> SD.extract_term_bounds(term; silent=true), terms)

    # find terms not strictly dominated by any other term
    keep_mask = map(enumerate(all_bounds)) do (i, bounds_i)
        !any(j -> i ≠ j && SD.strictly_dominates(
            all_bounds[j], bounds_i), eachindex(all_bounds))
    end

    kept_terms = terms[keep_mask]

    # safety check: never return empty formula
    return isempty(kept_terms) ? terms : kept_terms
end

# ---------------------------------------------------------------------------- #
#                              minimization core                               #
# ---------------------------------------------------------------------------- #

"""
    run_minimization(
        ::Val{:abc},
        extractor::LumenConfig,
        atoms::Vector{Vector{SL.Atom}}
    ) -> Vector{<:Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure}}

Minimize the DNF formula encoded by `atoms` using the ABC framework.

Delegates to `SD.abc_minimize` with the binary path from `extractor`, then
applies [`_refine_dnf`](@ref) to remove dominated terms.

# Arguments
- `extractor::LumenConfig`: Provides the ABC binary path and depth parameter.
- `atoms::Vector{Vector{SL.Atom}}`: Per-combination atom lists (one entry per
  input combination assigned to the target class).

# Returns
- Minimized and refined vector of conjunctive terms.
"""
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

"""
    run_minimization(
        ::Val{:mitespresso},
        extractor::LumenConfig,
        atoms::Vector{Vector{SL.Atom}}
    ) -> Vector{<:Union{SL.LeftmostConjunctiveForm{SL.Atom}, SyntaxStructure}}

Minimize the DNF formula encoded by `atoms` using the MIT Espresso minimizer.

Delegates to `SD.espresso_minimize` with the binary path from `extractor`, then
applies [`_refine_dnf`](@ref) to remove dominated terms.

# Arguments
- `extractor::LumenConfig`: Provides the Espresso binary path and
  depth parameter.
- `atoms::Vector{Vector{SL.Atom}}`: Per-combination atom lists.

# Returns
- Minimized and refined vector of conjunctive terms.
"""
function run_minimization(
    ::Val{:mitespresso},
    extractor::LumenConfig,
    atoms::Vector{Vector{SL.Atom}}
    # TODO mitespresso_kwargs...
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
"""
    lumen(config::LumenConfig, model::SM.AbstractModel) -> SM.DecisionSet

Core single-model entry point for the LUMEN algorithm.

Extracts a minimized [`DecisionSet`](@ref) from `model` using the parameters
encoded in `config`.

# Pipeline
1. Build [`ExtractRulesData`](@ref) from `config` and `model` (atom extraction,
   normalization to the canonical `<`/`≥` family, truth-table enumeration,
   per-class grouping).
2. For each class, call [`run_minimization`](@ref) on the derived atom vectors.
3. Filter out classes for which no formula could be produced.
4. Wrap the minimized formulas in `SM.Rule` objects and return a `DecisionSet`.

# Arguments
- `config::LumenConfig`: Algorithm configuration
  (minimization scheme, depth, etc.).
- `model::SM.AbstractModel`: A single decision-tree model.

# Returns
- `SM.DecisionSet`: The minimized rule set.

---

    lumen(config::LumenConfig, model::Vector{SM.AbstractModel}) -> LumenResult

Batch variant: applies `lumen(config, m)` to every model in the vector and
collects the results into a [`LumenResult`](@ref).

---

    lumen(model::SM.AbstractModel, args...; kwargs...) -> SM.DecisionSet

Convenience wrapper: constructs a `LumenConfig` from keyword arguments and
delegates to `lumen(config, model)`.

---

    lumen(model::Vector{SM.AbstractModel}, args...; kwargs...) -> LumenResult

Convenience wrapper for vector of models: constructs `LumenConfig` from keyword
arguments and maps over the vector.

# Examples
```julia
# Single model with default settings
ds = lumen(my_tree)

# Single model with custom minimization scheme
ds = lumen(my_tree; minimization_scheme=:mitespresso, depth=0.8)

# Explicit config object
config = LumenConfig(minimization_scheme=:abc, depth=0.7)
ds = lumen(config, my_tree)

# Batch processing
results = lumen(config, [tree1, tree2, tree3])
```

See also: [`LumenConfig`](@ref), [`LumenResult`](@ref),
[`ExtractRulesData`](@ref)
"""
function lumen(
    config::LumenConfig,
    model::SM.AbstractModel
)
    float_type = get_float_type(config)

    # extract conjuncts
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