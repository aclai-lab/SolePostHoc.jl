@inline evalenc(::Univariate) = 0x01
@inline evalenc(::Multivariate) = 0x02


function removeduals(values::Vector)
    newvalues = similar(values, (0,))
    for cond in values
        if !SoleData.hasdual(cond) || !(SoleData.dual(cond) in newvalues)
            push!(newvalues, cond)
        end
    end
    newvalues
end

# ---------------------------------------------------------------------------- #
#                                formula to pla                                #
# ---------------------------------------------------------------------------- #
"""
    formula_to_pla(
        formula::SoleLogics.Formula;
        allow_scalar_range_conditions::Bool=false, kwargs...
    ) -> (String, Vector{VariableValue})
    formula_to_pla(
        dnfformula::SoleLogics.DNF;
        allow_scalar_range_conditions::Bool=false, kwargs...
    ) -> (String, Vector{VariableValue})
    formula_to_pla(
        atoms::Vector{Vector{SoleLogics.Atom}};
        encoding::Symbol=:univariate,
        allow_scalar_range_conditions::Bool=false,
        offset::Union{Nothing,Vector{Vector{SoleLogics.Atom}}}=nothing,
        kwargs...
    ) -> (String, Vector{VariableValue})

Convert a logical formula into Programmable Logic Array (PLA) format.
See original module docstring for the full step-by-step description --
unchanged except for the NEW `offset` keyword documented below.

# `offset` (NEW, only on the `atoms::Vector{Vector{Atom}}` method)
If `nothing` (default): behavior IDENTICAL to before -- implicit `.type f`
PLA, Espresso computes the absolute complement. No change for any existing
caller that doesn't pass `offset`.

If given: `offset` is, like `atoms`, a `Vector{Vector{Atom}}` -- one cube
per row -- but represents cubes CONFIRMED OFF (not "unknown"). It is
encoded with `_encode_disjunct` using EXACTLY the same condition space
(`conditions`, `includes`, `excludes`, `feat_condindxss`) as the ON-set, so
PLA columns stay aligned between the two halves, then its rows are emitted
with output `"0"` instead of `"1"`. The `.type fr` header line is also
emitted, which is what tells Espresso "everything else is don't-care, not
complement".

Conditions that appear ONLY in `offset` (not in `atoms`) are still folded
into the global condition space (same mechanism already used for
`universe_conditions`), otherwise an off-set cube mentioning a
feature/threshold never seen in the on-set couldn't be encoded correctly.
"""
function write_pla(
    atoms::Vector{Vector{LumenAtom}},
    feat_idxs::Vector{R},
    encoding::Type{<:AbstractEncoding}=Univariate,
    # removewhitespaces::Bool=true,
    # pretty_op::Bool=false,
    # universe_conditions::Union{Nothing,Vector{<:SD.AbstractScalarCondition}}=nothing,
) where {R<:Unsigned}
    nfnames = length(feat_idxs)

    # sort!(conditions; by=SD._scalarcondition_sortby)
    # sort!(feat_idxs; by=syntaxstring)

    # conditions = SD.removeduals(conditions)

    # feat_condindxss = Vector{Vector{Int}}(undef, nfnames)
    # feat_condnames = Vector{Vector{String}}(undef, nfnames)

    # @inbounds for (i, feat) in enumerate(feat_idxs)
    #     feat_condindxs = findall(c->SD.feature(c) == feat, conditions)
    #     conds = filter(c->SD.feature(c) == feat, conditions)
    #     condname = SoleLogics.syntaxstring.(conds; removewhitespaces, pretty_op)

    #     feat_condindxss[i] = feat_condindxs
    #     feat_condnames[i] = condname
    # end

    # feat_nconds = length.(feat_condindxss)

    # includes, excludes = Vector{BitMatrix}(undef, nfnames),
    # Vector{BitMatrix}(undef, nfnames)
    # @inbounds for (i, feat_condindxs) in enumerate(feat_condindxss)
    #     includes[i] = BitMatrix([
    #         SD.includes(conditions[cond_i], conditions[cond_j]) for
    #         cond_i in feat_condindxs, cond_j in feat_condindxs
    #     ])
    #     excludes[i] = BitMatrix([
    #         SD.excludes(conditions[cond_j], conditions[cond_i]) for
    #         cond_i in feat_condindxs, cond_j in feat_condindxs
    #     ])
    # end

    # pla_header = if encoding == :multivariate
    #     _header(feat_nconds, feat_condnames; has_offset)
    # else
    #     _header(conditions, feat_condnames; has_offset)
    # end

    # # --- ON-set rows -----------------------------------------------------
    # conjuncts = _get_conjuncts(atoms)
    # pla_onset_rows = Vector{String}(undef, length(conjuncts))
    # Threads.@threads for i in eachindex(conjuncts)
    #     row = _encode_disjunct(
    #         conjuncts[i], feat_idxs, conditions, includes, excludes, feat_condindxss
    #     )
    #     pla_onset_rows[i] =
    #         encoding == :multivariate ? _onset_rows(feat_nconds, row) : _onset_rows(row)
    # end

    # # --- OFF-set rows (NEW, only if `offset` given) -----------------------
    # # Same conditions/includes/excludes/feat_condindxss as the ON-set: this
    # # reuse is exactly what guarantees column alignment between the two
    # # halves of the 
    # pla_offset_rows = String[]
    # if has_offset
    #     offset_conjuncts = _get_conjuncts(offset)
    #     pla_offset_rows = Vector{String}(undef, length(offset_conjuncts))
    #     Threads.@threads for i in eachindex(offset_conjuncts)
    #         row = _encode_disjunct(
    #             offset_conjuncts[i], feat_idxs, conditions, includes, excludes, feat_condindxss
    #         )
    #         pla_offset_rows[i] =
    #             encoding == :multivariate ? _offset_rows(feat_nconds, row) : _offset_rows(row)
    #     end
    # end

    # all_rows = vcat(pla_onset_rows, pla_offset_rows)

    # pla_content = join(
    #     [
    #         join(pla_header, "\n"),
    #         ".p $(length(all_rows))",
    #         join(all_rows, "\n"),
    #         ".e",
    #     ],
    #     "\n",
    # )

    # return pla_content, feat_idxs
end