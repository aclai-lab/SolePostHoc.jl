@inline evalst(::Fast) = 0x01
@inline evalst(::Balanced) = 0x02
@inline evalst(::Sop) = 0x03

# ---------------------------------------------------------------------------- #
#                                formula to pla                                #
# ---------------------------------------------------------------------------- #
# function formula_to_pla(
#     dnfformula::SL.DNF;
#     allow_scalar_range_conditions::Bool=false,
#     kwargs...
# )
#     dnfformula = SD.scalar_simplification(dnfformula; allow_scalar_range_conditions)
#     dnfformula = SL.dnf(dnfformula; profile=:nnf, allow_atom_flipping=true, kwargs...)

#     atoms_per_disjunct = Vector{Vector{SL.Atom}}([
#         collect(SL.atoms(d)) for d in SL.disjuncts(dnfformula)
#     ])

#     formula_to_pla(atoms_per_disjunct; allow_scalar_range_conditions, kwargs...)
# end

function formula_to_pla(
    atoms::Vector{LumenAtom{R,T}}
) where {R<:Unsigned,T<:AbstractFloat}
#     fnames = unique(SD.feature.(conditions))
#     nfnames = length(fnames)

#     sort!(conditions; by=SD._scalarcondition_sortby)
#     sort!(fnames; by=syntaxstring)

#     conditions = SD.removeduals(conditions)

#     feat_condindxss = Vector{Vector{Int}}(undef, nfnames)
#     feat_condnames = Vector{Vector{String}}(undef, nfnames)

#     @inbounds for (i, feat) in enumerate(fnames)
#         feat_condindxs = findall(c->SD.feature(c) == feat, conditions)
#         conds = filter(c->SD.feature(c) == feat, conditions)
#         condname = SoleLogics.syntaxstring.(conds; removewhitespaces, pretty_op)

#         feat_condindxss[i] = feat_condindxs
#         feat_condnames[i] = condname
#     end

#     feat_nconds = length.(feat_condindxss)

#     includes, excludes = Vector{BitMatrix}(undef, nfnames),
#     Vector{BitMatrix}(undef, nfnames)
#     @inbounds for (i, feat_condindxs) in enumerate(feat_condindxss)
#         includes[i] = BitMatrix([
#             SD.includes(conditions[cond_i], conditions[cond_j]) for
#             cond_i in feat_condindxs, cond_j in feat_condindxs
#         ])
#         excludes[i] = BitMatrix([
#             SD.excludes(conditions[cond_j], conditions[cond_i]) for
#             cond_i in feat_condindxs, cond_j in feat_condindxs
#         ])
#     end

#     pla_header = if encoding == :multivariate
#         _header(feat_nconds, feat_condnames; has_offset)
#     else
#         _header(conditions, feat_condnames; has_offset)
#     end

#     # --- ON-set rows -----------------------------------------------------
#     conjuncts = _get_conjuncts(atoms)
#     pla_onset_rows = Vector{String}(undef, length(conjuncts))
#     Threads.@threads for i in eachindex(conjuncts)
#         row = _encode_disjunct(
#             conjuncts[i], fnames, conditions, includes, excludes, feat_condindxss
#         )
#         pla_onset_rows[i] =
#             encoding == :multivariate ? _onset_rows(feat_nconds, row) : _onset_rows(row)
#     end

#     # --- OFF-set rows (NEW, only if `offset` given) -----------------------
#     # Same conditions/includes/excludes/feat_condindxss as the ON-set: this
#     # reuse is exactly what guarantees column alignment between the two
#     # halves of the 
#     pla_offset_rows = String[]
#     if has_offset
#         offset_conjuncts = _get_conjuncts(offset)
#         pla_offset_rows = Vector{String}(undef, length(offset_conjuncts))
#         Threads.@threads for i in eachindex(offset_conjuncts)
#             row = _encode_disjunct(
#                 offset_conjuncts[i], fnames, conditions, includes, excludes, feat_condindxss
#             )
#             pla_offset_rows[i] =
#                 encoding == :multivariate ? _offset_rows(feat_nconds, row) : _offset_rows(row)
#         end
#     end

#     all_rows = vcat(pla_onset_rows, pla_offset_rows)

#     pla_content = join(
#         [
#             join(pla_header, "\n"),
#             ".p $(length(all_rows))",
#             join(all_rows, "\n"),
#             ".e",
#         ],
#         "\n",
#     )

#     return pla_content, fnames
end

# # ---------------------------------------------------------------------------- #
# #                                pla to formula                                #
# # ---------------------------------------------------------------------------- #
# """
#     pla_to_formula(
#         pla::String,
#         fnames::Vector{<:VariableValue};
#         conditionstype::Type=SoleData.ScalarCondition,
#         conjunct::Bool=false,
#         float_type::Type=Float64
#     ) -> Union{SoleLogics.Formula, Vector{SyntaxStructure}}

# Convert a PLA string back into a formula. Unchanged by this patch: it
# already only reads rows starting with `['0','1','-','|']` and, of those,
# only the ones ending "1" ever become atoms in the ON-set fed back to the
# caller (`_onset_rows`/`_offset_rows` distinction is irrelevant here, since
# Espresso's OUTPUT PLA -- what this function parses -- is always the
# minimized single ON-set again, regardless of whether the INPUT PLA was
# `.type f` or `.type fr`).
# """
# function pla_to_formula(
#     pla::String,
#     fnames::Vector{<:VariableValue};
#     conditionstype::Type=SD.ScalarCondition,
#     conjunct::Bool=false,
#     float_type::Type=Float32
# )
#     lines = split(pla, '\n')
#     parsed_conditions = SoleLogics.Atom[]
#     binaries = String[]

#     for line in lines
#         startswith(line, ".ilb") &&
#             append!(parsed_conditions, _read_conditions(line, conditionstype, fnames; float_type))
#         startswith(line, ['0', '1', '-', '|']) && append!(binaries, [line[1:(end-2)]])
#     end

#     isempty(binaries) && return ⊤

#     disjuncts = Vector{SyntaxStructure}(undef, length(binaries))

#     Threads.@threads for i in eachindex(binaries)
#         binary = binaries[i]
#         disjuncts[i] = SD.scalar_simplification(
#             SL.LeftmostConjunctiveForm([
#                 SL.Literal(LiteralBool[value], parsed_conditions[idx]) for
#                 (idx, value) in enumerate(binary) if value ∈ ['1', '0']
#             ]);
#             allow_scalar_range_conditions=false
#         )
#     end

#     return conjunct ? LeftmostDisjunctiveForm(disjuncts) : disjuncts
# end

# ---------------------------------------------------------------------------- #
function abc_minimize(
    config::LumenShannonConfig{R,T},
    binary::String,
    atoms::Vector{LumenAtom{R,T}},
) where {R<:Unsigned,T<:AbstractFloat}
    # convert formula to pla string format
    # pla_string, fnames = formula_to_pla(atoms)
    formula_to_pla(atoms)

    # # create temporary files for input/output
    # mktempdir() do tmp
    #     inputfile = joinpath(tmp, "in.pla")
    #     outputfile = joinpath(tmp, "out.pla")

    #     write(inputfile, pla_string)

    #     abc_commands = if setup == 1
    #         "read $inputfile; strash; collapse; write $outputfile"
    #     elseif setup == 0
    #         "read $inputfile; strash; balance; rewrite; refactor; " *
    #         "balance; rewrite -z; collapse; sop; fx; strash; " *
    #         "balance; collapse; write $outputfile"
    #     else
    #         "read $inputfile; sop; strash; dc2; collapse; " *
    #         "strash; dc2; collapse; sop; write $outputfile"
    #     end

    #     # Execute ABC with error handling
    #     try
    #         run(`$binary -c $abc_commands`)
    #     catch e
    #         return LeftmostConjunctiveForm.(atoms)
    #     end

    #     minimized_pla_raw = read(outputfile, String)
    #     minimized_pla_raw = replace(minimized_pla_raw, ">=" => "≥")
    #     isempty(strip(minimized_pla_raw)) && return atoms

    #     minimized_pla = clean_abc_output(minimized_pla_raw)
    #     conditionstype = allow_scalar_range_conditions ?
    #                      SoleData.RangeScalarCondition :
    #                      SoleData.ScalarCondition

    #     return pla_to_formula(minimized_pla, fnames; conditionstype, float_type)
    # end
end

# ---------------------------------------------------------------------------- #
function run_minimization(
    ::Type{Abc},
    config::LumenShannonConfig{R,T},
    atoms::Vector{LumenAtom{R,T}}
) where {R<:Unsigned,T<:AbstractFloat}
    ABC_jll.abc() do binary
        minimized_formula = abc_minimize(
            config,
            binary,
            atoms
        )
    #     return _as_terms(refine_dnf(minimized_formula))
    end
end
