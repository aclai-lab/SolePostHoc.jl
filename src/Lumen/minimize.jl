@inline evalst(::Fast) = 0x01
@inline evalst(::Balanced) = 0x02
@inline evalst(::Sop) = 0x03

# function _write_rows(
#     io::IO, sp::ThresholdSpace, bs::BoxSet, L::PlaLayout,
#     row::Vector{UInt8}, outchar::UInt8, encoding::Symbol
# )
#     nfeat = nfeatures(sp)
#     nc = ncolumns(L)
#     @inbounds for i in 1:bs.n
#         lo = boxlo(bs, i); hi = boxhi(bs, i)

#         if encoding === :tight
#             for c in 1:nc
#                 j = L.colfeat[c]
#                 k = L.colthr[c]
#                 row[c] = k <= lo[j] - Int32(1) ? UInt8('1') :
#                          k >= hi[j]           ? UInt8('0') : UInt8('-')
#             end
#         else
#             fill!(row, UInt8('-'))
#             for j in 1:nfeat
#                 if lo[j] > 1
#                     c = _column(L, j, lo[j] - Int32(1))
#                     c != 0 && (row[c] = UInt8('1'))
#                 end
#                 if hi[j] <= Int32(nthresholds(sp, j))
#                     c = _column(L, j, hi[j])
#                     c != 0 && (row[c] = UInt8('0'))
#                 end
#             end
#         end

#         write(io, row)
#         write(io, UInt8(' '), outchar, UInt8('\n'))
#     end
#     return nothing
# end

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

function get_includes(a::LumenAtom, b::LumenAtom)
    (a.feat == b.feat) || return false
    return true
    # return issubset(tointervalset(b), tointervalset(a))
end

# function get_excludes(a::LumenAtom, b::LumenAtom)
#     (feature(a) == feature(b)) || return false
#     # @show tointervalset(a)
#     # @show tointervalset(b)
#     return isdisjoint(tointervalset(a), tointervalset(b))
# end

function formula_to_pla(
    buf::IOBuffer,
    atoms::Vector{LumenAtom{R,T}},
    cube::LumenCubeVec{R,T},
    feat_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    grouped = Dict(
    f => filter(a -> a.feat == f, atoms)
        for f in unique(getfield.(atoms, :feat)))

    natoms = length(atoms)
    nrows = length(cube)
    ngrp = length(grouped)

    print(buf, ".i ", natoms, "\n.o 1\n.ilb ")

    # for a in atoms
    #     print(buf, " [", a.feat, "]", evalop(a.op), a.thr)
    # end

    for i in 1:natoms
        print(buf, 'x', i, ' ')
    end

    print(buf, "\n.ob formula_output\n")
    print(buf, ".p ", nrows, "\n")

    # row = Vector{UInt8}(undef, natoms)
    # _write_rows(buf, sp, bs, L, row, UInt8('1'), encoding)
    # isnothing(offset_bs) || _write_rows(buf, sp, offset_bs, L, row, UInt8('0'), encoding)
    # print(buf, ".e\n")

    @show ngrp
    includes = Vector{BitMatrix}(undef, ngrp)
    excludes = Vector{BitMatrix}(undef, ngrp)
    @inbounds for (i, g) in enumerate(grouped)
        @show typeof(g)
        includes[i] = BitMatrix([
            get_includes(cond_i, cond_j) for
            cond_i in g, cond_j in g
        ])
        # excludes[i] = BitMatrix([
        #     get_excludes(conditions[cond_j], conditions[cond_i]) for
        #     cond_i in g, cond_j in g
        # ])
    end

    @show String(take!(buf))
    return nothing

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
    cube::LumenCubeVec{R,T},
    feat_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}

    # buf = IOBuffer(; sizehint = 64 + (ncolumns(L) + 4) * (nrows + 4))
    buf = IOBuffer()

    # convert formula to pla string format
    # pla_string, fnames = formula_to_pla(atoms)
    formula_to_pla(buf, atoms, cube, feat_idxs)

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
    atoms::Vector{LumenAtom{R,T}},
    cube::LumenCubeVec{R,T},
    feat_idxs::Vector{R}
) where {R<:Unsigned,T<:AbstractFloat}
    ABC_jll.abc() do binary
        minimized_formula = abc_minimize(
            config,
            binary,
            atoms,
            cube,
            feat_idxs
        )
    #     return _as_terms(refine_dnf(minimized_formula))
    end
end
