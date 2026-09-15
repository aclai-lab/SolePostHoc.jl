# ---------------------------------------------------------------------------- #
#                                formula to pla                                #
# ---------------------------------------------------------------------------- #

"""
    _colname(c::Integer) -> String

Label for PLA column `c`. Deliberately opaque and position-free: the minimisers
echo `.ilb` back in their output and may reorder or drop inputs, so the result
is decoded by NAME, never by position.
"""
@inline _colname(c::Integer) = string('x', c)

"""
    write_pla(io, sp, bs, L; offset_bs=nothing)

Write the PLA describing box set `bs` over layout `L`.

Each box becomes one product term. Two encodings are available, and which one
compresses better is a property of the backend, not of the mathematics — both
describe the same real region.

- `:lax` (default) pins only the box's two OWN boundary columns per feature —
  `'1'` at `V < thr[lo-1]`, `'0'` at `V ≥ thr[hi]` — and leaves every other
  column don't-care. That is sound because an assignment that disagrees with
  those two bounds on some other threshold of the same feature describes no
  real value at all, so the extra minterms swept in are unreachable.
- `:tight` pins every column the box determines: `'1'` below its upper bound,
  `'0'` at or above its lower one. The cube then covers exactly the box's own
  minterms and nothing else.

Passing `offset_bs` emits an explicit OFF-set and a `.type fr` header, which
tells the minimiser that everything in neither set is a genuine don't-care
rather than a zero. Combined with `:tight` that states the problem exactly:
this class ON, the other classes OFF, the impossible assignments free. Without
it the PLA is an ordinary ON-set-only `.type f` problem and the complement is
implied.
"""
function write_pla(
    io::IO,
    sp::ThresholdSpace,
    bs::BoxSet,
    L::PlaLayout;
    offset_bs::Union{Nothing,BoxSet}=nothing,
    encoding::Symbol=:lax
)
    nc = ncolumns(L)
    nrows = bs.n + (isnothing(offset_bs) ? 0 : offset_bs.n)

    print(io, ".i ", nc, "\n.o 1\n.ilb ")
    for c in 1:nc
        c > 1 && print(io, ' ')
        print(io, _colname(c))
    end
    print(io, "\n.ob lumen_out\n")
    isnothing(offset_bs) || print(io, ".type fr\n")
    print(io, ".p ", nrows, "\n")

    row = Vector{UInt8}(undef, nc)
    _write_rows(io, sp, bs, L, row, UInt8('1'), encoding)
    isnothing(offset_bs) || _write_rows(io, sp, offset_bs, L, row, UInt8('0'), encoding)
    print(io, ".e\n")
    return nothing
end

"""
    _write_rows(io, sp, bs, L, row, outchar)

Emit every box of `bs` as one PLA product term with output character `outchar`.
`row` is caller-owned scratch of length `ncolumns(L)`, reused across rows so the
whole emission allocates nothing per term.
"""
function _write_rows(
    io::IO, sp::ThresholdSpace, bs::BoxSet, L::PlaLayout,
    row::Vector{UInt8}, outchar::UInt8, encoding::Symbol
)
    nfeat = nfeatures(sp)
    nc = ncolumns(L)
    @inbounds for i in 1:bs.n
        lo = boxlo(bs, i); hi = boxhi(bs, i)

        if encoding === :tight
            for c in 1:nc
                j = L.colfeat[c]
                k = L.colthr[c]
                row[c] = k <= lo[j] - Int32(1) ? UInt8('1') :
                         k >= hi[j]           ? UInt8('0') : UInt8('-')
            end
        else
            fill!(row, UInt8('-'))
            for j in 1:nfeat
                if lo[j] > 1
                    c = _column(L, j, lo[j] - Int32(1))
                    c != 0 && (row[c] = UInt8('1'))
                end
                if hi[j] <= Int32(nthresholds(sp, j))
                    c = _column(L, j, hi[j])
                    c != 0 && (row[c] = UInt8('0'))
                end
            end
        end

        write(io, row)
        write(io, UInt8(' '), outchar, UInt8('\n'))
    end
    return nothing
end

"""
    pla_string(sp, bs, L; offset_bs=nothing) -> String

`write_pla` into a buffer sized from the problem, so the text is built once
without repeated string growth.
"""
function pla_string(
    sp::ThresholdSpace, bs::BoxSet, L::PlaLayout;
    offset_bs::Union{Nothing,BoxSet}=nothing,
    encoding::Symbol=:lax
)
    nrows = bs.n + (isnothing(offset_bs) ? 0 : offset_bs.n)
    buf = IOBuffer(; sizehint = 64 + (ncolumns(L) + 4) * (nrows + 4))
    write_pla(buf, sp, bs, L; offset_bs, encoding)
    return String(take!(buf))
end

# ---------------------------------------------------------------------------- #

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
function formula_to_pla(
    atoms::Vector{Vector{LumenAtom}};
    # encoding::Symbol=:univariate,
    # allow_scalar_range_conditions::Bool=false,
    # removewhitespaces::Bool=true,
    # pretty_op::Bool=false,
    # universe_conditions::Union{Nothing,Vector{<:SD.AbstractScalarCondition}}=nothing,
    # offset::Union{Nothing,Vector{Vector{SL.Atom}}}=nothing,
)
    # @assert encoding in [:univariate, :multivariate]

    # has_offset = !isnothing(offset) && !isempty(offset)

    # # extract domains: ON-set + (if present) OFF-set, so conditions
    # # mentioned ONLY in the offset still enter the shared condition space.
    # local_conditions = map(SL.value, reduce(vcat, atoms))
    # if has_offset
    #     offset_conditions = map(SL.value, reduce(vcat, offset))
    #     local_conditions = vcat(local_conditions, offset_conditions)
    # end

    # if isnothing(universe_conditions)
    #     conditions = unique(local_conditions)
    # else
    #     conditions = unique(vcat(collect(universe_conditions), local_conditions))
    # end

    # fnames = unique(SD.feature.(conditions))
    # nfnames = length(fnames)

    # sort!(conditions; by=SD._scalarcondition_sortby)
    # sort!(fnames; by=syntaxstring)

    # if allow_scalar_range_conditions
    #     original_conditions = conditions
    #     conditions = SD.scalartiling(conditions, fnames)
    #     @assert length(setdiff(original_conditions, conditions)) == 0
    # end

    # conditions = SD.removeduals(conditions)

    # feat_condindxss = Vector{Vector{Int}}(undef, nfnames)
    # feat_condnames = Vector{Vector{String}}(undef, nfnames)

    # @inbounds for (i, feat) in enumerate(fnames)
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
    #         conjuncts[i], fnames, conditions, includes, excludes, feat_condindxss
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
    #             offset_conjuncts[i], fnames, conditions, includes, excludes, feat_condindxss
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

    # return pla_content, fnames
end