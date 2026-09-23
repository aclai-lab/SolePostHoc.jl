# ---------------------------------------------------------------------------- #
#                        PLA round trip on Lumen structures                    #
#                                                                              #
#  Everything the minimiser sees is derived from `LumenAtom`s and the cube     #
#  views `_leaf_extract` hands out; nothing here touches SoleData conditions   #
#  or SoleLogics formulas. The only conversion to Sole types happens once, at  #
#  the very end of `lumen_shannon` (see `decision_set` in lumen_shannon.jl).   #
#                                                                              #
#  ENCODING. One PLA column per distinct (feature, threshold) pair that occurs #
#  in the cubes, columns of a feature grouped in a block and sorted by         #
#  threshold DESCENDING. Column value '1' means "x < thr" (or "x ≤ thr" for   #
#  the '≤' family), '0' its dual. For a fixed feature the consistent column    #
#  patterns are then the thermometer codes 1…10…0: once "x < t" is false for  #
#  some t it is false for every smaller t. Every other pattern is physically   #
#  impossible and is declared as DON'T CARE with one PLA row per adjacent      #
#  column pair ("0" at k, "1" at k+1). Declaring them matters: a box spanning  #
#  three or more regions is a single product term only if the minimiser may   #
#  cover impossible patterns for free. Without those rows Espresso would need  #
#  one term per pair of regions.                                               #
#                                                                              #
#  A cube (one region per feature) is written TIGHT: every column of every    #
#  feature is decided, i.e. one minterm per row. A product term read back is  #
#  decoded per block as "last '1' → upper bound, first '0' → lower bound", so  #
#  the minimiser's use of implied columns costs nothing and inconsistent       #
#  terms (a '0' before a '1', an empty interval) are simply dropped.           #
#                                                                              #
#  BACKENDS. Espresso reads `.type fd` natively and returns an irredundant     #
#  cover. ABC's PLA reader ignores the dc-set unless asked to fold it into the #
#  on-set (`read_pla -d`), which makes it cover the impossible patterns too;   #
#  the surplus terms either decode to empty intervals or are dominated by a    #
#  real one, and `refine!` removes both. There is no Espresso jll in General,  #
#  so the MIT binary comes from SoleData's artifact loader.                    #
# ---------------------------------------------------------------------------- #

@inline evalst(::Type{Fast}) = 0x01
@inline evalst(::Type{Balanced}) = 0x02
@inline evalst(::Type{Sop}) = 0x03
@inline evalst(s::AbstractMinimizeSetup) = evalst(typeof(s))

# '<' and '≤' bound a region from above, '≥' and '>' from below.
@inline isupper(op::UInt8) = op == 0x01 || op == 0x04
@inline upperop(op::UInt8) = isupper(op) ? op : dualop(op)

# ---------------------------------------------------------------------------- #
#                                  Lumen DNF                                   #
# ---------------------------------------------------------------------------- #
"""
    LumenDNF{R,T}

A disjunction of conjunctions of `LumenAtom`s, stored flat: term `i` is
`atoms[offsets[i]+1:offsets[i+1]]`. Indexing returns a view, so a term has
the same type as a `LumenCube` and can be fed back to `write_pla`.
An empty term is the tautology.
"""
struct LumenDNF{R<:Unsigned,T<:AbstractFloat}
    atoms::Vector{LumenAtom{R,T}}
    offsets::Vector{Int}
end

LumenDNF{R,T}() where {R<:Unsigned,T<:AbstractFloat} =
    LumenDNF{R,T}(LumenAtom{R,T}[], Int[0])

@inline nterms(d::LumenDNF) = length(d.offsets) - 1
Base.length(d::LumenDNF) = nterms(d)
Base.isempty(d::LumenDNF) = nterms(d) == 0
Base.eltype(::Type{LumenDNF{R,T}}) where {R,T} = LumenCube{R,T}
Base.@propagate_inbounds Base.getindex(d::LumenDNF, i::Int) =
    view(d.atoms, d.offsets[i]+1:d.offsets[i+1])
Base.iterate(d::LumenDNF, i::Int=1) = i > nterms(d) ? nothing : (d[i], i + 1)
Base.firstindex(::LumenDNF) = 1
Base.lastindex(d::LumenDNF) = nterms(d)
natoms(d::LumenDNF) = length(d.atoms)

function Base.push!(d::LumenDNF{R,T}, term::AbstractVector{LumenAtom{R,T}}) where {R,T}
    append!(d.atoms, term)
    push!(d.offsets, length(d.atoms))
    return d
end

function Base.append!(d::LumenDNF{R,T}, o::LumenDNF{R,T}) where {R,T}
    base = length(d.atoms)
    append!(d.atoms, o.atoms)
    @inbounds for i in 2:length(o.offsets)
        push!(d.offsets, base + o.offsets[i])
    end
    return d
end

function Base.empty!(d::LumenDNF)
    empty!(d.atoms)
    resize!(d.offsets, 1)
    d.offsets[1] = 0
    return d
end

function Base.show(io::IO, d::LumenDNF{R,T}) where {R,T}
    print(io, "LumenDNF{", R, ",", T, "}(", nterms(d), " terms, ", natoms(d), " atoms)")
end

# ---------------------------------------------------------------------------- #
#                                 PLA layout                                   #
# ---------------------------------------------------------------------------- #
"""
    PlaLayout{R,T}

Column layout of one PLA: column `c` stands for `feat[c] <upop[c]> thr[c]`
(`upop` is the family's upper-bound operator, `'<'` or `'≤'`); its dual is
the '0' polarity. Columns of a feature are contiguous, `bfirst[c]:blast[c]`,
and sorted by threshold descending. `col` maps `(feature, threshold)` back
to the column.
"""
struct PlaLayout{R<:Unsigned,T<:AbstractFloat}
    feat::Vector{R}
    thr::Vector{T}
    upop::Vector{UInt8}
    bfirst::Vector{Int}
    blast::Vector{Int}
    col::Dict{Tuple{R,T},Int}
end

@inline ncolumns(L::PlaLayout) = length(L.feat)
# number of don't-care rows: one per adjacent column pair inside a block
ndontcares(L::PlaLayout) = count(c -> L.blast[c] > c, 1:ncolumns(L))

"""
    PlaLayout(atoms)

Build the layout from the atoms that occur in the cubes. Dual atoms (same
feature and threshold, complementary operator) share one column, so the
input may or may not have been through `removeduals!`.
"""
function PlaLayout(
    atoms::AbstractVector{LumenAtom{R,T}}
) where {R<:Unsigned,T<:AbstractFloat}
    col = Dict{Tuple{R,T},Int}()
    feat = R[]
    thr = T[]
    upop = UInt8[]
    for a in atoms
        k = (a.feat, a.thr)
        haskey(col, k) && continue
        col[k] = 0
        push!(feat, a.feat); push!(thr, a.thr); push!(upop, upperop(a.op))
    end
    n = length(feat)
    # a key vector, not a `by` closure: the closure boxes on every comparison
    perm = sortperm([(feat[c], -thr[c]) for c in 1:n])
    feat = feat[perm]; thr = thr[perm]; upop = upop[perm]
    bfirst = Vector{Int}(undef, n)
    blast = Vector{Int}(undef, n)
    @inbounds begin
        s = 1
        for c in 1:n
            col[(feat[c], thr[c])] = c
            c > 1 && feat[c] != feat[c-1] && (s = c)
            bfirst[c] = s
        end
        e = n
        for c in n:-1:1
            c < n && feat[c] != feat[c+1] && (e = c)
            blast[c] = e
        end
    end
    return PlaLayout{R,T}(feat, thr, upop, bfirst, blast, col)
end

# ---------------------------------------------------------------------------- #
#                                  write pla                                   #
# ---------------------------------------------------------------------------- #
# Fill `row[off+1:off+nc]` with the tight encoding of `term`. An upper bound
# at column c implies '1' on every column of the block with a larger
# threshold, a lower bound implies '0' on every column with a smaller one; the
# two never overlap because a cube's bounds are adjacent thresholds.
@inline function _encode_row!(
    row::Vector{UInt8},
    off::Int,
    L::PlaLayout{R,T},
    term::AbstractVector{LumenAtom{R,T}}
) where {R,T}
    nc = ncolumns(L)
    @inbounds for k in 1:nc
        row[off+k] = UInt8('-')
    end
    @inbounds for a in term
        c = L.col[(a.feat, a.thr)]
        if isupper(a.op)
            for k in L.bfirst[c]:c
                row[off+k] = UInt8('1')
            end
        else
            for k in c:L.blast[c]
                row[off+k] = UInt8('0')
            end
        end
    end
    return row
end

"""
    write_pla(io::IO, L::PlaLayout, cubes)

Stream the on-set `cubes` (any indexable collection of atom vectors, e.g. a
`LumenDNF`) as a `.type fd` PLA on `L`'s columns, followed by the don't-care
rows for the impossible column patterns. One row of scratch is the only
buffer: the text never exists in memory as a whole.
"""
function write_pla(
    io::IO,
    L::PlaLayout{R,T},
    cubes
) where {R<:Unsigned,T<:AbstractFloat}
    nc = ncolumns(L)
    print(io, ".i ", nc, "\n.o 1\n.ilb")
    for c in 1:nc
        print(io, " x", c)
    end
    print(io, "\n.ob f\n.type fd\n.p ", length(cubes) + ndontcares(L), '\n')

    row = Vector{UInt8}(undef, nc + 3)
    @inbounds row[nc+1] = UInt8(' '); row[nc+2] = UInt8('1'); row[nc+3] = UInt8('\n')
    for term in cubes
        _encode_row!(row, 0, L, term)
        write(io, row)
    end

    @inbounds row[nc+2] = UInt8('-')
    @inbounds for c in 1:nc
        L.blast[c] > c || continue
        for k in 1:nc
            row[k] = UInt8('-')
        end
        row[c] = UInt8('0'); row[c+1] = UInt8('1')
        write(io, row)
    end
    print(io, ".e\n")
    return nothing
end

# ---------------------------------------------------------------------------- #
#                                   read pla                                   #
# ---------------------------------------------------------------------------- #
# ".ilb x3 x1 x2" (bytes i:e) -> [3, 1, 2]: output column k is layout column
# perm[k]. Labels are the ones `write_pla!` emitted, so anything else is an
# error rather than a silent misdecode.
function _read_ilb!(perm::Vector{Int}, pla::AbstractVector{UInt8}, i::Int, e::Int, nc::Int)
    empty!(perm)
    p = i + 4                       # past ".ilb"
    @inbounds while p ≤ e
        if pla[p] == UInt8(' ')
            p += 1
        elseif pla[p] == UInt8('x') && p < e && UInt8('0') ≤ pla[p+1] ≤ UInt8('9')
            c = 0; p += 1
            while p ≤ e && UInt8('0') ≤ pla[p] ≤ UInt8('9')
                c = 10c + (pla[p] - UInt8('0')); p += 1
            end
            1 ≤ c ≤ nc || throw(ArgumentError("PLA column label x$c out of range"))
            push!(perm, c)
        else
            throw(ArgumentError("unexpected PLA column label in: " *
                                String(pla[i:e])))
        end
    end
    return perm
end

# Decode one product term (already permuted into layout column order, one
# byte per column) into atoms appended to `out`; false if the term is empty
# as a region, i.e. it has a '0' before a '1' inside some block.
function _decode_row!(
    out::Vector{LumenAtom{R,T}},
    vals::Vector{UInt8},
    L::PlaLayout{R,T}
) where {R,T}
    n0 = length(out)
    nc = ncolumns(L)
    c = 1
    @inbounds while c ≤ nc
        e = L.blast[c]
        last1 = 0; first0 = 0
        for k in c:e
            v = vals[k]
            if v == UInt8('1')
                last1 = k
            elseif v == UInt8('0') && first0 == 0
                first0 = k
            end
        end
        if first0 != 0 && last1 != 0 && first0 < last1
            resize!(out, n0)
            return false
        end
        last1 != 0 && push!(out, LumenAtom{R,T}(L.feat[last1], L.thr[last1], L.upop[last1]))
        first0 != 0 && push!(out, LumenAtom{R,T}(L.feat[first0], L.thr[first0], dualop(L.upop[first0])))
        c = e + 1
    end
    return true
end

"""
    read_pla(pla::AbstractVector{UInt8}, L::PlaLayout) -> LumenDNF
    read_pla(pla::AbstractString, L::PlaLayout) -> LumenDNF

Parse a minimiser's output PLA back into terms on `L`'s columns. Only rows
with output `1` are terms; rows whose block pattern is inconsistent are
dropped. Handles a permuted `.ilb` and the `# comment` lines ABC emits. Works
on the raw bytes, one pass, no per-line strings.
"""
function read_pla(pla::AbstractVector{UInt8}, L::PlaLayout{R,T}) where {R<:Unsigned,T<:AbstractFloat}
    nc = ncolumns(L)
    perm = collect(1:nc)
    vals = Vector{UInt8}(undef, nc)
    dnf = LumenDNF{R,T}()
    term = LumenAtom{R,T}[]
    n = length(pla)
    i = firstindex(pla)

    @inbounds while i ≤ n
        j = findnext(==(UInt8('\n')), pla, i)
        e = isnothing(j) ? n : j - 1          # line is pla[i:e]
        next = e + 2
        if e < i
            i = next; continue
        end
        h = pla[i]
        if h == UInt8('.')
            e - i ≥ 3 && pla[i+1] == UInt8('i') && pla[i+2] == UInt8('l') &&
                pla[i+3] == UInt8('b') && _read_ilb!(perm, pla, i, e, nc)
            i = next; continue
        end
        if h == UInt8('0') || h == UInt8('1') || h == UInt8('-')
            sp = findnext(==(UInt8(' ')), pla, i)
            if !isnothing(sp) && sp ≤ e
                sp - i == length(perm) ||
                    throw(ArgumentError("PLA row has $(sp - i) inputs, expected $(length(perm))"))
                ob = sp + 1
                while ob ≤ e && pla[ob] == UInt8(' ')
                    ob += 1
                end
                if ob ≤ e && pla[ob] == UInt8('1')
                    for k in 1:nc
                        vals[perm[k]] = pla[i+k-1]
                    end
                    empty!(term)
                    _decode_row!(term, vals, L) && push!(dnf, term)
                end
            end
        end
        i = next
    end
    return dnf
end

read_pla(pla::AbstractString, L::PlaLayout) = read_pla(codeunits(pla), L)

# ---------------------------------------------------------------------------- #
#                             dominance refinement                             #
# ---------------------------------------------------------------------------- #
# a ⊇ b as regions: every bound of `a` is matched in `b` by a bound on the
# same feature and side that is at least as tight. Terms hold at most one
# bound per (feature, side), so this is a plain double loop.
function dominates(
    a::AbstractVector{LumenAtom{R,T}},
    b::AbstractVector{LumenAtom{R,T}}
) where {R,T}
    @inbounds for x in a
        up = isupper(x.op)
        ok = false
        for y in b
            y.feat == x.feat && isupper(y.op) == up || continue
            ok = up ? y.thr ≤ x.thr : y.thr ≥ x.thr
            ok && break
        end
        ok || return false
    end
    return true
end

# (upper, lower) bound signature of a term, one bit per feature modulo 64.
# `a ⊇ b` needs every (feature, side) of `a` to be present in `b`, so
# `mask(a) ⊆ mask(b)` is a necessary condition that costs two ANDs.
function _sidemask(term::AbstractVector{LumenAtom{R,T}}) where {R,T}
    up = zero(UInt64); lo = zero(UInt64)
    @inbounds for a in term
        bit = one(UInt64) << ((a.feat - one(R)) & 63)
        isupper(a.op) ? (up |= bit) : (lo |= bit)
    end
    return up, lo
end

"""
    refine!(d::LumenDNF) -> d

Drop every term contained in another term of `d` (duplicates included).
Coverage is unchanged; only redundant terms go.
"""
function refine!(d::LumenDNF{R,T}) where {R,T}
    n = nterms(d)
    n ≤ 1 && return d
    masks = [_sidemask(d[i]) for i in 1:n]
    # a term can only be dominated by one with no more bounds than itself
    order = sortperm(1:n; by=i -> d.offsets[i+1] - d.offsets[i])
    keep = trues(n)
    @inbounds for oi in 1:n
        i = order[oi]
        keep[i] || continue
        ti = d[i]
        ui, li = masks[i]
        for oj in oi+1:n
            j = order[oj]
            keep[j] || continue
            uj, lj = masks[j]
            (ui & ~uj == 0 && li & ~lj == 0) || continue
            dominates(ti, d[j]) && (keep[j] = false)
        end
    end
    all(keep) && return d

    atoms = LumenAtom{R,T}[]
    offsets = Int[0]
    @inbounds for i in 1:n
        keep[i] || continue
        append!(atoms, d[i])
        push!(offsets, length(atoms))
    end
    resize!(d.atoms, length(atoms)); copyto!(d.atoms, atoms)
    resize!(d.offsets, length(offsets)); copyto!(d.offsets, offsets)
    return d
end

# ---------------------------------------------------------------------------- #
#                            adjacent box merging                              #
# ---------------------------------------------------------------------------- #
# A term is a box: one interval per feature it mentions. Two boxes that agree
# on every feature but one, and whose intervals on that feature touch or
# overlap, have a box as their union, so replacing them by it changes nothing
# in what the rule covers. This is what the plain union of two cofactors leaves
# behind (siblings differing only in the split feature's complementary bound)
# and what a minimiser working on one leaf could never see.

# (feature, lower, upper, upper-bound op); a missing bound is ±Inf
const _Box{R,T} = Vector{Tuple{R,T,T,UInt8}}

function _tobox(term::AbstractVector{LumenAtom{R,T}}) where {R,T}
    box = _Box{R,T}()
    for a in term
        i = findfirst(e -> e[1] == a.feat, box)
        f, lo, hi, up = isnothing(i) ? (a.feat, typemin(T), typemax(T), upperop(a.op)) : box[i]
        isupper(a.op) ? (hi = min(hi, a.thr)) : (lo = max(lo, a.thr))
        isnothing(i) ? push!(box, (f, lo, hi, up)) : (box[i] = (f, lo, hi, up))
    end
    return sort!(box; by=first)
end

# hash of the box with entry `skip` left out
function _boxhash(box::_Box, skip::Int)
    h = hash(length(box) - 1)
    @inbounds for (i, e) in enumerate(box)
        i == skip && continue
        h = hash(e, h)
    end
    return h
end

function _sameexcept(a::_Box, ia::Int, b::_Box, ib::Int)
    # same feature set, same position of the skipped feature
    (length(a) == length(b) && ia == ib) || return false
    @inbounds for i in eachindex(a)
        i == ia && continue
        a[i] == b[i] || return false
    end
    return true
end

# one merging pass over `boxes`; dead boxes are emptied. Returns whether
# anything merged.
function _mergepass!(boxes::Vector{_Box{R,T}}) where {R,T}
    changed = false
    feats = unique!(sort!([e[1] for b in boxes for e in b]))
    buckets = Dict{UInt64,Vector{Int}}()
    for f in feats
        empty!(buckets)
        for (t, b) in enumerate(boxes)
            isempty(b) && continue
            i = findfirst(e -> e[1] == f, b)
            isnothing(i) && continue
            push!(get!(Vector{Int}, buckets, _boxhash(b, i)), t)
        end
        for ts in values(buckets)
            length(ts) ≥ 2 || continue
            # sweep along the feature by lower bound, absorbing neighbours
            sort!(ts; by=t -> boxes[t][findfirst(e -> e[1] == f, boxes[t])][2])
            for x in 1:length(ts)-1
                a = boxes[ts[x]]; isempty(a) && continue
                ia = findfirst(e -> e[1] == f, a)
                for y in x+1:length(ts)
                    b = boxes[ts[y]]; isempty(b) && continue
                    ib = findfirst(e -> e[1] == f, b)
                    _sameexcept(a, ia, b, ib) || continue
                    fa, loa, hia, up = a[ia]; _, lob, hib, _ = b[ib]
                    (hia ≥ lob && hib ≥ loa) || continue
                    lo = min(loa, lob); hi = max(hia, hib)
                    if lo == typemin(T) && hi == typemax(T)
                        deleteat!(a, ia)
                    else
                        a[ia] = (fa, lo, hi, up)
                    end
                    empty!(b)
                    changed = true
                    # `a` grew along `f` only, so its key in this bucket still
                    # holds and the sweep goes on; unless `f` became unbounded
                    # and left the box, in which case nothing else here can
                    # agree with it on `f`
                    ia = findfirst(e -> e[1] == f, a)
                    isnothing(ia) && break
                end
            end
        end
    end
    return changed
end

"""
    compact!(d::LumenDNF) -> d

Drop dominated terms and merge adjacent boxes, to a fixed point. Coverage
is unchanged; only the description gets shorter.
"""
function compact!(d::LumenDNF{R,T}) where {R,T}
    refine!(d)
    nterms(d) ≤ 1 && return d
    boxes = [_tobox(d[i]) for i in 1:nterms(d)]
    while _mergepass!(boxes)
    end
    empty!(d)
    for b in boxes
        isempty(b) && continue
        n0 = length(d.atoms)
        for (f, lo, hi, up) in b
            hi != typemax(T) && push!(d.atoms, LumenAtom{R,T}(f, hi, up))
            lo != typemin(T) && push!(d.atoms, LumenAtom{R,T}(f, lo, dualop(up)))
        end
        push!(d.offsets, length(d.atoms))
    end
    # an unbounded box is the tautology, and it stands alone
    if any(i -> d.offsets[i+1] == d.offsets[i], 1:nterms(d))
        empty!(d); push!(d.offsets, 0)
    end
    return refine!(d)
end

# ---------------------------------------------------------------------------- #
#                                  minimizer                                   #
# ---------------------------------------------------------------------------- #
"""
    Minimizer{MS<:AbstractMinimization}

A resolved backend: the command to run, the setup level and one scratch
directory for the whole run. Build it once per extraction
(`Minimizer(Abc, config)`, `Minimizer(MitEspresso, config)`), reuse it for
every leaf and class, and `close` it when done; `lumen_shannon` does so in a
`finally`, and Julia removes the directory at exit regardless.

Problems reach the backend as files, streamed row by row: the PLA text is
never held in memory, and a file redirect costs no pipe buffer.
"""
mutable struct Minimizer{MS<:AbstractMinimization}
    cmd::Cmd
    setup::UInt8
    dir::String
end

Minimizer(::Type{Abc}, config::LumenShannonConfig) =
    Minimizer{Abc}(ABC_jll.abc(), evalst(config.minimization_setup), mktempdir(; prefix="lumen_"))

function Minimizer(::Type{MitEspresso}, config::LumenShannonConfig)
    binary = joinpath(SD.load(SD.MITESPRESSOLoader()), "espresso")
    isfile(binary) || error("espresso binary not found at $binary")
    Minimizer{MitEspresso}(Cmd([binary]), evalst(config.minimization_setup), mktempdir(; prefix="lumen_"))
end

Base.close(m::Minimizer) = (rm(m.dir; force=true, recursive=true); nothing)

# a fresh input file name; classes are minimised concurrently, so one fixed
# name per session would not do
_tempfile(m::Minimizer) = tempname(m.dir; cleanup=false)

# ABC command scripts, one per setup level. `read_pla -d` folds the dc-set
# into the on-set (ABC has no other way to use it); `collapse` yields the SOP.
function _abc_script(setup::UInt8, inp::String, outp::String)
    body = setup == evalst(Fast) ?
        "strash; collapse" :
    setup == evalst(Balanced) ?
        "strash; balance; rewrite; refactor; balance; rewrite -z; " *
        "collapse; sop; fx; strash; balance; collapse" :
        "sop; strash; dc2; collapse; strash; dc2; collapse; sop"
    return "read_pla -d $inp; $body; write_pla $outp"
end

# Espresso flags per setup level: heuristic, strong heuristic, exact.
function _espresso_flags(setup::UInt8)
    setup == evalst(Fast) ? String[] :
    setup == evalst(Balanced) ? ["-estrong"] : ["-Dexact"]
end

# Every stream is a file, never a pipe: Julia backs each captured pipe with a
# 122 KiB buffer per spawn, which on a run of thousands of small leaves is
# the dominant allocation. The backend's diagnostics are read only on failure.
function _minimize(m::Minimizer{Abc}, L::PlaLayout, cubes)
    base = _tempfile(m)
    inp = base * ".in.pla"; outp = base * ".out.pla"; errp = base * ".err"
    try
        open(io -> write_pla(io, L, cubes), inp, "w")
        cmd = `$(m.cmd) -c $(_abc_script(m.setup, inp, outp))`
        ok = success(pipeline(cmd; stdout=devnull, stderr=errp))
        (ok && isfile(outp)) ||
            error("ABC failed:\n", isfile(errp) ? read(errp, String) : "")
        return refine!(read_pla(read(outp), L))
    finally
        rm(inp; force=true); rm(outp; force=true); rm(errp; force=true)
    end
end

function _minimize(m::Minimizer{MitEspresso}, L::PlaLayout, cubes)
    base = _tempfile(m)
    inp = base * ".in.pla"; outp = base * ".out.pla"; errp = base * ".err"
    try
        open(io -> write_pla(io, L, cubes), inp, "w")
        cmd = `$(m.cmd) $(_espresso_flags(m.setup)) $inp`
        ok = success(pipeline(cmd; stdout=outp, stderr=errp))
        ok || error("espresso failed:\n", isfile(errp) ? read(errp, String) : "")
        return read_pla(read(outp), L)
    finally
        rm(inp; force=true); rm(outp; force=true); rm(errp; force=true)
    end
end

# ---------------------------------------------------------------------------- #
#                               run minimization                               #
# ---------------------------------------------------------------------------- #
"""
    run_minimization(m::Minimizer, atoms, cubes) -> LumenDNF
    run_minimization(MS::Type{<:AbstractMinimization}, config, atoms, cubes)

Minimise the disjunction of `cubes` (atom vectors, e.g. one class' output of
`_leaf_extract`) over the columns spanned by `atoms`, which must contain
every (feature, threshold) pair the cubes use. Zero or one cube never reaches
the backend.
"""
function run_minimization(
    m::Minimizer,
    atoms::AbstractVector{LumenAtom{R,T}},
    cubes
) where {R<:Unsigned,T<:AbstractFloat}
    n = length(cubes)
    n == 0 && return LumenDNF{R,T}()
    n == 1 && return push!(LumenDNF{R,T}(), first(cubes))
    return _minimize(m, PlaLayout(atoms), cubes)
end

function run_minimization(
    ::Type{MS},
    config::LumenShannonConfig,
    atoms::AbstractVector{LumenAtom{R,T}},
    cubes
) where {MS<:AbstractMinimization,R<:Unsigned,T<:AbstractFloat}
    m = Minimizer(MS, config)
    try
        return run_minimization(m, atoms, cubes)
    finally
        close(m)
    end
end
