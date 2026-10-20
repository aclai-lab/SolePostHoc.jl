using Test
using SolePostHoc
using SoleModels, SoleData, SoleLogics
using CategoricalArrays, DataFrames, Random

const LU = SolePostHoc.Lumen
const SM = SoleModels

# ---------------------------------------------------------------------------- #
#                                   helpers                                    #
# ---------------------------------------------------------------------------- #
const R = UInt32
const T = Float64
const CLS = categorical(["a", "b", "c"])

A(f, thr, op) = LU.LumenAtom{R,T}(R(f), T(thr), LU.evalop(op))

atom(fn, f, op, t) = Atom(ScalarCondition(VariableValue(f, fn[f]), op, t))
leaf(c) = ConstantModel(CLS[c])

forest(trees, nfeat) = DecisionEnsemble(
    trees,
    (; featurenames=[Symbol("V", i) for i in 1:nfeat],
       supporting_labels=String.(levels(CLS))),
)

# three trees on two features, every operator family used at least once
function small_forest()
    fn = [:V1, :V2]
    t1 = Branch(atom(fn, 1, <, 2.0), leaf(1),
        Branch(atom(fn, 2, ≥, 1.0), leaf(2), leaf(3)))
    t2 = Branch(atom(fn, 2, <, 3.0),
        Branch(atom(fn, 1, ≤, 4.0), leaf(1), leaf(2)), leaf(3))
    t3 = Branch(atom(fn, 1, >, 1.0), leaf(2), leaf(1))
    return forest([t1, t2, t3], 2)
end

# a random tree whose root is always a split; thresholds are drawn from a
# small pool so that trees share cut points, as real forests do
function randtree(rng, nfeat, depth, pool; root=true)
    if !root && (depth == 0 || rand(rng) < 0.2)
        return leaf(rand(rng, 1:length(CLS)))
    end
    fn = [Symbol("V", i) for i in 1:nfeat]
    f = rand(rng, 1:nfeat)
    op = rand(rng, (<, ≤, >, ≥))
    return Branch(atom(fn, f, op, rand(rng, pool)),
        randtree(rng, nfeat, depth - 1, pool; root=false),
        randtree(rng, nfeat, depth - 1, pool; root=false))
end

random_forest(rng, ntrees, nfeat, depth, pool) =
    forest([randtree(rng, nfeat, depth, pool) for _ in 1:ntrees], nfeat)

# the thresholds of a model's splits, per feature
function model_thresholds(model, nfeat)
    thrs = [Float64[] for _ in 1:nfeat]
    walk(n) = if n isa Branch
        c = SoleLogics.value(SM.antecedent(n))
        push!(thrs[SoleData.i_variable(SoleData.feature(c))], SoleData.threshold(c))
        walk(SM.posconsequent(n)); walk(SM.negconsequent(n))
    end
    foreach(walk, SM.models(model))
    return thrs
end

# probe points: every threshold, its float neighbours, and values far off,
# combined at random so that every region and every boundary gets hit
function probe_points(rng, model, nfeat, n)
    cand = map(model_thresholds(model, nfeat)) do ts
        isempty(ts) && return [0.0]
        unique!(vcat(ts, prevfloat.(ts), nextfloat.(ts),
            minimum(ts) - 1, maximum(ts) + 1, (ts .+ 0.25)))
    end
    X = DataFrame([Symbol("V", j) => Float64[] for j in 1:nfeat])
    for _ in 1:n
        push!(X, [rand(rng, c) for c in cand])
    end
    return X
end

# class of row `i` under a decision set: the first rule that fires
function ds_predict(ds, L, i)
    for r in SM.rules(ds)
        SoleLogics.check(SM.antecedent(r), L, i) &&
            return string(SM.outcome(SM.consequent(r)))
    end
    return nothing
end

# number of rows on which the decision set disagrees with the model
function disagreements(ds, model, X)
    L = scalarlogiset(X; allow_propositional=true)
    pm = string.(SM.apply(model, L; suppress_parity_warning=true))
    pd = [ds_predict(ds, L, i) for i in 1:size(X, 1)]
    return count(pm .!= something.(pd, "none"))
end

# does a set of terms cover point `x` (one value per feature)?
covers(term, x) = all(a -> LU.evalop(a.op, x[a.feat], a.thr), term)
covers(d::LU.LumenDNF, x) = any(t -> covers(t, x), d)

atomset(term) = Set(collect(term))
termsets(d) = Set(atomset(t) for t in d)

# the PLA text `write_pla` produces
plastring(L, cubes) = sprint(io -> LU.write_pla(io, L, cubes))

config64(; kw...) = LumenConfig(; float_resolution=Float64, kw...)

# ---------------------------------------------------------------------------- #
#                                    tests                                     #
# ---------------------------------------------------------------------------- #
@testset "Lumen Shannon" begin

@testset "operators" begin
    for (op, sym) in ((<, '<'), (>, '>'), (≥, '≥'), (≤, '≤'))
        b = LU.evalop(op)
        @test LU.evalop(b) == sym
        @test LU._opfun(b) === op
        @test LU.dualop(LU.dualop(b)) == b
        for x in (0.5, 1.0, 1.5)
            @test LU.evalop(b, x, 1.0) == op(x, 1.0)
            # an operator and its dual split ℝ in two
            @test LU.evalop(b, x, 1.0) != LU.evalop(LU.dualop(b), x, 1.0)
        end
    end
    @test LU.isupper(LU.evalop(<)) && LU.isupper(LU.evalop(≤))
    @test !LU.isupper(LU.evalop(>)) && !LU.isupper(LU.evalop(≥))
    @test LU.upperop(LU.evalop(≥)) == LU.evalop(<)
    @test LU.upperop(LU.evalop(>)) == LU.evalop(≤)
    @test LU.upperop(LU.evalop(<)) == LU.evalop(<)
end

@testset "LumenAtom" begin
    a = A(2, 1.5, <)
    @test LU.dual(a) == A(2, 1.5, ≥)
    @test LU.dual(LU.dual(a)) == a
    @test hash(a) == hash(A(2, 1.5, <))
    @test a != A(2, 1.5, ≤)

    e = LU.LumenAtom{R,T}()
    @test LU.isempty_atom(e)
    @test !LU.isempty_atom(a)
    @test LU.isempty_atom(LU.LumenAtom())
    @test LU.LumenAtom() isa LU.LumenAtom{UInt32,Float64}

    atoms = [A(1, 1.0, <), A(2, 2.0, ≥), A(1, 3.0, ≤)]
    @test LU.get_features(atoms) == R[1, 2, 1]
    @test LU.get_thresholds(atoms) == [1.0, 2.0, 3.0]
    @test LU.get_thresholds(atoms, R(1)) == [1.0, 3.0]
    @test LU.get_op(atoms) == LU.evalop.([<, ≥, ≤])
    @test LU.featidxs(atoms, R(1)) == [1, 3]
    @test sort(atoms) == [A(1, 1.0, <), A(1, 3.0, ≤), A(2, 2.0, ≥)]
end

@testset "LumenEnsemble" begin
    model = small_forest()
    cfg = config64()
    e = LumenEnsemble(cfg, model)

    @test length(e.roots) == 3
    @test length(e) == 13
    @test count(LU.isleaf, e.nodes) == 8
    @test e.roots == R[1, 6, 11]
    # children always come after their parent
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && continue
        @test n.left > i && n.right > i
        @test !LU.isempty_atom(e.atoms[i])
    end
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && @test LU.isempty_atom(e.atoms[i])
    end
    @test e.atoms[1] == A(1, 2.0, <)
    @test LU._label.(e.classes) == ["a", "b", "c"]

    atoms = LU.get_atoms(e)
    @test length(atoms) == 5
    @test !any(LU.isempty_atom, atoms)

    # BFS order, ceil(depth * n) atoms per tree
    @test Set(LU.bfs_atoms(e, 1.0)) == Set(atoms)
    @test LU.bfs_atoms(e, 0.5) == [A(1, 2.0, <), A(2, 3.0, <), A(1, 1.0, >)]

    @test occursin("3 trees", sprint(show, e))

    # the compiled forest votes exactly like the Sole model, ties included
    rng = Xoshiro(7)
    X = probe_points(rng, model, 2, 200)
    tbl = Matrix{T}(X)
    preds = LU.apply(e, view(tbl, :, :), 3)
    L = scalarlogiset(X; allow_propositional=true)
    ref = SM.apply(model, L; suppress_parity_warning=true)
    @test CLS[preds] == ref

    buf = zeros(R, size(tbl, 1) + 2)
    @test LU.apply!(buf, e, view(tbl, :, :), 3) === buf
    @test buf[1:end-2] == preds

    @test_throws MethodError LumenEnsemble(cfg, leaf(1))
end

@testset "ThresholdSpace" begin
    e = LumenEnsemble(config64(), small_forest())
    ts = ThresholdSpace(e, R[1, 2], R[1, 2, 3], 1.0)

    # '>' and '≤' are normalised to '<' at nextfloat(thr); thresholds are
    # stored descending, then the boundary value below the smallest one
    @test ts.nlev == R[4, 3]
    @test ts.thrs_offset == R[1, 5]
    @test ts.thrs == [nextfloat(4.0), 2.0, nextfloat(1.0), 1.0,
                      3.0, 1.0, prevfloat(1.0)]
    @test all(==(LU.evalop(<)), ts.op_families)
    @test occursin("2 features", sprint(show, ts))

    # duplicate cut points collapse: x ≥ 2 and x < 2 are the same cut
    fn = [:V1]
    m = forest([Branch(atom(fn, 1, <, 2.0), leaf(1), leaf(2)),
                Branch(atom(fn, 1, ≥, 2.0), leaf(2), leaf(1))], 1)
    ts1 = ThresholdSpace(LumenEnsemble(config64(), m), R[1], R[1, 2, 3], 1.0)
    @test ts1.nlev == R[2]

    # a feature no tree reads gets one level and a NaN sample
    ts3 = ThresholdSpace(e, R[1, 2, 3], R[1, 2, 3], 1.0)
    @test ts3.nlev == R[4, 3, 1]
    @test isnan(ts3.thrs[end])

    # depth < 1 keeps only the BFS-first atoms of each tree
    tsd = ThresholdSpace(e, R[1, 2], R[1, 2, 3], 0.5)
    @test tsd.nlev == R[3, 2]

    @test LU._thrs_boundary([3.0, 1.0], LU.evalop(<)) == prevfloat(1.0)
    @test LU._thrs_boundary([1.0, 3.0], LU.evalop(>)) == nextfloat(3.0)
    @test LU._thrs_boundary([[3.0, 1.0], [1.0, 2.0]], LU.evalop.([<, ≥])) ==
        [prevfloat(1.0), nextfloat(2.0)]
end

@testset "gather_atoms" begin
    e = LumenEnsemble(config64(), small_forest())
    ts = ThresholdSpace(e, R[1, 2], R[1, 2, 3], 1.0)
    out = Vector{LU.LumenAtom{R,T}}(undef, 4)

    # first level: no upper bound
    n = LU.gather_atoms(view(out, :), 0, ts, R[1, 1])
    @test n == 2
    @test out[1:2] == [A(1, nextfloat(4.0), ≥), A(2, 3.0, ≥)]

    # inner level on V1, last level on V2: no lower bound on V2
    n = LU.gather_atoms(view(out, :), 0, ts, R[2, 3])
    @test n == 3
    @test out[1:3] == [A(1, nextfloat(4.0), <), A(1, 2.0, ≥), A(2, 1.0, <)]

    # the cube of a grid point contains that point's sample value
    for l1 in 1:4, l2 in 1:3
        n = LU.gather_atoms(view(out, :), 0, ts, R[l1, l2])
        x = [ts.thrs[ts.thrs_offset[1]+l1-1], ts.thrs[ts.thrs_offset[2]+l2-1]]
        @test covers(view(out, 1:n), x)
    end

    @test LU.removeduals!([A(1, 1.0, <), A(1, 1.0, ≥), A(1, 2.0, ≥)]) ==
        [A(1, 1.0, <), A(1, 2.0, ≥)]
    out2 = LU.LumenAtom{R,T}[]
    pool = [A(1, 2.0, ≥), A(1, 1.0, <), A(1, 2.0, <), A(1, 1.0, <)]
    @test LU._alphabet!(out2, pool) == [A(1, 1.0, <), A(1, 2.0, <)]
    @test length(pool) == 4
end

@testset "LumenDNF" begin
    d = LU.LumenDNF{R,T}()
    @test isempty(d) && length(d) == 0 && LU.natoms(d) == 0

    push!(d, [A(1, 1.0, <), A(2, 2.0, ≥)])
    push!(d, [A(1, 3.0, ≥)])
    @test LU.nterms(d) == 2 && LU.natoms(d) == 3
    @test collect(d[1]) == [A(1, 1.0, <), A(2, 2.0, ≥)]
    @test d[2] isa LU.LumenCube{R,T}
    @test collect.(collect(d)) == [collect(d[1]), collect(d[2])]
    @test eachindex(d) == 1:2 && lastindex(d) == 2

    o = push!(LU.LumenDNF{R,T}(), [A(2, 5.0, <)])
    append!(d, o)
    @test LU.nterms(d) == 3
    @test collect(d[3]) == [A(2, 5.0, <)]

    empty!(d)
    @test isempty(d) && d.offsets == [0]
    @test occursin("0 terms", sprint(show, d))
end

@testset "PLA round trip" begin
    atoms = [A(1, 1.0, <), A(1, 2.0, <), A(1, 3.0, <), A(2, 5.0, ≥)]
    L = LU.PlaLayout(atoms)
    @test LU.ncolumns(L) == 4
    @test LU.ndontcares(L) == 2
    # one feature's columns are contiguous, thresholds descending
    @test L.thr == [3.0, 2.0, 1.0, 5.0]
    @test L.bfirst == [1, 1, 1, 4] && L.blast == [3, 3, 3, 4]
    @test L.upop == LU.evalop.([<, <, <, <])
    # dual atoms share one column
    @test LU.ncolumns(LU.PlaLayout([A(1, 1.0, <), A(1, 1.0, ≥)])) == 1

    cubes = LU.LumenDNF{R,T}()
    push!(cubes, [A(1, 3.0, <), A(1, 2.0, ≥), A(2, 5.0, <)])
    push!(cubes, [A(1, 1.0, ≥)])
    push!(cubes, [A(1, 1.0, <), A(2, 5.0, ≥)])

    pla = plastring(L, cubes)
    @test occursin(".type fd", pla)
    @test occursin(".p 5", pla)  # 3 cubes + 2 don't-care rows
    lines = split(pla, '\n')
    # cube 1: '1' up to its upper bound (3), '0' from its lower bound (2) on
    @test "1001 1" in lines
    @test "--0- 1" in lines   # cube 2: lower bound only
    # one don't-care row per adjacent column pair of V1
    @test "01-- -" in lines && "-01- -" in lines

    back = LU.read_pla(pla, L)
    @test termsets(back) == termsets(cubes)
    @test termsets(LU.read_pla(codeunits(pla), L)) == termsets(cubes)

    # ABC-style permuted inputs and comment lines
    L2 = LU.PlaLayout([A(1, 2.0, <), A(1, 1.0, <)])
    perm = ".i 2\n.o 1\n# comment\n.ilb x2 x1\n.ob f\n.p 1\n01 1\n.e\n"
    @test termsets(LU.read_pla(perm, L2)) ==
        Set([Set([A(1, 2.0, <), A(1, 1.0, ≥)])])
    # without the permutation the same row is an empty region and is dropped
    plain = ".i 2\n.o 1\n.ilb x1 x2\n.ob f\n.p 1\n01 1\n.e\n"
    @test isempty(LU.read_pla(plain, L2))
    # output '0' rows are not terms; an all '-' row is the tautology
    @test isempty(LU.read_pla(".i 2\n.o 1\n10 0\n.e\n", L2))
    taut = LU.read_pla(".i 2\n.o 1\n-- 1\n.e\n", L2)
    @test LU.nterms(taut) == 1 && isempty(taut[1])

    @test_throws ArgumentError LU.read_pla(".ilb x1 x9\n", L2)
    @test_throws ArgumentError LU.read_pla(".ilb y1 x2\n", L2)
    @test_throws ArgumentError LU.read_pla(".i 2\n101 1\n", L2)
end

@testset "refine!" begin
    @test LU.dominates([A(1, 2.0, <)], [A(1, 1.0, <)])
    @test !LU.dominates([A(1, 1.0, <)], [A(1, 2.0, <)])
    @test LU.dominates(LU.LumenAtom{R,T}[], [A(1, 1.0, <)])
    @test !LU.dominates([A(2, 1.0, <)], [A(1, 1.0, <)])
    @test LU.dominates([A(1, 1.0, ≥)], [A(1, 2.0, ≥), A(2, 0.0, <)])

    u, l = LU._sidemask([A(1, 1.0, <), A(3, 1.0, ≥)])
    @test u == 0x1 && l == 0x4

    # duplicates and contained terms go
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 2.0, <)])
    push!(d, [A(1, 1.0, <)])
    push!(d, [A(1, 2.0, <)])
    push!(d, [A(2, 0.0, ≥)])
    push!(d, [A(2, 1.0, ≥), A(1, 0.5, <)])
    LU.refine!(d)
    @test termsets(d) == Set([Set([A(1, 2.0, <)]), Set([A(2, 0.0, ≥)])])
    @test LU.nterms(d) == 2

    # a contained term placed BEFORE its container, with the same number of
    # atoms, survives: dominance is only checked forward in size order
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 1.0, <)])
    push!(d, [A(1, 2.0, <)])
    LU.refine!(d)
    @test_broken LU.nterms(d) == 1
end

@testset "compact!" begin
    # siblings differing only in the split feature's bound merge into one
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 2.0, <), A(2, 1.0, <)])
    push!(d, [A(1, 2.0, ≥), A(1, 4.0, <), A(2, 1.0, <)])
    LU.compact!(d)
    @test termsets(d) == Set([Set([A(1, 4.0, <), A(2, 1.0, <)])])

    # a chain of three adjacent boxes collapses to one
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 1.0, <)])
    push!(d, [A(1, 1.0, ≥), A(1, 2.0, <)])
    push!(d, [A(1, 2.0, ≥), A(1, 3.0, <)])
    LU.compact!(d)
    @test termsets(d) == Set([Set([A(1, 3.0, <)])])

    # disjoint boxes stay apart
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 1.0, <)])
    push!(d, [A(1, 2.0, ≥)])
    LU.compact!(d)
    @test LU.nterms(d) == 2

    # covered points never change
    rng = Xoshiro(3)
    for _ in 1:200
        d = LU.LumenDNF{R,T}()
        for _ in 1:rand(rng, 2:8)
            term = LU.LumenAtom{R,T}[]
            for f in 1:2
                lo, hi = sort(rand(rng, 0:4, 2))
                lo > 0 && push!(term, A(f, lo, ≥))
                hi < 4 && hi > lo && push!(term, A(f, hi, <))
            end
            push!(d, term)
        end
        pts = [(x, y) for x in -0.5:0.5:4.5, y in -0.5:0.5:4.5]
        before = [covers(d, p) for p in pts]
        # a union covering everything is the tautology case, broken below
        all(before) && continue
        LU.compact!(d)
        @test [covers(d, p) for p in pts] == before
    end

    # two halves of the axis make the tautology, a single empty term
    d = LU.LumenDNF{R,T}()
    push!(d, [A(1, 2.0, <)])
    push!(d, [A(1, 2.0, ≥)])
    LU.compact!(d)
    @test_broken LU.nterms(d) == 1 && isempty(d[1])
end

@testset "minimizer settings" begin
    @test LU.evalst(Fast) == 0x01
    @test LU.evalst(Balanced) == 0x02
    @test LU.evalst(Sop) == 0x03
    @test LU.evalst(Fast()) == 0x01
    @test LU._espresso_flags(0x01) == String[]
    @test LU._espresso_flags(0x02) == ["-estrong"]
    @test LU._espresso_flags(0x03) == ["-Dexact"]
    for s in (0x01, 0x02, 0x03)
        sc = LU._abc_script(s, "in.pla", "out.pla")
        @test startswith(sc, "read_pla -d in.pla;")
        @test endswith(sc, "write_pla out.pla")
    end

    m = LU.Minimizer(Abc, config64())
    @test isdir(m.dir)
    @test startswith(LU._tempfile(m), m.dir)
    close(m)
    @test !isdir(m.dir)
end

@testset "run_minimization (ABC)" begin
    cfg = config64()
    alphabet = [A(1, 1.0, <), A(1, 2.0, <), A(1, 3.0, <), A(2, 5.0, <)]

    # zero and one cube never reach the backend
    @test isempty(LU.run_minimization(Abc, cfg, alphabet, LU.LumenDNF{R,T}()))
    single = push!(LU.LumenDNF{R,T}(), [A(1, 1.0, <)])
    @test termsets(LU.run_minimization(Abc, cfg, alphabet, single)) == termsets(single)

    # the four regions of V1 below 3, all with V2 < 5: one term
    cubes = LU.LumenDNF{R,T}()
    push!(cubes, [A(1, 1.0, <), A(2, 5.0, <)])
    push!(cubes, [A(1, 1.0, ≥), A(1, 2.0, <), A(2, 5.0, <)])
    push!(cubes, [A(1, 2.0, ≥), A(1, 3.0, <), A(2, 5.0, <)])
    for setup in (Fast, Balanced, Sop)
        c = config64(; p_minimization_setup=setup)
        d = LU.run_minimization(Abc, c, alphabet, cubes)
        @test termsets(d) == Set([Set([A(1, 3.0, <), A(2, 5.0, <)])])
    end

    # coverage is exactly preserved on random problems
    rng = Xoshiro(11)
    thr = [1.0, 2.0, 3.0]
    alpha = [A(f, t, <) for f in 1:2 for t in thr]
    for _ in 1:10
        cubes = LU.LumenDNF{R,T}()
        for l1 in 1:4, l2 in 1:4
            rand(rng) < 0.5 || continue
            term = LU.LumenAtom{R,T}[]
            for (f, l) in ((1, l1), (2, l2))
                l > 1 && push!(term, A(f, thr[l-1], ≥))
                l < 4 && push!(term, A(f, thr[l], <))
            end
            push!(cubes, term)
        end
        d = LU.run_minimization(Abc, cfg, alpha, cubes)
        pts = [(x, y) for x in 0.5:0.5:3.5, y in 0.5:0.5:3.5]
        @test [covers(d, p) for p in pts] == [covers(cubes, p) for p in pts]
        @test LU.nterms(d) ≤ max(LU.nterms(cubes), 1)
    end
end

@testset "_extract: leaf size does not change coverage" begin
    model = small_forest()
    e = LumenEnsemble(config64(), model)
    ts = ThresholdSpace(e, R[1, 2], R[1, 2, 3], 1.0)
    lo = ones(R, 2)
    hi = copy(ts.nlev)
    pts = [(x, y) for x in -1:0.25:6, y in -1:0.25:5]

    ref = LU._extract(config64(; M=4096), ts, e, lo, hi)
    @test lo == ones(R, 2) && hi == ts.nlev  # inputs are never mutated
    for M in (1, 2, 3, 5, 12)
        out = LU._extract(config64(; M), ts, e, lo, hi)
        for c in 1:3
            @test [covers(out[c], p) for p in pts] == [covers(ref[c], p) for p in pts]
        end
    end
    # every grid point is covered by exactly one class
    for p in pts
        @test count(c -> covers(ref[c], p), 1:3) == 1
    end
end

@testset "lumen_shannon on a small forest" begin
    model = small_forest()
    rng = Xoshiro(1)
    X = probe_points(rng, model, 2, 400)

    for M in (1, 4, 4096), setup in (Fast, Balanced, Sop)
        cfg = config64(; M, p_minimization_setup=setup)
        ds = lumen_shannon(cfg, model)
        @test ds isa SM.DecisionSet
        @test disagreements(ds, model, X) == 0
    end

    # "c" never wins a vote, so it gets no rule; rules follow
    # `supporting_labels` order
    ds = lumen_shannon(config64(), model)
    @test [string(SM.outcome(SM.consequent(r))) for r in SM.rules(ds)] == ["a", "b"]

    # the model's info is left alone
    labels = copy(SM.info(model, :supporting_labels))
    fnames = copy(SM.info(model, :featurenames))
    m2 = DecisionEnsemble(SM.models(model),
        (; featurenames=[:V1, :V2, :V1], supporting_labels=["b", "a", "b", "c"]))
    lumen_shannon(config64(), model)
    @test SM.info(model, :supporting_labels) == labels
    @test SM.info(model, :featurenames) == fnames
    lumen_shannon(config64(), m2)
    @test_broken SM.info(m2, :supporting_labels) == ["b", "a", "b", "c"]
end

@testset "lumen_shannon on random forests" begin
    rng = Xoshiro(2024)
    pool = collect(0.5:0.5:4.0)
    for (ntrees, nfeat, depth, M) in (
        (1, 2, 3, 4096), (3, 3, 3, 16), (4, 3, 3, 8),  # 4 trees: many ties
        (5, 4, 4, 64), (7, 2, 4, 2),
    )
        model = random_forest(rng, ntrees, nfeat, depth, pool)
        X = probe_points(rng, model, nfeat, 300)
        ds = lumen_shannon(config64(; M), model)
        @test disagreements(ds, model, X) == 0
    end
end

@testset "lumen_shannon depth < 1" begin
    model = small_forest()
    ds = lumen_shannon(config64(; depth=0.5), model)
    @test ds isa SM.DecisionSet
    @test !isempty(SM.rules(ds))
end

@testset "assign and decision_set" begin
    labels, codes = LU.assign(R, ["b", "a", "b"])
    @test labels == ["b", "a"] && codes == R[1, 2, 1]
    labels, codes = LU.assign(R, ["b", "a", "b"], true)
    @test labels == ["a", "b"] && codes == R[2, 1, 2]
    @test labels[codes] == ["b", "a", "b"]

    d1 = push!(LU.LumenDNF{R,T}(), [A(1, 2.0, <)])
    push!(d1, [A(1, 4.0, ≥), A(2, 1.0, <)])
    d2 = LU.LumenDNF{R,T}()
    d3 = push!(LU.LumenDNF{R,T}(), LU.LumenAtom{R,T}[])
    ds = LU.decision_set([d1, d2, d3], [:x, :y], ["p", "q", "r"], [3, 1, 2])
    rs = SM.rules(ds)
    @test length(rs) == 2  # empty `d2` gets no rule
    @test SM.outcome(SM.consequent(rs[1])) == "r"
    @test SM.antecedent(rs[1]) == SoleLogics.⊤
    @test SM.outcome(SM.consequent(rs[2])) == "p"
    @test SM.antecedent(rs[2]) isa SoleLogics.LeftmostDisjunctiveForm

    a = LU._toatom(A(2, 1.5, ≥), [:x, :y])
    c = SoleLogics.value(a)
    @test SoleData.threshold(c) == 1.5
    @test SoleData.test_operator(c) === (≥)
    @test SoleData.i_variable(SoleData.feature(c)) == 2

    m = small_forest()
    @test LU._ruleorder(m, ["a", "b", "c"]) == [1, 2, 3]
    m2 = DecisionEnsemble(SM.models(m), (; supporting_labels=["c", "a", "c"]))
    @test LU._ruleorder(m2, ["a", "b", "c"]) == [3, 1, 2]
    @test LU._featurenames(m, LumenEnsemble(config64(), m)) == [:V1, :V2]
    m3 = DecisionEnsemble(SM.models(m), (; supporting_labels=["a"]))
    @test LU._featurenames(m3, LumenEnsemble(config64(), m3)) == [:V1, :V2]
end

@testset "known issues" begin
    model = small_forest()

    # the MIT Espresso backend reads `config.minimization_setup`, a field
    # `LumenConfig` does not have
    cfg = LumenConfig(; minimization_scheme=:mitespresso,
        p_minimization_scheme=MitEspresso, float_resolution=Float64)
    @test_broken try
        lumen_shannon(cfg, model) isa SM.DecisionSet
    catch
        false
    end

    # the one-off `_leaf_extract` passes its arguments in the wrong order
    e = LumenEnsemble(config64(), model)
    ts = ThresholdSpace(e, R[1, 2], R[1, 2, 3], 1.0)
    @test_broken try
        LU._leaf_extract(config64(), ts, e, ones(R, 2), copy(ts.nlev)) isa Tuple
    catch
        false
    end

    # with the default Float32 thresholds, `x ≤ 4.0` becomes
    # `x < nextfloat(4f0)`, which also admits Float64 values just above 4
    ds32 = lumen_shannon(LumenConfig(), model)
    X = DataFrame(V1=[nextfloat(4.0)], V2=[0.0])
    @test_broken disagreements(ds32, model, X) == 0

    # a class predicted everywhere, split into two leaves (`M = 1`), loses
    # its rule: `compact!` merges the halves into an empty box and drops it
    fn = [:V1]
    m1 = forest([Branch(atom(fn, 1, <, 2.0), leaf(1), leaf(1))], 1)
    ds = lumen_shannon(config64(; M=1), m1)
    @test_broken length(SM.rules(ds)) == 1
end

end
