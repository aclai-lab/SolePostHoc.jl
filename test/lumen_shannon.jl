using Test
using SolePostHoc
const LU = SolePostHoc.Lumen

using SoleModels
using ModalDecisionLists

using CategoricalArrays
using DecisionTree
using RDatasets
using StableRNGs
using XGBoost

# ---------------------------------------------------------------------------------------- #
#                                      configuration                                       #
# ---------------------------------------------------------------------------------------- #
rng = StableRNG(42)
const R = UInt32
const T = Float32

iris_df = dataset("datasets", "iris")
y = iris_df.Species
X = T.(select(iris_df, Not(:Species)))

# ---------------------------------------------------------------------------------------- #
#                                         utils                                            #
# ---------------------------------------------------------------------------------------- #
lumenatom(f, thr, op) = LU.LumenAtom{R,T}(R(f), T(thr), LU.evalop(op))

function decision_tree(
    X::DataFrame,
    y::CategoricalArray,
    max_depth::Int,
    rng
)
    n_subfeatures = 0

    model = DecisionTree.build_tree(
        y, Matrix(X),
        n_subfeatures,
        max_depth;
        rng
    )

    return SoleModels.solemodel(DecisionTree.wrap(model.node, (featurenames=names(X),)))
end

function random_forest(
    X::DataFrame,
    y::CategoricalArray,
    n_trees::Int,
    max_depth::Int,
    rng
)
    n_subfeatures = -1
    partial_sampling = 0.7

    model = DecisionTree.build_forest(
        y, Matrix(X),
        n_subfeatures,
        n_trees,
        partial_sampling,
        max_depth;
        rng
    )

    return SoleModels.solemodel(model; featurenames=names(X))
end

function xg_boost(
    X::DataFrame,
    y::CategoricalArray,
    num_round::Int,
    max_depth::Int
)
    yl = levelcode.(y) .- 1 # CategoricalArrays: 1-based codes → 0-based
    model = xgboost((Matrix(X), yl);
        num_round,
        max_depth,
        objective = "multi:softmax",
        num_class = length(levels(y)),
    )

    return SoleModels.solemodel(
        XGBoost.trees(model),
        Matrix(X), y;
        featurenames = names(X),
        classlabels  = levels(y),
    )
end

function decision_list(
    X::DataFrame,
    y::CategoricalArray,
    num_models::Int,
    max_rulebase_length::Int, # 3
    max_rule_length::Int, # 3
    rng
)
    logiset = PropositionalLogiset(X)

    model = ModalDecisionLists.build_ensemble(
        logiset,
        y,
        num_models;
        model_wrapper=ModalDecisionLists.sequentialcovering,
        max_rulebase_length,
        max_rule_length,
        loss_function=ModalDecisionLists.LossFunctions.LaplaceMetric(),
        featurenames=propertynames(X),
        min_rule_coverage=max(1, round(Int, size(X, 1) * 0.05)),
        discretizedomain=false,
        rng
    )

    return solemodel(model; featurenames=names(X))
end

# ---------------------------------------------------------------------------------------- #
#                                  operators and atoms                                     #
# ---------------------------------------------------------------------------------------- #
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
    a = lumenatom(2, 1.5, <)
    @test LU.dual(a) == lumenatom(2, 1.5, ≥)
    @test LU.dual(LU.dual(a)) == a
    @test hash(a) == hash(lumenatom(2, 1.5, <))
    @test a != lumenatom(2, 1.5, ≤)

    e = LU.LumenAtom{R,T}()
    @test LU.isempty_atom(e)
    @test !LU.isempty_atom(a)
    @test LU.isempty_atom(LU.LumenAtom())
    @test LU.LumenAtom() isa LU.LumenAtom{UInt32,Float64}

    atoms = [lumenatom(1, 1.0, <), lumenatom(2, 2.0, ≥), lumenatom(1, 3.0, ≤)]
    @test LU.get_features(atoms) == R[1, 2, 1]
    @test LU.get_thresholds(atoms) == [1.0, 2.0, 3.0]
    @test LU.get_thresholds(atoms, R(1)) == [1.0, 3.0]
    @test LU.get_op(atoms) == LU.evalop.([<, ≥, ≤])
    @test LU.featidxs(atoms, R(1)) == [1, 3]
    @test sort(atoms) == [lumenatom(1, 1.0, <), lumenatom(1, 3.0, ≤), lumenatom(2, 2.0, ≥)]
end

# ---------------------------------------------------------------------------------------- #
#                                        ensembles                                         #
# ---------------------------------------------------------------------------------------- #
@testset "LumenEnsemble: Decision Tree" begin
    model = decision_tree(X, y, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    @test length(e.roots) == 1
    # children always come after their parent
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && continue
        @test n.left > i && n.right > i
        @test !LU.isempty_atom(e.atoms[i])
    end
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && @test LU.isempty_atom(e.atoms[i])
    end

    atoms = @test_nowarn LU.get_atoms(e)
    @test !any(LU.isempty_atom, atoms)
end

@testset "LumenEnsemble: Random Forest" begin
    model = random_forest(X, y, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    @test length(e.roots) == 3
    # children always come after their parent
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && continue
        @test n.left > i && n.right > i
        @test !LU.isempty_atom(e.atoms[i])
    end
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && @test LU.isempty_atom(e.atoms[i])
    end

    atoms = @test_nowarn LU.get_atoms(e)
    @test !any(LU.isempty_atom, atoms)
end

@testset "LumenEnsemble: XGBoost" begin
    model = xg_boost(X, y, 3, 3)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    @test length(e.roots) == 9
    # children always come after their parent
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && continue
        @test n.left > i && n.right > i
        @test !LU.isempty_atom(e.atoms[i])
    end
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && @test LU.isempty_atom(e.atoms[i])
    end

    atoms = @test_nowarn LU.get_atoms(e)
    @test !any(LU.isempty_atom, atoms)
end

@testset "LumenEnsemble: Decision List" begin
    model = decision_list(X, y, 3, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    @test length(e.roots) == 3
    # children always come after their parent
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && continue
        @test n.left > i && n.right > i
        @test !LU.isempty_atom(e.atoms[i])
    end
    for (i, n) in enumerate(e.nodes)
        LU.isleaf(n) && @test LU.isempty_atom(e.atoms[i])
    end

    atoms = @test_nowarn LU.get_atoms(e)
    @test !any(LU.isempty_atom, atoms)
end

@testset "ThresholdSpace" begin
    feat_idxs, class_idxs = R[1, 2, 3, 4], R[1, 2, 3]

    model = decision_tree(X, y, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    @test_nowarn thrs.nlev
    @test_nowarn thrs.thrs_offset
    @test_nowarn thrs.thrs
    @test all(==(LU.evalop(<)), thrs.op_families)

    model = random_forest(X, y, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    @test_nowarn thrs.nlev
    @test_nowarn thrs.thrs_offset
    @test_nowarn thrs.thrs
    @test all(==(LU.evalop(<)), thrs.op_families)

    model = xg_boost(X, y, 3, 3)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    @test_nowarn thrs.nlev
    @test_nowarn thrs.thrs_offset
    @test_nowarn thrs.thrs
    @test all(==(LU.evalop(<)), thrs.op_families)

    model = decision_list(X, y, 3, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)

    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    @test_nowarn thrs.nlev
    @test_nowarn thrs.thrs_offset
    @test_nowarn thrs.thrs
    @test all(==(LU.evalop(<)), thrs.op_families)
end

@testset "Extraction" begin
    feat_idxs, class_idxs = R[1, 2, 3, 4], R[1, 2, 3]

    model = decision_tree(X, y, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))

    lumendnf = LU._extract(cfg, thrs, e, lo, hi)
    @test lumendnf isa Vector{<:LumenDNF}
    @test length(lumendnf) == 3

    model = random_forest(X, y, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))

    lumendnf = LU._extract(cfg, thrs, e, lo, hi)
    @test lumendnf isa Vector{<:LumenDNF}
    @test length(lumendnf) == 3

    model = xg_boost(X, y, 3, 3)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))

    lumendnf = LU._extract(cfg, thrs, e, lo, hi)
    @test lumendnf isa Vector{<:LumenDNF}
    @test length(lumendnf) == 3

    model = decision_list(X, y, 3, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))

    lumendnf = LU._extract(cfg, thrs, e, lo, hi)
    @test lumendnf isa Vector{<:LumenDNF}
    @test length(lumendnf) == 3
end

@testset "LumenDNF" begin
    d = LU.LumenDNF{R,T}()
    @test isempty(d) && length(d) == 0 && LU.natoms(d) == 0

    push!(d, [lumenatom(1, 1.0, <), lumenatom(2, 2.0, ≥)])
    push!(d, [lumenatom(1, 3.0, ≥)])
    @test LU.nterms(d) == 2 && LU.natoms(d) == 3
    @test collect(d[1]) == [lumenatom(1, 1.0, <), lumenatom(2, 2.0, ≥)]
    @test d[2] isa LU.LumenCube{R,T}
    @test collect.(collect(d)) == [collect(d[1]), collect(d[2])]
    @test eachindex(d) == 1:2 && lastindex(d) == 2

    o = push!(LU.LumenDNF{R,T}(), [lumenatom(2, 5.0, <)])
    append!(d, o)
    @test LU.nterms(d) == 3
    @test collect(d[3]) == [lumenatom(2, 5.0, <)]

    empty!(d)
    @test isempty(d) && d.offsets == [0]
    @test occursin("0 terms", sprint(show, d))
end

@testset "PLA round trip" begin
    # the PLA text `write_pla` produces
    plastring(L, cubes) = sprint(io -> LU.write_pla(io, L, cubes))
    termsets(d) = Set(collect(t) for t in d)

    atoms = [
        lumenatom(1, 1.0, <),
        lumenatom(1, 2.0, <),
        lumenatom(1, 3.0, <),
        lumenatom(2, 5.0, ≥)
    ]
    L = LU.PlaLayout(atoms)
    @test LU.ncolumns(L) == 4
    @test LU.ndontcares(L) == 2
    # one feature's columns are contiguous, thresholds descending
    @test L.thr == [3.0, 2.0, 1.0, 5.0]
    @test L.bfirst == [1, 1, 1, 4] && L.blast == [3, 3, 3, 4]
    @test L.upop == LU.evalop.([<, <, <, <])
    # dual atoms share one column
    @test LU.ncolumns(LU.PlaLayout([lumenatom(1, 1.0, <), lumenatom(1, 1.0, ≥)])) == 1

    cubes = LU.LumenDNF{R,T}()
    push!(cubes, [lumenatom(1, 3.0, <), lumenatom(1, 2.0, ≥), lumenatom(2, 5.0, <)])
    push!(cubes, [lumenatom(1, 1.0, ≥)])
    push!(cubes, [lumenatom(1, 1.0, <), lumenatom(2, 5.0, ≥)])

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

@testset "compact!" begin
    feat_idxs, class_idxs = R[1, 2, 3, 4], R[1, 2, 3]

    model = decision_tree(X, y, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))
    terms = LU._extract(cfg, thrs, e, lo, hi)

    @test_nowarn foreach(compact!, terms)

    model = random_forest(X, y, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))
    terms = LU._extract(cfg, thrs, e, lo, hi)

    @test_nowarn foreach(compact!, terms)

    model = xg_boost(X, y, 3, 3)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))
    terms = LU._extract(cfg, thrs, e, lo, hi)

    @test_nowarn foreach(compact!, terms)

    model = decision_list(X, y, 3, 3, 3, rng)
    cfg = LumenConfig(; float_resolution=T)
    e = LumenEnsemble(cfg, model)
    thrs = ThresholdSpace(e, feat_idxs, class_idxs, cfg.depth)
    hi = thrs.nlev
    lo = ones(R, length(hi))
    terms = LU._extract(cfg, thrs, e, lo, hi)

    @test_nowarn foreach(compact!, terms)
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
end

# ---------------------------------------------------------------------------------------- #
#                                      lumen shannon                                       #
# ---------------------------------------------------------------------------------------- #
@testset "Lumen Shannon" begin

end

end
