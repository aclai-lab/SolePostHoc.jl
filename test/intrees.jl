using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

Xc, yc = @load_iris
Xc = DataFrame(Xc)

train_ratio = 0.7
rng=Xoshiro(1)
ttpairs = MLJ.MLJBase.train_test_pairs(Holdout(; shuffle=true), 1:length(yc), yc)
train = ttpairs[1][1]
test = ttpairs[1][2]

DTModel = MLJ.@load RandomForestClassifier pkg=DecisionTree verbosity=0
model = DTModel(n_trees=2)
mach = machine(model, Xc, yc)

MLJ.fit!(mach, rows=train, verbosity=0)
featurenames = MLJ.report(mach).features
classlabels = mach.fitresult[2][sortperm((mach).fitresult[3])]
solem_rf = solemodel(MLJ.fitted_params(mach).forest; featurenames, classlabels)
logiset = scalarlogiset(Xc[test, :], allow_propositional=true)
apply!(solem_rf, logiset, yc[test])

# ---------------------------------------------------------------------------- #
#                           intrees rules extraction                           #
# ---------------------------------------------------------------------------- #
config = InTreesConfig(pruning=PruningConfig(
    prune_rules=false), rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(pruning=PruningConfig(
    prune_rules=true), rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

# ---------------------------------------------------------------------------- #
config = InTreesConfig(pruning=PruningConfig(
    prune_rules=true,
    percentage_degradation=false),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(pruning=PruningConfig(
    prune_rules=true,
    percentage_degradation=true),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

# ---------------------------------------------------------------------------- #
config = InTreesConfig(rule_selection=nothing, rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(rule_selection=CBC(), rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

# ---------------------------------------------------------------------------- #
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=false),
    rule_selection=nothing,
    post_process=nothing,
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=false),
    rule_selection=CBC(),
    post_process=nothing,
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

config = InTreesConfig(
    pruning=PruningConfig(prune_rules=true),
    rule_selection=nothing,
    post_process=nothing,
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=true),
    rule_selection=CBC(),
    post_process=nothing,
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

# ---------------------------------------------------------------------------- #
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=false),
    rule_selection=nothing,
    post_process=STEL(),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=false),
    rule_selection=CBC(),
    post_process=STEL(),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

config = InTreesConfig(
    pruning=PruningConfig(prune_rules=true),
    rule_selection=nothing,
    post_process=STEL(),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])
    
config = InTreesConfig(
    pruning=PruningConfig(prune_rules=true),
    rule_selection=CBC(),
    post_process=STEL(),
    rng=Xoshiro(42))
@test_nowarn extracted_rules =
    RuleExtraction.extractrules(config, solem_rf, Xc[test, :], yc[test])

# ---------------------------------------------------------------------------- #
config = InTreesConfig(
    pruning=PruningConfig(
        prune_rules=true,
        percentage_degradation=true),
    rule_selection=CBC(),
    post_process=STEL(min_coverage=0.33), # 3 classes
    rng=Xoshiro(42)
)
@test_nowarn extracted_rules = RuleExtraction.extractrules(
    config,
    solem_rf,
    Xc[test, :],
    yc[test]
)

@test_nowarn extracted_rules = RuleExtraction.extractrules(
    InTreesRuleExtractor,
    solem_rf,
    Xc[test, :],
    yc[test];
    pruning=PruningConfig(
        prune_rules=true,
        percentage_degradation=true),
    rule_selection=CBC(),
    post_process=STEL(min_coverage=0.33), # 3 classes
    rng=Xoshiro(42)
)

@test_throws MethodError InTreesConfig(invalid=true)
@test_throws MethodError RuleExtraction.extractrules(
    solem_rf,
    Xc[test, :],
    yc[test];
    invalid=true,
)
@test_throws MethodError RuleExtraction.extractrules(
    config,
    solem_rf,
    Xc[test, :],
    yc[test];
    invalid=true,
)