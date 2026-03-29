using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

using SoleData.Artifacts

# fill your Artifacts.toml file;
# fillartifacts()


# Loader lists
abcloader = ABCLoader()
mitloader = MITESPRESSOLoader()

Xc, yc = @load_iris
Xc = DataFrame(Xc)

train_ratio = 0.7
rng=Xoshiro(42)
ttpairs = MLJ.MLJBase.train_test_pairs(Holdout(; shuffle = true, rng), 1:length(yc), yc)
train = ttpairs[1][1]
test = ttpairs[1][2]

DTModel = MLJ.@load RandomForestClassifier pkg=DecisionTree verbosity=0
model = DTModel(n_trees = 2, rng=rng)
mach = machine(model, Xc, yc)

MLJ.fit!(mach, rows = train, verbosity = 0)
featurenames = MLJ.report(mach).features
classlabels = mach.fitresult[2][sortperm((mach).fitresult[3])]
solem_rf = solemodel(MLJ.fitted_params(mach).forest; featurenames, classlabels)
logiset = scalarlogiset(Xc[test, :], allow_propositional = true)
apply!(solem_rf, logiset, yc[test])

# ---------------------------------------------------------------------------- #
#                           lumen rules extraction                             #
# ---------------------------------------------------------------------------- #
extractor = LumenRuleExtractor()

extracted_rules = RuleExtraction.extractrules(
    extractor,
    solem_rf;
    minimization_scheme = :abc,
);
