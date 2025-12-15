using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

using SoleData.Artifacts

# Loader lists
abcloader = ABCLoader()
mitloader = MITESPRESSOLoader()

Xc, yc = @load_iris
Xc = DataFrame(Xc)
train_ratio = 0.7
rng=Xoshiro(1)
ttpairs = MLJ.MLJBase.train_test_pairs(Holdout(; shuffle = true), 1:length(yc), yc)
train = ttpairs[1][1]
test = ttpairs[1][2]

DTModel = @load DecisionTreeClassifier pkg=DecisionTree verbosity=0
model = DTModel()
mach = machine(model, Xc, yc)
MLJ.fit!(mach, rows = train, verbosity = 0)
featurenames = MLJ.report(mach).features
classlabels = sort(MLJ.report(mach).classes_seen)
solem_dt = solemodel(MLJ.fitted_params(mach).tree; featurenames, classlabels)

DTModel = @load RandomForestClassifier pkg=DecisionTree verbosity=0
model = DTModel(n_trees = 2)
mach = machine(model, Xc, yc)
MLJ.fit!(mach, rows = train, verbosity = 0)
featurenames = MLJ.report(mach).features
classlabels = mach.fitresult[2][sortperm((mach).fitresult[3])]
solem_rf = solemodel(MLJ.fitted_params(mach).forest; featurenames, classlabels)

extractor = LumenRuleExtractor()
extracted_rules = RuleExtraction.modalextractrules(extractor, solem_dt);
extracted_rules = RuleExtraction.modalextractrules(extractor, solem_rf);