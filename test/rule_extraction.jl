using Test
using SoleModels, SolePostHoc
using MLJ
using DataFrames, Random

using SoleData.Artifacts

# fill your Artifacts.toml file;
# @test_nowarn fillartifacts()
fillartifacts()


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

DTModel = MLJ.@load DecisionTreeClassifier pkg=DecisionTree verbosity=0
model = DTModel()
mach = machine(model, Xc, yc)

MLJ.fit!(mach, rows = train, verbosity = 0)
featurenames = MLJ.report(mach).features
classlabels = sort(MLJ.report(mach).classes_seen)
solem_dt = solemodel(MLJ.fitted_params(mach).tree; featurenames, classlabels)
logiset = scalarlogiset(Xc[test, :], allow_propositional = true)
apply!(solem_dt, logiset, yc[test])

DTModel = MLJ.@load RandomForestClassifier pkg=DecisionTree verbosity=0
model = DTModel(n_trees = 2)
mach = machine(model, Xc, yc)

MLJ.fit!(mach, rows = train, verbosity = 0)
featurenames = MLJ.report(mach).features
classlabels = mach.fitresult[2][sortperm((mach).fitresult[3])]
solem_rf = solemodel(MLJ.fitted_params(mach).forest; featurenames, classlabels)
logiset = scalarlogiset(Xc[test, :], allow_propositional = true)
apply!(solem_rf, logiset, yc[test])

# ---------------------------------------------------------------------------- #
#                          intrees rules extraction                           #
# ---------------------------------------------------------------------------- #
extractor = InTreesRuleExtractor()

extracted_rules =
    RuleExtraction.modalextractrules(extractor, solem_dt, Xc[test, :], yc[test])
extracted_rules = RuleExtraction.modalextractrules(
    extractor,
    solem_dt,
    Xc[test, :],
    yc[test];
    min_coverage = 1.0,
)
@test_throws MethodError RuleExtraction.modalextractrules(
    extractor,
    solem_dt,
    Xc[test, :],
    yc[test];
    invalid = true,
)

# ---------------------------------------------------------------------------- #
#                           lumen rules extraction                             #
# ---------------------------------------------------------------------------- #
extractor = LumenRuleExtractor()

extracted_rules = RuleExtraction.modalextractrules(extractor, solem_dt);
extracted_rules = RuleExtraction.modalextractrules(
    extractor,
    solem_dt;
    minimization_scheme = :mitespresso,
);

extracted_rules = RuleExtraction.modalextractrules(extractor, solem_rf);
extracted_rules = RuleExtraction.modalextractrules(
    extractor,
    solem_rf;
    minimization_scheme = :mitespresso,
);

# ---------------------------------------------------------------------------- #
#                          batrees rules extraction                            #
# ---------------------------------------------------------------------------- #
extractor=BATreesRuleExtractor()

extracted_rules =
    RuleExtraction.modalextractrules(extractor, solem_rf; dataset_name = "Sole_Analysis")
extracted_rules = RuleExtraction.modalextractrules(
    extractor,
    solem_rf;
    dataset_name = "Sole_Analysis",
    num_trees = 5,
)

# ---------------------------------------------------------------------------- #
#                         rulecosi rules extraction                            #
# ---------------------------------------------------------------------------- #
    # extractor=RULECOSIPLUSRuleExtractor()
    # 
    # extracted_rules =
    #     RuleExtraction.modalextractrules(extractor, solem_rf, Xc[test, :], yc[test])
    # 
# ---------------------------------------------------------------------------- #
#                           refne rules extraction                             #
# ---------------------------------------------------------------------------- #
extractor=REFNERuleExtractor()

Xmin = map(minimum, eachcol(Xc[test, :]))
Xmax = map(maximum, eachcol(Xc[test, :]))
extracted_rules = RuleExtraction.modalextractrules(extractor, solem_rf, Xmin, Xmax; L = 2)

# ---------------------------------------------------------------------------- #
#                          trepan rules extraction                             #
# ---------------------------------------------------------------------------- #
extractor=TREPANRuleExtractor()

extracted_rules = RuleExtraction.modalextractrules(extractor, solem_rf, Xc[test, :])
