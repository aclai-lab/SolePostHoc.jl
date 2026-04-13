# ---------------------------------------------------------------------------- #
#                               test per Perry                                 #
# ---------------------------------------------------------------------------- #
# SoleData: branch devPaso
# ModalDecisionLists: branch devPaso

using DataTreatments
const DT = DataTreatments
using SoleXplorer
const SX = SoleXplorer

using ModalDecisionLists
using SolePostHoc

using MLJ
using DataFrames, Random

Xc, yc = @load_iris
Xc = DataFrame(Xc)

iris = DT.load_dataset(Xc, yc)

function model_wrapper(X, y, w; featurenames, rng, iteration, kwargs...)
    irepstar(X, y; featurenames, max_k=1, rng)
end

# natopsloader = SX.NatopsLoader()
# Xts, yts = SX.load(natopsloader)

# natops = DT.load_dataset(
#     Xts, yts,
#     TreatmentGroup(
#         dims=1,
#         aggrfunc=SX.aggregate(
#             # win=SX.adaptivewindow(nwindows=3, overlap=0.3),
#             # features=(maximum, mean)
#             win=SX.wholewindow(),
#             features=mean
#         )
#     )
# )
# ---------------------------------------------------------------------------- #
#                                random forest                                 #
# ---------------------------------------------------------------------------- #
@btime SX.solexplorer(
    iris,
    model=SX.RandomForestClassifier(n_trees=50,),
    # extractor=SX.LumenRuleExtractor(minimization_scheme=:abc),
    extractor=SX.LumenRuleExtractor(),
    seed=42
)
#   210.685 s (727569802 allocations: 403.57 GiB)
# ModelSet{DataSet{RandomForestClassifier, Int64}}:
#   Dataset: DataSet{RandomForestClassifier, Int64}
#   Models:  1 symbolic models
#   Rules: 3 extracted rules per model
#   Measures:
#     Accuracy() = 0.9555555555555556
#     Kappa() = 0.9317147192716237

# ---------------------------------------------------------------------------- #
#                           random forest float32                              #
# ---------------------------------------------------------------------------- #
@btime SX.solexplorer(
    iris,
    model=SX.RandomForestClassifier(n_trees=50,),
    # extractor=SX.LumenRuleExtractor(minimization_scheme=:abc),
    extractor=SX.LumenRuleExtractor(float_type=Float32,),
    seed=42
)
#   96.081 s (325824597 allocations: 197.66 GiB)
# ModelSet{DataSet{RandomForestClassifier, Int64}}:
#   Dataset: DataSet{RandomForestClassifier, Int64}
#   Models:  1 symbolic models
#   Rules: 3 extracted rules per model
#   Measures:
#     Accuracy() = 0.9555555555555556
#     Kappa() = 0.9317147192716237

# ---------------------------------------------------------------------------- #
#                                   xgboost                                    #
# ---------------------------------------------------------------------------- #
@btime SX.solexplorer(
    iris,
    model=SX.XGBoostClassifier(num_round=100,),
    # extractor=SX.LumenRuleExtractor(minimization_scheme=:abc),
    extractor=SX.LumenRuleExtractor(),
    seed=42
)
#   271.988 ms (1839800 allocations: 158.14 MiB)
# ModelSet{DataSet{XGBoostClassifier, Int64}}:
#   Dataset: DataSet{XGBoostClassifier, Int64}
#   Models:  1 symbolic models
#   Rules: 3 extracted rules per model
#   Measures:
#     Accuracy() = 0.9555555555555556
#     Kappa() = 0.9310872894333843

# ---------------------------------------------------------------------------- #
#                               xgboost float32                                #
# ---------------------------------------------------------------------------- #
@btime SX.solexplorer(
    iris,
    model=SX.XGBoostClassifier(num_round=100,),
    # extractor=SX.LumenRuleExtractor(minimization_scheme=:abc),
    extractor=SX.LumenRuleExtractor(float_type=Float32,),
    seed=42
)
#   284.528 ms (1840659 allocations: 133.86 MiB)
# ModelSet{DataSet{XGBoostClassifier, Int64}}:
#   Dataset: DataSet{XGBoostClassifier, Int64}
#   Models:  1 symbolic models
#   Rules: 3 extracted rules per model
#   Measures:
#     Accuracy() = 0.9555555555555556
#     Kappa() = 0.9310872894333843

# ---------------------------------------------------------------------------- #
#                       random decision list ensemble                          #
# ---------------------------------------------------------------------------- #
rng = Xoshiro(42)
Xc, yc = @load_iris
Xc = DataFrame(Xc)

featurenames = Symbol.(names(Xc))
logiset = scalarlogiset(Xc; featurenames, allow_propositional=true)
ensemble_model = build_ensemble(logiset, yc, 50; featurenames, model_wrapper)

@btime lumen(ensemble_model);
# 422.738 s (567032367 allocations: 97.72 GiB)

# ---------------------------------------------------------------------------- #
#                   random decision list ensemble float32                      #
# ---------------------------------------------------------------------------- #
rng = Xoshiro(42)
Xc, yc = @load_iris
Xc = DataFrame(Xc)

featurenames = Symbol.(names(Xc))
logiset = scalarlogiset(Xc; featurenames, allow_propositional=true)
ensemble_model = build_ensemble(logiset, yc, 50; featurenames, model_wrapper)

@btime lumen(ensemble_model; float_type=Float32);
# 366.494 s (538308022 allocations: 52.31 GiB)

