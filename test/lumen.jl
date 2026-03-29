using Test
using SoleXplorer
const SX = SoleXplorer

using DataTreatments
const DT = DataTreatments

using MLJ
using DataFrames, Random

Xc, yc = @load_iris
Xc = DataFrame(Xc)

Xr, yr = @load_boston
Xr = DataFrame(Xr)

natopsloader = SX.NatopsLoader()
Xts, yts = SX.load(natopsloader)

# ---------------------------------------------------------------------------- #
#                           How to use SoleXplorer                             #
# ---------------------------------------------------------------------------- #
# Import dataset with SoleXplorer/DataTreatments combo
ds = SX.setup_dataset(
    Xc,
    yc,
    # DataTreatments directives
    # dims=0 -> use only tabular portion of the dataset
    SX.TreatmentGroup(dims=0);
    # Choose model and parameters
    model=SX.RandomForestClassifier(
        max_depth = -1,
        n_trees=20
    ),
    # Training directives
    resampling=SX.Holdout(fraction_train=0.7, shuffle=true),
    seed=42
)
# Train model with SoleXplorer
solem = SX.train_test(ds)


# _solexplorer(ds, solem; extractor, measures)
