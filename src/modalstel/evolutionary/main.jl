using Evolutionary
import Evolutionary: default_values
using Evolutionary: funargnum, tournament, twowaycomp, default_options

using DrWatson
using Hyperopt

using DataFrames
using StatsBase
using DataStructures
using Random
using StatsBase
using MLJ
using JLD2
using ProgressMeter
# using Revise
using Plots
using Metaheuristics
using ThreadSafeDicts
using BenchmarkTools

using SoleLogics
using SoleLogics: Operator, ntokens, atoms
using SoleLogics: dual, LeftmostLinearForm, op, Atom
using SoleModels
using SoleModels: Label, AbstractModel, AbstractConditionalAlphabet
using SoleModels: DecisionForest, trees, bestguess, apply!, apply, info
using SoleModels: alphabet, meandelaydl, LeafModel
using ModalDecisionTrees
using ModalDecisionTrees: translate

include("../../lib.jl")
include("../../datasets/dataset-utils.jl")
include("../../datasets/apply-multimodal-modes.jl")

# Include utility functions and data used in the flow
include("../utilities.jl")
include("data.jl")

# Include functions to implement evolutionary programming
include("evolutionary.jl")
include("crossover.jl")
include("mutations.jl")
include("experiment.jl")

############################################################################################
###################################### MAIN ################################################
############################################################################################

tasks = begin
    ks = collect(keys(inner))
    if length(ARGS) == 0
        ks
    else
        availabletasks = filter(a -> a ∈ ["T1", "T2", "T3", "S1", "S2", "S3"], ARGS)
        if length(availabletasks) == 0
            error("At least one of the passed arguments is unknown, ARGS: $(ARGS)")
        else
            availabletasks
        end
    end
end

checkpoint_stdout("Passed tasks: $(tasks)")


#return [(_error(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train))]
#return [(_error(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train) * absnrulescomplexity(indiv; Y=Y_train))]
#return [(_error(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train) * absnrulescomplexity(indiv; Y=Y_train) * nsymbolscomplexity(indiv))]
#return [(_error(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train) * absnrulescomplexity(indiv; Y=Y_train) * nsymbolscomplexity(indiv) * meandelay(indiv; X=X_train, memostruct = memostruct_train))]
#return [(_error(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train)), (absnrulescomplexity(indiv; Y=Y_train))]
Evolutionary_metricfuns = [_error, absnrulescomplexity] # Note: must be > 0.0
#meanmetric = map(bdl->(_error(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val) * absnrulescomplexity(bdl; Y=Y_val)),_best_dls)
#meanmetric = map(bdl->(_error(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val) * absnrulescomplexity(bdl; Y=Y_val) * nsymbolscomplexity(bdl)),_best_dls)
#meanmetric = map(bdl->(_error(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val) * absnrulescomplexity(bdl; Y=Y_val) * nsymbolscomplexity(bdl) * meandelay(bdl; X=X_val, memostruct = memostruct_val)),_best_dls)
#meanmetric = map(bdl->mean([(_error(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val)), (absnrulescomplexity(bdl; Y=Y_val))]),_best_dls)
#meanmetric = map(bdl->(_error(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val)),_best_dls)
Hyperopt_metricfun = _error

ntrees = 10

for task in tasks
    checkpoint_stdout("Running Task $(task)")

    isprop, isspatial, models, dataset, nametask = inner[task]

    # Path to save models and corresponding csv file
    father_dir = "/home/lele7x/results9/rule-extraction/evolutionary/"
    results_dir = father_dir * "models_cache/$(nametask)"
    csvpath = father_dir * "csv_cache/$(nametask).csv"
    touch(csvpath)
    efg = open(csvpath, "w")

    # Computing dataset
    X, Y = compute_dataset(dataset, isspatial, isprop)

    # Extracting alphabet and classes
    println("\nExtracting Alphabet in ...")
    _alphabet = @time SoleModels.alphabet(X)
    println("Extracting Classes in ...")
    _classes = @time sort(unique(Y))
    memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i = 1:ninstances(X)]

    # Global variables
    best_hyperparam = nothing
    tablesforest = nothing

    @showprogress "Computing Forests..." for (i, m) in enumerate(models)

        modelpath, nontest_ids, test_ids = m

        rng = MersenneTwister(1)
        println("Computing Testing Dataset in ...")
        X_test, Y_test = @time slicedataset((X, Y), test_ids; return_view = true)
        println("Computing Training Dataset in ...")
        X_nontest, Y_nontest = @time slicedataset((X, Y), nontest_ids; return_view = true)

        memostruct_test = @view memostruct[test_ids]
        memostruct_nontest = @view memostruct[nontest_ids]

        checkpoint_stdout(
            "\nExtracting and translating Forest $(i)/$(length(models)) in ...",
        )
        model = @time begin
            m = JLD2.load(modelpath)
            m = m["model_pruned"]
            ModalDecisionTrees.translate(m)
        end

        # Printing
        println("\nTraining Indices: $(nontest_ids)")
        println("Testing Indices: $(test_ids)")
        println()
        @show model

        if i == 1
            println()
            checkpoint_stdout("Training to get the best parameterization...")
            best_hyperparam, best_best_metric, best_best_dls_cv, best_best_niterations =
                scan_evolutionary_rule_extraction(
                    model,
                    X_nontest,
                    Y_nontest;
                    exec_metricfuns = [Evolutionary_metricfuns],
                    exec_ntrees = [ntrees],
                    Hyperopt_metricfun = Hyperopt_metricfun,
                    dataset_slices = 2,
                    memostruct = memostruct_nontest,
                    #
                    _alphabet = _alphabet,
                    _classes = _classes,
                    rng = rng,
                )

            println("\nObtained best hyper parameterization model with path: $(modelpath)")
            println("Best Hyper Parameterization: $(best_hyperparam)")
            println("Best obtained metric: $(best_best_metric)")
            println("Best number of iterations: $(best_best_niterations)")
            println("Best Population: $(best_best_dls_cv)")
        end

        @assert length(best_hyperparam) == 8
        checkpoint_stdout("Running Forest $(i)/$(length(models))")
        _, best_best_metric, best_best_dls, best_best_niterations =
            scan_evolutionary_rule_extraction(
                model,
                X,
                Y;
                exec_populationSize = [best_hyperparam[1]],
                exec_crossoverRate = [best_hyperparam[2]],
                exec_mutationRate = [best_hyperparam[3]],
                exec_selection = [best_hyperparam[4]],
                exec_crossover = [best_hyperparam[5]],
                exec_mutation = [best_hyperparam[6]],
                exec_metricfuns = [best_hyperparam[7]],
                exec_ntrees = [best_hyperparam[8]],
                #
                Hyperopt_metricfun = Hyperopt_metricfun,
                Hyperopt_niterations = 1,
                traintesting = [(X_nontest, Y_nontest), (X_test, Y_test)],
                dataset_slices = [(nontest_ids, test_ids),],
                memostruct = memostruct,
                #
                _alphabet = _alphabet,
                _classes = _classes,
                rng = rng,
            )

        # Computing metrics for initial DecisionForest
        moriginalforest = finalmetrics(model, X_test, Y_test; memostruct = memostruct_test)

        # Computing average metrics for DecisionTrees in DecisionForest
        maveragetrees = finalmetrics(
            trees(model),
            X_test,
            Y_test;
            memostruct = memostruct_test,
            originalmodel = model,
        )

        # Computing metrics for DecisionList with maximum kappa
        kappas = map(
            dli->(kappa(dli, X_test, Y_test; memostruct = memostruct_test)),
            best_best_dls,
        )
        idkappamax = argmax(kappas)
        kmax = best_best_dls[idkappamax]
        kmaxmetrics = finalmetrics(kmax, X_test, Y_test; memostruct = memostruct_test)

        # Computing metrics for DecisionList with minimum delay
        meandelays =
            map(dli->meandelay(dli, X_test; memostruct = memostruct_test), best_best_dls)
        iddelaymin = argmin(meandelays)
        dmin = best_best_dls[iddelaymin]
        dminmetrics = finalmetrics(dmin, X_test, Y_test; memostruct = memostruct_test)

        # Computing metrics for DecisionList with minimum symbols number
        symbolsnumbers = symbolnumber.(best_best_dls)
        idsymbolmin = argmin(symbolsnumbers)
        smin = best_best_dls[idsymbolmin]
        sminmetrics = finalmetrics(smin, X_test, Y_test; memostruct = memostruct_test)

        # Computing metrics for DecisionList with minimum rules number
        rulesnumbers = nrules.(best_best_dls)
        idrulesmin = argmin(rulesnumbers)
        rmin = best_best_dls[idrulesmin]
        rminmetrics = finalmetrics(rmin, X_test, Y_test; memostruct = memostruct_test)

        println("\nSaving resulting model in ...")
        submodelspath = @time [
            begin
                _hash = get_hash_sha256(m)
                mpath = results_dir * "/model_" * _hash * ".jld2"
                JLD2.save_object(mpath, m)
                mpath
            end for m in [kmax, dmin, smin, rmin]
        ]

        # Table construction
        itable = DataFrame(
            Whoami = ["DF", "ADT", "DK", "DD", "DS", "DR"],
            Path = [modelpath, modelpath, submodelspath...],
            Kappa = [
                moriginalforest[2],
                maveragetrees[2],
                kappasmax[1],
                delaysmin[1],
                symbolsmin[1],
                rulesmin[1],
            ],
            Accuracy = [
                (100.0 - moriginalforest[3]),
                (100.0 - maveragetrees[3]),
                (100.0 - kappasmax[2]),
                (100.0 - delaysmin[2]),
                (100.0 - symbolsmin[2]),
                (100.0 - rulesmin[2]),
            ],
            Error = [
                moriginalforest[3],
                maveragetrees[3],
                kappasmax[2],
                delaysmin[2],
                symbolsmin[2],
                rulesmin[2],
            ],
            NRules = [
                NaN,
                maveragetrees[4],
                kappasmax[3],
                delaysmin[3],
                symbolsmin[3],
                rulesmin[3],
            ],
            NSymbols = [
                NaN,
                maveragetrees[5],
                kappasmax[4],
                delaysmin[4],
                symbolsmin[4],
                rulesmin[4],
            ],
            Delay = [
                NaN,
                maveragetrees[6],
                kappasmax[5],
                delaysmin[5],
                symbolsmin[5],
                rulesmin[5],
            ],
            AbsNRules = [
                NaN,
                maveragetrees[7],
                kappasmax[6],
                delaysmin[6],
                symbolsmin[6],
                rulesmin[6],
            ],
            NLeaves = [moriginalforest[1], maveragetrees[1], NaN, NaN, NaN, NaN],
            NIterations = [best_best_niterations, NaN, NaN, NaN, NaN, NaN],
        )
        push!(tablesforest, itable)

        println("Results for Forest $(i): ")
        println(itable)
        println()
        println("Kappa Max Individual: ")
        println(kmax)
        println()
        println("Delay Min Individual: ")
        println(dmin)
        println()
        println("Symbol Min Individual: ")
        println(smin)
        println()
        println("NRules Min Individual: ")
        println(rmin)
        println()
    end

    # Computing average final results
    avgcolumns = []
    stdcolumns = []

    for ncolumn = 3:11
        cols = [filter(a -> !isnan(a), t[:, ncolumn]) for t in tablesforest]
        push!(avgcolumns, mean(cols, dims = 1))
        push!(stdcolumns, std(cols))
    end

    finaltable = DataFrame(
        Whoami = [
            "DF",
            "DF std",
            "ADT",
            "ADT std",
            "DK",
            "DK std",
            "DD",
            "DD std",
            "DS",
            "DS std",
            "DR",
            "DR std",
        ],
        Kappa = combining(avgcolumns[1], stdcolumns[1]),
        Accuracy = combining(avgcolumns[2], stdcolumns[2]),
        Error = combining(avgcolumns[3], stdcolumns[3]),
        NRules = [NaN, NaN, combining(avgcolumns[4], stdcolumns[4])...],
        NSymbols = [NaN, NaN, combining(avgcolumns[5], stdcolumns[5])...],
        Delay = [NaN, NaN, combining(avgcolumns[6], stdcolumns[6])...],
        AbsNRules = [NaN, NaN, combining(avgcolumns[7], stdcolumns[7])...],
        NLeaves = [
            combining(avgcolumns[8], stdcolumns[8])...,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
        ],
        NIterations = [
            combining(avgcolumns[9], stdcolumns[9])...,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
        ],
    )

    println("Final Average Table for task $(task)")
    println(finaltable)
end
