using DrWatson

using CSV
using DataFrames
using StatsBase
using DataStructures
using Random
using StatsBase
using MLJ
using JLD
using JLD2
using ProgressMeter
# using Revise
using ThreadSafeDicts
using BenchmarkTools

using SoleData
using SoleLogics
using SoleLogics: Operator, ntokens, atoms
using SoleLogics: dual, LeftmostLinearForm, op, Atom
using SoleModels
using SoleModels: Label, AbstractModel, AbstractConditionalAlphabet
using SoleModels: DecisionForest, trees, bestguess, apply!, apply, info
using SoleModels: alphabet, meandelaydl
using SolePostHoc
using SolePostHoc.RuleExtraction: intrees, bellatrex
using ModalDecisionTrees
using ModalDecisionTrees: translate, trees

include("../../lib.jl")
include("../../scanner.jl")
include("../../datasets/dataset-utils.jl")
include("../../datasets/apply-multimodal-modes.jl")

include("../utilities.jl")
include("data.jl")

############################################################################################
################################# Window terminal arguments ################################
############################################################################################

tasks = begin
    ks = collect(keys(inner))
    if length(ARGS) == 0
        ks
    else
        filter(a -> a ∈ ks, ARGS)
    end
end

checkpoint_stdout("Chosen tasks: $(tasks)")

algorithm = intrees
#algorithm = bellatrex

for task in tasks
    checkpoint_stdout("Running Task $(task)")

    isprop, isspatial, models, dataset =  inner[task]

    # Path to save models and corresponding csv file
    father_dir =
        "/home/lele7x/results9/rule-extraction/granular-computing/$(task)/$(algorithm)/"
    results_dir = father_dir * "models_cache/"
    csvpath = father_dir * "$(task)-$(algorithm).csv"
    touch(csvpath)
    efg = open(csvpath, "w")

    # Computing dataset
    X, Y = compute_dataset(dataset, isspatial, isprop)
    #memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i in 1:ninstances(X)]

    # Global variables
    tablesforest = []

    @showprogress for (i,model) in enumerate(models)
        modelpath, nontest_ids, test_ids = model

        rng = MersenneTwister(1)
        println("Computing Testing Dataset in ...")
        X_test,Y_test = @time slicedataset((X,Y), test_ids; return_view = true)
        println("Computing Training Dataset in ...")
        X_nontest,Y_nontest = @time slicedataset((X,Y), nontest_ids; return_view = true)

        #memostruct_test = @view memostruct[test_ids]
        #memostruct_nontest = @view memostruct[nontest_ids]

        checkpoint_stdout("\nExtracting and translating Forest $(i)/$(length(models)) in ...")
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

        finaldl = algorithm(model,X_nontest,Y_nontest; rng=rng)
        #=
        finaldl = algorithm(
            model,X_test,Y_test;
            exec_ntrees=[0.2,0.5,0.8],
            exec_ndims=[2,5,nothing],
            exec_nclusters=[1,2,3],
            rng=rng,
        )
        push!(tablesforest,(MLJ.Accuracy()(finaldl, Y_test))*100)
        @show tablesforest[i]
        =#

        GC.gc()

        println("\nSaving resulting model in ...")
        submodelpath = @time begin
            _hash = get_hash_sha256(finaldl)
            modelpath = results_dir * "model_" * _hash * ".jld2"
            JLD2.save_object(modelpath, finaldl)
            modelpath
        end

        # Computing metrics for initial DecisionForest
        moriginalforest = finalmetrics(model,X_test,Y_test)

        # Computing average metrics for DecisionTrees in DecisionForest
        maveragetrees = finalmetrics(
            trees(model),X_test,Y_test;
            originalmodel=model,
        )

        # Computing metrics for final DecisionList
        mfinaldl = finalmetrics(finaldl,X_test,Y_test)

        # Table construction
        itable = DataFrame(
            Whoami = ["DF","ADT","DL"],
            Path = [
                modelpath,
                modelpath,
                submodelpath,
            ],
            Kappa = [
                moriginalforest[2],
                maveragetrees[2],
                mfinaldl[1],
            ],
            Accuracy = [
                (100.0 - moriginalforest[3]),
                (100.0 - maveragetrees[3]),
                (100.0 - mfinaldl[2]),
            ],
            Error = [
                moriginalforest[3],
                maveragetrees[3],
                mfinaldl[2],
            ],
            NRules = [
                NaN,
                maveragetrees[4],
                mfinaldl[3],
            ],
            NSymbols = [
                NaN,
                maveragetrees[5],
                mfinaldl[4],
            ],
            Delay = [
                NaN,
                maveragetrees[6],
                mfinaldl[5],
            ],
            AbsNRules = [
                NaN,
                maveragetrees[7],
                mfinaldl[6],
            ],
            NLeaves = [moriginalforest[1], maveragetrees[1], NaN],
        )
        push!(tablesforest,itable)
        CSV.write(csvpath, itable, append=true)

        println("Results for Forest $(i): ")
        println(itable)
        println()
        println("Final Decision List: ")
        println(finaldl)
        println()

        GC.gc()
    end

    # Computing average final results
    avgcolumns = []
    stdcolumns = []

    for ncolumn in 3:10
        cols = [filter(a -> !isnan(a), t[:,ncolumn]) for t in tablesforest]
        push!(avgcolumns, first(mean(cols, dims=1)))
        push!(stdcolumns, std(cols))
    end

    finaltable = DataFrame(
        Whoami = [
            "DF","DF std",
            "ADT","ADT std",
            "DL","DL std",
        ],
        Path = [NaN,NaN,NaN,NaN,NaN,NaN],
        Kappa = combining(avgcolumns[1],stdcolumns[1]),
        Accuracy = combining(avgcolumns[2],stdcolumns[2]),
        Error = combining(avgcolumns[3],stdcolumns[3]),
        NRules = [NaN,NaN,combining(avgcolumns[4],stdcolumns[4])...],
        NSymbols = [NaN,NaN,combining(avgcolumns[5],stdcolumns[5])...],
        Delay = [NaN,NaN,combining(avgcolumns[6],stdcolumns[6])...],
        AbsNRules = [NaN,NaN,combining(avgcolumns[7],stdcolumns[7])...],
        NLeaves = [combining(avgcolumns[8],stdcolumns[8])...,NaN,NaN],
    )
    CSV.write(csvpath, finaltable, append=true)

    checkpoint_stdout(
        "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Final Results " *
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    )
    for (i,t) in enumerate(tablesforest)
        println("Results for Forest $(i): ")
        println(t)
    end

    println()
    println("Final Average Table for task $(task)")
    println(finaltable)
    #@show mean(tablesforest)
end
