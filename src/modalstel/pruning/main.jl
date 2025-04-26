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
using ModalDecisionTrees
using ModalDecisionTrees: translate, trees

include("../../lib.jl")
include("../../datasets/dataset-utils.jl")
include("../../datasets/apply-multimodal-modes.jl")

include("../utilities.jl")
include("data.jl")
include("sampling.jl")

tasks = begin
    ks = collect(keys(inner))
    if length(ARGS) == 0
        [
            [(k,:rantic) for k in ks]...,
            [(k,:increment) for k in ks]...,
            [(k,:genetic) for k in ks]...,
        ]
    else
        # Tasks
        availabletasks = filter(a -> a ∈ ["T1","B1"], ARGS)
        @assert length(availabletasks) != 0 "At least one of the passed arguments is " *
            "unknown, ARGS: $(ARGS)"

        availablepolicies = filter(p -> p ∈ ["rantic","increment","genetic"], ARGS)

        if length(availablepolicies) == 0
            [
                [(k,:rantic) for k in ks]...,
                [(k,:increment) for k in ks]...,
                [(k,:genetic) for k in ks]...,
            ]
        else
            tasks = []
            for p in availablepolicies
                for t in availabletasks
                    push!(tasks, (t,:p))
                end
            end
            tasks
        end
    end
end

checkpoint_stdout("Chosen tasks: $(tasks)")

# Global variables
#tablesforest = nothing

for (task,policy) in tasks
    checkpoint_stdout("Running Task $(task) with Policy $(policy)")

    rng = MersenneTwister(1)

    isprop, isspatial, models, dataset, nametask =  inner["$(task-policy)"]

    # Path to save models and corresponding csv file
    father_dir = "/home/lele7x/results9/rule-extraction/pruning/"
    results_dir = father_dir * "models_cache/$(nametask)"
    csvpath = father_dir * "csv_cache/$(nametask).csv"
    plotpath = "plots_cache/$(nametask)/"
    touch(csvpath)
    efg = open(csvpath, "w")

    # Computing dataset
    X, Y = compute_dataset(dataset, isspatial, isprop)

    # Extracting alphabet and classes
    println("\nExtracting Alphabet in ...")
    _alphabet = @time SoleData.alphabet(X)
    println("Extracting Classes in ...")
    _classes = @time sort(unique(Y))
    memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i in 1:ninstances(X)]

    @showprogress "Computing Forests..." for (i,m) in enumerate(models)

        modelpath, nontest_ids, test_ids = m

        println("Computing Testing Dataset in ...")
        X_test,Y_test = @time slicedataset((X,Y), test_ids; return_view = true)
        println("Computing Training Dataset in ...")
        X_nontest,Y_nontest = @time slicedataset((X,Y), nontest_ids; return_view = true)

        memostruct_test = @view memostruct[test_ids]
        memostruct_nontest = @view memostruct[nontest_ids]

        checkpoint_stdout("\nExtracting and translating Forest in ...")
        model = @time begin
            m = JLD2.load_object(modelpath)
        end

        println("\nTraining Indices: $(nontest_ids)")
        println("Testing Indices: $(test_ids)")
        println()
        @show model

        println("\nSelecting Subforest in ...")
        subforest = @time samplingensemble(
            model,
            policy,
            X_nontest,
            Y_nontest;
            step=10,
            nforests=1000,
            npopulation=10,
            bounds=[[0,1],[100,10]],
            rng=rng,
            fileplots=plotpath,
        )

        println("\nSaving resulting model in ...")
        submodelpath = @time begin
            _hash = get_hash_sha256(subforest)
            modelpath = results_dir * "/model_" * _hash * ".jld2"
            JLD2.save_object(modelpath, subforest)
            modelpath
        end

        # Metrics: kappa, accuracy, sensitivity, specificity

        # Original DForest
        cmorigforest = finalmetrics(model,X_test,Y_test; suppress_parity_warning=true)
        # Pruning DForest
        cmsubforest = finalmetrics(subforest,X_test,Y_test; suppress_parity_warning=true)
        # Average DTree
        cmavgsubforest = begin
            cmalls = nothing

            for tree in ModalDecisionTrees.trees(subforest)
                currentcm = finalmetrics(tree; X=X_test, Y=Y_test, suppress_parity_warning=true)

                cmalls = isnothing(cmalls) ? currentcm' : [cmalls; currentcm']
            end

            mean(cmalls, dims=1)
        end

        df = DataFrame(
            Task = ["$(data)" for _ in 1:3],
            Type = ["$(isstatic ? "static" : "modal")" for _ in 1:3],
            Forest = ["Original", "Subforest", "Avgforest"],
            Path = [modelpath, submodelpath, "Nothing"],
            NTrees = [
                length(ModalDecisionTrees.trees(model)),
                length(ModalDecisionTrees.trees(subforest)),
                length(ModalDecisionTrees.trees(subforest)),
            ],
            Kappa = [cmorigforest[1], cmsubforest[1], cmavgsubforest[1]],
            Accuracy = [cmorigforest[2], cmsubforest[2], cmavgsubforest[2]],
            Sensitivity = [cmorigforest[3], cmsubforest[3], cmavgsubforest[3]],
            Specificity = [cmorigforest[4], cmsubforest[4], cmavgsubforest[4]],
        )
        #push!(tablesforest,df)

        CSV.write(csvpath, df, append=true)
        println("\nResults for task $(data) for policy $(policy):")
        println(df)
    end
end
