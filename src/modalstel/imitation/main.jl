using DrWatson

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
using Metrics # Heavy dep? TODO maybe

using SoleLogics
using SoleLogics: Operator, ntokens, atoms
using SoleLogics: dual, LeftmostLinearForm, op, Atom
using SoleModels
using SoleModels: Label, AbstractModel, AbstractConditionalAlphabet
using SoleModels: displaymodel
using SoleModels: DecisionForest, trees, bestguess, apply!, apply, info
using SoleModels: alphabet, meandelaydl
using ModalDecisionTrees
using ModalDecisionTrees: translate

include("../../lib.jl")
include("../../datasets/dataset-utils.jl")
include("../../datasets/apply-multimodal-modes.jl")
include("../utilities.jl")
include("../evolutionary/data.jl")

############################################################################################
############################################################################################
############################################################################################

function surrogatetree(
    model::Union{AbstractModel,DecisionForest},
    X,
    Y::AbstractVector{<:Label},
    preds::AbstractVector{<:Label};
    method::Symbol = :modal_decision_tree,
    #memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i in 1:ninstances(X)],
    kwargs...,
)
    @assert method == :modal_decision_tree

    treey, treepreds = begin
        if method == :modal_decision_tree
            XMFMD = MultiLogiset(X)
            treey =
                ModalDecisionTrees.translate(ModalDecisionTrees.build_tree(XMFMD, Y))
            treepreds = ModalDecisionTrees.translate(
                ModalDecisionTrees.build_tree(XMFMD, preds),
            ) # TODO create alias

            (treey, treepreds)
        else
            error("Unknown surrogate method specified: $(method)")
        end
    end

    return (treey, treepreds)
end

#=
# Regression Problem
function r2(
    m::AbstractModel;
    X,
    Y::AbstractVector{<:Label},
    memostruct,
    kwargs...
)
    Y_pred = apply(m, X; check_kwargs = (; use_memo = memostruct), kwargs...)
    return Metrics.r2_score(Y_pred, Y)
end
=#

############################################################################################
############################################################################################
############################################################################################

tasks, stat, temp = begin
    ks = keys(inner)
    if length(ARGS) == 0
        (ks, true, true)
    else
        availabletasks = filter(a -> a ∈ ["T1", "T2", "T3", "S1", "S2", "S3"], ARGS)
        if length(availabletasks) == 0
            error("At least one of the passed arguments is unknown, ARGS: $(ARGS)")
        elseif length(findall(ARGS .== "static")) == 0 &&
               length(findall(ARGS .== "temporal")) == 0
            (availabletasks, true, true)
        elseif length(findall(ARGS .== "static")) == 1 &&
               length(findall(ARGS .== "temporal")) == 1
            (availabletasks, true, true)
        elseif length(findall(ARGS .== "static")) == 1
            (availabletasks, true, false)
        elseif length(findall(ARGS .== "temporal")) == 1
            (availabletasks, false, true)
        end
    end
end
tasks = begin
    if stat && temp
        [[(t, true) for t in tasks]..., [(t, false) for t in tasks]...]
    elseif stat
        [(t, true) for t in tasks]
    else
        temp
        [(t, false) for t in tasks]
    end
end
checkpoint_stdout("Chosen tasks: $(tasks)")

for (data, isstatic) in tasks
    checkpoint_stdout("Running Task $(data) $(isstatic ? "static" : "temporal")")
    rng = Random.MersenneTwister(1)
    isspatial, staticmodels, temporalmodels, dataset = inner[data]
    models = isstatic ? staticmodels : temporalmodels

    X, Y = compute_dataset(dataset, isspatial, isstatic)
    memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i = 1:ninstances(X)]

    accforests = []
    acctreesy = []
    acctreespreds = []
    fidforests = []
    fidtreesy = []
    fidtreespreds = []
    #r2forests = []
    #r2trees = []

    @showprogress "Computing Forests..." for (i, m) in enumerate(models)

        checkpoint_stdout("Running Forest number $(i)")

        modelpath, nontest_ids, test_ids = m

        println("Computing Testing Dataset in ...")
        X_test, Y_test = @time slicedataset((X, Y), test_ids; return_view = true)
        println("Computing Training Dataset in ...")
        X_nontest, Y_nontest = @time slicedataset((X, Y), nontest_ids; return_view = true)

        memostruct_test = @view memostruct[test_ids]
        memostruct_nontest = @view memostruct[nontest_ids]

        checkpoint_stdout("Extracting and translating Forest $(i)/$(length(models)) in ...")
        model = @time begin
            m = JLD2.load(modelpath)
            m = m["model_pruned"]
            ModalDecisionTrees.translate(m)
        end

        println("\nTraining Indices: $(nontest_ids)")
        println("Testing Indices: $(test_ids)")
        println()
        @show model

        preds_train = apply(
            model,
            X_nontest;
            check_kwargs = (; use_memo = memostruct_nontest),
            suppress_parity_warning = true,
        )
        preds_test = apply(
            model,
            X_test;
            check_kwargs = (; use_memo = memostruct_test),
            suppress_parity_warning = true,
        )

        println("\nComputing Surrogate Tree in ...")
        treey, treepreds = @time surrogatetree(
            model,
            X_nontest,
            Y_nontest,
            preds_train;
            method = :modal_decision_tree,
        )
        println(displaymodel(treey; header = false))
        println(displaymodel(treepreds; header = false))

        #accuracy
        accf = accuracy(
            model;
            X = X_test,
            Y = Y_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(accforests, accf)
        accy = accuracy(
            treey;
            X = X_test,
            Y = Y_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(acctreesy, accy)
        accpreds = accuracy(
            treepreds;
            X = X_test,
            Y = Y_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(acctreespreds, accpreds)
        #fidelity
        fidf = accuracy(
            model;
            X = X_test,
            Y = preds_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(fidforests, fidf)
        fidy = accuracy(
            treey;
            X = X_test,
            Y = preds_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(fidtreesy, fidy)
        fidpreds = accuracy(
            treepreds;
            X = X_test,
            Y = preds_test,
            memostruct = memostruct_test,
            suppress_parity_warning = true,
        )
        push!(fidtreespreds, fidpreds)

        #r2f = r2(model; X=X_test, Y=Y_test, memostruct=memostruct_test, suppress_parity_warning=true)
        #push!(r2forests, r2f)
        #r2t = r2(tree; X=X_test, Y=Y_test, memostruct=memostruct_test, suppress_parity_warning=true)
        #push!(r2trees, r2t)

        println("Results for Forest $(i):")
        df = DataFrame(
            Who = ["Forest", "Tree", "Surrogate Tree"],
            Accuracy = [accf, accy, accpreds],
            Fidelity = [fidf, fidy, fidpreds],
            #R2 = [r2f, r2t],
        )
        println(df)
    end

    meanaccforests = mean(accforests)
    stdaccforests = std(accforests)
    meanacctreesy = mean(acctreesy)
    stdacctreesy = std(acctreesy)
    meanacctreespreds = mean(acctreespreds)
    stdacctreespreds = std(acctreespreds)

    meanfidforests = mean(fidforests)
    stdfidforests = std(fidforests)
    meanfidtreesy = mean(fidtreesy)
    stdfidtreesy = std(fidtreesy)
    meanfidtreespreds = mean(fidtreespreds)
    stdfidtreespreds = std(fidtreespreds)

    #meanr2forests = mean(r2forests)
    #stdr2forests = std(r2forests)
    #meanr2trees = mean(r2trees)
    #stdr2trees = std(r2trees)

    println("Average Results for Task $(data) $(isstatic ? "static" : "temporal"):")
    df = DataFrame(
        Who = [
            "Forest mean",
            "Forest std",
            "Tree mean",
            "Tree std",
            "STree mean",
            "STree std",
        ],
        Accuracy = [
            meanaccforests,
            stdaccforests,
            meanacctreesy,
            stdacctreesy,
            meanacctreespreds,
            stdacctreespreds,
        ],
        Fidelity = [
            meanfidforests,
            stdfidforests,
            meanfidtreesy,
            stdfidtreesy,
            meanfidtreespreds,
            stdfidtreespreds,
        ],
        #R2 = [meanr2forests, stdr2forests, meanr2trees, stdr2trees],
    )
    println(df)
end
