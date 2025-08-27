############################################################################################
############################################################################################
############################################################################################

# Debug
ho = nothing

function scan_evolutionary_rule_extraction(
    # Starting model
    model::DecisionForest,
    X,
    Y::AbstractVector{<:Label};
    # Search space parameters
    exec_populationSize = [nothing],
    exec_crossoverRate  = collect(0.75:0.05:0.95),
    exec_mutationRate   = collect(0.05:0.05:0.2),
    exec_selection      = [Evolutionary.tournament(2, select=twowaycomp)],
    exec_crossover      = [GENERATE_CROSSOVER_FUNCTION()], # SBX
    exec_mutation       = [
        GENERATE_MUTATION_FUNCTION(nmaxmutations = nmaxmutations)
        for nmaxmutations in [1]
    ],
    exec_ntrees        = [100],
    exec_metricfuns::AbstractVector{<:AbstractVector{<:Function}}, # Evolutionary metrics
    # Hyper-optimization
    Hyperopt_metricfun::Function,
    Hyperopt_niterations::Union{Integer,AbstractFloat} = 0.5,
    # Other
    traintesting::Union{Vector{<:Tuple{Any,AbstractVector{<:Label}}},Nothing} = nothing,
    dataset_slices::Union{Integer,Vector{<:Tuple{Vector,Vector}}} = 10,
    memostruct = [ThreadSafeDict{SyntaxTree,Vector{worldtype(X)}}() for i in 1:ninstances(X)],
    _alphabet::AbstractConditionalAlphabet = SoleModels.alphabet(X),
    _classes::AbstractVector{<:Label} = sort(unique(Y)),
    rng = default_rng(),
)
    # Balanced cross validation
    dataset_slices = begin
        if dataset_slices isa Integer
            @assert isnothing(traintesting) "Can't pass training and testing dataset if" *
                " dataset_slices isa Integer: $(dataset_slices)"
            balanced_cv_dataset_slices(Y, dataset_slices, rng; strict=false)
        else
            dataset_slices
        end
    end

    @assert isnothing(traintesting) || (!isnothing(traintesting) && length(traintesting) == 2 &&
        length(dataset_slices) == 1) "With traintesting, dataset_slices must be " *
        "[(Training Dataset, Training Labels), (Testing Dataset, Testing Labels)]" *
        "\ndataset_slices = $(dataset_slices)"

    println("Dataset slices for internal cross validation")
    @show dataset_slices

    println("\nComputing Reference to Alphabet in ...")
    _alphabet = @time Ref(_alphabet)
    println("Computing Reference to Classes in ...")
    _classes = @time Ref(_classes)

    cube_dim_lengths = length.([
        exec_populationSize,
        exec_crossoverRate,
        exec_mutationRate,
        exec_selection,
        exec_crossover,
        exec_mutation,
        exec_metricfuns,
        exec_ntrees,
    ])
    Hyperopt_tot_niterations = prod(cube_dim_lengths)

    if Hyperopt_niterations isa AbstractFloat
        Hyperopt_niterations = round(Int, Hyperopt_niterations * Hyperopt_tot_niterations)
    end

    println("\nSearch space (size = $(join(cube_dim_lengths, "×")) = $(Hyperopt_tot_niterations))")
    println(" - # Hyper-Iterations: $(Hyperopt_niterations)")
    println(" - Population Size: $(exec_populationSize)")
    println(" - Crossover Rate: $(exec_crossoverRate)")
    println(" - Mutation Rate: $(exec_mutationRate)")
    println(" - Selection: $(exec_selection)")
    println(" - Crossover: $(exec_crossover)")
    println(" - Mutation: $(exec_mutation)")
    println(" - Metrics: $(exec_metricfuns)")
    println(" - # Trees: $(exec_ntrees)")

    @assert ! (Hyperopt_niterations == 1 && Hyperopt_tot_niterations != 1 &&
        !isnothing(dataset_slices)) "Hyperopt_niterations = $(Hyperopt_niterations), " *
        "cube_dim_lengths = $(cube_dim_lengths) (Hyperopt_tot_niterations = $(Hyperopt_tot_niterations))"

    @assert Hyperopt_niterations <= Hyperopt_tot_niterations "$(Hyperopt_niterations) " *
        "<= $(Hyperopt_tot_niterations) must hold for Hyper-optimization to make sense!"

    best_best_dls = nothing
    best_best_metric = nothing
    best_best_niterations = nothing

    # https://github.com/baggepinnen/Hyperopt.jl
    global ho
    #Hyperband(R=50, η=3, inner=RandomSampler(rng)),
    ho = @hyperopt for i = Hyperopt_niterations,
        sampler = RandomSampler(rng),
        populationSize = exec_populationSize,
        crossoverRate  = exec_crossoverRate,
        mutationRate   = exec_mutationRate,
        selection      = exec_selection,
        crossover      = exec_crossover,
        mutation       = exec_mutation,
        metricfuns     = exec_metricfuns,
        ntrees         = exec_ntrees

        best_dls = []
        best_niterations = []
        best_metric = [begin

            train_ids, val_ids = dataset_slice

            X_train, Y_train, X_val, Y_val = begin
                if isnothing(traintesting)
                    println("\nComputing Training Dataset in ...")
                    X_train,Y_train = @time slicedataset((X,Y), train_ids; return_view = true)
                    println("Computing Testing Dataset in ...")
                    X_val,Y_val = @time slicedataset((X,Y), val_ids; return_view = true)

                    (X_train, Y_train, X_val, Y_val)
                else
                    # Note: should use train_ids and val_ids here...
                    (first(traintesting)..., last(traintesting)...)
                end
            end

            memostruct_train = @view memostruct[train_ids]
            memostruct_val = @view memostruct[val_ids]

            function evaluate(indiv::DLIndividual)
                return map(f->f(indiv; X=X_train, Y=Y_train, memostruct = memostruct_train), metricfuns)
            end
            function evaluate(indivs::AbstractVector{<:DLIndividual})
                return map(indiv->evaluate(indiv), indivs)
            end

            @show countmap(Y_train)
            @show countmap(Y_val)

            dls = begin
                tmpdls = extract_decision_lists(model, (X_train, Y_train))
                ntmpdls = length(tmpdls)
                if ntmpdls > ntrees
                    length(tmpdls) % 2 == 0 ? tmpdls : [tmpdls..., rand(rng, tmpdls)]
                elseif ntmpdls < ntrees
                    [tmpdls..., StatsBase.sample(rng, tmpdls, (ntrees-ntmpdls))...]
                else
                    tmpdls
                end
            end
            println("Numero di individui: $(length(dls))")
            _alldls = Ref(dls)
            indivs = map(dl->DLIndividual(dl, _alphabet, _classes, _alldls), dls)
            constrs = Evolutionary.NoConstraints()
            F = zeros(length(metricfuns)) # TODO check fill(-Inf, length(metricfuns))
            # https://wildart.github.io/Evolutionary.jl/stable/moea/
            alg = Evolutionary.NSGA2(;
                populationSize = isnothing(populationSize) ? length(dls) : populationSize,
                crossoverRate = crossoverRate,
                mutationRate = mutationRate,
                selection = selection,
                crossover = crossover,
                mutation = mutation,
            )
            opts = Evolutionary.Options(rng=rng, iterations=500, parallelization=:thread)

            println("\nEvolutionary Computation in ...")
            # l = ReentrantLock()
            # res = Evolutionary.optimize(indiv->(lock(l); evaluate(indiv); unlock(l)), F, constrs, alg, indivs, opts)
            res = @time Evolutionary.optimize(indiv->evaluate(indiv), F, constrs, alg, indivs, opts)

            _best_dls = Evolutionary.minimizer(res)
            push!(best_niterations, Evolutionary.iterations(res))
            append!(best_dls, _best_dls)

            meanmetric = map(bdl->(Hyperopt_metricfun(bdl; X=X_val, Y=Y_val, memostruct = memostruct_val)),_best_dls)
            meanmetric = mean(meanmetric)

            meanmetric
        end for (i_fold, dataset_slice) in enumerate(dataset_slices)]

        best_metric = mean(best_metric)
        if isnothing(best_best_metric) || best_best_metric > best_metric
            best_best_metric = best_metric
            best_best_dls = best_dls
            best_best_niterations = mean(best_niterations)
        end

        println("Best Average Metric for iteration $(i): $(best_metric)")
        cost = best_metric
    end
    best_hyperparam = ho.minimizer

    println("Best Average Metric Obtained: $(ho.minimum)")

    best_hyperparam, best_best_metric, best_best_dls, best_best_niterations
end

############################################################################################
############################################################################################
############################################################################################
