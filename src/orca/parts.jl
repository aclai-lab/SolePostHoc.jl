"""
    parts

Bitvector-driven forest compression pipeline.

This module consumes a `BitVector` whose length depends on which the pipeline
stages are active. If all three stages are enabled the vector has the usual
`3 * n_trees` layout:

  [ b_1 … b_N | b_(N+1) … b_(2N) | b_(2N+1) … b_(3N) ]
    presence     depth / pruning     alphabet

When only a subset of stages is requested the vector is compacted: only the
segments for the active stages are present, in ascending order of stage number.
The helper `_segment` resolves the correct slice for each stage at runtime.

# Pipeline
1. `first_part`  – select trees   (first N bits)
2. `second_part` – prune trees    (second N bits)
3. `third_part`  – modify trees' alphabet (third N bits)

# Exports
- `first_part`   – tree selection based on presence mask
- `second_part`  – tree pruning based on depth mask
- `third_part`   – tree alphabet manipulation based on alphabet mask
- `core`         – evolutionary optimisation of the bitvector
"""




# ─────────────────────────────────────────────────────────────────────────────
# INTERNAL HELPER — segment resolver
# ─────────────────────────────────────────────────────────────────────────────

"""
    _segment(stage, n_trees, active_parts) → UnitRange{Int}

Return the slice of the compacted bitvector that corresponds to `stage`.
`active_parts` must be a sorted collection of integers from {1, 2, 3}.

Example: active_parts = [1, 3], n_trees = 4
  stage 1 → 1:4
  stage 3 → 5:8
"""
function _segment(stage::Int, n_trees::Int, active_parts)
    idx = findfirst(==(stage), active_parts)
    idx === nothing && error("Stage $stage is not in active_parts = $active_parts")
    start = (idx - 1) * n_trees + 1
    stop  = idx * n_trees
    return start:stop
end


# ─────────────────────────────────────────────────────────────────────────────
# STAGE 1 — TREE SELECTION
# ─────────────────────────────────────────────────────────────────────────────

"""
    first_part(n_trees, original_f, original_vec, active_parts=[1,2,3])
              → (selected_forest, presence_mask)

Select a subset of trees from `original_f` using the first `n_trees` bits of
`original_vec` as a presence mask.

A bit value of `1` includes the corresponding tree; `0` drops it entirely.

# Arguments
- `n_trees::Int`      : Number of trees in the forest.
- `original_f`        : Full `DecisionEnsemble` from `learn_and_convert`.
- `original_vec`      : Bitvector of length `length(active_parts) * n_trees`.
- `active_parts`      : Which stages are active (used to locate the correct segment).

# Returns
- `selected_forest`   : New `DecisionEnsemble` with only the present trees.
- `presence_mask`     : Boolean vector of length `n_trees` (first N bits).
"""
function first_part(n_trees, original_f, original_vec, active_parts=[1,2,3])

    # Extract the presence-segment bits: each bit indicates whether the
    # corresponding tree should be included (1) or excluded (0) from the new forest
    presence_mask = original_vec[_segment(1, n_trees, active_parts)]

    # Keep only the models (SoleTrees) whose presence bit is 1
    remaining_trees = original_f.models[presence_mask]

    # Build a new forest with only the selected trees
    selected_forest = DecisionEnsemble(remaining_trees)

    return selected_forest, presence_mask
end

# ─────────────────────────────────────────────────────────────────────────────
# STAGE 2 — TREE PRUNING
# ─────────────────────────────────────────────────────────────────────────────

"""
    second_part(n_trees, presence_mask, original_vec,
                active_parts=[1,2,3]; silent=true)
              → (pruned_forest, target_depth, depth_mask)

Prune the surviving trees using the depth-bit segment of `original_vec`.
When stage 1 was skipped `presence_mask` should be `trues(n_trees)`.

*Reference trees* have both presence bit = 1 and depth bit = 1. Their maximum
depth defines the `target_depth` budget for the whole forest. All other present
trees are pruned down to `target_depth` if they exceed it.

# Arguments
- `n_trees::Int`      : Number of trees in the forest.
- `presence_mask`     : Boolean vector returned by `first_part`.
                        (needed for raw node access and depth measurement).
- `original_vec`      : Bitvector of length `length(active_parts) * n_trees`.
- `active_parts`      : Which stages are active (used to locate the correct segment).

# Returns
- `pruned_forest`     : Final `DecisionEnsemble`, each tree ≤ `target_depth` deep.
- `target_depth`      : The depth ceiling that was enforced.
- `depth_mask`        : Boolean vector of length `n_trees` (depth segment bits).
"""

function second_part(n_trees, presence_mask, original_f, original_vec, active_parts=[1,2,3]; silent=true)

    # ── Step 1: extract the depth-segment bits (depth / pruning flags) ────────
    depth_mask = original_vec[_segment(2, n_trees, active_parts)]

    # ── Step 2: determine the target depth ───────────────────────────────────
    # "Reference" trees are those with presence = 1 AND depth_flag = 1.
    # They represent the unconstrained trees and define the allowed depth budget.

    interesting_trees_indeces = findall(presence_mask .& depth_mask)

    if isempty(interesting_trees_indeces)
        # Edge case: no reference tree exists; default to a minimal depth of 1.
        # This is a fallback and may produce a poor forest — a different
        # default strategy could be explored here.

        present_indices = findall(presence_mask)

        if isempty(present_indices)
            target_depth = 42
        else
            # Measure the actual depth of each reference tree in the original forest
            ref_depths = [tree_depth(original_f.models[i]) for i in present_indices]
            target_depth = maximum(ref_depths)
        end
    else
        ref_depths = [tree_depth(original_f.models[i]) for i in interesting_trees_indeces]

        # Use the deepest reference tree as the pruning ceiling

        target_depth = maximum(ref_depths)
        silent || println("target_depth = ", target_depth)
    end

    # ── Step 3: build the final list of trees ────────────────────────────────
    final_models = AbstractModel[]

    for i in 1:n_trees

        # Skip trees that were excluded in the presence step
        if presence_mask[i] == false
            continue
        end

        # At this point presence_mask[i] is guaranteed to be true.

        if depth_mask[i] == false
            # This tree is present but NOT a reference tree: it must be pruned
            # if it exceeds the target depth

            curr_depth = tree_depth(original_f.models[i])
            if curr_depth > target_depth
                
                # The tree is too deep → prune it.
                # `prune_tree` starts at depth 1 (root) and cuts at `target_depth`.
                pruned_tree = prune_tree(original_f.models[i], 1, target_depth)
                # Pushing the pruned tree into the final pruned forest
                
                push!(final_models, pruned_tree)
            else
               
                # The tree is already within the depth budget, no pruning needed.
                push!(final_models, original_f.models[i])
            end

        else
            # presence = 1, depth_flag = 1 → reference tree, keep it intact.
           
            push!(final_models, original_f.models[i])
        end
    end

    
    
    # Fix del convert: gestisce mix di Node e Leaf
    # ── Step 4: wrap the final list into a new DecisionEnsemble ──────────────
    pruned_forest = DecisionEnsemble(final_models)

    return pruned_forest, target_depth, depth_mask
end


# ─────────────────────────────────────────────────────────────────────────────
# THIRD PART OF THE BITVECTOR: ALPHABET MODIFICATION
# ─────────────────────────────────────────────────────────────────────────────



"""
    third_part(n_trees, presence_mask, pruned_forest, original_vec,
               active_parts=[1,2,3])
             → (sole_final_forest, alphabet_table, alphabet_mask)

Modify the alphabet of the trees using the alphabet-bit segment of `original_vec`.
When stage 1 was skipped `presence_mask` should be `trues(n_trees)`.

The modification steps are the following:

1) Creation of `alphabet_table` by collecting each unique feature in each node
   of each tree in the forest with its alphabet bit = 1, and mapping it to the list
   of values that feature operates with in any node in any tree.

2) After creating `alphabet_table` using the trees with alphabet bit = 1,
   every node in each tree with alphabet bit = 0 is analized.
   Let's considerate a certain node N talking about a feature F_k(a). There are 3 scenarios:

    2.1) In the alphabet table there is NOT `F_k` -> the tree gets pruned and the node `N`
         becomes a leaf.

    2.2) In the alphabet table there is `F_k` but in its corrisponding list there is not
         the value `a` -> in the node `N`, the value `a` gets replaced by the value inside
         the `F_k` feature's list in the alphabet table that's closest to `a`.

    2.3) In the alphabet table there is `F_k` and in its corrisponding list there is
         `a` -> nothing must be done and the node will remain unchanged.
---------------------------------------------------------------------
# Arguments
- `n_trees::Int`      : Number of trees in the forest.
- `presence_mask`     : Boolean vector returned by `first_part`.
- `pruned_forest`     : `DecisionEnsemble` returned by `second_part`.
- `original_vec`      : Bitvector of length `length(active_parts) * n_trees`.
- `active_parts`      : Which stages are active (used to locate the correct segment).

# Returns
- `sole_final_forest` : Final forest in the form of `DecisionEnsemble`,
                        with its trees' alphabet modified.
- `alphabet_table`    : A dictionary that maps each unique feature to a
                        list of all the values they operates with.
- `alphabet_mask`     : Boolean vector of length n_trees.
"""
function third_part(n_trees, presence_mask, pruned_forest, original_vec,
                    active_parts=[1,2,3])

    # I create the mask I'll use to determine if i have to
    # modify the current tree
    alphabet_mask = original_vec[_segment(3, n_trees, active_parts)]

    #=
     In the first column of this table I'll put every feature
     that appears in any node of any tree with alphabet bit = 1.

     Then each feature's row will contain the values that
     feature operates with in any node of any tree.

     After being filled with values this table will be used
     while visiting the trees with alphabet bit = 0.
    =#
    alphabet_table = Dict{Any, Vector{Float64}}()

    #=
     Assumption: n_trees is strictly greater than the number
     of trees in pruned_forest.

     This assumption prevents me from just using n_trees as index.
     To solve this problem active_tree_index will be the key.

     This special index will be the one used to access the trees in
     the forest, while n_trees will be used to iterate the for loop.

     This ensures that there will be no out of bound error.
    =#
    active_tree_index = 1

    # First for loop, where the goal is to fill the alphabet
    # table with all the features and their operating values.
    for i in 1:n_trees

        # I check if the currently analized tree is in the forest
        if presence_mask[i] == false

            # If the current tree is not in the forest, i just Skip
            # this iteration of the cycle.
            #
            # Incrementing i but not active_tree_index

            continue
        end

        # I select the current tree using active_tree_index
        current_tree = pruned_forest.models[active_tree_index]
        # I am interested in only the trees with the alphabet bit = 1
        if alphabet_mask[i] == true

            # I use the collect_alphabet function to fill the table
            # with the features and operating values of the current
            # tree's nodes.
            collect_alphabet!(current_tree, alphabet_table)
        end

        # I increment active_tree_index.
        #
        # Note that when I increment this variable
        # only when the presence bit = 1
        active_tree_index += 1
    end

    # After filling the table with all the values
    # I remove any duplicate and sort everything.
    for (feat, vals) in alphabet_table

        unique_vals = unique(vals)
        sort!(unique_vals)
        alphabet_table[feat] = unique_vals
    end

    #------------------------------------------------

    # Now that I have my alphabet table, it's time
    # to modify each tree with alphabet bit = 0.
    # To do that i need to reset the special index.
    active_tree_index = 1

    # This will be the forest this function will return
    # I will put in here every modified tree and every
    # tree that i used to gather the alphabet
    final_forest = AbstractModel[]

    for i in 1:n_trees

        # As before, i ignore any tree that has
        # presence bit = 0
        if presence_mask[i] == false

            continue
        end

        # I select the current tree using active_tree_index as before
        current_tree = pruned_forest.models[active_tree_index]

        # I want to modify only the trees with alphabet bit = 0
        if alphabet_mask[i] == false

            #=
             Using the modify_alphabet function I visit the tree
             node by node, checking if the current node's feature
             is in the alphabet table.

             If the feature is NOT present I prune the tree and
             transform the current node into a leaf.

             If the feature is present I check if the value it's
             operating with is present in that feature's row in
             the alpahbet table.

             If the value is present I simply do nothing.

             If the value is not present i replace the node's
             value with the value that is closest in the feature's
             row, in the alphabet table
            =#
            modified_tree = modify_alphabet(current_tree, alphabet_table)

            push!(final_forest, modified_tree)
        else

            # If the tree has not been modified I insert it
            # in the final forest as is.
            push!(final_forest, current_tree)
        end

        # I increment the active_tree_index only with trees with
        # presence bit = 1
        active_tree_index += 1

    end

    # I transform the standard Julia ensemble into a solemodel
    sole_final_forest = DecisionEnsemble(final_forest)

    # I return the forest
    # and the alphabet table and alphabet mask as extra information
    return sole_final_forest, alphabet_table, alphabet_mask
end


# ─────────────────────────────────────────────────────────────────────────────
# STAGE 4 — EVOLUTIONARY OPTIMISATION OF THE BITVECTOR
# ─────────────────────────────────────────────────────────────────────────────

"""
    evaluate_bitvector(bitvec, n_trees, original_f, f_val, l_val;
                       active_parts=[1,2,3]) → (error_rate, n_ones)

Apply only the stages listed in `active_parts` and measure the classification
error on the validation set `(f_val, l_val)`.

The bitvector length must be `length(active_parts) * n_trees`.

Stages that are skipped receive a default pass-through treatment:
- Stage 1 skipped → all trees are kept (`presence_mask = trues(n_trees)`)
- Stage 2 skipped → no pruning is applied
- Stage 3 skipped → alphabet is left unchanged

# Returns
- `error_rate` : fraction of misclassified validation samples  (0.0 = perfect).
- `n_ones`     : number of `1` bits in `bitvec`  (lower = more compressed).
"""
function evaluate_bitvector(bitvec::BitVector, n_trees::Int,
                             original_f, f_val::Matrix{Float64}, 
                             l_val::Vector{String};
                             active_parts::Vector{Int} = [1,2,3])

    # ── Stage 1: tree selection ───────────────────────────────────────────────
    if 1 in active_parts
        selected_forest, presence_mask = first_part(n_trees, original_f, bitvec, active_parts)

        # Edge-case: if no tree was selected the forest is useless → maximum penalty
        if isempty(selected_forest.models)
            return 1.0, sum(bitvec)
        end
    else
        # Stage 1 skipped: keep all trees
        selected_forest = original_f
        presence_mask   = trues(n_trees)
    end

    # ── Stage 2: tree pruning ─────────────────────────────────────────────────
    if 2 in active_parts
        pruned_forest, _target_depth, _depth_mask =
            second_part(n_trees, presence_mask, original_f, bitvec, active_parts)
    else
        # Stage 2 skipped: build a forest from the kept trees (no pruning)
        kept_indices  = findall(presence_mask)
        pruned_forest = DecisionEnsemble(original_f.models[kept_indices])
    end

    # ── Stage 3: alphabet modification ───────────────────────────────────────
    if 3 in active_parts
        final_forest, _alphabet_table, _alphabet_mask =
            third_part(n_trees, presence_mask, pruned_forest, bitvec, active_parts)
    else
        # Stage 3 skipped: use as-is 
        final_forest = pruned_forest
    end

    # ── Accuracy evaluation on the validation set ────────────────────────────
    # SoleModels.apply requires a DataFrame, not a raw Matrix{Float64}.
    df_val      = DataFrame(f_val, :auto)
    predictions = string.(SoleModels.apply(final_forest, df_val))
    n_correct   = sum(predictions .== l_val)
    error_rate  = 1.0 - (n_correct / length(l_val))

    return error_rate, sum(bitvec)
end


"""
    core(n_trees, original_f, f_val, l_val;
         population_size, n_generations, penalty_weight,
         active_parts=[1,2,3]) → (best_bitvector, best_cost)

Run a Genetic Algorithm to find the `BitVector` of length
`length(active_parts) * n_trees` that simultaneously minimises the
classification error on `(f_val, l_val)` and the number of active bits (`1`s),
running only the stages listed in `active_parts`.

The fitness to minimise is:

    fitness = error_rate + penalty_weight * (n_ones / n_bits)

`penalty_weight` trades compression against accuracy:
- higher → aggressively prune (may hurt accuracy)
- lower  → preserve accuracy (less compression)

# Arguments
- `n_trees`         : Number of trees in the original forest.
- `original_f`      : `DecisionEnsemble` from `learn_and_convert`.
- `f_val`           : Validation feature matrix  (n_samples × n_features).
- `l_val`           : Validation label vector    (n_samples,).
- `population_size` : Number of individuals per generation  (default: 50).
- `n_generations`   : Number of GA iterations               (default: 100).
- `penalty_weight`  : Weight of the compression penalty     (default: 0.3).
- `active_parts`    : Which stages to run — any non-empty subset of {1,2,3} (default: [1,2,3]).

# Returns
- `best_bitvector` : The winning `BitVector` found by the GA.
- `best_cost`      : Fitness value of `best_bitvector`.
"""
function core(
        n_trees::Int,
        original_f,
        f_val::Matrix{Float64},
        l_val::Vector{String};
        population_size::Int      = 50,
        n_generations::Int        = 100,
        penalty_weight::Float64   = 0.3,
        active_parts::Vector{Int} = [1,2,3])

    # The bitvector only contains segments for the active stages
    n_bits = length(active_parts) * n_trees

    # ── Fitness closure (captures all data needed for evaluation) ─────────────
    #=
      We minimise:
        fitness = error_rate + penalty_weight * (n_ones / n_bits)

      Normalising n_ones by n_bits keeps both terms in [0, 1], making
      penalty_weight directly interpretable as the relative importance
      of compression vs. accuracy.
    =#
    function fitness_fn(bitvec::BitVector)
        error_rate, n_ones = evaluate_bitvector(
            bitvec, n_trees, original_f, f_val, l_val;
            active_parts = active_parts)
        compression_penalty = penalty_weight * (n_ones / n_bits)
        return error_rate + compression_penalty
    end

    # ── Initial population: random bitvectors ────────────────────────────────
    initial_population = [BitVector(rand(Bool, n_bits)) for _ in 1:population_size]

    # ── GA configuration ─────────────────────────────────────────────────────
    algo = GA(
        populationSize = population_size,
        selection      = tournament(3),      # tournament selection, k=3
        crossover      = UX,                 # uniform crossover — works well for BitVectors
        mutationRate   = 1.0 / n_bits        # expected one bit flip per individual
    )

    opts = Evolutionary.Options(
        iterations = n_generations,
        show_trace = true,
        show_every = max(1, n_generations ÷ 10)   # print ~10 progress lines
    )

    # ── Run optimisation ─────────────────────────────────────────────────────
    result = Evolutionary.optimize(fitness_fn, initial_population, algo, opts)

    best_bitvector = Evolutionary.minimizer(result)
    best_cost      = Evolutionary.minimum(result)

    return best_bitvector, best_cost
end
