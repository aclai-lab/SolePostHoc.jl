# ---------------------------------------------------------------------------- #
#                                    Cbc                                       #
# ---------------------------------------------------------------------------- #
"""
    cbc(ruleset, X, y, config::CBC) -> Vector{<:Rule}

Select a compact subset of `ruleset` via the Complexity-guided Boolean
Combination (CBC) procedure: rule-coverage checkmasks are used as features to
fit a random forest against `y`, and rules whose normalized impurity
importance exceeds `get_cbc_threshold(config)` are kept, sorted by
importance (desc), error (asc), and complexity (asc).
"""
function cbc(ruleset, X, y, config::AbstractConfig)
    n_rules = length(ruleset)
    metrics = Matrix{Float64}(undef, n_rules, 3)
    checkmasks = Vector{BitVector}(undef, n_rules)

    Threads.@threads for i in eachindex(ruleset)
        eval_result = rulemetrics(ruleset[i], X, y)
        checkmasks[i] = eval_result[:checkmask,]
        metrics[i, 1] = eval_result[:coverage]
        metrics[i, 2] = eval_result[:error]
        metrics[i, 3] = eval_result[get_complexity_metric(config)]
    end

    # build random forest for feature importance
    rf = DT.build_forest(
        y, hcat(checkmasks...),
        get_nsubfeatures(config),
        get_ntrees(config),
        get_partial_sampling(config),
        get_max_depth(config);
        rng=get_rng(config))
    importance = DT.impurity_importance(rf)
    importances = importance ./ maximum(importance)

    # select features with sufficient importance
    selected_idxs = findall(importances .> get_threshold(config))
    isempty(selected_idxs) && return ruleset

    # combine metrics with importance and original indices
    combined = hcat(
        metrics[selected_idxs, :],
        importances[selected_idxs],
        selected_idxs
    )

    # sort by importance (desc), error (asc), complexity (asc)
    sorted = sortslices(combined, dims=1, by=x->(-x[4], x[2], x[3]))

    # extract final indices, limiting if max_rules is set
    max_rules = get_max_rules(config)
    n_selected = max_rules > 0 ? min(max_rules, size(sorted, 1)) : size(sorted, 1)
    final_idxs = Int64.(sorted[1:n_selected, 5])

    return ruleset[final_idxs]
end