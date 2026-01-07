module RULECOSIPLUS

using PyCall
using Conda
using DataFrames
using CSV
using Random
using JSON
using SoleModels      # Assicurati che i moduli SoleModels e SoleLogics siano configurati
using SoleLogics
import Dates: now
using DataStructures
using CategoricalArrays

include("apiRuleCosi.jl")


export rulecosiplus

if !@isdefined rulecosi
    const rulecosi = PyNULL()
end
if !@isdefined sklearn
    const sklearn = PyNULL()
end

function __init__()
    # First ensure pip interop is enabled
    Conda.pip_interop(true)
    
    # Try to import rulecosi, if it fails, install it via pip
    try
        copy!(rulecosi, pyimport("rulecosi"))
    catch
        @info "Installing rulecosi via pip..."
        Conda.pip("install", "git+https://github.com/jobregon1212/rulecosi.git")
        copy!(rulecosi, pyimport("rulecosi"))
    end

    # Try to import sklearn, if it fails, install it via pip
    try
        copy!(sklearn, pyimport("sklearn.ensemble"))
    catch
        @info "Installing scikit-learn via pip..."
        PyCall.Conda.pip("install", "scikit-learn", PyCall.Conda.ROOTENV)
        copy!(sklearn, pyimport("sklearn.ensemble"))
    end

    py"""
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.tree._tree import Tree
    from collections import Counter
    import io, sys

    def build_sklearn_tree(tree_state):
        n_nodes = tree_state["node_count"]
        n_features = tree_state["n_features"]
        n_classes = tree_state["n_classes"]
        t = Tree(n_features, np.array([n_classes], dtype=np.intp), 1)
        dt = np.dtype([
            ('left_child','<i8'),
            ('right_child','<i8'),
            ('feature','<i8'),
            ('threshold','<f8'),
            ('impurity','<f8'),
            ('n_node_samples','<i8'),
            ('weighted_n_node_samples','<f8'),
            ('missing_go_to_left','u1')
        ])
        cl = tree_state["children_left"]
        cr = tree_state["children_right"]
        ft = tree_state["feature"]
        th = tree_state["threshold"]
        v_orig = tree_state["value"]
        v = np.ascontiguousarray(v_orig, dtype=np.float64)
        nodes_list = []
        for i in range(n_nodes):
            lc = cl[i]
            rc = cr[i]
            fe = ft[i]
            tr = th[i]
            counts = v[i,0,:]
            tot = int(np.sum(counts))
            imp = 0.0
            nodes_list.append((lc, rc, fe, tr, imp, tot, float(tot), 0))
        structured_nodes = np.array(nodes_list, dtype=dt)
        st = {"max_depth": 20, "node_count": n_nodes, "nodes": structured_nodes, "values": v}
        t.__setstate__(st)
        return t

    class RealDecisionTreeClassifier(DecisionTreeClassifier):
        def __init__(self, tree_state, classes):
            super().__init__()
            self.tree_ = build_sklearn_tree(tree_state)
            self.n_features_ = tree_state["n_features"]
            self.n_outputs_ = 1
            self.n_classes_ = tree_state["n_classes"]
            self.classes_ = np.array(classes)
            self.fitted_ = True
        def fit(self, X, y=None):
            return self

    class RealRandomForestClassifier(RandomForestClassifier):
        def __init__(self, forest_json):
            super().__init__(n_estimators=0)
            self.classes_ = np.array(forest_json["classes"])
            self.estimators_ = []
            for tree_state in forest_json["trees"]:
                self.estimators_.append(RealDecisionTreeClassifier(tree_state, forest_json["classes"]))
            self.fitted_ = True
        def fit(self, X, y=None):
            return self
        def predict(self, X):
            arr = np.array([est.predict(X) for est in self.estimators_])
            final = []
            for i in range(X.shape[0]):
                c = Counter(arr[:, i]).most_common(1)[0][0]
                final.append(c)
            return np.array(final)

    def build_sklearn_model_from_julia(dict_model):
        return RealRandomForestClassifier(dict_model)

    def get_simplified_rules(rc, heuristics_digits=4, condition_digits=1):
        buffer = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buffer
        try:
            rc.simplified_ruleset_.print_rules(heuristics_digits=heuristics_digits, condition_digits=condition_digits)
        finally:
            sys.stdout = old_stdout
        rules_str = buffer.getvalue()
        return rules_str.strip().splitlines()
    """
end


##############################
# 1) Struct & BFS/DFS
##############################

mutable struct SkNode
    left::Int
    right::Int
    feature::Int
    threshold::Float64
    counts::Vector{Float64}
end

const TreeType = Union{String,CategoricalArrays.CategoricalValue{String,UInt32}}

"""
build_sklearnlike_arrays(branch, n_classes, class_to_idx)

Performs post-order DFS on a Branch{String} (or ConstantModel{String}),
creating children_left, children_right, feature, threshold, value, n_nodes.
"""
function build_sklearnlike_arrays(branch, n_classes::Int, class_to_idx::Dict{String,Int})
    nodes = SkNode[]
    function dfs(b)
        i = length(nodes)
        push!(nodes, SkNode(-1, -1, -2, -2.0, fill(0.0, n_classes)))
        if b isa ConstantModel{<:TreeType}
            cidx = haskey(class_to_idx, b.outcome) ? class_to_idx[b.outcome] : 0
            nodes[i+1].counts[cidx+1] = 10.0
        elseif b isa Branch{<:TreeType}
            thr = b.antecedent.value.threshold
            fx = b.antecedent.value.metacond.feature.i_variable - 1
            left_i = dfs(b.posconsequent)
            right_i = dfs(b.negconsequent)
            nodes[i+1].feature = fx
            nodes[i+1].threshold = thr
            nodes[i+1].left = left_i
            nodes[i+1].right = right_i
            for c = 1:n_classes
                nodes[i+1].counts[c] =
                    nodes[left_i+1].counts[c] + nodes[right_i+1].counts[c]
            end
        else
            error("Unknown node type: $(typeof(b))")
        end
        return i
    end
    dfs(branch)
    n_nodes = length(nodes)
    children_left = Vector{Int}(undef, n_nodes)
    children_right = Vector{Int}(undef, n_nodes)
    feature = Vector{Int}(undef, n_nodes)
    threshold = Vector{Float64}(undef, n_nodes)
    value = Array{Float64,3}(undef, n_nodes, 1, n_classes)
    for i = 1:n_nodes
        children_left[i] = nodes[i].left
        children_right[i] = nodes[i].right
        feature[i] = nodes[i].feature
        threshold[i] = nodes[i].threshold
        for c = 1:n_classes
            value[i, 1, c] = nodes[i].counts[c]
        end
    end
    return children_left, children_right, feature, threshold, value, n_nodes
end

function serialize_branch_sklearn(
    b,
    classes::Vector{String},
    class_to_idx::Dict{String,Int},
)
    n_classes = length(classes)
    cl, cr, ft, th, val, nn = build_sklearnlike_arrays(b, n_classes, class_to_idx)
    posfeats = ft[ft .>= 0]
    maxfeat = isempty(posfeats) ? 0 : maximum(posfeats)
    n_features = maxfeat + 1
    return Dict(
        "children_left" => cl,
        "children_right" => cr,
        "feature" => ft,
        "threshold" => th,
        "value" => val,
        "node_count" => nn,
        "n_features" => n_features,
        "n_classes" => n_classes,
    )
end

"""
serialize_julia_ensemble(ensemble, classes)

Creates a dictionary with:
  "classes" => classes,
  "trees"   => [<tree_state1>, <tree_state2>, ...]
"""
function serialize_julia_ensemble(ensemble, classes::Vector{String})
    class_to_idx = Dict{String,Int}()
    for (i, cl) in enumerate(classes)
        class_to_idx[cl] = i - 1
    end
    trees = [serialize_branch_sklearn(b, classes, class_to_idx) for b in ensemble.models]
    return Dict("classes" => classes, "trees" => trees)
end

##############################
# 2) BLOCK PYTHON
##############################



##############################
# 3) FUN process_rules_decision_list
##############################
"""
process_rules_decision_list(rules::Vector{String})

Processes raw rules from RuleCOSI into a vector of tuples (condition, outcome).
The last rule with empty condition is treated as the default rule.
"""
function process_rules_decision_list(rules::Vector{String})
    structured = Tuple{String,String}[]
    default_outcome = nothing

    for line in rules
        if occursin("cov", line)
            continue
        end

        parts = split(line, ":")
        if length(parts) < 2
            continue
        end

        rule_part = strip(parts[2])
        parts2 = split(rule_part, "→")
        if length(parts2) < 2
            continue
        end

        cond = strip(parts2[1], ['(', ')', ' '])
        # Fix parentheses and operators
        cond = replace(cond, "˄" => " ∧ ")
        # Ensure balanced parentheses
        if !startswith(cond, "(") && occursin("∧", cond)
            cond = "(" * cond * ")"
        end

        out = strip(parts2[2])
        out = replace(out, "[" => "")
        out = replace(out, "]" => "")
        out = replace(out, "'" => "")

        if isempty(cond) || cond == "( )"
            default_outcome = out
            continue
        end

        push!(structured, (cond, out))
    end

    return structured, default_outcome
end

function build_rule_list(
    decision_list::Vector{Tuple{String,String}},
    default_outcome::String,
)
    rules = Rule[]

    for (cond_str, out) in decision_list
        try
            # Ensure proper spacing around operators
            cond_str = replace(cond_str, ">" => " > ")
            cond_str = replace(cond_str, "≤" => " ≤ ")
            cond_str = replace(cond_str, "∧" => " ∧ ")

            φ = SoleLogics.parseformula(
                cond_str;
                atom_parser = a -> Atom(
                    parsecondition(
                        ScalarCondition,
                        a;
                        featuretype = VariableValue,
                        featvaltype = Real,
                    ),
                ),
            )

            push!(rules, Rule(φ, out))
        catch e
            println("Warning: Failed to parse rule: $cond_str")
            println("Error: ", e)
            continue
        end
    end

    return rules, default_outcome
end

function convertToDecisionList(raw_rules::Vector{String})
    # Process the raw rules and extract the default outcome
    structured_rules, default_outcome = process_rules_decision_list(raw_rules)

    # Build the rule list and get the default outcome
    rules, default = build_rule_list(structured_rules, default_outcome)
    #rules = dnf(rules)
    # Create the DecisionList
    return DecisionList(rules, default)
end

##############################
# 4) Build ds on dnf
##############################

# Use the Rule type from SoleModels to build a DecisionSet.
struct MyRule
    formula::Formula   # The parsed symbolic formula
    outcome::String
end

##############################
# 5) FUNCTION __init__
##############################
"""
    rulecosiplus(ensemble::Any, X_train::Any, y_train::Any)

Extract interpretable rules from decision tree ensembles using the RuleCOSI+ algorithm.

This function implements the RuleCOSI+ methodology for rule extraction from trained ensemble
classifiers, producing a simplified and interpretable rule-based model. The method combines
and simplifies rules extracted from individual trees in the ensemble to create a more
compact and understandable decision list.

# Reference
Obregon, J. (2022). RuleCOSI+: Rule extraction for interpreting classification tree ensembles.
*Information Fusion*, 89, 355-381.
Available at: https://www.sciencedirect.com/science/article/pii/S1566253522001129

# Arguments
- `ensemble::Any`: A trained ensemble classifier (e.g., Random Forest, Gradient Boosting)
  that will be serialized and converted to a compatible format for rule extraction.
- `X_train::Any`: Training feature data. Can be a DataFrame or Matrix. If DataFrame,
  column names will be preserved in the extracted rules; otherwise, generic names (V1, V2, ...)
  will be generated.
- `y_train::Any`: Training target labels corresponding to `X_train`. Will be converted to
  string format for processing.

# Returns
- `DecisionList`: A simplified decision list containing the extracted and combined rules
  from the ensemble, suitable for interpretable classification.

# Details
The function performs the following steps:
1. Converts input data to appropriate matrix format
2. Generates or extracts feature column names
3. Serializes the Julia ensemble to a Python-compatible format
4. Builds an sklearn-compatible model using the serialized ensemble
5. Applies RuleCOSI+ algorithm with the following default parameters:
   - `metric="fi"`: Optimization metric for rule combination
   - `n_estimators=100`: Number of estimators considered
   - `tree_max_depth=100`: Maximum depth of trees
   - `conf_threshold=0.25` (α): Confidence threshold for rule filtering
   - `cov_threshold=0.1` (β): Coverage threshold for rule filtering
   - `verbose=2`: Detailed output during processing
6. Extracts and converts rules to a decision list format

# Configuration
The algorithm uses fixed parameters optimized for interpretability:
- Confidence threshold (α) = 0.25: Rules below this confidence are discarded
- Coverage threshold (β) = 0.1: Rules covering fewer samples are excluded
- Maximum rules = max(20, n_classes × 5): Adaptive limit based on problem complexity

# Example
```julia
# Assuming you have a trained ensemble and training data
ensemble = ... # your trained ensemble
X_train = ... # training features
y_train = ... # training labels

# Extract interpretable rules
decision_list = rulecosiplus(ensemble, X_train, y_train)
```

# Notes
- The function prints diagnostic information including the number of trees and dataset statistics
- Raw rules are displayed before conversion to decision list format
- Requires Python interoperability and the RuleCOSI implementation
- The resulting decision list provides an interpretable alternative to the original ensemble
"""
function rulecosiplus(ensemble::Any, X_train::Any, y_train::Any; silent::Bool = true)

    # Convert training data to matrix format if needed
    if !isa(X_train, Matrix)
        X_train_matrix = Matrix(X_train)
    else
        X_train_matrix = X_train
    end

    # Generate column names
    if isa(X_train, DataFrame)
        column_names = names(X_train)
    else
        num_features = size(X_train_matrix, 2)
        column_names = ["V$i" for i = 1:num_features]
    end

    # Serialize the ensemble
    classes = map(String, collect(unique(y_train)))
    dict_model = serialize_julia_ensemble(ensemble, classes)

    # Build sklearn model
    builder = py"build_sklearn_model_from_julia"
    base_ensemble = builder(dict_model)

    silent || println("======================")
    silent || println("Ensemble serialize:", dict_model)
    silent || println("======================")

    num_estimators = pycall(pybuiltin("len"), Int, base_ensemble["estimators_"])
    silent || println("number of trees in base_ensemble:", num_estimators)

    n_samples = size(X_train_matrix, 1)
    n_classes = length(unique(y_train))

    silent || println("Dataset: $n_samples campioni, $n_classes classi")

    #=
            RuleCOSIExtractor

        Extracts, combines, and simplifies rules from decision tree ensembles,
        producing a single interpretable rule-based model.
        Supports various types of binary ensembles (e.g., RandomForest, GradientBoosting, etc.)
        and returns rules suitable for classification.

        Parameters
        ----------
        base_ensemble : BaseEnsemble or None, optional (default=None)
            Instance of an already trained or to be trained ensemble. Supported:
            - sklearn.ensemble.RandomForestClassifier
            - sklearn.ensemble.BaggingClassifier
            - sklearn.ensemble.GradientBoostingClassifier
            - xgboost.XGBClassifier
            - catboost.CatBoostClassifier
            - lightgbm.LGBMClassifier
            If None, a default ensemble (e.g., GradientBoostingClassifier) will be created (and trained).

        metric : str, optional (default="f1")
            Metric used to optimize the combination/simplification of rules.
            Common options:
            - "f1"
            - "roc_auc"
            - "accuracy"

        n_estimators : int, optional (default=5)
            Number of estimators (trees) to train if `base_ensemble` is not already trained.

        tree_max_depth : int, optional (default=3)
            Maximum depth of the trees within the ensemble.
            A higher value allows more complex trees and more rules.

        conf_threshold : float, optional (default=0.5)
            Confidence threshold below which rules are discarded during the combination phase.
            The higher the threshold, the fewer the final rules.

        cov_threshold : float, optional (default=0.0)
            Coverage threshold: excludes rules that cover fewer samples than indicated.
            With 0.0, only rules that cover no samples are discarded.

        c : float, optional (default=0.25)
            Statistical confidence level used for rule error estimation (Laplace correction or similar).
            Helps to avoid overfitting in rules.

        percent_training : float or None, optional (default=None)
            Percentage of the training dataset to use for rule combination/selection.
            If None, uses all data.

        early_stop : float, optional (default=0)
            If > 0, allows the combination algorithm to stop if no improvements in the metric
            indicated in `metric` are observed for a certain number of iterations (equal to `n_estimators * early_stop`).

        rule_order : str, optional (default="supp")
            Initial sorting criterion for rules (in the iterative combination phase).
            Typical options:
            - "cov" (sort by coverage),
            - "conf" (sort by confidence),
            - "supp" (sort by support).

        sort_by_class : bool or list or None, optional (default=None)
            If True, sorts the final rule set by class (in lexicographic order).
            If a list, the class order is as specified by the user.
            If None, does not enforce class sorting.

        column_names : list of str or None, optional (default=None)
            Names of the dataset columns/features. If provided, the extracted rules will use
            these names instead of numerical indices.

        random_state : int, RandomState or None, optional (default=None)
            Seed to control the reproducibility of some stochastic processes
            (e.g., training an ensemble if not already trained).

        verbose : int, optional (default=0)
            Verbosity level of the output:
            - 0 = no output
            - 1 = general prints
            - 2 = detailed information of each combination iteration

        df : pandas.DataFrame or None, optional (default=None)
            If you directly provide a DataFrame (e.g., of features), the class
            can use it to extract `column_names` or other information.
            If None, it is assumed that the features (X) and target vector (y) will be
            passed during `fit` as separate arrays/series.

        Attributes
        ----------
        X_ : numpy.ndarray of shape (n_samples, n_features)
            Training features used during :meth:`fit`.

        y_ : numpy.ndarray of shape (n_samples,)
            Training labels/targets used during :meth:`fit`.

        classes_ : numpy.ndarray of shape (n_classes,)
            Classes encountered in :meth:`fit`.

        original_rulesets_ : list of RuleSet, shape (n_estimators,)
            Set of rules extracted directly from the base ensemble,
            before combination/simplification.

        simplified_ruleset_ : RuleSet
            Final rule set, combined and simplified.

        n_combinations_ : int
            Number of iterations or "steps" of combination performed by the algorithm.

        combination_time_ : float
            Time spent in the extraction, combination, and simplification phase.

        ensemble_training_time_ : float
            Time spent training the ensemble (if not already trained).
            If the ensemble was already ready, remains 0.

        References
        ----------
        .. [1] Obregon, J., Kim, A., & Jung, J. Y.,
            "RuleCOSI: Combination and simplification of production rules from boosted
            decision trees for imbalanced classification", 2019.

        Examples
        -------
        >>> from sklearn.datasets import make_classification
        >>> X, y = make_classification(n_samples=1000, n_features=4, random_state=0)
        >>> extractor = RuleCOSIExtractor(n_estimators=10, metric='roc_auc')
        >>> extractor.fit(X, y)
        RuleCOSIExtractor(n_estimators=10, metric='roc_auc')
        >>> y_pred = extractor.predict(X)
        >>> extractor.score(X, y)
        0.95
    =#
    rc = rulecosi.RuleCOSIClassifier(
        base_ensemble = base_ensemble,
        metric = "fi",
        n_estimators = 100,
        tree_max_depth = 100,
        conf_threshold = 0.25, # α
        cov_threshold = 0.1, # β
        random_state = 3,
        column_names = column_names,
        verbose = 2,
    )

    @time rc.fit(X_train_matrix, Vector(String.(y_train)))


    max_rules = max(20, n_classes * 5)

    raw_rules = py"get_simplified_rules"(rc, 4, 1)
    silent || println("Raw rules:")

    if silent == false
        for r in raw_rules
            println(r)
        end
    end

    dl = convertToDecisionList(raw_rules)

    return dl
end

end

# Required modules (adjust paths as needed)
#include("apiRuleCosi.jl")
