module RULECOSIPLUS
__precompile__()

using PyCall
using Conda
using DataFrames
using CSV
using Random
using JSON
using SoleModels      # Assicurati che i moduli SoleModels e SoleLogics siano configurati
using SoleLogics

include("apiRuleCosi.jl")


export rulecosiplus

if !@isdefined rulecosi
    const rulecosi = PyNULL()
end
if !@isdefined sklearn
    const sklearn = PyNULL()
end

##############################
# 1) Struttura e BFS/DFS
##############################

mutable struct SkNode
    left::Int
    right::Int
    feature::Int
    threshold::Float64
    counts::Vector{Float64}
end

"""
build_sklearnlike_arrays(branch, n_classes, class_to_idx)

Esegue DFS post-order su un Branch{String} (o ConstantModel{String}),
creando children_left, children_right, feature, threshold, value, n_nodes.
"""
function build_sklearnlike_arrays(branch, n_classes::Int, class_to_idx::Dict{String,Int})
    nodes = SkNode[]
    function dfs(b)
        i = length(nodes)
        push!(nodes, SkNode(-1, -1, -2, -2.0, fill(0.0, n_classes)))
        if b isa ConstantModel{String}
            cidx = haskey(class_to_idx, b.outcome) ? class_to_idx[b.outcome] : 0
            nodes[i+1].counts[cidx+1] = 10.0
        elseif b isa Branch{String}
            thr = b.antecedent.value.threshold
            fx  = b.antecedent.value.metacond.feature.i_variable - 1
            left_i  = dfs(b.posconsequent)
            right_i = dfs(b.negconsequent)
            nodes[i+1].feature   = fx
            nodes[i+1].threshold = thr
            nodes[i+1].left      = left_i
            nodes[i+1].right     = right_i
            for c in 1:n_classes
                nodes[i+1].counts[c] = nodes[left_i+1].counts[c] + nodes[right_i+1].counts[c]
            end
        else
            error("Tipo di nodo sconosciuto: $(typeof(b))")
        end
        return i
    end
    dfs(branch)
    n_nodes = length(nodes)
    children_left  = Vector{Int}(undef, n_nodes)
    children_right = Vector{Int}(undef, n_nodes)
    feature        = Vector{Int}(undef, n_nodes)
    threshold      = Vector{Float64}(undef, n_nodes)
    value = Array{Float64,3}(undef, n_nodes, 1, n_classes)
    for i in 1:n_nodes
        children_left[i]  = nodes[i].left
        children_right[i] = nodes[i].right
        feature[i]        = nodes[i].feature
        threshold[i]      = nodes[i].threshold
        for c in 1:n_classes
            value[i,1,c] = nodes[i].counts[c]
        end
    end
    return children_left, children_right, feature, threshold, value, n_nodes
end

function serialize_branch_sklearn(b, classes::Vector{String}, class_to_idx::Dict{String,Int})
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
        "n_classes" => n_classes
    )
end

"""
serialize_julia_ensemble(ensemble, classes)

Crea un dizionario con:
  "classes" => classes,
  "trees"   => [<tree_state1>, <tree_state2>, ...]
"""
function serialize_julia_ensemble(ensemble, classes::Vector{String})
    class_to_idx = Dict{String,Int}()
    for (i, cl) in enumerate(classes)
        class_to_idx[cl] = i-1
    end
    trees = [serialize_branch_sklearn(b, classes, class_to_idx) for b in ensemble.models]
    return Dict("classes" => classes, "trees" => trees)
end

##############################
# 2) BLOCCO PYTHON
##############################

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

##############################
# 3) FUNZIONE process_rules_decision_list: costruisce la decision list nel formato desiderato
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
        
        cond = strip(parts2[1], ['(',')',' '])
        # Fix parentheses and operators
        cond = replace(cond, "˄" => " ∧ ")
        # Ensure balanced parentheses
        if !startswith(cond, "(") && occursin("∧", cond)
            cond = "(" * cond * ")"
        end
        
        out = strip(parts2[2])
        out = replace(out, "["=>"")
        out = replace(out, "]"=>"")
        out = replace(out, "'"=>"")
        
        if isempty(cond) || cond == "( )"
            default_outcome = out
            continue
        end
        
        push!(structured, (cond, out))
    end
    
    return structured, default_outcome
end

function build_rule_list(decision_list::Vector{Tuple{String,String}}, default_outcome::String)
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
                        featvaltype = Real
                    )
                )
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
    
    # Create the DecisionList
    return DecisionList(
        rules,
        default,
        (source = "RuleCOSI", timestamp = now())
    )
end

##############################
# 4) COSTRUZIONE DEL DECISION SET IN DNF
##############################

# Utilizziamo il tipo Rule di SoleModels per costruire un DecisionSet.
struct MyRule
    formula::Formula   # La formula parsata simbolicamente
    outcome::String
end

function build_dnf_rules(decision_list::Vector{Tuple{String,String}})
    minimized_rules = Rule[]
    for (out, cond_str) in decision_list
        cond_str = replace(cond_str, "˄" => "∧")
        φ = SoleLogics.parseformula(
            cond_str;
            atom_parser = a -> Atom(
                parsecondition(
                    ScalarCondition,
                    a;
                    featuretype = VariableValue,
                    featvaltype = Real
                )
            )
        )
        # Se SoleModels.symbolize è disponibile, usala; altrimenti, usa φ direttamente.
        #sym_φ = SoleModels.symbolize(φ)
        push!(minimized_rules, Rule(φ, out))
    end
    return minimized_rules
end

function convertApi(decision_list::Vector{Tuple{String,String}})
    minimized_rules = build_dnf_rules(decision_list)
    ds = DecisionSet(minimized_rules)
    return ds
end

##############################
# 5) FUNZIONE __init__
##############################

function rulecosiplus(ensemble::Any)
    Conda.pip_interop(true, PyCall.Conda.ROOTENV)
    PyCall.Conda.pip("install", "git+https://github.com/jobregon1212/rulecosi.git", PyCall.Conda.ROOTENV)
    PyCall.Conda.pip("install", "scikit-learn", PyCall.Conda.ROOTENV)
    
    copy!(rulecosi, pyimport("rulecosi"))
    copy!(sklearn, pyimport("sklearn.ensemble"))
    
    current_dir = dirname(@__FILE__)
    csv_path = joinpath(current_dir, "data", "wisconsin.csv")
    data = CSV.read(csv_path, DataFrame)
    X = select(data, Not(:Class))
    y = data.Class
    
    unique_classes = unique(y)
    classes = collect(string.(unique_classes))
    
    Random.seed!(1212)
    idx = shuffle(1:nrow(data))
    ntrain = floor(Int, 0.9*nrow(data))
    train_idx = idx[1:ntrain]
    test_idx = idx[ntrain+1:end]
    X_train = Matrix(X[train_idx, :])
    y_train = y[train_idx]
    X_test  = Matrix(X[test_idx, :])
    y_test  = y[test_idx]
    
    classes = map(String, collect(unique(y)))
    dict_model = serialize_julia_ensemble(ensemble, classes)
    
    builder = py"build_sklearn_model_from_julia"
    base_ensemble = builder(dict_model)
    
    println("======================")
    println("Ensemble serializzato", dict_model)
    println("======================")
    
    num_estimators = pycall(pybuiltin("len"), Int, base_ensemble["estimators_"])
    println("Numero di alberi caricati in base_ensemble:", num_estimators)
    
#=
            RuleCOSIExtractor

        Estrae, combina e semplifica regole da ensemble di alberi di decisione,
        producendo un singolo modello basato su regole interpretabili.
        Supporta diversi tipi di ensemble binari (es. RandomForest, GradientBoosting, ecc.)
        e restituisce regole adatte a essere utilizzate per la classificazione.

        Parametri
        ----------
        base_ensemble : BaseEnsemble o None, opzionale (default=None)
            Istanza di un ensemble già addestrato o da addestrare. Sono supportati:
            - sklearn.ensemble.RandomForestClassifier
            - sklearn.ensemble.BaggingClassifier
            - sklearn.ensemble.GradientBoostingClassifier
            - xgboost.XGBClassifier
            - catboost.CatBoostClassifier
            - lightgbm.LGBMClassifier
            Se None, verrà creato (e addestrato) un ensemble di default (ad es. GradientBoostingClassifier).

        metric : str, opzionale (default="f1")
            Metrica utilizzata per ottimizzare la combinazione/simplificazione delle regole.
            Opzioni comuni:
            - "f1"
            - "roc_auc"
            - "accuracy"

        n_estimators : int, opzionale (default=5)
            Numero di stimatori (alberi) da addestrare se `base_ensemble` non è già addestrato.

        tree_max_depth : int, opzionale (default=3)
            Profondità massima degli alberi all’interno dell’ensemble.  
            Un valore più alto consente alberi più complessi e regole più numerose.

        conf_threshold : float, opzionale (default=0.5)
            Soglia di confidenza (o accuratezza) sotto la quale le regole vengono
            scartate durante la fase di combinazione. Più è alta la soglia,
            minore sarà il numero di regole finali.

        cov_threshold : float, opzionale (default=0.0)
            Soglia di copertura (coverage): esclude le regole che coprono meno
            campioni di quanto indicato. Con 0.0 si scartano solo le regole che
            non coprono alcun campione.

        c : float, opzionale (default=0.25)
            Livello di confidenza statistica usato per la stima dell’errore
            della regola (correzione di Laplace o simile). Serve a evitare
            overfitting nelle regole.

        percent_training : float o None, opzionale (default=None)
            Percentuale del dataset di training da usare per la
            combinazione/selezione delle regole. Se None, usa tutti i dati.

        early_stop : float, opzionale (default=0)
            Se > 0, consente di interrompere l’algoritmo di combinazione se
            per un certo numero di iterazioni (pari a `n_estimators * early_stop`)
            non si osservano miglioramenti nella metrica indicata in `metric`.

        rule_order : str, opzionale (default="supp")
            Criterio di ordinamento iniziale delle regole (nella fase iterativa di combinazione).
            Opzioni tipiche: 
            - "cov" (ordina per coverage), 
            - "conf" (ordina per confidenza), 
            - "supp" (ordina per supporto).

        sort_by_class : bool o list o None, opzionale (default=None)
            Se True, ordina il set di regole finali in base alla classe (in ordine lessicografico).
            Se è una lista, l’ordine delle classi è quello specificato dall’utente.
            Se None, non forza un ordinamento per classe.

        column_names : list di str o None, opzionale (default=None)
            Nomi delle colonne/feature del dataset. Se forniti, le regole estratte useranno
            questi nomi invece di indici numerici.

        random_state : int, RandomState o None, opzionale (default=None)
            Semenza (seed) per controllare la riproducibilità di alcuni processi stocastici
            (es. l’addestramento di un ensemble se non è già addestrato).

        verbose : int, opzionale (default=0)
            Livello di verbosità dell’output:
            - 0 = nessun output
            - 1 = stampe generali
            - 2 = informazioni dettagliate di ogni iterazione di combinazione

        df : pandas.DataFrame o None, opzionale (default=None)
            Se fornisci direttamente un DataFrame (ad es. delle feature), la classe
            può utilizzarlo per estrarre `column_names` o altre informazioni.
            Se None, si presume che le feature (X) e il vettore obiettivo (y) verranno
            passati in fase di `fit` come array/serie separate.

        Attributi
        ---------
        X_ : numpy.ndarray di shape (n_samples, n_features)
            Feature di training usate durante :meth:`fit`.

        y_ : numpy.ndarray di shape (n_samples,)
            Label/target di training usate durante :meth:`fit`.

        classes_ : numpy.ndarray di shape (n_classes,)
            Le classi incontrate in :meth:`fit`.

        original_rulesets_ : list di RuleSet, shape (n_estimators,)
            Insieme delle regole estratte direttamente dall’ensemble base,
            prima della combinazione/simplificazione.

        simplified_ruleset_ : RuleSet
            Set di regole finale, combinato e semplificato.

        n_combinations_ : int
            Numero di iterazioni o “passi” di combinazione effettuati dall’algoritmo.

        combination_time_ : float
            Tempo speso nella fase di estrazione, combinazione e semplificazione.

        ensemble_training_time_ : float
            Tempo speso per addestrare l’ensemble (se non era già addestrato).
            Se l’ensemble era già pronto, rimane 0.

        Riferimenti
        -----------
        .. [1] Obregon, J., Kim, A., & Jung, J. Y., 
            "RuleCOSI: Combination and simplification of production rules from boosted
            decision trees for imbalanced classification", 2019.

        Esempi
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
        base_ensemble=base_ensemble,
        metric="fi",
        n_estimators=100,
        tree_max_depth=100,
        conf_threshold=0.1,
        cov_threshold=0.0,
        random_state=3,
        column_names=names(X)
    )
    
    @time rc.fit(X_train, y_train)
    
    raw_rules = py"get_simplified_rules"(rc, 4, 1)
    println("Regole grezze:")
    for r in raw_rules
        println(r)
    end
    
    dl = convertToDecisionList(raw_rules)

    return dl
    #ds = convertApi(dl)
    #return ds
end

end