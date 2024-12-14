module Lumen

using Revise
using Logging
using Dates
using DataStructures
using SoleModels
using AbstractTrees
using SoleData
using SoleData: MultivariateScalarAlphabet, UnivariateScalarAlphabet

using ModalDecisionTrees
using SoleLogics
using BenchmarkTools, StatProfilerHTML
using Base: intersect
using Base.Threads: @threads
using Base.Threads: Atomic, atomic_add!
using Profile
using ConcurrentCollections
using ProgressMeter

"""
    lumen(model, tempo_inizio, vertical = 1.0, orizontal = 1.0, ott_mode = false, 
          controllo = false, minimization_scheme = :espresso; 
          minimization_kwargs = (;), filteralphabetcallback! = identity, kwargs...)

Logic-driven Unified Minimal Extractor of Notions (LUMEN): A function that extracts and minimizes 
logical rules from a decision tree ensemble model into DNF (Disjunctive Normal Form) formulas.

# Arguments
- `model`: The decision tree ensemble model to analyze
- `tempo_inizio`: Start time for performance measurement
- `vertical::Real`: Vertical coverage parameter (α) for rule extraction (0.0 < vertical ≤ 1.0)
- `orizontal::Real`: Horizontal coverage parameter (β) for rule extraction (0.0 < orizontal ≤ 1.0)
- `ott_mode::Bool`: Flag to enable optimized processing mode
- `controllo::Bool`: Flag to enable validation mode for comparing results
- `minimization_scheme::Symbol`: DNF minimization algorithm to use (e.g., :espresso)
- `minimization_kwargs::NamedTuple`: Additional arguments for the minimization algorithm
- `filteralphabetcallback!`: Callback function for alphabet processing (default: identity)
- `kwargs...`: Additional keyword arguments

# Processing Steps
1. Rule Extraction:
   - Extracts logical rules from the decision trees in the model
   - Validates and processes input parameters

2. Alphabet Construction:
   - Extracts atoms from the rules
   - Constructs and processes the logical alphabet
   - Applies custom alphabet filtering if provided

3. Truth Table Generation:
   - Creates truth combinations based on the extracted atoms
   - Handles both standard and optimized processing modes
   - Generates L-class path forest for rule representation

4. DNF Minimization:
   - Applies the specified minimization scheme to simplify the DNF formulas
   - Verifies the correctness of the simplification
   - Reports the reduction in term count

5. Result Validation (when controllo = true):
   - Compares results between different processing modes
   - Validates optimization correctness

# Returns
No explicit return value, but produces:
- Simplified DNF formulas for each class
- Performance metrics and validation results
- Statistical reports (when enabled)

# Notes
- Parameters `vertical` and `orizontal` must be in range (0.0, 1.0]
- Setting both parameters to 1.0 enforces strong rule extraction
- The function aligns with Algorithm 1 from the reference implementation
- Performance statistics are generated based on processing time and rule reduction

# Example
```julia
model = load_decision_tree_model()
start_time = time()
lumen(model, start_time, 1.0, 1.0, false, false, :espresso)
```
"""
function lumen(
  model,
  modelJ, # attualmente truth_combinations usa model non di sole ? 
  tempo_inizio,
  vertical::Real = 1.0,
  orizontal::Real = 1.0,
  ott_mode::Bool = false,
  controllo::Bool = false,
  minimization_scheme::Symbol = :espresso;
  minimization_kwargs::NamedTuple = (;),
  filteralphabetcallback! = identity,
  kwargs...
)
  if vertical <= 0.0 || vertical > 1.0 || orizontal <= 0.0 || orizontal > 1.0
      @warn "Inserito parametri non validi"
      @warn "Verranno settati entrambi a 1"
      vertical = 1.0 # inserito
      orizontal = 1.0 # TODO MARCO METTILA, NON TI SCORDARE!! 
  end

  spa() && println(
      "\n\n$COLORED_TITLE$TITLE\n PARTE 2.a ESTRAZIONE DELLE REGOLE DAGLI ALBERI \n$TITLE$RESET",
  )

  all_rules = vcat([listrules(tree) for tree in model.models]...)
  spa() && println(all_rules)

  spa() && println(
      "\n\n$COLORED_TITLE$TITLE\n PARTE 2.b ESTRAZIONE DEGLI ATOMI DALLE REGOLE \n$TITLE$RESET",
  )

  num_all_atoms, my_atoms, my_alphabet = begin
    all_atoms = [atom for rule in all_rules for atom in atoms(antecedent(rule))] # all_atoms = collect(atoms(SoleModels.alphabet(model, false)))
    num_all_atoms = length(all_atoms)
    my_atoms = unique(all_atoms)

    spa() && println(my_atoms)

    spa() && println(
        "\n\n$COLORED_TITLE$TITLE\n PARTE 2.c ESTRAZIONE DELL'ALFABETO pt2 :< \n$TITLE$RESET",
    )

    # Get number of features from the maximum feature index in atoms
    n_features = maximum(atom.value.metacond.feature.i_variable for atom in my_atoms)
    my_alphabet = Lumen.process_alphabet(my_atoms, n_features)
    #=
        
        # @show my_alphabet
        filteralphabetcallback!(my_alphabet)
        # @show my_alphabet

        all(x -> (x == (<)), SoleData.test_operator.(subalphabets(my_alphabet))) ||
            error("Atoms with test operators other from < are not supported.")

    =#
    spa() && println(my_alphabet)
    num_all_atoms, my_atoms, my_alphabet
  end

  

  if (controllo == false)
      spa() && println("\n\n$COLORED_TITLE$TITLE\n PARTE 3 TABELLA \n$TITLE$RESET")

      if (ott_mode == true)
          results, label_count =
              Lumen.truth_combinations_ott(modelJ, my_alphabet, my_atoms, vertical, apply_forest)
      else
          results, label_count =
              Lumen.truth_combinations(modelJ, my_alphabet, my_atoms, vertical, apply_forest)
      end

      spa() && println(
          "\n\n$COLORED_TITLE$TITLE\n PARTE 3 GENERAZIONE L-CLASS PATH FOREST & CONSEGUENTE SEMPLIFICAZIONE CON METODI AD HOC\n$TITLE$RESET",
      )

      num_atoms = length(my_atoms)

      thresholds_by_feature = Dict(
          subalpha.featcondition[1].feature.i_variable =>
              sort(subalpha.featcondition[2]) for subalpha in my_alphabet.subalphabets
      )
      atoms_by_feature = Dict{Int,Vector{Tuple{Float64,Bool}}}()
      for atom in my_atoms
          feat = atom.value.metacond.feature.i_variable
          threshold = atom.value.threshold
          push!(
              get!(Vector{Tuple{Float64,Bool}}, atoms_by_feature, feat),
              (threshold, true),
          )
      end
      for (_, atom_list) in atoms_by_feature
          sort!(atom_list, by = first)
      end

      combined_results =
          Lumen.concat_results(results, num_atoms, thresholds_by_feature, atoms_by_feature)

      start_time = 0
      for (result, formula) in combined_results
          spa() && println("Risultato: $result")
          spa() && stampa_dnf(stdout, formula)
          spa() && println()

          @info "Iniziando la semplificazione per il risultato $result"
          start_time = time()

          formula_semplificata = @timed Lumen.minimizza_dnf(
              Val(minimization_scheme),
              formula;
              minimization_kwargs...,
          )
          try
              @info "Semplificazione completata in $(formula_semplificata.time) secondi"
              spa() && println("$COLORED_INFO**************⬆️**************$RESET")

              spa() && println("Termini originali: ", length(formula.combinations))
              spa() && println(
                  "Termini dopo la semplificazione: ",
                  length(formula_semplificata.value.combinations),
              )
              spa() && println()

              Lumen.printIO_custom_or_formula(stdout, formula_semplificata.value)

              # Verifica della semplificazione
              is_congruent = Lumen.verify_simplification(formula, formula_semplificata.value)
              spa() && println(
                  "\n\n$COLORED_INFO$TITLE\n PARTE 3.a Semplificazione valida ?\n$TITLE$RESET",
              )
              if is_congruent
                  @info "La semplificazione è stata verificata e risulta corretta."
              else
                  @warn "ATTENZIONE: La semplificazione non è congruente con la formula originale!"
              end

          catch e
              @error "Errore durante la semplificazione: $e"
          end
          spa() && println(
              "----------------------------------------------------------------------------",
          )
      end
      tempo_fine = time()

      spa() &&
          println("\n\n$COLORED_TITLE$TITLE\n PARTE 4 DOCUMENTO I DATI \n$TITLE$RESET")

      tempo_esecuzione = tempo_fine - tempo_inizio

      spa() && begin
          if (ott_mode == true)
              #nome_file_report = "report_statistico_Parallelo_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS")).txt"
              #genera_report_statistiche_ott(nome_file_report, all_rules, num_all_atoms, num_atoms, results, label_count, combined_results, tempo_esecuzione, model)
          else
              #nome_file_report = "report_statistico_Seriale_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS")).txt"
              #genera_report_statistiche(nome_file_report, all_rules, num_all_atoms, num_atoms, results, label_count, combined_results, tempo_esecuzione, model)
          end
      end
  end

  if (controllo == true)
      spa() && println(
          "\n\n$COLORED_INFO$TITLE\n PARTE 2.d È un ottimizzazione valida?\n$TITLE$RESET",
      )

      are_results_equal =
          Lumen.compare_truth_combinations(modelJ, my_alphabet, my_atoms, vertical; kwargs...)

      if are_results_equal
          @info "\nL'ottimizzazione è valida: i risultati sono identici."
      else
          @warn "\nATTENZIONE: L'ottimizzazione potrebbe non essere valida. Ci sono differenze nei risultati."
      end
  end
end

##
# Types
##

# TritVector TODO Choose and implement one 
include("types/trit-vector.jl")
include("types/balanced-trit-vector.jl")
include("types/balanced-ternary-vector.jl")

# ConstVariable e myOwnTypes
include("types/types.jl")

##
# Utils
##

# Gestione dei report su file
include("utils/report.jl")

# IO-Algoritmic-Utils - IO algoritmico del progetto
include("utils/IO.jl")

# Minor-Algoritmic-Utils - core algoritmico del progetto
include("utils/minor.jl")

# Algoritmic-Utils - core algoritmico del progetto
include("utils/core.jl")

# Algoritmic-Optimization-Utils - core algoritmico del progetto se avviato in ott_mode
include("utils/coreOttMode.jl")

include("utils/minimization.jl")

include("deprecate.jl")

end
