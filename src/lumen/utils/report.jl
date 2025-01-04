
"""
Generates a statistical report for the execution of a set of rules and atoms.

The report includes the following sections:
1. Rules and atoms:
   - Total number of rules
   - Number of propositions before and after deduplication
   - Reduction percentage of propositions
   - Number of atoms

2. Combinations:
   - Total number of valid combinations

3. Label distribution:
   - Count and percentage of each label

4. Formula simplification:
   - Original and simplified formula lengths
   - Reduction percentage

5. Performance:
   - Total execution time
   - Average time per combination

6. Complexity:
   - Number of trees in the forest

7. Simplified formulas:
   - Original and simplified formulas
   - Reduction percentage

Args:
    nome_file (str): The name of the file to write the report to.
    all_rules (list): The list of all rules.
    ntotatoms (int): The number of propositions before deduplication.
    nuniqatoms (int): The number of propositions after deduplication.
    results (dict): The dictionary of results.
    label_count (dict): The dictionary of label counts.
    combined_results (dict): The dictionary of combined results.
    elapsed_time (float): The total execution time.
    model (object): The model object.
"""
function genera_report_statistiche(
    nome_file,
    all_rules,
    ntotatoms,
    nuniqatoms,
    results,
    label_count,
    combined_results,
    elapsed_time,
    model,
)
    open(nome_file, "w") do file
        println(file, "=== REPORT STATISTICO DELL'ESECUZIONE ===")
        println(file, "Data e ora: ", Dates.now())
        println(file, "====================================")

        # Statistiche sulle regole e gli atomi
        println(file, "\n1. REGOLE E ATOMI")
        println(file, "  Numero totale di regole: ", length(all_rules))
        println(file, "  Numero di proposizioni prima dell'uniq: ", ntotatoms)
        println(file, "  Numero di proposizioni dopo l'uniq: ", nuniqatoms)
        println(
            file,
            "  Riduzione delle proposizioni: ",
            round((1 - nuniqatoms / ntotatoms) * 100, digits = 2),
            "%",
        )
        println(file, "  Numero di atomi: ", nuniqatoms)

        # Statistiche sulle combinazioni
        println(file, "\n2. COMBINAZIONI")
        num_combinazioni_totali = sum(length(combs) for (_, combs) in results)
        println(file, "  Numero totale di combinazioni valide: ", num_combinazioni_totali)

        # Statistiche sulle etichette
        println(file, "\n3. DISTRIBUZIONE DELLE ETICHETTE")
        for (label, count) in sort(collect(label_count), by = x -> x[2], rev = true)
            percentage = (count / num_combinazioni_totali) * 100
            println(file, "  $label: $count (", round(percentage, digits = 2), "%)")
        end

        # Statistiche sulla semplificazione
        println(file, "\n4. SEMPLIFICAZIONE DELLE FORMULE")
        for (result, formula) in combined_results
            formula_semplificata = minimizza_dnf(formula)
            riduzione =
                (
                    1 -
                    nterms(formula_semplificata) /
                    nterms(formula)
                ) * 100
            println(file, "  Etichetta $result:")
            println(file, "    Termini originali: ", nterms(formula))
            println(
                file,
                "    Termini dopo la semplificazione: ",
                nterms(formula_semplificata),
            )
            println(file, "    Riduzione: ", round(riduzione, digits = 2), "%")
        end

        # Prestazioni
        println(file, "\n5. PRESTAZIONI")
        println(
            file,
            "  Tempo totale di esecuzione: ",
            round(elapsed_time, digits = 2),
            " secondi",
        )
        println(
            file,
            "  Tempo medio per combinazione: ",
            round(elapsed_time / num_combinazioni_totali * 1000, digits = 2),
            " millisecondi",
        )

        # Complessità
        println(file, "\n6. COMPLESSITÀ")
        println(file, "  Numero di alberi nella foresta: ", length(model.trees))

        # Aggiunta della stampa delle formule semplificate
        println(file, "\n7. FORMULE SEMPLIFICATE")
        for (result, formula) in combined_results
            println(file, "  Etichetta $result:")
            formula_semplificata = minimizza_dnf(formula)
            println(file, "    Formula originale:")
            stampa_dnf(file, formula, 10)
            println(file, "    Formula semplificata:")
            stampa_dnf(file, formula_semplificata, 10)
            println(
                file,
                "    Riduzione: ",
                round(
                    (
                        1 -
                        nterms(formula_semplificata) /
                        nterms(formula)
                    ) * 100,
                    digits = 2,
                ),
                "%",
            )
            println(file)
        end

        println(file, "\n====================================")
        println(file, "Fine del report")

    end
    println("Report generato con successo: $nome_file")
end
