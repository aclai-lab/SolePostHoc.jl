
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
    file_name,
    all_rules,
    ntotatoms,
    nuniqatoms,
    results,
    label_count,
    combined_results,
    elapsed_time,
    model,
)
    open(file_name, "w") do file
        println(file, "=== STATISTICAL REPORT OF EXECUTION ===")
        println(file, "Date and time: ", Dates.now())
        println(file, "====================================")

        # Statistics on rules and atoms
        println(file, "\n1. RULES AND ATOMS")
        println(file, "  Total number of rules: ", length(all_rules))
        println(file, "  Number of propositions before deduplication: ", ntotatoms)
        println(file, "  Number of propositions after deduplication: ", nuniqatoms)
        println(
            file,
            "  Reduction of propositions: ",
            round((1 - nuniqatoms / ntotatoms) * 100, digits = 2),
            "%",
        )
        println(file, "  Number of atoms: ", nuniqatoms)

        # Statistics on combinations
        println(file, "\n2. COMBINATIONS")
        total_combinations = sum(length(combs) for (_, combs) in results)
        println(file, "  Total number of valid combinations: ", total_combinations)

        # Statistics on labels
        println(file, "\n3. LABEL DISTRIBUTION")
        for (label, count) in sort(collect(label_count), by = x -> x[2], rev = true)
            percentage = (count / total_combinations) * 100
            println(file, "  $label: $count (", round(percentage, digits = 2), "%)")
        end

        # Statistics on simplification
        println(file, "\n4. FORMULA SIMPLIFICATION")
        for (result, formula) in combined_results
            simplified_formula = minimize_dnf(formula)
            reduction = (1 - nterms(simplified_formula) / nterms(formula)) * 100
            println(file, "  Label $result:")
            println(file, "    Original terms: ", nterms(formula))
            println(file, "    Terms after simplification: ", nterms(simplified_formula))
            println(file, "    Reduction: ", round(reduction, digits = 2), "%")
        end

        # Performance
        println(file, "\n5. PERFORMANCE")
        println(
            file,
            "  Total execution time: ",
            round(elapsed_time, digits = 2),
            " seconds",
        )
        println(
            file,
            "  Average time per combination: ",
            round(elapsed_time / total_combinations * 1000, digits = 2),
            " milliseconds",
        )

        # Complexity
        println(file, "\n6. COMPLEXITY")
        println(file, "  Number of trees in the forest: ", length(model.trees))

        # Adding the print of simplified formulas
        println(file, "\n7. SIMPLIFIED FORMULAS")
        for (result, formula) in combined_results
            println(file, "  Label $result:")
            simplified_formula = minimize_dnf(formula)
            println(file, "    Original formula:")
            print_dnf(file, formula, 10)
            println(file, "    Simplified formula:")
            print_dnf(file, simplified_formula, 10)
            println(
                file,
                "    Reduction: ",
                round((1 - nterms(simplified_formula) / nterms(formula)) * 100, digits = 2),
                "%",
            )
            println(file)
        end

        println(file, "\n====================================")
        println(file, "End of report")

    end
    println("Report successfully generated: $file_name")
end
