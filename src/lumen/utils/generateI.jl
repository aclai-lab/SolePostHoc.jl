"""
CODICE GENERATO CON IA 
"""

"""
Represents a single constraint with a variable, operator, and value.
"""
struct Constraint
    variable::String
    operator::String
    value::Float64
end

"""
Represents a collection of constraints that define a branch in a logical expression.
"""
struct Branch
    constraints::Vector{Constraint}
end

"""
Checks if a given value satisfies a constraint.

This function takes a `Constraint` and a `Float64` value, and returns `true` if the value satisfies the constraint, and `false` otherwise.

The function supports two types of constraints:
- Greater than or equal to (`>=`)
- Less than (`<`)

If the constraint operator is not one of these two, the function will return `false`.
"""
function is_satisfied(constraint::Constraint, value::Float64)
    if constraint.operator == ">="
        return value >= constraint.value
    elseif constraint.operator == "<"
        return value < constraint.value
    end
    return false
end

function is_satisfied(branch::Branch, values::Dict{String,Float64})
    all(
        is_satisfied(constraint, values[constraint.variable]) for
        constraint in branch.constraints
    )
end

"""
Parses a formula string into a vector of `Branch` objects, where each `Branch` represents a collection of constraints.

The formula string is expected to be in the format:

OR[1]: SoleLogics.SyntaxBranch: V4 ≥ 1.7000000000000002 ∧ V2 ≥ 2.6500000000000004 ∧ V3 ≥ 5.05


The function will extract the constraints from each line, split them by the `∧` operator, and create a `Constraint` object for each constraint. These `Constraint` objects are then used to create a `Branch` object, which is added to the returned vector of `Branch` objects.
"""
function parse_formula(formula_lines::Vector{<:AbstractString})  # Modificata questa riga
    branches = Branch[]

    for line in formula_lines
        isempty(strip(line)) && continue

        # Rimuove "OR[n]:" e "SoleLogics.SyntaxBranch:"
        constraints_str = split(line, ":", limit = 3)[end] |> strip

        # Divide in singoli vincoli
        constraint_parts = split(constraints_str, "∧")

        constraints = Constraint[]
        for part in strip.(constraint_parts)
            if contains(part, ">=")
                var, value = split(part, ">=")
                push!(constraints, Constraint(strip(var), ">=", parse(Float64, value)))
            elseif contains(part, "<")
                var, value = split(part, "<")
                push!(constraints, Constraint(strip(var), "<", parse(Float64, value)))
            end
        end

        push!(branches, Branch(constraints))
    end

    return branches
end

"""
Finds a valid instance of the constraints defined in the given `branches` vector.

The function iterates through the constraints in the first branch, and initializes a dictionary of variable values that satisfy the constraints. It then adjusts the variable values to ensure that all constraints in the branch are satisfied.

The function returns the dictionary of variable values that satisfy all constraints in the first branch.
"""
function find_valid_instance(branches::Vector{Branch})
    # Prova il primo branch
    branch = branches[1]

    # Dizionario per i valori delle variabili
    values = Dict{String,Float64}()

    # Per ogni vincolo nel branch
    for constraint in branch.constraints
        var = constraint.variable

        if constraint.operator == ">="
            # Se è >=, usa un valore leggermente più grande
            values[var] = constraint.value + 0.1
        elseif constraint.operator == "<"
            # Se è <, usa un valore leggermente più piccolo
            values[var] = constraint.value - 0.1
        end
    end

    # Aggiusta i valori per soddisfare tutti i vincoli
    for constraint in branch.constraints
        var = constraint.variable
        current_value = values[var]

        # Trova tutti i vincoli per questa variabile
        var_constraints = filter(c -> c.variable == var, branch.constraints)

        # Trova il minimo valore che soddisfa tutti i vincoli
        min_value = current_value
        for c in var_constraints
            if c.operator == ">="
                min_value = max(min_value, c.value + 0.01)
            elseif c.operator == "<"
                min_value = min(min_value, c.value - 0.01)
            end
        end

        values[var] = min_value
    end

    return values
end

"""
Runs the `find_valid_instance` function to find a valid instance of the constraints defined in the given `branches` vector.

This function is an example of how to use the `find_valid_instance` function to find a valid instance of the constraints. It first creates an example input formula, then parses the formula into a vector of `Branch` objects, and finally calls `find_valid_instance` to find a valid instance of the constraints. The valid instance is then printed to the console.
"""
function run_generate_istance()
    # Esempio di input
    formula = """
  OR[1]: SoleLogics.SyntaxBranch: V4 ≥ 1.7000000000000002 ∧ V2 ≥ 2.6500000000000004 ∧ V3 ≥ 5.05
    """

    # Converti le sottostringhe in stringhe complete
    formula_lines = String.(split(strip(formula), "\n"))
    branches = parse_formula(formula_lines)
    instance = find_valid_instance(branches)

    println("Istanza valida trovata:")
    for (var, value) in sort(collect(instance))
        println("$var = $value")
    end
end
