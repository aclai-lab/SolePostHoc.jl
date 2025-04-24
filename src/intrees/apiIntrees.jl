# Custom Atom parser for SoleLogics.parseformula.
# Since conditions are now written as scalar conditions, we use the parser for ScalarCondition here.
atom_parser = function (a::String)
    #println("Parsing atom: ", a)
    return Atom(parsecondition(SoleData.ScalarCondition, a;
        featuretype=SoleData.VariableValue,
        featvaltype=Real))
end

# Function to convert a DNF formula back to SyntaxBranch representation
function dnf_to_syntaxbranch(dnf_formula)
    if dnf_formula isa LeftmostLinearForm
        # This is a disjunction of conjunctions (typical DNF structure)
        disjuncts = dnf_formula.grandchildren

        if length(disjuncts) == 1
            # If there's only one disjunct, we just need the conjunction
            return conjunction_to_syntaxbranch(disjuncts[1])
        elseif length(disjuncts) == 2
            # Binary disjunction of two conjunctions
            child1 = conjunction_to_syntaxbranch(disjuncts[1])
            child2 = conjunction_to_syntaxbranch(disjuncts[2])
            return SyntaxBranch(NamedConnective{:∨}(), (child1, child2))
        else
            # For more than two disjuncts, create nested binary disjunctions
            current_branch = conjunction_to_syntaxbranch(disjuncts[1])

            # Combine with each remaining disjunct using binary disjunction
            for i in 2:length(disjuncts)
                next_disjunct = conjunction_to_syntaxbranch(disjuncts[i])
                current_branch = SyntaxBranch(NamedConnective{:∨}(), (current_branch, next_disjunct))
            end

            return current_branch
        end
    elseif dnf_formula isa LeftmostConjunctiveForm
        # This is just a conjunction
        return conjunction_to_syntaxbranch(dnf_formula)
    elseif dnf_formula isa Literal
        # This is a single literal
        return literal_to_syntaxbranch(dnf_formula)
    elseif dnf_formula isa Atom
        # This is a single atom
        return dnf_formula
    else
        error("Unexpected formula type: $(typeof(dnf_formula))")
    end
end

# Convert a conjunction (LeftmostConjunctiveForm) to SyntaxBranch
function conjunction_to_syntaxbranch(conjunction)
    literals = conjunction.grandchildren

    if length(literals) == 1
        # If there's only one literal, we don't need a conjunction
        return literal_to_syntaxbranch(literals[1])
    elseif length(literals) == 2
        # Create a conjunction of two literals (binary operation)
        child1 = literal_to_syntaxbranch(literals[1])
        child2 = literal_to_syntaxbranch(literals[2])
        return SyntaxBranch(NamedConnective{:∧}(), (child1, child2))
    else
        # For more than two literals, create nested binary conjunctions
        # Convert first literal to SyntaxBranch
        current_branch = literal_to_syntaxbranch(literals[1])

        # Combine with each remaining literal using binary conjunction
        for i in 2:length(literals)
            next_literal = literal_to_syntaxbranch(literals[i])
            current_branch = SyntaxBranch(NamedConnective{:∧}(), (current_branch, next_literal))
        end

        return current_branch
    end
end

# Convert a Literal to SyntaxBranch

#### already defined in apiREFNESole.jl

function literal_to_syntaxbranch(literal)
    if literal.ispos
        # Positive literal, just return the atom
        return literal.atom
    else
        # Negative literal, wrap the atom in a negation
        return SyntaxBranch(NamedConnective{:¬}(), (literal.atom,))
    end
end

