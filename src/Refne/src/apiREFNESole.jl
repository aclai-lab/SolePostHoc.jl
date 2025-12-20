using DecisionTree
using SoleModels
using SoleLogics

# Custom Atom parser for SoleLogics.parseformula.
# Since conditions are now written as scalar conditions, we use the parser for ScalarCondition here.
atom_parser = function (a::String)
    # println("Parsing atom: ", a)
    return Atom(
        parsecondition(
            SoleData.ScalarCondition,
            a;
            featuretype = SoleData.VariableValue,
            featvaltype = Real,
        ),
    )
end
