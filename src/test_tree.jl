################################################
        #              p
        #      ┌───────┴─────────────┐
        #      │                     r
        #      q                 ┌───┴───┐
        #      │                 s      "yes"
        #  ┌───┴───┐         ┌───┴───┐
        # "yes"   "no"      "yes"   "no"
##################################################

using SoleLogics

formula_p = Formula(FNode(Letter("p")))
formula_q = Formula(FNode(Letter("q")))
formula_r = Formula(FNode(Letter("r")))
formula_s = Formula(FNode(Letter("s")))

branch_q = Branch(formula_q,("yes","no"),(;))
branch_s = Branch(formula_s,("yes","no"),(;))
#branch_r = Branch(formula_r,(branch_s,"yes"),(;))

#dt_q = DecisionTree(branch_r,(;))


#Possibile path
path_all = [formula_p,formula_q,formula_s,formula_r,"yes"]
path_2 = [formula_p,formula_q,"yes"]
path_1 = [formula_p,"yes"]
path_0 = ["yes"]
