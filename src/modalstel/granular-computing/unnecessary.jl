import Base: convert

function unnecessary(m::Union{AbstractModel,AbstractSyntaxStructure})
    error("Provide a method for $(typeof(m))")
end

unnecessary(m::Vector{<:Rule}) = unnecessary.(m)
unnecessary(m::Rule) = Rule(unnecessary(antecedent(m)),consequent(m),info(m))
# unneccessary(m::SyntaxTree)
unnecessary(m::LeftmostLinearForm) = children(m)[end]
unnecessary(m::Literal) = m

############################################################################################
############################################################################################
############################################################################################

function convert(
    tree::DTree,
    info = (;),
)
    new_root = convert(ModalDecisionTrees.root(tree))

    info = merge(info, SoleModels.info(new_root))
    info = merge(info, (;))

    return DecisionTree(new_root, info)
end

function convert(
    tree::DTInternal,
    info = (;),
)
    f = formula(ModalDecisionTrees.decision(tree))
    p = MultiFormula(i_modality(tree), SyntaxTree(get_atom(f)))
    ◊ = DiamondRelationalConnective{typeof(relation(f))}()

    info = merge(info, (;
        supporting_labels = ModalDecisionTrees.supp_labels(node),
    ))

    return SoleModels.Branch(◊p, convert(left(node), (;)), convert(right(node), (;)), info)
end

function convert(
    tree::DTLeaf,
    info = (;),
)
    info = merge(info, (;
        supporting_labels      = ModalDecisionTrees.supp_labels(tree),
        supporting_predictions = ModalDecisionTrees.predictions(tree),
    ))

    return SoleModels.ConstantModel(ModalDecisionTrees.prediction(tree), info)
end
