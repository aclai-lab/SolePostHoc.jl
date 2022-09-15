using SolePostHocAnalysis
using Documenter

DocMeta.setdocmeta!(SolePostHocAnalysis, :DocTestSetup, :(using SolePostHocAnalysis); recursive=true)

makedocs(;
    modules=[SolePostHocAnalysis],
    authors="Eduard I. STAN, Giovanni PAGLIARINI",
    repo="https://github.com/aclai-lab/SolePostHocAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="SolePostHocAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aclai-lab.github.io/SolePostHocAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aclai-lab/SolePostHocAnalysis.jl",
)
