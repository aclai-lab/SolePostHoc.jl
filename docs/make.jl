using SolePostHoc
using Documenter

DocMeta.setdocmeta!(SolePostHoc, :DocTestSetup, :(using SolePostHoc); recursive=true)

makedocs(;
    modules=[SolePostHoc],
    authors="Eduard I. STAN, Giovanni PAGLIARINI",
    repo="https://github.com/aclai-lab/SolePostHoc.jl/blob/{commit}{path}#{line}",
    sitename="SolePostHoc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aclai-lab.github.io/SolePostHoc.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aclai-lab/SolePostHoc.jl",
)
