using SolePostHoc
using Documenter

DocMeta.setdocmeta!(SolePostHoc, :DocTestSetup, :(using SolePostHoc); recursive=true)

makedocs(;
    modules=[SolePostHoc],
    authors="Michele Ghiotti, Giovanni Pagliarini, Eduard I. Stan",
    repo=Documenter.Remotes.GitHub("aclai-lab", "SolePostHoc.jl"),
    sitename="SolePostHoc.jl",
    format=Documenter.HTML(;
        size_threshold = 4000000,
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
    target = "build",
    branch = "gh-pages",
    versions = ["main" => "main", "stable" => "v^", "v#.#", "dev" => "dev"],
)
