using SolePostHoc

using Documenter

DocMeta.setdocmeta!(SolePostHoc, :DocTestSetup, :(using SolePostHoc); recursive = true)

makedocs(;
    modules = [SolePostHoc],
    authors = "Michele Ghiotti, Giovanni Pagliarini, Eduard I. Stan, Marco Perrotta",
    repo = Documenter.Remotes.GitHub("aclai-lab", "SolePostHoc.jl"),
    sitename = "SolePostHoc.jl",
    checkdocs = :none,
    format = Documenter.HTML(;
        size_threshold = 4000000,
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://aclai-lab.github.io/SolePostHoc.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Extract with algorithms" => "extract-algorithms.md",

        # "Hands on" => "hands-on.md",
        "Contributing" => "contributing.md",
    ],
    warnonly = :true,
)

deploydocs(;
    repo = "github.com/aclai-lab/SolePostHoc.jl",
    devbranch = "main",
    target = "build",
    branch = "gh-pages",
    versions = ["main" => "main", "stable" => "v^", "v#.#", "dev" => "dev"],
)
