using Documenter
using DocumenterVitepress
using Diplodocus
using DiplodocusCollisions
using DiplodocusTransport
using DiplodocusPlots


# build local docs but don't deploy
makedocs(;
    modules = [Diplodocus],
    repo = Remotes.GitHub("cneverett","DiplodocusTransport.jl#dev"),
    authors = "Christopher Everett",
    sitename = "Diplodocus.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/cneverett/Diplodocus.jl",
        devurl = "dev",
        devbranch = "main",
        #deploy_url = "https://cneverett.github.io/Diplodocus.jl/",
        # for liveServer, COMMENT OUT BEFORE DEPLOYING
        #md_outpath = ".",
        #build_vitepress = false,
    ),
    draft = false,
    source = "src",
    build = "build",
    pages=[
        "Overview" => [
            "Overview" => "Overview/overview.md"
            "Installation" => "Overview/installation.md"
            "Tutorials" => [
                "Hard Spheres" => "Overview/Examples/hardsphere.md"
            ]
        ],
        "Collisions" => [
            "Overview" => "DiplodocusCollisions/overview.md",
            "Cross Sections" => "DiplodocusCollisions/cross sections.md",
        ],
        "Transport" => [
            "Overview" => "DiplodocusTransport/overview.md"
        ],
        "Plots" => [
            "Overview" => "DiplodocusPlots/overview.md"
        ],
    ],
    #linkcheck = true,
    doctest = false,
    warnonly = [:missing_docs],
    clean = true,
    checkdocs = :none
    # for liveServer, COMMENT OUT BEFORE DEPLOYING
    #clean = false,
)

deploydocs(;
    repo= "github.com/cneverett/Diplodocus.jl",    
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,

    #devurl = "dev",
    #versions = ["stable" => "v^", "v#.#"]
)
