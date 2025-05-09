using Diplodocus
using Documenter
using DocumenterVitepress


# build local docs but don't deploy
makedocs(;
    modules = Module[],
    repo = Remotes.GitHub("cneverett","Diplodocus"),
    authors = "Christopher Everett",
    sitename = "Diplodocus",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/cneverett/Diplodocus",
        devurl = "dev",
        devbranch = "main",
        deploy_url = "https://cneverett.github.io/Diplodocus",
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
            "Examples" => [
                "Examples" => "Overview/Examples/examples.md"
            ]
        ],
        "Collisions" => [
            "Overview" => "DiplodocusCollisions/overview.md"
        ],
        "Transport" => [
            "Overview" => "DiplodocusTransport/overview.md"
        ],
        "Plots" => [
            "Overview" => "DiplodocusPlots/overview.md"
        ],
    ],
    #warnonly = true,
    # for liveServer, COMMENT OUT BEFORE DEPLOYING
    #clean = false,
)

deploydocs(;
    repo= "github.com/cneverett/Diplodocus",    
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,

    #devurl = "dev",
    #versions = ["stable" => "v^", "v#.#"]
)
