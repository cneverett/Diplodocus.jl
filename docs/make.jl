using Documenter
using DocumenterVitepress
using DocumenterCitations
using Diplodocus
using DiplodocusCollisions
using DiplodocusTransport
#using DiplodocusPlots

push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"/docs/")

bib = CitationBibliography(
    joinpath(@__DIR__,"src","refs.bib");
    style = :authoryear
)

# build local docs but don't deploy
makedocs(;
    modules = [Diplodocus,DiplodocusCollisions],
    plugins = [bib],
    repo = Remotes.GitHub("cneverett","Diplodocus.jl"),
    authors = "Christopher Everett",
    sitename = "Diplodocus.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/cneverett/Diplodocus.jl",
        devurl = "dev",
        devbranch = "main",
        deploy_url = "https://cneverett.github.io/Diplodocus.jl/",
        # for liveServer, COMMENT OUT BEFORE DEPLOYING
        #md_outpath = ".",
        #build_vitepress = false,
    ),
    draft = false,
    source = "src",
    build = "build",
    pages=[
        "Get Started" => [
            "Overview" => "Get Started/overview.md",
            "Installation" => "Get Started/installation.md",
            "Conventions" => "Get Started/conventions.md",
        ],
        "Tutorials" => [
            "Overview" => "Tutorials/overview.md",
            "1. Hard Spheres" => "Tutorials/hardsphere.md",
            "2. Radiation Reaction" => "Tutorials/radreact.md",
            "3. Synchrotron" => "Tutorials/synchrotron.md",
            "4. Synchrotron Self-Compton" => "Tutorials/SSC.md",
        ],
        "Collisions" => [
            "Overview" => "DiplodocusCollisions/overview.md",
            "Implemented Collisions" => "DiplodocusCollisions/implemented collisions.md",
            "Adding New Collisions" => "DiplodocusCollisions/adding new collisions.md",
        ],
        "Transport" => [
            "Overview" => "DiplodocusTransport/overview.md",
            "Implemented Coordinates and Forces" => "DiplodocusTransport/implemented terms.md",
        ],
        "Plots" => [
            "Overview" => "DiplodocusPlots/overview.md"
        ],
    ],
    #linkcheck = true,
    doctest = false,
    warnonly = [:missing_docs],
    #clean = true,
    checkdocs = :none
    # for liveServer, COMMENT OUT BEFORE DEPLOYING
    #clean = false,
)

DocumenterVitepress.deploydocs(;
    repo= "github.com/cneverett/Diplodocus.jl",    
    target = joinpath(@__DIR__, "build"),
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,

    #devurl = "dev",
    versions = ["stable" => "v^", "v#.#","dev"=>"dev"]
)
