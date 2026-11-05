# Inside make.jl
push!(LOAD_PATH,"../src/")
using Documenter
using DiplodocusCollisions

makedocs(
        sitename = "DiplodocusCollisions.jl",
        modules  = [DiplodocusCollisions],
        authors = "Christopher Everett",
        pages=[
                "Overview" => "index.md",
                "Getting Started" => "quickstart.md",
                "Grids, Particles and Cross Sections" => "crosssections.md",
                "Synchrotron Emissions" => "synchrotron.md",
                "Internal Functions" => "internalfunctions.md"
               ],
        checkdocs = :export
        #format = Documenter.HTML(
        #    mathengine = MathJax3(Dict(
        #        :loader => Dict("load" => ["[tex]/physics"]),
        #        :tex => Dict(
        #            "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        #            "tags" => "ams",
        #            "packages" => ["base", "ams", "autoload", "physics"],
        #        ),
        #    ))
        #),
)

deploydocs(
    repo="github.com/cneverett/DiplodocusCollisions.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"]
)