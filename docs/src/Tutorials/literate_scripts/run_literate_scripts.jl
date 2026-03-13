using Literate

# Tutorial 1

INPUT = joinpath(@__DIR__,"Tutorial1")
OUTPUT_MD = joinpath(@__DIR__,"..","Tutorial1")
OUTPUT_JL = joinpath(@__DIR__,"..","..","..","..","Tutorials","Tutorial1")

Literate.markdown(joinpath(INPUT,"tutorial1_collisions.jl"), OUTPUT_MD,credit=false,codefence="````julia" => "````")
Literate.script(joinpath(INPUT,"tutorial1_collisions.jl"), OUTPUT_JL,credit=false)

Literate.markdown(joinpath(INPUT,"tutorial1_transport.jl"), OUTPUT_MD,credit=false,codefence="````julia" => "````")
Literate.script(joinpath(INPUT,"tutorial1_transport.jl"), OUTPUT_JL,credit=false)

Literate.markdown(joinpath(INPUT,"tutorial1_plots.jl"), OUTPUT_MD,credit=false,codefence="````julia" => "````")
Literate.script(joinpath(INPUT,"tutorial1_plots.jl"), OUTPUT_JL,credit=false)