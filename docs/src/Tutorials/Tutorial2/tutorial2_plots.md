```@meta
EditURL = "../literate_scripts/Tutorial2/tutorial2_plots.jl"
```

# Tutorial 2c: Plotting Results

---

Now that we've run the simulation, we can use a variety of pre-built plotting functions from `DiplodocusPlots` to plot and analyse the results.

## Load the Simulation Output

First thing to do is load the simulation output using the same `fileLocation` and `fileName` as before:

````julia
using Diplodocus

fileName = "tutorial2_output.jld2"
fileLocation = joinpath(pwd(),"Data")
(PhaseSpace,sol) = SolutionFileLoad(fileLocation,fileName)
````

## Particle Distribution

