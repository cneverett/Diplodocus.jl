# Tutorials

Several tutorials have been designed to familiarise users with how to set up and run successful simulations using the *Diplodocus.jl* code.

- [Tutorial 1: Evolution of a Population of Hard Spheres](@ref)
- [Tutorial 2: Cooling of Electrons via Radiation Reaction](@ref)
- [Tutorial 3: Synchrotron Emissions](@ref)
- [Tutorial 4: Synchrotron Self Compton with Pair Production](@ref)

::: warning

    These tutorials only cover systems which are homogenous in space. Tutorials covering spatial variation are in development

:::

# Tutorial 0: A Brief Overview of Simulating with Diplodocus.jl

Simulating the evolution of an astrophysical system with Diplodocus.jl typically consists of three phases:
1. Building collision matrices for any particle interactions
2. Setting up and running a simulation of particle transport
3. Plotting/analysing the results

## Building Collision Matrices

One of the advantages of Diplodocus.jl is that all particle interactions are pre-computed in the form of collision matrices ahead of the main simulation of particle transport and then re-used in subsequent simulations. This pre-computation is the job of the sub-package DiplodocusCollisions.jl (loaded either as part of Diplodocus.jl or on its own). 

DiplodocusCollisions.jl currently supports a range of binary and emissive leptonic interactions, details of which you can find here [Implemented Collisions](@ref) and users may add their own interactions following the guidelines here [Adding New Collisions](@ref).

The simple templates of how to build collision matrices for binary or emissive interactions can be found in the two scripts `src/tutorials/Tutorial0_templates/BinaryCollisionMatrix_template.jl` and `src/tutorials/Tutorial0_templates/EmissionCollisionMatrix_template.jl` on the Diplodocus.jl GitHub.

Building such matrices takes the following standard steps:

1. Define the names of the particles involved (and interaction type for emissive processes) following the conventions found in [Particles, Grids and Units](@ref)
2. Define the momentum space domain for each particle and its discretisation onto a grid
3. Tell the Monte-Carlo sampler how many incoming and outgoing points to sample
4. Run the Monte-Carlo integration
5. Check the resulting collision matrices and re-run the sampling if needed to improve accuracy

This process will generate a new file in the `.jld2` format that contains the collision matrices for that interaction. 

::: warning

The `fileName` generated in this process is used internally when running simulations to ensure the correct data is loaded, and allows for more Monte-Carlo sampling to be run to improve accuracy without starting over, it is not recommended to changed the file name of any collision matrices generated.

:::

## Setting Up and Running a Simulation

Running simulations is the job of the DiplodocusTransport.jl sub-package (loaded either as part of Diplodocus.jl or on its own). 

Running a simulation follows these standard steps:

1. Define the coordinates and grids for time and space and any boundary conditions
2. Define the particles in the simulation and their momentum domain and grids
3. Define which binary and emissive interactions to include between particles, and any external forces that will be applied to these particles
4. Set up the initial conditions for the particles
5. Choose a solver (currently only forward Euler)
6. Run the simulation

Due to the generality of the DIPLODOCUS framework, this setup may require more steps than other codes but at the same time gives more flexibility in the physics that a user may want to include. A setup template can be found in `src/tutorials/Tutorial0_templates/Transport_template.jl`, but it is recommended that new users work through the following tutorials in sequence as they provide much more detail.

## Plotting Results

Once any collision matrices are generated and simulations are run, DiplodocusPlots.jl (loaded either as part of Diplodocus.jl or on its own) hosts a suite of plotting functions for analysing results. It is built upon the [Makie.jl](https://docs.makie.org/stable/) plotting ecosystem making it easy to edit plots for users to write their own plotting scripts.