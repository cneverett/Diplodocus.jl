# Getting Started

## Installation 
BoltzmannCollisionIntegral.jl is available to download from the Julia package manager. Inside a Julia session, enter the package manager with `]`, then run the command

```julia
pkg> add BoltzmannCollisionIntegral
```
finally load the package by running

```julia
using BoltzmannCollisionIntegral
```

## Integrating 
BoltzmannCollisionIntegral.jl contains the following modules:
- Binary Interactions 12->34
- Synchrotron emission 2->12, where 2 is a charged particle and 1 is a photon   

### Quick Start for Binary Interactions Module

!!! note
    An example script `Run_Integration.jl` for setting up and running the evaluation of the discrete binary collision integral can be found under `src/Common/` of the package. It is recommended that you copy this script to your working directory, edit the relevant fields and then running, either using the command 
    ```julia-repl
    include("Run_Integration.jl)
    ```
    in a julia-repl session, or by running the script line by line in your favourite code editor.

The example script `Run_Integration.jl` operates as follows:
- Define the names of the 4 particles involved in the interaction (``12\to34``) as the strings `name1` `name2` `name3` `name4`
    - These should of the form of three letters, which abbreviate the particles full name (see [Particles](@ref) for list of currently implemented particles).
    - They should be ordered to match a currently [Implemented Interactions](@ref)
- Define the momentum space discretisation for particles 1,2,3 and 4. This includes the upper and lower bounds of momentum magnitude (`p_low...` and `p_up...`), the grid type (`p_grid...`, see [Grids](@ref) for options) and the number of divisions (bins) (`p_num...`).
- Define the discretisation for the angular momentum space (i.e. u=cos(theta) space) for the particle species 1,2,3 and 4. This includes grid type (`u_grid...`) and number of bins (`u_num...`). The upper and lower bounds are by default assumed to be $u=\pm1$. 
- Define the number of Monte-Carlo samples to perform (as a rule of thumb, on a modern CPU it takes approximately 200ns per sample).
    - `numTiter` for the number of random sets of $\{\vec{p}_1,\vec{p}_2\}$ to sample. 
    - `numSiter` for the number of random $\{\vec{p}_3\}$ states to sample per $\{\vec{p}_1,\vec{p}_2\}$.
- If multithreading then define `nThreads` that will be used. This generates `nThreads` workers that perform evaluation in parallel, utilising `locks` to prevent data races. (see [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to set up multi-threading in Julia)
- Define the `fileLocation` where the output file ([JLD2](https://github.com/JuliaIO/JLD2.jl)) named `fileName` is to be written.
- Evaluate the emission and absorption spectrum using the [`SpectraEvaluateSerial`](@ref) function for serial and [`SpectraEvaluateMultiThread`](@ref) for multithread. Once run, these functions will save the results to the output file.

### Quick Start for Synchrotron Module

!!! note
    An example script `Run_Integration_Sync.jl` for setting up and running the evaluation of the discrete binary collision integral can be found under `src/Synchrotron/Common/` of the package. It is recommended that you copy this script to your working directory, edit the relevant fields and then running, either using the command 
    ```julia-repl
    include("Run_Integration_Sync.jl)
    ```
    in a julia-repl session, or by running the script line by line in your favourite code editor.

The example script `Run_Integration_Sync.jl` operates as follows:
- Define the name of the emitting charged particle as the strings `name2`
    - This should of the form of three letters, which abbreviate the particles full name (see [Particles](@ref) for list of currently implemented particles).
- Define the momentum space discretisation for the emitting particle and photons. This includes the upper and lower bounds of momentum magnitude (`p_low...` and `p_up...`), the grid type (`p_grid...`, see [Grids](@ref) for options) and the number of divisions (bins) (`p_num...`).
- Define the discretisation for the angular momentum space (i.e. u=cos(theta) space) for the emitting particles and photons. This includes grid type (`u_grid...`) and number of bins (`u_num...`). The upper and lower bounds are by default assumed to be $u=\pm1$.
- Define the number of Monte-Carlo samples to perform (as a rule of thumb, on a modern CPU it takes approximately 200ns per sample).
    - `numTiter` for the number of random sets of $\{\vec{p}_2\}$ to sample, i.e. emitting particle states. 
    - `numSiter` for the number of random $\{\vec{p}_2\}$ states to sample per $\{\vec{p}_2\}$, i.e. number of emitted photons to sample per emitting particle state.
- If multithreading then define `nThreads` that will be used. This generates `nThreads` workers that perform evaluation in parallel, utilising `locks` to prevent data races. (see [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to set up multi-threading in Julia)
- Define the `fileLocation` where the output file ([JLD2](https://github.com/JuliaIO/JLD2.jl)) named `fileName` is to be written.
- Evaluate the emission and absorption spectrum using the [`SyncEvaluateSerial`](@ref) function for serial and [`SyncEvaluateMultiThread`](@ref) for multithread. Once run, these functions will save the results to the output file.


## Output Files

### Output for Binary Interactions
The data stored in an output file can be loaded back into the workspace as a tuple using the [fload_All](@ref) function.
        
```julia-repl
    (Run_Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,t3MinMax,t4MinMax,SConv3,SConv4,TConv) = fload_All(fileLocation,fileName)
```

This returns a `Tuple` containing various arrays:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
- `Stot3` : A 6D matrix totalling all the emission spectrum values sampled for ``12\to34`` interaction.
- `Stot4` : A 6D matrix totalling all the emission spectrum values sampled for ``12\to43`` interaction.
- `Ttot` : A 4D matrix totalling all the absorption spectrum values sampled.
- `Stal3` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->34 interaction.
- `Stal4` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->43 interaction.
- `Ttal` : A 4D matrix of tallies of the number of absorption spectrum values sampled.
- `SMatrix3` : A 6D matrix of the emission spectrum for ``12\to34`` interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for ``12\to43`` interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for ``12\to34`` interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for ``21\to34`` interaction i.e. by permutation of TMatrix1 and correct application of phase space factors if species 1 != species 2.
- `p3Max` : The maximum value of the momentum space variable p3 sampled for each bin. (Useful for correcting numerical diffusion)
- `t3MinMax` : The minimum and maximum values of the momentum space variable t3 sampled for each bin. (Useful for correcting numerical diffusion)
- `p4Max` : The maximum value of the momentum space variable p4 sampled for each bin. (Useful for correcting numerical diffusion)
- `t4MinMax` : The minimum and maximum values of the momentum space variable t4 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv3` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for ``12\to34`` interaction.
- `SConv4` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for ``12\to43`` interaction.
- `TConv` : A 4D matrix of the convergence of the absorption spectrum compared to the previous run with given `Run_Parameters`.

Conservation of particle number and energy can be checked using the [DoesConserve](@ref) function.
- The key statistic is `ratioN` and `ratioE` which dictate the ratio of particle number and energy before and after the interaction and should be close to 1.

### Output for Synchrotron
The data stored in an output file can be loaded back into the workspace as a tuple using the `fload_All_Sync` function.
        
```julia-repl
    (Run_Parameters,SAtot,SAtal,SMatrix,#=pMax,tMinMax,=#,SConv) = fload_All_Sync(fileLocation,fileName)
```

Returns a tuple of the data stored in the file. The fields are as follows:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
- `Stot` : A 4D matrix totalling all the synchrotron emission spectrum values
- `Stal` : A 4D matrix of tallies of the number of synchrotron emission spectrum values sampled
- `SMatrix` : A 4D matrix of the synchrotron emission spectrum.
- `pMax` : The maximum value of the momentum space variable p1 (photon mommentum) sampled for each bin. (Useful for correcting numerical diffusion)
- `tMinMax` : The minimum and maximum values of the momentum space variable t1 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv` : A 4D matrix of the convergence of the synchrotron emission spectrum compared to the previous run with given `Run_Parameters`.


