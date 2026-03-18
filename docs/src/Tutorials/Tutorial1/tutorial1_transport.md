```@meta
EditURL = "../literate_scripts/Tutorial1/tutorial1_transport.jl"
```

# Tutorial 1b: Running the Simulation

---

With the generation of a good set of collision matrices, we are now ready to run the simulation. This simulation will be run using the functions contained within the `DiplodocusTransport` sub-package.

## Backend Setup

Before we can run any simulation we need to set some backend parameters, those being the `Precision` of the simulation (either `Float64` or `Float32`) and the `Backend`. This `Backend` tells the solver what hardware to behind the scenes. Here we choose to run our simulation on a CPU by specifying `CPUBackend()`. If we had a Nvidea GPU we could use this instead by setting it to `CUDABackend()`.

````julia
using Diplodocus

Precision::DataType = Float32
Backend::BackendType = CPUBackend()
````

## Characteristic Scales

Next we define any simulation specific scales in SI units. This will override the default scales, which are specified in ... LINK TO SCALES ...

````julia
CHAR_number_density::Float64 = 1.0 # m^-3
CHAR_time = 1.0 / DiplodocusTransport.CONST_σT / DiplodocusTransport.CONST_c / CHAR_number_density # s
````

Here we have defined the characteristic number density to be to ``n=1 \mathrm{m}^{-3}``. We then set the characteristic time to be ``t=1/(n*σT*c)`` where ``σT`` is the Thompson scattering cross section (i.e. the cross section used for the elastic scattering of hard spheres), and ``c`` is the speed of light. As we are going to be examining the relativistic collisions of these spheres, this time corresponds to the characteristic scattering time.

From here on out, all values are assumed to be in term of the characteristic scales of the system.

## Phase-Space Setup

### Spatial Grids

We need to define the spatial geometry on which we are going to simulate. In this tutorial we are going to assume our spheres are homogenous in space so we only need a single spatial cell. This cell we will take as a unit cube with `Cartesian` coordinates (the `Cartesian` function be default specifies periodic boundary conditions so we don't need to worry about that here).

````julia
space_coords::CoordinateType = Cartesian()

x_up::Float64 = 0.5
x_low::Float64 = -0.5
x_grid::String = "u"
x_num::Int64 = 1

y_up::Float64 = 0.5
y_low::Float64 = -0.5
y_grid::String = "u"
y_num::Int64 = 1

z_up::Float64 = 0.5
z_low::Float64 = -0.5
z_grid::String = "u"
z_num::Int64 = 1
````

No matter the type of spatial coordinates, they are always mapped to ``x,y,z``, for other coordinate options, see [Implemented Coordinates and Forces](@ref).

### Momentum Grids

Next we similarly define the momentum space grids and the list of particles to be included in the simulation. For momentum space, Diplodocus uses a local spherical coordinate basis, therefore, the momentum coordinates should be `Spherical`. As for the grids for each coordinate, we should choose the same values as were used to generate the collision matrices.

````julia
name_list::Vector{String} = ["Sph",];

momentum_coords = Spherical() # px = p, py = u, pz = phi

px_up_list::Vector{Float64} = [4.0,]
px_low_list::Vector{Float64} = [-5.0,]
px_grid_list::Vector{String} = ["l",]
px_num_list::Vector{Int64} = [72,]

py_up_list::Vector{Float64} = [1.0,]
py_low_list::Vector{Float64} = [-1.0,]
py_grid_list::Vector{String} = ["u",]
py_num_list::Vector{Int64} = [9,]

pz_up_list::Vector{Float64} = [2.0*pi,]
pz_low_list::Vector{Float64} = [0.0,]
pz_grid_list::Vector{String} = ["u",]
pz_num_list::Vector{Int64} = [1,]
````

::: warning

Momentum-space grids must match those used in generating the collision matrices

:::

### Build the Phase Space

Last part of the phase-space setup is to actually build the phase space using what we have defined above

````julia
PhaseSpace = PhaseSpaceStruct()
````

## Interactions and Forces

With the phase space now set up we can define what interactions and forces take place on this phase space.

### Defining what Interactions to Include

In our case we only have a single binary interaction to consider. For this we need to define two things: the `Binary_list`, a `Vector` of `BinaryStruct`s that define the binary interactions to include; and a `Binary_Domain` that defines where these interactions will take place.

::: tip

For larger simulations with large spatial grids, you can use the `Binary_Domain` to specify to tell the solver only to consider binary interactions in certain areas and speed up simulations.

:::

````julia
Binary_list::Vector{BinaryStruct} = [BinaryStruct("Sph","Sph","Sph","Sph")]
Binary_Domain = CollisionDomain(PhaseSpace)
````

Here the `CollisionDomain` is over all phase space, i.e. the single spatial cell of our simulation.

Even though we don't have any forces or emissive interactions we still need to specify this by defining them as empty `Vector`s:

````julia
Emission_list::Vector{EmiStruct} = []
Forces_list::Vector{ForceType} = []
````

### Building the Interaction and Force Matrices

Now that we hae defined what interactions and forces to include, we need to build the matrices that the solver will actually use to dictate the time stepping of the simulation.

First we define the `DataDirectory` where any saved collision matrices are saved and will be loaded from.

````julia
DataDirectory = joinpath(pwd(),"Data")
````

Then we build three sets of matrices: `BinM` for the binary interactions, `EmiM` for the emissive interactions and `FluxM` for phase-space transport.

````julia
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=false); # jl
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=false) # md

EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true); # jl
EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true) # md

FluxM = BuildFluxMatrices(PhaseSpace,Forces_list); # jl
FluxM = BuildFluxMatrices(PhaseSpace,Forces_list) # md
````

Here we have used the argument `Bin_corrected=true` to use the number and energy density correct binary collision matrices, and the argument `Bin_sparse=false` to define the binary collision matrix as a dense matrix.
The size of all these matrices is determined by the resolution of the phase-space grids and the sparsity of the interactions. The total size will be printed to the standard output once these matrices are built.

::: tip

If there are lots of particles involved in a simulation the binary collision matrix (typically the largest matrix) may become sparse enough to warrant using `Bin_sparse=true` to save memory and speed up the solver. This is only works if the sparsity is grater than around 70%.

:::

## Initial and Injection Conditions

With the phase space and interactions all set up, we can now define the initial and injection conditions of our simulation. To do this we first initialise two state vectors `Initial` and `Injection` using:

````julia
Initial = Initialise_Initial_Condition(PhaseSpace)
Injection = Initialise_Injection_Condition(PhaseSpace)
````

These two state vectors are currently empty, but we can use several functions to add to them a variety of initial and injection conditions ... LINK TO INIT/INJ OPTIONS ... . For this tutorial, we don't need to inject any particles, so there is nothing we need to do to `Injection`.
As for `Initial`, in this tutorial we will consider the spheres to be initially distributed with momenta ``p`` between ``10^{3}m_\text{Ele}c`` and ``10^{3.1}m_\text{Ele}c``, angles ``u`` between ``-1/9`` and ``1/9`` and angles ``h`` between ``0`` and ``2\pi``. Therefore initially **non-thermal and anisotropic**, but with zero bulk velocity. This population with have a number density ``n=1 \mathrm{m}^{-3}``, and can be generated using

````julia
Initial_Constant!(Initial,PhaseSpace,"Sph",pmin=10.0,pmax=13.0,umin=-0.11,umax=0.11,hmin=0.0,hmax=2.0,num_Init=1.0)
````

::: tip

The function `Initial_Constant!` modifies `Initial` but does not reset it, therefore if you call this function multiple times or additionally call other other initial condition functions, they will all add to each other, allowing a wide range of initial conditions to be applied simultaneously.

:::

::: warning

The way bin location are calculated rounds up, i.e. if a value is on a grid boundary then it is placed in the next bin.

:::

## Running the Solver

### Defining the Time-Stepping and Scheme

Before running the solver we need to tell it how to it will actually step through time, when to save results, and where to save the resulting file.

The time stepping is defined by two parameters: the initial timestep size `dt_initial`, and a `Vector` of times at which to save the simulation `t_save`:

````julia
dt_initial::Float64 = 1.0
t_save::Vector{Precision} = range(0.0,1e3,length=1001)
````

Here the initial timestep will be `1.0` characteristic times, and the simulation will run from ``t=0.0`` to ``t=10^3`` characteristic times, saving the results every timestep (1001, steps including the ``t=0.0``).

We will use a simple forward Euler method as the `scheme`, defined by `ForwardEulerStruct`:

````julia
scheme = ForwardEulerStruct(PhaseSpace,Initial,Injection,BinM,EmiM,FluxM)
````

File saving only requires us to define a `fileName` and `fileLocation` where the output will be saved. The output is saved in the JLD2 file format, a Julia native version of HDF5.

````julia
fileName = "tutorial1_output.jld2"
fileLocation = joinpath(pwd(),"Data")
````

### Actually Running the Solver

Everything is now all set up and we can run the solver using the `Solve` function:

````julia
sol = Solve(scheme,dt_initial,t_save;fileName=fileName,fileLocation=fileLocation,progress=true,Verbose=2)
````

If all went well, the output will be saved in `FileLocation` and we can now move on to plotting and examining the results!

