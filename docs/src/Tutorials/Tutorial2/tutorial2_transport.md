```@meta
EditURL = "../literate_scripts/Tutorial2/tutorial2_transport.jl"
```

# Tutorial 2c: Running the Simulation

---

With the generation of the Compton and Synchrotron matrices, we are now ready to run the simulation. This simulation will be run using the functions contained within the `DiplodocusTransport` sub-package.

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
CHAR_length::Float64 = 1e14 # m = 1e16 cm
````

Here we have defined the characteristic length to be the radius of our spherical region.

::: tip

The characteristic speed is always taken to be the speed of light, so the characteristic time for this simulation is then defined as ``t=CHAR_length/CONST_c`` where ``CONST_c`` is the speed of light in SI units.

:::

## Phase-Space Setup

### Spatial Grids

We need to define the spatial geometry on which we are going to simulate. In this tutorial we are considering a single spherical region. This region will therefore have `Spherical` coordinates. By default, the radial boundary condition for `Spherical` coordinates is taken to be `Closed` i.e. no particle may leave. We can change this by modifying the keyword parameter `xp_BC=Escape()`, when defining the `space_coords`. As out particles are isotropic we must use `Escape` which ignores the particles direction and just uses their momentum to determine whether they escape or not, as otherwise if we used `Open` the angle averaged direction of propagation is zero and no particles would leave the zone.

````julia
space_coords::CoordinateType = Spherical(xp_BC=Escape())

x_up::Float64 = 1.0
x_low::Float64 = 0.0
x_grid::String = "u"
x_num::Int64 = 1

y_up::Float64 = pi
y_low::Float64 = 0.0
y_grid::String = "u"
y_num::Int64 = 1

z_up::Float64 = 2pi
z_low::Float64 = 0.0
z_grid::String = "u"
z_num::Int64 = 1
````

No matter the type of spatial coordinates, they are always mapped to ``x,y,z`` internally, so for spherical coordinates ``x=r``, ``y=\theta``, and ``z=\psi`` (for more detail see [Implemented Coordinates and Forces](@ref)).

### Momentum Grids

Next we similarly define the momentum space grids and the list of particles to be included in the simulation. For momentum space, Diplodocus uses a local spherical coordinate basis, therefore, the momentum coordinates should be `Spherical`. As for the grids for each coordinate, we should choose the same values as were used to generate the collision matrices.

````julia
name_list::Vector{String} = ["Ele","Pho"];

momentum_coords = Spherical() # px = p, py = u, pz = phi

px_up_list::Vector{Float64} = [9.0,9.0]
px_low_list::Vector{Float64} = [-3.0,-16.0]
px_grid_list::Vector{String} = ["l","l"]
px_num_list::Vector{Int64} = [288,600]

py_up_list::Vector{Float64} = [1.0,1.0]
py_low_list::Vector{Float64} = [-1.0,-1.0]
py_grid_list::Vector{String} = ["u","u"]
py_num_list::Vector{Int64} = [1,1]

pz_up_list::Vector{Float64} = [2.0*pi,2.0*pi]
pz_low_list::Vector{Float64} = [0.0,0.0]
pz_grid_list::Vector{String} = ["u","u"]
pz_num_list::Vector{Int64} = [1,1]
````

::: warning

Momentum-space grids must match those used in generating the collision matrices

:::

### Define the Electromagnetic field

::: warning

This process may change in future version.

:::

We can tell the code what the magnetic field is as a function of position according to the local-orthonormal frame that is used for momentum-space. This ortho-normal frame is always defined with respect to the magnetic field direction to ensure accurate synchrotron emissions (see [Implemented Coordinates and Forces](@ref)). In our region this is just a constant field of ``B=1e-3\,\mathrm{T}``:

````julia
ElectroMagneticField::ElectroMagneticFieldStruct = ElectroMagneticField_Constant(parameters=[1e-3,1e0])
````

We will later tell the solver to treat this field as isotropic when calculating synchrotron emission.

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
Binary_list::Vector{BinaryStruct} = [BinaryStruct("Ele","Pho","Ele","Pho")]
Binary_Domain = CollisionDomain(PhaseSpace)
````

Here the `CollisionDomain` is over all phase space, i.e. the single spatial cell of our simulation.

We now need to tell the solver what emissive interaction to include. For step has a few more options. First we must list out the values of the external parameter with which we sampled the synchrotron emission matrices, here `B_sample`. Then we define the domain overwhich to include the synchrotron `Sync_Domain` (`CollisionDomain() defaults to everywhere) and then build an `Emission_list` that is a `Vector` of `EmiStruct`s:

````julia
B_sampled::Vector{Float64} = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1];
Sync_Domain = CollisionDomain(PhaseSpace);
Emission_list::Vector{EmiStruct} = [EmiStruct(name1="Ele",name2="Ele",name3="Pho",EmiName="Sync",Ext_sampled=B_sampled,mode=Iso(),Force=true,Domain=Sync_Domain)]
````

The `EmiStruct` takes a `mode` argument, which we set to `Iso()` to tell the code to treat the magnetic field as isotropic when calculating synchrotron emissions (this is technically unneeded here as our momentum domains are isotropic, but can be used to replicate turbulent, tangled, magnetic fields with anisotropic momentum grids). We also set `Force=true` to tell the code to apply the corresponding synchrotron cooling force (radiation reaction)to the electrons. This will also correct the synchrotron emission spectra to ensure that energy is conserved between cooling and emission.

Last we specify that there are no other forces to include:

````julia
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
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=true); # jl
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=true) # md

EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true); # jl
EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true) # md

FluxM = BuildFluxMatrices(PhaseSpace,Forces_list); # jl
FluxM = BuildFluxMatrices(PhaseSpace,Forces_list) # md
````

Here we have used the argument `Bin_corrected=true` to use the number and energy density correct binary collision matrices, and the argument `Bin_sparse=true` to define the binary collision matrix as a sparse matrix.
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

These two state vectors are currently empty, but we can use several functions to add to them a variety of initial and injection conditions ... LINK TO INIT/INJ OPTIONS ... .
For this tutorial, we don't need to set any initial electrons but we do need to set some to be injected.
This injected population is of the form of a power-law, so we can use `Injection_PowerLaw!`:

````julia
Injection_PowerLawExpDecay!(Injection,PhaseSpace,"Ele",pmin=10^1.9,pmax=10^2.1,index=2.0,num_Inj=1e3,rate_Inj=1.0)
````

This will inject a power-law distribution of electrons with minimum momentum of ``p=10^{1.9}``, maximum momentum of ``p=10^{2.1}``, and a power-law index of `2`. The total number of particles injected is set by `num_Inj` and the rate at which they are injected is set by `rate_Inj`, both given in characteristic units.

::: tip

The function `Injection_PowerLaw!` modifies `Injection` but does not reset it, therefore if you call this function multiple times or additionally call other other injection condition functions, they will all add to each other, allowing a wide range of injection conditions to be applied simultaneously.

:::

## Running the Solver

### Defining the Time-Stepping and Scheme

Before running the solver we need to tell it how to it will actually step through time, when to save results, and where to save the resulting file.

The time stepping is defined by two parameters: the initial timestep size `dt_initial`, and a `Vector` of times at which to save the simulation `t_save`:

````julia
dt_initial::Float64 = 0.1
t_save::Vector{Precision} = range(0.0,1e1,length=101)
````

Here the initial timestep will be `1.0` characteristic times, and the simulation will run from ``t=0`` to ``t=10`` characteristic times, saving the results every timestep (101, steps including the ``t=0``).

We will use a simple forward Euler method as the `scheme`, defined by `ForwardEulerStruct`:

````julia
scheme = ForwardEulerStruct(PhaseSpace,Initial,Injection,BinM,EmiM,FluxM)
````

File saving only requires us to define a `fileName` and `fileLocation` where the output will be saved. The output is saved in the JLD2 file format, a Julia native version of HDF5.

````julia
fileName = "tutorial2_output.jld2"
fileLocation = joinpath(pwd(),"Data")
````

### Actually Running the Solver

Everything is now all set up and we can run the solver using the `Solve` function:

````julia
sol = Solve(scheme,dt_initial,t_save;fileName=fileName,fileLocation=fileLocation,progress=true,Verbose=2)
````

If all went well, the output will be saved in `FileLocation` and we can now move on to plotting and examining the results!

