# # Tutorial 2a: Generating the Compton Matrices
#
# ---
#
# Before we examine the *inverse-Compton catastrophe* we need to pre-compute some Compton scattering collision matrices. This pre-computation is conducted by the functions contained within the `DiplodocusCollisions` sub-package of `Diplodocus`, which generate collision matrices that describe the rate at which some incoming state of particles will be scattered into some outgoing state.
#
# This process is required whenever you want to include a *binary interaction* in your simulation, that is any interaction that involves two incoming particles and two outgoing particles, e.g. two incoming spheres that collide to create two outgoing spheres.
#
# # Generating the Collision Matrices
#
# ## Particle Names
#
# For Compton scattering, we are considering an incoming electron and photon scattering to produce a scattered electron and scattered photon. In Diplodocus, each particle species is designated a unique three letter `string` to identify it. These unique identifiers are then used internally to define the properties of the particles (please refer to [Particles, Grids and Units](@ref) for more details). In binary interactions ``12\rightarrow34``, we follow the convention that the pairs of names (12) and (34) should be in alphabetical order. Therefore for the Compton scattering between electrons and photons is declared as (see [Implemented Collisions](@ref)):

using Diplodocus

name1 = "Ele"
name2 = "Pho"
name3 = "Ele"
name4 = "Pho"

# ## Momentum Space Grids
# 
# Next is to define how momentum space is to be discretised for each particle species. This includes the upper and lower bounds of momentum magnitude ``p``, upper and lower bounds for polar angle cosine ``u`` and azimuthal angle ``h`` are not requrired as they are assume to be bounded by ``[-1,1]`` and ``[0, 2\pi]``), grid types (see [Particles, Grids and Units](@ref)) and the number of grid bins. These must be of the format `p_low_name`, `p_low_name`, `p_grid_name`, `p_num_name`, `u_grid_name`, `u_num_name`, `h_grid_name` and  `h_num_name` where `name` is the abbreviated three letter name of the particle species.
# 
# For this tutorial we use a generous range of momenta such that the collision matrices we will generate could be re-used for a wider range of application. The upper bound in momenta for both electrons and photons will be will ``10^{9}m_\text{Ele}c``, with the lower bound for electrons being ``10^{-3}m_\text{Ele}c`` and ``10^{-16}m_\text{Ele}c`` for photons. We will a logarithmic (base 10) grid with 24 bins per decade for both electron and photon momenta and assume isotropy so we only include 1 bin in both ``u`` and ``\phi`` (`h`).

    p_low_Ele = -3.0 
    p_up_Ele = 9.0
    p_grid_Ele = "l"
    p_num_Ele = 288
    u_grid_Ele = "u"
    u_num_Ele = 1
    h_grid_Ele = "u"
    h_num_Ele = 1

    p_low_Pho = -16.0
    p_up_Pho = 9.0
    p_grid_Pho = "l"
    p_num_Pho = 600
    u_grid_Pho = "u"
    u_num_Pho = 1
    h_grid_Pho = "u"
    h_num_Pho = 1

# ## Monte-Carlo Setup
#
# The elements of collision matrices are generated using a Monte-Carlo method as they are integrals of high dimension. We must define how many times each momentum bin is to be sampled by this Monte-Carlo approach. This is defined by three values: `numLoss` defines the number of incoming states (i.e. sets of ``(\vec{p}_1,\vec{p}_2)`` vectors, which in this case define the momentum vectors for the two incoming spheres about to collide) to sample per incoming state momentum bin, `numGain` defines the number of outgoing states (i.e. ``\vec{p}_3`` and ``\vec{p}_4)``, defining the momentum of the spheres after colliding) to sample for each outgoing state bin in momentum space, and `numThreads` defines the number of threads to use (see [the documentation for multi-threading in Julia](https://docs.julialang.org/en/v1/manual/multi-threading/)).
#
# ::: tip
# 
# The larger `numLoss` and `numGain` are, the more  accurate the integration will be, but this will take more time. 
#
# :::

numLoss = 16
numGain = 16
numThreads = Threads.nthreads()

# ::: warning
#
# By default this script uses the maximum number of available threads `Threads.nthreads()`, if you would prefer the scripts to run with fewer threads then modify the variable `numThreads` to the desired amount. 
#
# :::
#
# Further, to improve the Monte-Carlo sampling, we can deploy importance sampling by defining a non-zero `scale` factor. This `scale` factor weights the sampling of outgoing state towards the centre-of-momentum direction. Useful for high energy scatterings where particles are preferentially scattered towards that direction. As the user, you can set this scale factor to take a range of values `a:b:c` where `a` is its minimum value, `b` is the step size and `c` its maximum value (a multiple of `b`).

scale = 0.0:0.1:1.0

# ::: tip
#
# For the first round of integration it is advised to run with a scale of zero i.e. `scale = 0.0:0.1:0.0` as weighting may not be needed for accurate integration. If inaccuracy is found, the scale should be increased gradually. 
#
# :::

# ## Where to Save the Collision Matrices
#
# We now need to define where to store the generated collision matrices. This can be defined by specifying `fileLocation` as the path to the desired folder. If left unspecified, this location will default to `joinpath(pwd(),"Data")` i.e. a folder named `Data` found in the current directory. If the specified `fileLocation` does not exist, it will be created when this script is run.

fileLocation = joinpath(pwd(),"Data")

# ## Final Setup and Running
# 
# We can bundle all the above user defined parameters into a tuple called `Setup`. This is generated by the `UserBinaryParameters()` function, which reads the above parameters from `Main`. It also generates the `fileName` under which the collision matrices will be saved and stored.

(Setup,fileName) = UserBinaryParameters()

# ::: warning
# 
# It is not advised to edit `fileName` as its format is used internally by `DiplodocusCollisions` and `DiplodocusTransport` to identify the appropriate collision matrices to load.
#
# :::
#
# Finally, we can run the Monte-Carlo integration using the function `BinaryInteractionIntegration`:

BinaryInteractionIntegration(Setup)

# ::: tip
#
# Even on a good system, integration typically takes several hours if not days. This is not because the MC sampler isn't optimised (it only takes around 100ns per sample) but more because there are a lot of bins and points within each bin to sample. But once it's done the matrices can be used over and over again. So don't worry, and while you wait perhaps sit back, relax and read that book  that's been on your shelf for years just waiting to be read.
#
# :::
#
# # Checking the Collision Matrices
# 
# The evaluation of the collision matrices relies on Monte-Carlo sampling, therefore accuracy isn't guaranteed. There are several functions within Diplodocus aimed to assist with determining if the integration accuracy is sufficient for use. 
#
# ## Plotting the Spectra
#
# The most visual approach is to use the plotting function `InteractiveBinaryGainLossPlot` from `DiplodocusPlots`. This creates an interactive plot in which the user can look at different incoming states and view the resulting outgoing spectra. To use this function, first load the file containing the collision matrices using `BinaryFileLoad_Matrix`:

Output = BinaryFileLoad_Matrix(fileLocation,fileName);
InteractiveBinaryGainLossPlot(Output)

# This will display the plot in a separate window where you can interact with it. There are sliders for the the discrete bins of the two incoming particles momentum `p1` and `p2`, and their propagation directions, in spherical polars `u` (``u=\cos\theta``) and `h` (``phi``). There are then also sliders for the outgoing particles directions with the spectra of outgoing states as a function of momentum being plotted.
#
# ::: tip
#
# If a specific range of incoming states seems to be particularly poorly sampled, you can focus on that regime by specifying the input parameters `p1loc_low`, `p1loc_up`, `p2loc_low`, `p2loc_up`, these can be set to any `Int64` values that specify range of incoming momentum bins to sample from `low` to `up`. The integrator will then not sample incoming states from outside this region.
#
# :::
#
# Another way to judge the accuracy of integration is through the `DoesConserve` function. This prints to the terminal a series of statistics about the integration:

DoesConserve(Output)

# It includes, the total gain and loss of particles as well as the total gain and loss of energy (note as all particles are identical the gain and loss of particles 2 and 4 are neglected as they are identical to 1 and 3).  Nevertheless, the ratio of particles and energy in and out `ratioN` and `ratioE` respectively should be around 1.0 for good integration. These values give a sense of overall convergence but also provided are statistics for individual incoming to outgoing states. For a single incoming state (one set of ``p1,u1,h1,p2,u2,h2`` bins) an error in particle number and energy can be calculated by comparing the loss matrix element for that state to the sum of all outgoing gain matrix elements that correspond to that incoming state. `mean error in N` then gives the mean of that error and `std of error in N` gives the standard deviation.
#
# ## Corrected Spectra
#  
# There are two ways to improve the accuracy of the MC sampling: first is to just increase the number of sampled points by re-running the integration as described above. The second is to apply a correction to the outgoing spectra to enforce conservation of number and energy. This latter method guarantees conservation to numerical precision and is already built into `DiplodocusCollisions`. In brief, it applies two corrective scalings, one to the entire spectrum and one to the highest momentum bins. This ensures that the total gain of outgoing particles and their energy matches that of what is lost from the incoming states, which maintaining the shape of the spectra as much as possible.
#
# This corrective step is applied as part of the integration and the resulting collision matrices can be accessed by setting `corrected=true` as an option of `BinaryFileLoad_Matrix` e.g. 

Output = BinaryFileLoad_Matrix(fileLocation,fileName,corrected=true); 
DoesConserve(Output)

# Now the conservation statistics for particle number and energy are accurate to machine precision, allowing much longer simulations to be run without worrying about conservation issues.
#
# Next we need to do a similar process to generate the emission matrices for synchrotron.