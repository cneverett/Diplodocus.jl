# define constants with which to import data 

using Diplodocus

# Define space and time domain 

    # time domain

        t_up::Float64 = 1e-11 # seconds * (σT*c)
        t_low::Float64 = 0.0 # seconds * (σT*c)
        t_num::Int64 = 100000
        t_grid::String = "u"

        t_up::Float64 = -10e0 # seconds * (σT*c)
        t_low::Float64 = -18e0 # seconds * (σT*c)
        t_num::Int64 = 1000
        t_grid::String = "l"

        time = TimeStruct(t_up,t_low,t_num,t_grid)

    # space domain

        space_coords = Cylindrical(0e0,0.0,0e0)

        # cartesian x, cylindrical r, spherical r
        x_up::Float64 = 1.0
        x_low::Float64 = 0f0
        x_grid::String = "u"
        x_num::Int64 = 1
        # cartesian y, cylindrical phi, spherical u=cos(theta) 
        y_up::Float64 = 2.0*pi
        y_low::Float64 = 0.0
        y_grid::String = "u"
        y_num::Int64 = 1
        # cartesian z, cylindrical z, spherical phi
        z_up::Float64 = 1.0
        z_low::Float64 = 0.0
        z_grid::String = "u"
        z_num::Int64 = 1

        space = SpaceStruct(space_coords,x_up,x_low,x_grid,x_num,y_up,y_low,y_grid,y_num,z_up,z_low,z_grid,z_num)

# define particle names to include in the simulation as a vector i.e. [name1, name2, ...]

    name_list::Vector{String} = ["Ele",];

# define momentum space domain for each particle as a vector

    momentum_coords = Spherical()

    # cartesian x, cylindrical r, spherical r
        px_up_list::Vector{Float64} = [7.0,];
        px_low_list::Vector{Float64} = [-7.0,];
        px_grid_list::Vector{String} = ["l",];
        px_num_list::Vector{Int64} = [112,];
    # cartesian y, cylindrical phi, spherical u=cos(theta) 
        py_up_list::Vector{Float64} = [1.0,];
        py_low_list::Vector{Float64} = [-1.0,];
        py_grid_list::Vector{String} = ["u",];
        py_num_list::Vector{Int64} = [16,];
    # cartesian z, cylindrical z, spherical phi
        pz_up_list::Vector{Float64} = [2.0*pi,];
        pz_low_list::Vector{Float64} = [0.0,];
        pz_grid_list::Vector{String} = ["u",];
        pz_num_list::Vector{Int64} = [1,];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

# tuple of vectors, each vector representing [1, 2, 3, 4] particles in the interactions to include in the simulation i.e. [[interaction1], [interaction2], ...]

    #interaction_list_Binary::Vector{Vector{String}} = [["Ele","Pho","Ele","Pho"],];
    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [];

# define forces to include in the simulation as a vector i.e. [force1, force2, ...]

    forces::Vector{ForceType} = [SyncRadReact(Axi()),];

# build phase space struct to pass to solver

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd() * "/examples/Synchrotron/Data/"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);

    FluxM = BuildFluxMatrices(PhaseSpace);

# Set initial conditions

    # Maxwell-Juttner distribution of electrons
    initial_low = Initial_MaxwellJuttner(PhaseSpace,"Ele",1e5,1e6);
    initial_high = Initial_MaxwellJuttner(PhaseSpace,"Ele",1e15,1e6);


# Solver Setup and Time Stepping

    low = ArrayPartition(initial_low,);
    high = ArrayPartition(initial_high,);

# Run Transport Solver

    scheme = EulerStruct(low,PhaseSpace,BigM,FluxM,false)

    fileName = "RadReact_low.jld2";
    fileLocation = pwd() * "/examples/Synchrotron/Data/";

    sol = Solve(low,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);

# Plot results 

    (PhaseSpace, sol) = Diplodocus.DiplodocusTransport.fload(fileLocation,fileName);

# Get scaling for plots (initial values)

    (numInit_list,engInit_list,tempInit_list) = Diplodocus.DiplodocusTransport.getInitialScaling(sol,PhaseSpace);
    
# Plotting

    Diplodocus.DiplodocusTransport.AllPlots_Ani(sol,PhaseSpace,numInit_list,engInit_list,tempInit_list,"test1.mp4";fps=24,istart=1,istop=100,iframe=nothing,step=5)

    #DT.PDistributionPlot_AllSpecies(sol,num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,dt*100)


    fig = MomentumDistributionPlot(sol,"Ele",PhaseSpace,uDis=false,step=1,logt=false,plot_limits=((-8.0,8.0),(-10,10.0)))
    Diplodocus.DiplodocusTransport.save(fileLocation*"/LowTemperatureMomentumDistribution.png",fig)

    

    DT.FracNumPlot(sol,"Pho",PhaseSpace,fig=nothing)
    Diplodocus.DiplodocusTransport.FracNumPlot(sol,"Ele",PhaseSpace,fig=nothing)

##### Observer plotting



tes = DT.ObserverFlux(PhaseSpace,sol,[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],1.0);



DT.ObserverFluxPlot(PhaseSpace,sol,[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],1.0,plot_limits=((-15,2),(-8.0,10.0)))