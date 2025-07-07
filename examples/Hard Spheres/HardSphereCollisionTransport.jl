using Diplodocus

# ==== Define domains of time and space ====== #

    t_up::Float64 = 2.5e3 # seconds * (σT*c)
    t_low::Float64 = 0.0 # seconds * (σT*c)
    t_num::Int64 = 250000
    t_grid::String = "u"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Cylindrical()  # x = r, y = phi, z = z

    x_up::Float64 = 1.0
    x_low::Float64 = 0f0
    x_grid::String = "u"
    x_num::Int64 = 1

    y_up::Float64 = 2.0*pi
    y_low::Float64 = 0.0
    y_grid::String = "u"
    y_num::Int64 = 1

    z_up::Float64 = 1.0
    z_low::Float64 = 0.0
    z_grid::String = "u"
    z_num::Int64 = 1

    space = SpaceStruct(space_coords,x_up,x_low,x_grid,x_num,y_up,y_low,y_grid,y_num,z_up,z_low,z_grid,z_num)

# ==== Define particles and their momentum space domains ====== #

    name_list::Vector{String} = ["Sph",];

    momentum_coords = Spherical() # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [4.0,];
    px_low_list::Vector{Float64} = [-5.0,];
    px_grid_list::Vector{String} = ["l",];
    px_num_list::Vector{Int64} = [72,];

    py_up_list::Vector{Float64} = [1.0,];
    py_low_list::Vector{Float64} = [-1.0,];
    py_grid_list::Vector{String} = ["u",];
    py_num_list::Vector{Int64} = [8,];

    pz_up_list::Vector{Float64} = [2.0*pi,];
    pz_low_list::Vector{Float64} = [0.0,];
    pz_grid_list::Vector{String} = ["u",];
    pz_num_list::Vector{Int64} = [1,];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

    Binary_list::Vector{BinaryStruct} = [BinaryStruct("Sph","Sph","Sph","Sph")];
    Emi_list::Vector{EmiStruct} = [];
    Forces::Vector{ForceType} = [];

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# ==== Build Interaction and Flux Matrices ====== #

    DataDirectory = pwd()*"\\examples\\Hard Spheres\\Data"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=false);
    FluxM = BuildFluxMatrices(PhaseSpace);

# ===== Set initial conditions ================== #

    Initial_Sph = Initial_Constant(PhaseSpace,"Sph",3.0,3.1,-0.24,0.25,0.0,2.0,1f0);
    f_init = ArrayPartition(Initial_Sph,);

# Run BoltzmannEquationSolver

    scheme = EulerStruct(f_init,PhaseSpace,BigM,FluxM,false)

    fileName = "HardSphere.jld2";
    fileLocation = pwd()*"\\examples\\Hard Spheres\\Data";

    sol = Solve(f_init,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);

# Plot results 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);

# Get scaling for plots (initial values)

    (numInit_list,engInit_list,tempInit_list) = Diplodocus.DiplodocusTransport.getInitialScaling(sol,PhaseSpace);
    
# Plotting

    Diplodocus.DiplodocusTransport.AllPlots_Ani(sol,PhaseSpace,numInit_list,engInit_list,tempInit_list,"test1.mp4";fps=12,istart=1,istop=1000,iframe=nothing,step=5)

    #DT.PDistributionPlot_AllSpecies(sol,num_species,name_list,meanp_list,dp_list,du_list,tempInit_list,mass_list,p_num_list,u_num_list,dt*100)

    NumberDensityPlot(sol,PhaseSpace,theme=DP.DiplodocusLight())
    fig = FracNumberDensityPlot(sol,PhaseSpace,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("HardSphereFracNumLight.png",fig)

    EnergyDensityPlot(sol,PhaseSpace)
    FracEnergyDensityPlot(sol,PhaseSpace)

    MomentumDistributionPlot(sol,"Sph",PhaseSpace,step=1,uDis=false,logt=false,thermal=true,order=1)


    fig = MomentumDistributionPlot(sol,"Sph",PhaseSpace,uDis=false,step=100,logt=false,plot_limits=(-4.0,5.0,-10.0,0.0))
    #DT.save("Amsterdam/SyncAmstPlots/PhoUavg.png",fig)


    

    Diplodocus.DiplodocusTransport.FracNumPlot(sol,"Sph",PhaseSpace,fig=nothing)
    DT.FracNumPlot(sol,"Ele",PhaseSpace,fig=nothing)