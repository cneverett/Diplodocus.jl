using Diplodocus

# ==== Define domains of time and space ====== #

    #t_up::Float64 = 1e-11 # seconds * (σT*c)
    #t_low::Float64 = 0.0 # seconds * (σT*c)
    #t_num::Int64 = 100000
    #t_grid::String = "u"

    t_up::Float64 = -12e0 # seconds * (σT*c)
    t_low::Float64 = -20e0 # seconds * (σT*c)
    t_num::Int64 = 8000
    t_grid::String = "l"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Cylindrical(0.0,0.0,0.0)# x = r, y = phi, z = z

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

    name_list::Vector{String} = ["Ele","Pho"];

    momentum_coords = Spherical() # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [7.0,0.0];
    px_low_list::Vector{Float64} = [-3.0,-20.0];
    px_grid_list::Vector{String} = ["l","l"];
    px_num_list::Vector{Int64} = [80,80];

    py_up_list::Vector{Float64} = [1.0,1.0];
    py_low_list::Vector{Float64} = [-1.0,-1.0];
    py_grid_list::Vector{String} = ["u","u"];
    py_num_list::Vector{Int64} = [15,15];

    pz_up_list::Vector{Float64} = [2.0*pi,2.0*pi];
    pz_low_list::Vector{Float64} = [0.0,0.0];
    pz_grid_list::Vector{String} = ["u","u"];
    pz_num_list::Vector{Int64} = [1,1];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

# ==== Define Interactions  ====== #

    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[1e-4],Iso())];
    Forces::Vector{ForceType} = [SyncRadReact(Iso(),1e-4),];

    # build phase space
    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd()*"\\examples\\Synchrotron\\Data"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

# Set initial conditions

    Initial_Ele = Initial_PowerLaw(PhaseSpace,"Ele",1e3,1e6,-1.0,1.0,0.0,2.0,2.0,1e6);
    Initial_Pho = Initial_Constant(PhaseSpace,"Pho",1,80,1,15,1,1,0.0);
    Initial = ArrayPartition(Initial_Ele,Initial_Pho);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "SyncTest.jld2";
    fileLocation = pwd()*"\\examples\\Synchrotron\\Data";

    sol = Solve(Initial,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);

    MomentumAndPolarAngleDistributionPlot(sol,"Ele",PhaseSpace,(1,10,50),order=1)
    MomentumAndPolarAngleDistributionPlot(sol,"Pho",PhaseSpace,(1,10,50),order=1)

    MomentumDistributionPlot(sol,"Ele",PhaseSpace,step=10,order=2,plot_limits=(-3,7,0,10),wide=true)
    MomentumDistributionPlot(sol,"Pho",PhaseSpace,step=10,order=2,plot_limits=(-21,1,-3,10),wide=true)

    MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,step=80,order=2,wide=true,TimeUnits=SIUnitsTime)

    NumberDensityPlot(sol,PhaseSpace,species="Pho",theme=DiplodocusDark(),title=nothing)
    NumberDensityPlot(sol,PhaseSpace,species="Ele",theme=DiplodocusDark(),title=nothing)
    
    EnergyDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),title=nothing)


# ==== Saving plots for tutorial /paper ==== #

    SyncPDisPlotDark = MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,step=80,order=2,wide=true,TimeUnits=SIUnitsTime,theme=DiplodocusDark())
    SyncPDisPlotLight = MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,step=80,order=2,wide=true,TimeUnits=SIUnitsTime,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.pdf",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.pdf",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.pdf",SyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.svg",SyncPDisPlotLight)

    # ==== AM3 comparison plots ==== # 

    AM3ElePDisPlotDark = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test.jld2",SIUnitsTime(1e-12),SIUnitsTime(1e-20),"l",wide=true,plot_limits=(-21,8,-0.5345678329467773, 9.465432167053223),theme=DiplodocusDark())
    AM3ElePDisPlotLight = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test.jld2",SIUnitsTime(1e-12),SIUnitsTime(1e-20),"l",wide=true,plot_limits=(-21,9,-0.5345678329467773, 9.465432167053223),theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotDark.pdf",AM3ElePDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotDark.svg",AM3ElePDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotLight.pdf",AM3ElePDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotLight.svg",AM3ElePDisPlotLight)
