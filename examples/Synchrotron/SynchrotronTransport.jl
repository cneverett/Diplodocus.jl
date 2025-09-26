using Diplodocus

# ==== Define domains of time and space ====== #

    B = 1e-4; # Magnetic field strength in Tesla

    t_up::Float64 = log10(SIToCodeUnitsTime(1e8)) # seconds * (σT*c)
    t_low::Float64 = log10(SIToCodeUnitsTime(1e0)) # seconds * (σT*c)
    t_num::Int64 = 400
    t_grid::String = "l"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Spherical()# x = r, y = theta, z = phi

    x_up::Float64 = 1.0
    x_low::Float64 = 0f0
    x_grid::String = "u"
    x_num::Int64 = 1 

    y_up::Float64 = pi
    y_low::Float64 = 0.0
    y_grid::String = "u"
    y_num::Int64 = 1

    z_up::Float64 = 2.0*pi
    z_low::Float64 = 0.0
    z_grid::String = "u"
    z_num::Int64 = 1

    space = SpaceStruct(space_coords,x_up,x_low,x_grid,x_num,y_up,y_low,y_grid,y_num,z_up,z_low,z_grid,z_num)

# ==== Define particles and their momentum space domains ====== #

    name_list::Vector{String} = ["Ele","Pho"];

    momentum_coords = Spherical() # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [7.0,7.0];
    px_low_list::Vector{Float64} = [-3.0,-15.0];
    px_grid_list::Vector{String} = ["l","l"];
    px_num_list::Vector{Int64} = [80,88];

    py_up_list::Vector{Float64} = [1.0,1.0];
    py_low_list::Vector{Float64} = [-1.0,-1.0];
    py_grid_list::Vector{String} = ["u","u"];
    py_num_list::Vector{Int64} = [9,9];

    pz_up_list::Vector{Float64} = [2.0*pi,2.0*pi];
    pz_low_list::Vector{Float64} = [0.0,0.0];
    pz_grid_list::Vector{String} = ["u","u"];
    pz_num_list::Vector{Int64} = [1,1];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

# ==== Define Interactions  ====== #

    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[B],Iso())];
    Forces::Vector{ForceType} = [CoordinateForce(),SyncRadReact(Iso(),B),];

    # build phase space
    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd()*"\\examples\\Data"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

# Set initial conditions

    Initial = Initialise_Initial_Condition(PhaseSpace);
    Initial_PowerLaw!(Initial,PhaseSpace,"Ele",pmin=1e3,pmax=1e6,umin=-1.0,umax=1.0,hmin=0.0,hmax=2.0,index=2.0,num_Init=1e6);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "Sync.jld2";
    fileLocation = pwd()*"\\examples\\Data";

    sol = Solve(Initial,scheme;save_steps=5,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);

    MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,Static(),step=10,order=2,wide=true,plot_limits=((-15.0,7.0),(2.5,9.5)),TimeUnits=CodeToSIUnitsTime)

    MomentumDistributionPlot(sol,["Ele","Pho"],PhaseSpace,Animated(),thermal=false,order=2,plot_limits=((-15.0,7.0),(2.5,9.5)),wide=true,TimeUnits=CodeToSIUnitsTime,filename="SyncPDisAnimated.mp4")

    EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSIUnitsTime)

# ==== Saving plots for tutorial /paper ==== #

    #=

    SyncPDisPlotDark = MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),plot_limits=((-15.0,7.0),(2.5,9.5)))
    SyncPDisPlotLight = MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),plot_limits=((-15.0,7.0),(2.5,9.5)))
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.pdf",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.svg",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.pdf",SyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.svg",SyncPDisPlotLight)

    SyncEngDenPlotDark = EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark())
    SyncEngDenPlotLight = EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("SyncEngDenPlotDark.pdf",SyncEngDenPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncEngDenPlotDark.svg",SyncEngDenPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncEngDenPlotLight.pdf",SyncEngDenPlotLight)
    Diplodocus.DiplodocusPlots.save("SyncEngDenPlotLight.svg",SyncEngDenPlotLight)

    # ==== AM3 comparison plots ==== # 

    AM3SyncPDisPlotDark = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",1e8,1e0,"l",wide=true,plot_limits=((-15.0,7.0),(2.5,9.5)),theme=DiplodocusDark())
    AM3SyncPDisPlotLight = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",1e8,1e0,"l",wide=true,plot_limits=((-15.0,7.0),(2.5,9.5)),theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("AM3SyncPDisPlotDark.pdf",AM3SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3SyncPDisPlotDark.svg",AM3SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3SyncPDisPlotLight.pdf",AM3SyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3SyncPDisPlotLight.svg",AM3SyncPDisPlotLight)

    AM3DIPSyncPDisPlotDark = Diplodocus.DiplodocusPlots.AM3_DIP_Combo_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",sol,PhaseSpace,1e8,1e0,"l",plot_limits=((-15.0,7.0),(2.5,9.5)),theme=DiplodocusDark(),yticks=2:2:16)
    AM3DIPSyncPDisPlotLight = Diplodocus.DiplodocusPlots.AM3_DIP_Combo_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",sol,PhaseSpace,1e8,1e0,"l",plot_limits=((-15.0,7.0),(2.5,9.5)),theme=DiplodocusLight(),yticks=2:2:16)
    Diplodocus.DiplodocusPlots.save("AM3DIPSyncPDisPlotDark.pdf",AM3DIPSyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3DIPSyncPDisPlotDark.svg",AM3DIPSyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3DIPSyncPDisPlotLight.pdf",AM3DIPSyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3DIPSyncPDisPlotLight.svg",AM3DIPSyncPDisPlotLight)
    
    =#

    