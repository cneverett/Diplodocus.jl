# define constants with which to import data 

using Diplodocus

# ==== Define domains of time and space ====== #

    B = 1e-4; # Magnetic field strength in Tesla

    t_up::Float64 = SyncToCodeUnitsTime(1.0,B=B) # seconds * (σT*c)
    t_low::Float64 = SyncToCodeUnitsTime(0.0,B=B)# seconds * (σT*c)
    t_num::Int64 = 1000
    t_grid::String = "u"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Cylindrical() # x = r, y = phi, z = z

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

    name_list::Vector{String} = ["Ele",];

    momentum_coords = Spherical()  # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [5.0,];
    px_low_list::Vector{Float64} = [-10.0,];
    px_grid_list::Vector{String} = ["l",];
    px_num_list::Vector{Int64} = [480,];

    py_up_list::Vector{Float64} = [1.0,];
    py_low_list::Vector{Float64} = [-1.0,];
    py_grid_list::Vector{String} = ["u",];
    py_num_list::Vector{Int64} = [65,];

    pz_up_list::Vector{Float64} = [2.0*pi,];
    pz_low_list::Vector{Float64} = [0.0,];
    pz_grid_list::Vector{String} = ["u",];
    pz_num_list::Vector{Int64} = [1,];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [];
    Forces::Vector{ForceType} = [CoordinateForce(),SyncRadReact(Ani(),B),];

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# ==== Build Interaction and Flux Matrices ====== #

    DataDirectory = pwd() * "/examples/Data/"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

# ===== Set Initial Conditions ================== #

    function pth_to_T(pth)
        m = 9.10938356e-31
        c = Float64(299792458)
        kb = 1.380649e-23
        return m*c^2*pth^2/kb
    end
    Initial = Initialise_Initial_Condition(PhaseSpace);
    Initial_MaxwellJuttner!(Initial,PhaseSpace,"Ele",T=pth_to_T(2.0),umin=-1.0,umax=1.0,hmin=0.0,hmax=2.0,num_Init=1e0);

# ===== Run the Solver ================== #

    fileLocation = pwd() * "/examples/Data/";
    fileName = "RadReact.jld2"

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    sol = Solve(Initial,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);

    MomentumDistributionPlot(sol,["Ele"],PhaseSpace,Static(),step=33,order=-2,paraperp=true,plot_limits=((-3.4,1.9),(-5.9,-0.1)),TimeUnits=CodeToSyncUnitsTime)

    MomentumAndPolarAngleDistributionPlot(sol,"Ele",PhaseSpace,Static(),(1,52,102),order=-2,TimeUnits=CodeToSyncUnitsTime)
    
    FracNumberDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSyncUnitsTime)
    EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSyncUnitsTime)

    MomentumComboAnimation(sol,["Ele"],PhaseSpace;plot_limits_momentum=((-3.4,1.9),(-5.9,-0.1)),order=-2,thermal=false,paraperp=true,initial=false,filename="RadReactMomentumComboAnimation.mp4",TimeUnits=CodeToSyncUnitsTime)

# ====== Saving Plots for tutorial/paper ======= #

    #=
    
    PDisPlotDark = MomentumDistributionPlot(sol,["Ele"],PhaseSpace,Static(),step=33,thermal=false,paraperp=true,order=-2,plot_limits=((-3.4,1.9),(-5.9,-0.1)),theme=DiplodocusDark(),TimeUnits=CodeToSyncUnitsTime)
    PDisPlotLight = MomentumDistributionPlot(sol,["Ele"],PhaseSpace,Static(),step=33,thermal=false,paraperp=true,order=-2,plot_limits=((-3.4,1.9),(-5.9,-0.1)),theme=DiplodocusLight(),TimeUnits=CodeToSyncUnitsTime)

    Diplodocus.DiplodocusPlots.save("RadReactPDisPlotDark.svg",PDisPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactPDisPlotDark.pdf",PDisPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactPDisPlotLight.svg",PDisPlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactPDisPlotLight.pdf",PDisPlotLight)

    PAndUDisPlotDark = MomentumAndPolarAngleDistributionPlot(sol,"Ele",PhaseSpace,Static(),(1,52,102),order=-2,TimeUnits=CodeToSyncUnitsTime,theme=DiplodocusDark())
    PAndUDisPlotLight = MomentumAndPolarAngleDistributionPlot(sol,"Ele",PhaseSpace,Static(),(1,52,102),order=-2,TimeUnits=CodeToSyncUnitsTime,theme=DiplodocusLight())

    Diplodocus.DiplodocusPlots.save("RadReactPAndUDisPlotDark.svg",PAndUDisPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactPAndUDisPlotDark.pdf",PAndUDisPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactPAndUDisPlotLight.svg",PAndUDisPlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactPAndUDisPlotLight.pdf",PAndUDisPlotLight)

    MomentumComboAnimation(sol,["Ele"],PhaseSpace;plot_limits_momentum=((-3.4,1.9),(-5.9,-0.1)),order=-2,thermal=false,paraperp=true,initial=false,filename="RadReactMomentumComboAnimation.mp4",TimeUnits=CodeToSyncUnitsTime)

    FracNumPlotDark = FracNumberDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),TimeUnits=CodeToSyncUnitsTime)
    EngPlotDark = EnergyDensityPlot(sol,species="Ele",PhaseSpace,theme=DiplodocusDark(),TimeUnits=CodeToSyncUnitsTime)

    FracNumPlotLight = FracNumberDensityPlot(sol,species="Ele",PhaseSpace,theme=DiplodocusLight(),TimeUnits=CodeToSyncUnitsTime)
    EngPlotLight = EnergyDensityPlot(sol,species="Ele",PhaseSpace,theme=DiplodocusLight(),TimeUnits=CodeToSyncUnitsTime)

    Diplodocus.DiplodocusPlots.save("RadReactFracNumPlotDark.svg",FracNumPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactFracNumPlotDark.pdf",FracNumPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactEngPlotDark.svg",EngPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactEngPlotDark.pdf",EngPlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactFracNumPlotLight.svg",FracNumPlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactFracNumPlotLight.pdf",FracNumPlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactEngPlotLight.svg",EngPlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactEngPlotLight.pdf",EngPlotLight)

    TimePlotLight = Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol.f[1],1;paraperp=true,p_timescale=true,wide=false,legend=false,horz_lines=[1.0,0.667,0.33],plot_limits=((-4,2),(-2,2)),TimeUnits=CodeToSyncUnitsTime,theme=DiplodocusLight())
    TimePlotDark = Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol.f[1],1;paraperp=true,p_timescale=true,wide=false,legend=false,horz_lines=[1.0,0.667,0.33],plot_limits=((-4,2),(-2,2)),TimeUnits=CodeToSyncUnitsTime,theme=DiplodocusDark())

    Diplodocus.DiplodocusPlots.save("RadReactTimePlotLight.svg",TimePlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactTimePlotLight.pdf",TimePlotLight)
    Diplodocus.DiplodocusPlots.save("RadReactTimePlotDark.svg",TimePlotDark)
    Diplodocus.DiplodocusPlots.save("RadReactTimePlotDark.pdf",TimePlotDark)

    =#

