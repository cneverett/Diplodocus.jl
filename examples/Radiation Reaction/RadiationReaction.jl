# define constants with which to import data 

using Diplodocus

# ==== Define domains of time and space ====== #

    tau_to_t::Float64 = 3.644e-13

    t_up::Float64 = 3.0tau_to_t # seconds * (σT*c)
    t_low::Float64 = 0.0tau_to_t # seconds * (σT*c)
    t_num::Int64 = 15000
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
    px_num_list::Vector{Int64} = [1920,];

    py_up_list::Vector{Float64} = [0.001,];
    py_low_list::Vector{Float64} = [-0.001,];
    py_grid_list::Vector{String} = ["u",];
    py_num_list::Vector{Int64} = [1,];

    pz_up_list::Vector{Float64} = [2.0*pi,];
    pz_low_list::Vector{Float64} = [0.0,];
    pz_grid_list::Vector{String} = ["u",];
    pz_num_list::Vector{Int64} = [1,];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [];
    Forces::Vector{ForceType} = [SyncRadReact(Ani(),1e-4),];

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# ==== Build Interaction and Flux Matrices ====== #

    DataDirectory = pwd() * "/examples/Radiation Reaction/Data/"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

# ===== Set Initial Conditions ================== #

    function pth_to_T(pth)
        return 5.9e9*pth^2
    end

    Initial_low = Initial_MaxwellJuttner(PhaseSpace,"Ele",pth_to_T(0.5),1,1,1,1,1e0);
    Initial_med = Initial_MaxwellJuttner(PhaseSpace,"Ele",pth_to_T(1.0),1,1,1,1,1e0);
    Initial_high = Initial_MaxwellJuttner(PhaseSpace,"Ele",pth_to_T(1.5),1,1,1,1,1e0);
    low = ArrayPartition(Initial_low,);
    med = ArrayPartition(Initial_med,);
    high = ArrayPartition(Initial_high,);

    Initial_test = Initial_MaxwellJuttner(PhaseSpace,"Ele",pth_to_T(0.5),1,1,1,1,1e6);
    test = ArrayPartition(Initial_test,);

# ===== Run the Solver ================== #

    fileLocation = pwd() * "/examples/Radiation Reaction/Data/";
    fileName_low = "RadReact_low.jld2";
    fileName_med = "RadReact_med.jld2";
    fileName_high = "RadReact_high.jld2";

    scheme_low = EulerStruct(low,PhaseSpace,BigM,FluxM,false)
    sol_low = Solve(low,scheme_low;save_steps=1000,progress=true,fileName=fileName_low,fileLocation=fileLocation);

    scheme_med = EulerStruct(med,PhaseSpace,BigM,FluxM,false)
    sol_med = Solve(med,scheme_med;save_steps=1000,progress=true,fileName=fileName_med,fileLocation=fileLocation);

    scheme_high = EulerStruct(high,PhaseSpace,BigM,FluxM,false)
    sol_high = Solve(high,scheme_high;save_steps=1000,progress=true,fileName=fileName_high,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace_low, sol_low) = SolutionFileLoad(fileLocation,fileName_low);
    (PhaseSpace_med, sol_med) = SolutionFileLoad(fileLocation,fileName_med);
    (PhaseSpace_high, sol_high) = SolutionFileLoad(fileLocation,fileName_high);

    MomentumAndPolarAngleDistributionPlot(sol_low,"Ele",PhaseSpace_low,(1.0tau_to_t,2.0tau_to_t,3.0tau_to_t),order=1)
    MomentumAndPolarAngleDistributionPlot(sol_med,"Ele",PhaseSpace_med,(0.0tau_to_t,0.5tau_to_t,1.0tau_to_t),order=1)
    MomentumAndPolarAngleDistributionPlot(sol_high,"Ele",PhaseSpace_high,(0.0tau_to_t,0.5tau_to_t,1.0tau_to_t),order=1)

    MomentumDistributionPlot(sol_low,"Ele",PhaseSpace_low,step=50,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)))
    MomentumDistributionPlot(sol_med,"Ele",PhaseSpace_med,step=50,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)))
    MomentumDistributionPlot(sol_high,"Ele",PhaseSpace_high,step=50,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)))

    FracNumberDensityPlot(sol_low,PhaseSpace_low)
    EnergyDensityPlot(sol_low,PhaseSpace_low)

    FracNumberDensityPlot(sol_med,PhaseSpace_med)
    EnergyDensityPlot(sol_med,PhaseSpace_med)

    FracNumberDensityPlot(sol_high,PhaseSpace_high)
    EnergyDensityPlot(sol_high,PhaseSpace_high)
    
# ====== Saving Plots for tutorial/paper ======= #

    #=PDisLowPlotDark = MomentumDistributionPlot(sol_low,"Ele",PhaseSpace_low,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusDark())
    PDisMedPlotDark = MomentumDistributionPlot(sol_med,"Ele",PhaseSpace_med,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusDark())
    PDisHighPlotDark = MomentumDistributionPlot(sol_high,"Ele",PhaseSpace_high,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusDark())

    PDisLowPlotLight = MomentumDistributionPlot(sol_low,"Ele",PhaseSpace_low,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusLight())
    PDisMedPlotLight = MomentumDistributionPlot(sol_med,"Ele",PhaseSpace_med,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusLight())
    PDisHighPlotLight = MomentumDistributionPlot(sol_high,"Ele",PhaseSpace_high,step=50,thermal=false,perp=true,order=-2,plot_limits=((-5.0,1.0),(-4.0,5.0)),theme=DiplodocusLight())

    Diplodocus.DiplodocusPlots.save("PDisLowPlotDark.svg",PDisLowPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisLowPlotDark.pdf",PDisLowPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisLowPlotLight.svg",PDisLowPlotLight)
    Diplodocus.DiplodocusPlots.save("PDisLowPlotLight.pdf",PDisLowPlotLight)

    Diplodocus.DiplodocusPlots.save("PDisMedPlotDark.svg",PDisMedPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisMedPlotDark.pdf",PDisMedPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisMedPlotLight.svg",PDisMedPlotLight)
    Diplodocus.DiplodocusPlots.save("PDisMedPlotLight.pdf",PDisMedPlotLight)

    Diplodocus.DiplodocusPlots.save("PDisHighPlotDark.svg",PDisHighPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisHighPlotDark.pdf",PDisHighPlotDark)
    Diplodocus.DiplodocusPlots.save("PDisHighPlotLight.svg",PDisHighPlotLight)
    Diplodocus.DiplodocusPlots.save("PDisHighPlotLight.pdf",PDisHighPlotLight)

    FracNumPlotLowDark = FracNumberDensityPlot(sol_low,PhaseSpace_low,theme=DiplodocusDark())
    EngPlotLowDark = EnergyDensityPlot(sol_low,PhaseSpace_low,theme=DiplodocusDark())

    FracNumPlotMedDark = FracNumberDensityPlot(sol_med,PhaseSpace_med,theme=DiplodocusDark())
    EngPlotMedDark = EnergyDensityPlot(sol_med,PhaseSpace_med,theme=DiplodocusDark())

    FracNumPlotHighDark = FracNumberDensityPlot(sol_high,PhaseSpace_high,theme=DiplodocusDark())
    EngPlotHighDark = EnergyDensityPlot(sol_high,PhaseSpace_high,theme=DiplodocusDark())

    FracNumPlotLowLight = FracNumberDensityPlot(sol_low,PhaseSpace_low,theme=DiplodocusLight())
    EngPlotLowLight = EnergyDensityPlot(sol_low,PhaseSpace_low,theme=DiplodocusLight())

    FracNumPlotMedLight = FracNumberDensityPlot(sol_med,PhaseSpace_med,theme=DiplodocusLight())
    EngPlotMedLight = EnergyDensityPlot(sol_med,PhaseSpace_med,theme=DiplodocusLight())

    FracNumPlotHighLight = FracNumberDensityPlot(sol_high,PhaseSpace_high,theme=DiplodocusLight())
    EngPlotHighLight = EnergyDensityPlot(sol_high,PhaseSpace_high,theme=DiplodocusLight())

    Diplodocus.DiplodocusPlots.save("FracNumPlotLowDark.svg",FracNumPlotLowDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotLowDark.pdf",FracNumPlotLowDark)
    Diplodocus.DiplodocusPlots.save("EngPlotLowDark.svg",EngPlotLowDark)
    Diplodocus.DiplodocusPlots.save("EngPlotLowDark.pdf",EngPlotLowDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotMedDark.svg",FracNumPlotMedDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotMedDark.pdf",FracNumPlotMedDark)
    Diplodocus.DiplodocusPlots.save("EngPlotMedDark.svg",EngPlotMedDark)
    Diplodocus.DiplodocusPlots.save("EngPlotMedDark.pdf",EngPlotMedDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotHighDark.svg",FracNumPlotHighDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotHighDark.pdf",FracNumPlotHighDark)
    Diplodocus.DiplodocusPlots.save("EngPlotHighDark.svg",EngPlotHighDark)
    Diplodocus.DiplodocusPlots.save("EngPlotHighDark.pdf",EngPlotHighDark)
    Diplodocus.DiplodocusPlots.save("FracNumPlotLowLight.svg",FracNumPlotLowLight)
    Diplodocus.DiplodocusPlots.save("FracNumPlotLowLight.pdf",FracNumPlotLowLight)
    Diplodocus.DiplodocusPlots.save("EngPlotLowLight.svg",EngPlotLowLight)
    Diplodocus.DiplodocusPlots.save("EngPlotLowLight.pdf",EngPlotLowLight)
    Diplodocus.DiplodocusPlots.save("FracNumPlotMedLight.svg",FracNumPlotMedLight)
    Diplodocus.DiplodocusPlots.save("FracNumPlotMedLight.pdf",FracNumPlotMedLight)
    Diplodocus.DiplodocusPlots.save("EngPlotMedLight.svg",EngPlotMedLight)
    Diplodocus.DiplodocusPlots.save("EngPlotMedLight.pdf",EngPlotMedLight)
    Diplodocus.DiplodocusPlots.save("FracNumPlotHighLight.svg",FracNumPlotHighLight)
    Diplodocus.DiplodocusPlots.save("FracNumPlotHighLight.pdf",FracNumPlotHighLight)
    Diplodocus.DiplodocusPlots.save("EngPlotHighLight.svg",EngPlotHighLight)
    Diplodocus.DiplodocusPlots.save("EngPlotHighLight.pdf",EngPlotHighLight)=#



