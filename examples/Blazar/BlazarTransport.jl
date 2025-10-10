using Diplodocus

# ==== Define domains of time and space ====== #

    B = 1e-4

    t_up::Float64 = SIToCodeUnitsTime(1e5) #log10(SIToCodeUnitsTime(1e6))
    t_low::Float64 = SIToCodeUnitsTime(0e0) # log10(SIToCodeUnitsTime(1e0))
    t_num::Int64 = 1000
    t_grid::String = "u"

    #t_up::Float64 = SIToCodeUnitsTime(1e2) # -12e0 seconds * (σT*c)
    #t_low::Float64 = SIToCodeUnitsTime(0e0) #-20e0 seconds * (σT*c)
    #t_num::Int64 = 2000
    #t_grid::String = "u"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Cylindrical(0.0,0.5*1/3,0.0)# x = r, y = phi, z = z

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

    px_up_list::Vector{Float64} = [7.0,7.0];
    px_low_list::Vector{Float64} = [-3.0,-15.0];
    px_grid_list::Vector{String} = ["l","l"];
    px_num_list::Vector{Int64} = [80,88];

    py_up_list::Vector{Float64} = [1.0,1.0];
    py_low_list::Vector{Float64} = [-1.0,-1.0];
    py_grid_list::Vector{String} = ["u","u"];
    py_num_list::Vector{Int64} = [9,9];

    pz_up_list::Vector{Float64} = [2.0pi,2.0pi];
    pz_low_list::Vector{Float64} = [0.0,0.0];
    pz_grid_list::Vector{String} = ["u","u"];
    pz_num_list::Vector{Int64} = [1,1]
    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

# ==== Define Interactions  ====== #

    Binary_list::Vector{BinaryStruct} = [BinaryStruct("Ele","Pho","Ele","Pho"),];
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[B],Ani()),];
    Forces::Vector{ForceType} = [CoordinateForce(),SyncRadReact(Ani(),B),];

    # build phase space
    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd()*"\\Data"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

    GC.gc()

# Set initial conditions

    Initial = Initialise_Initial_Condition(PhaseSpace);
    Initial_UnBoostedPowerLaw!(Initial,PhaseSpace,"Ele",pmin=1e3,pmax=1e6,Gamma=2.0,index=2.0,num_Init=1e6);
    #Initial_PowerLaw!(Initial,PhaseSpace,"Ele",pmin=1e3,pmax=1e6,umin=-1.0,umax=1.0,hmin=0.0,hmax=2.0,index=2.0,num_Init=1e6);
    #Initial_BlackBody!(Initial,PhaseSpace,"Pho",T=3.0,umin=-1.0,umax=1.0,hmin=0.0,hmax=2.0,num_Init=1e15);

    N = Diplodocus.DiplodocusTransport.FourFlow(Initial[1:80*9],80,9,1,PhaseSpace.Grids.pxr_list[1],PhaseSpace.Grids.pyr_list[1],PhaseSpace.Grids.pzr_list[1],1.0)
    U = Diplodocus.DiplodocusTransport.HydroFourVelocity(N)

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "Blazar_new_Ani_5div3.jld2";
    fileLocation = pwd()*"\\examples\\Data";

    #VSCodeServer.Profile.clear()
    #@profview sol = Solve(Initial,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);
    sol = Solve(Initial,scheme;save_steps=20,progress=true,fileName=fileName,fileLocation=fileLocation);

    (PhaseSpace,sol) = SolutionFileLoad(fileLocation,fileName);

    Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol.f[62],62;wide=true,plot_limits=((-15,8),(0,20)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),p_timescale=false)

    MomentumDistributionPlot(sol,["Pho","Ele"],PhaseSpace,Static(),step=5,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,plot_limits=((-15,7),(-10,11)))

# ===== Load and Plot Results ================== # 

    fileLocation = pwd()*"\\examples\\Data";
    fileName = "Blazar_new_Iso.jld2";
    (PS_Iso, sol_Iso) = SolutionFileLoad(fileLocation,fileName);
    fileName = "Blazar_new_Ani_0.jld2";
    (PS_Ani0, sol_Ani0) = SolutionFileLoad(fileLocation,fileName);
    fileName = "Blazar_new_Ani_5.jld2";
    (PS_Ani5, sol_Ani5) = SolutionFileLoad(fileLocation,fileName);
    fileName = "Blazar_new_Ani_5div3.jld2";
    (PS_Ani5div3, sol_Ani5div3) = SolutionFileLoad(fileLocation,fileName);
    fileName = "Blazar_new_Ani_10div3.jld2";
    (PS_Ani10div3, sol_Ani10div3) = SolutionFileLoad(fileLocation,fileName);

    SSCBFieldDark = Diplodocus.DiplodocusPlots.BFieldObserverPlot([sol_Ani0,sol_Ani5div3,sol_Ani10div3,sol_Ani5,sol_Iso],[PS_Ani0,PS_Ani5div3,PS_Ani10div3,PS_Ani5,PS_Iso],52,0.1,plot_limits=((-16.0,7.0),(-21.5,-12.5)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),R=1e14,Z=1e21,ObserverDistance=1e26)
    SSCBFieldLight = Diplodocus.DiplodocusPlots.BFieldObserverPlot([sol_Ani0,sol_Ani5div3,sol_Ani10div3,sol_Ani5,sol_Iso],[PS_Ani0,PS_Ani5div3,PS_Ani10div3,PS_Ani5,PS_Iso],52,0.1,plot_limits=((-16.0,6.0),(-21.5,-12.5)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),R=1e14,Z=1e21,ObserverDistance=1e26)

    Diplodocus.DiplodocusPlots.save("SSCBFieldDark.svg",SSCBFieldDark)
    Diplodocus.DiplodocusPlots.save("SSCBFieldDark.pdf",SSCBFieldDark)
    Diplodocus.DiplodocusPlots.save("SSCBFieldLight.svg",SSCBFieldLight)
    Diplodocus.DiplodocusPlots.save("SSCBFieldLight.pdf",SSCBFieldLight)