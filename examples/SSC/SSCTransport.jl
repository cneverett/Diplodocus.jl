using Diplodocus

# ==== Define domains of time and space ====== #

    B = 1e-4

    t_up::Float64 = SIToCodeUnitsTime(1e6)#log10(SIToCodeUnitsTime(1e6)) # -12e0 seconds * (σT*c)
    t_low::Float64 = SIToCodeUnitsTime(0e0) #log10(SIToCodeUnitsTime(1e0)) #-20e0 seconds * (σT*c)
    t_num::Int64 = 10000#12000#30000
    t_grid::String = "u"#"l"

    #t_up::Float64 = SIToCodeUnitsTime(1e2) # -12e0 seconds * (σT*c)
    #t_low::Float64 = SIToCodeUnitsTime(0e0) #-20e0 seconds * (σT*c)
    #t_num::Int64 = 2000
    #t_grid::String = "u"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Spherical()# x = r, y = phi, z = z

    x_up::Float64 = 1.0
    x_low::Float64 = 0f0
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

    space = SpaceStruct(space_coords,x_up,x_low,x_grid,x_num,y_up,y_low,y_grid,y_num,z_up,z_low,z_grid,z_num)

# ==== Define particles and their momentum space domains ====== #

    name_list::Vector{String} = ["Ele","Pos","Pho"];

    momentum_coords = Spherical() # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [7.0,7.0,7.0];
    px_low_list::Vector{Float64} = [-3.0,-3.0,-15.0];
    px_grid_list::Vector{String} = ["l","l","l"];
    px_num_list::Vector{Int64} = [80,80,88];

    py_up_list::Vector{Float64} = [1.0,1.0,1.0];
    py_low_list::Vector{Float64} = [-1.0,-1.0,-1.0];
    py_grid_list::Vector{String} = ["u","u","u"];
    py_num_list::Vector{Int64} = [9,9,9];

    pz_up_list::Vector{Float64} = [2.0pi,2.0pi,2.0pi];
    pz_low_list::Vector{Float64} = [0.0,0.0,0.0];
    pz_grid_list::Vector{String} = ["u","u","u"];
    pz_num_list::Vector{Int64} = [1,1,1];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

# ==== Define Interactions  ====== #

    Binary_list::Vector{BinaryStruct} = [BinaryStruct("Ele","Pho","Ele","Pho"),BinaryStruct("Pos","Pho","Pos","Pho")#= ,BinaryStruct("Pho","Pho","Ele","Pos") =#];
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[B],Iso()),EmiStruct("Pos","Pos","Pho","Sync",[B],Iso())];
    Forces::Vector{ForceType} = [SyncRadReact(Iso(),B),];

    # build phase space
    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd()*"\\examples\\Data"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true,Bin_corrected=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

    GC.gc()

# Set initial conditions

    Initial = Initialise_Initial_Condition(PhaseSpace);
    Initial_PowerLaw!(Initial,PhaseSpace,"Ele",pmin=1e3,pmax=1e6,umin=-1.0,umax=1.0,hmin=0.0,hmax=2.0,index=2.0,num_Init=10e6);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "SSC_new_6_UniformT_no_pair_prod.jld2";
    fileLocation = pwd()*"\\examples\\Data";

    sol = Solve(Initial,scheme;save_steps=10,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);

    MomentumAndPolarAngleDistributionPlot(sol,"Ele",PhaseSpace,Static(),(1,32,10000),order=1,TimeUnits=CodeToSIUnitsTime)
    MomentumAndPolarAngleDistributionPlot(sol,"Pos",PhaseSpace,Static(),(1,32,52),order=1,TimeUnits=CodeToSIUnitsTime)
    MomentumAndPolarAngleDistributionPlot(sol,"Pho",PhaseSpace,Static(),(1,32,52),order=1,TimeUnits=CodeToSIUnitsTime)

    MomentumDistributionPlot(sol,["Pho","Ele","Pos"],PhaseSpace,Static(),step=1,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,plot_limits=((-14,7),(3.5,10.5)))
    MomentumDistributionPlot(sol,["Ele"],PhaseSpace,Static(),step=10,order=2,wide=false,TimeUnits=CodeToSIUnitsTime)
    MomentumDistributionPlot(sol,["Pos"],PhaseSpace,Static(),step=1,order=2,wide=false,TimeUnits=CodeToSIUnitsTime)
    MomentumDistributionPlot(sol,["Pho"],PhaseSpace,Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime)

    NumberDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),title=nothing)
    NumberDensityPlot(sol,PhaseSpace,species="Pho",theme=DiplodocusDark(),title=nothing)
    NumberDensityPlot(sol,PhaseSpace,species="Ele",theme=DiplodocusDark(),title=nothing)
    FracNumberDensityPlot(sol,PhaseSpace,species="Ele",theme=DiplodocusDark(),title=nothing)
    FracNumberDensityPlot(sol,PhaseSpace,species="Pho",theme=DiplodocusDark(),title=nothing)
    FracNumberDensityPlot(sol,PhaseSpace,species="All",theme=DiplodocusDark(),title=nothing)

    EnergyDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),title=nothing,logt=true,TimeUnits=CodeToSIUnitsTime)
    FracEnergyDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),title=nothing)
    FracEnergyDensityPlot(sol,PhaseSpace,theme=DiplodocusDark(),title=nothing,only_all=true)

    MomentumDistributionPlot(sol_Iso,["Ele","Pho"],PhaseSpace,Animated(),thermal=false,order=2,plot_limits=((-15,7),(2,10)),wide=true,TimeUnits=CodeToSIUnitsTime,filename="SyncElePhoAnimated.mp4")

    DataDirectory = pwd()*"\\examples\\Data"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true,Bin_corrected=true);
    FluxM = BuildFluxMatrices(PhaseSpace);
    scheme = EulerStruct(sol.f[1],PhaseSpace,BigM,FluxM,false)

    CodeToSIUnitsTime.([PhaseSpace.Grids.dt[t*save_steps-2]])
    t = 10000
    save_steps = 10
    sol.t[t]
    PhaseSpace.Grids.tr[(t-2)*save_steps]
    Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol,[1,2,3,12,85,1002,10002],["Pho","Ele","Pos"];wide=true,plot_limits=((-15,8),(0.0,18.0)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),u_avg=true,p_timescale=true,plot_dt=true,logt=true)

    Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,Initial,1;wide=true,paraperp=true,plot_limits=((-15,8),(0,20)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark())
    #Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,1e12 .* ones(Float64,size(sol.f[52])),52;wide=true,plot_limits=((-15,8),(0,20)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark())

    dstate = zeros(eltype(Initial),size(Initial))
    timescale = zeros(eltype(Initial),size(Initial))

    t_idx = 450
    dt0 = PhaseSpace.Grids.tr[2] - PhaseSpace.Grids.tr[1]
    dt = PhaseSpace.Grids.tr[t_idx+1] - PhaseSpace.Grids.tr[t_idx]
    t = PhaseSpace.Grids.tr[t_idx]

    scheme(dstate,sol.f[Int(t_idx/10)],dt0,dt,t)

    test = scheme.df ./ sol.f[Int(t_idx/10)]
    
    replace!(test,Inf32=>0.0,-Inf32=>0.0,NaN=>0.0)
    maximum(abs.(test))

# ====== Plot Observer angle dependence ======= #

    (PhaseSpace, sol_Ani) = SolutionFileLoad(fileLocation,fileName);

    ObserverFluxPlot(PhaseSpace,sol,118,[0.01,0.1,0.25,0.5],1.0,TimeUnits=CodeToSIUnitsTime,plot_limits=(-15.5,7.5,-0.5,12.5))
    ObserverFluxPlot(PhaseSpace,sol,102,[0.05,0.1,0.15,0.2,0.25],1.0,TimeUnits=CodeToSIUnitsTime,plot_limits=(-15.5,7.5,-0.5,8.5))

    Diplodocus.DiplodocusPlots.AngleDistributionPlot(sol,["Pho"],PhaseSpace,Static(),118;theme=DiplodocusDark(),order=2,TimeUnits=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=true,legend=true,angle_step=1)

    Diplodocus.DiplodocusTransport.StressEnergyTensor(copy(Location_Species_To_StateVector(sol.f[52],PhaseSpace,species_index=3)),88,9,1,PhaseSpace.Grids.pxr_list[3],PhaseSpace.Grids.pyr_list[3],PhaseSpace.Grids.pzr_list[3],0.0)
    Diplodocus.DiplodocusTransport.StressEnergyTensor(copy(Location_Species_To_StateVector(sol.f[52],PhaseSpace,species_index=1)),80,9,1,PhaseSpace.Grids.pxr_list[1],PhaseSpace.Grids.pyr_list[1],PhaseSpace.Grids.pzr_list[1],0.0)

# ==== Saving plots for tutorial /paper ==== #

    
    SSCPDisPlotDark = MomentumDistributionPlot(sol,["Pho","Ele","Pos"],PhaseSpace,Diplodocus.DiplodocusPlots.Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),plot_limits=((-15,7),(3,11)))
    SSCPDisPlotLight = MomentumDistributionPlot(sol,["Pho","Ele","Pos"],PhaseSpace,Diplodocus.DiplodocusPlots.Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),plot_limits=((-15,7),(3,11)))
    Diplodocus.DiplodocusPlots.save("SSCPDisPlotDark.pdf",SSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SSCPDisPlotDark.svg",SSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SSCPDisPlotLight.pdf",SSCPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("SSCPDisPlotLight.svg",SSCPDisPlotLight)

    SSCEngDenPlotDark = EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),logt=true)
    SSCEngDenPlotLight = EnergyDensityPlot(sol,PhaseSpace,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),logt=true)
    Diplodocus.DiplodocusPlots.save("SSCEngDenPlotDark.pdf",SSCEngDenPlotDark)
    Diplodocus.DiplodocusPlots.save("SSCEngDenPlotDark.svg",SSCEngDenPlotDark)
    Diplodocus.DiplodocusPlots.save("SSCEngDenPlotLight.pdf",SSCEngDenPlotLight)
    Diplodocus.DiplodocusPlots.save("SSCEngDenPlotLight.svg",SSCEngDenPlotLight)

    # ==== AM3 comparison plots ==== # 

    AM3SSCPDisPlotDark = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/ssc_test.jld2",1e6,1e0,"l",wide=true,plot_limits=((-15,7),(3,11)),theme=DiplodocusDark())
    AM3SSCPDisPlotLight = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/ssc_test.jld2",1e6,1e0,"l",wide=true,plot_limits=((-15,7),(3,11)),theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("AM3SSCPDisPlotDark.pdf",AM3SSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3SSCPDisPlotDark.svg",AM3SSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3SSCPDisPlotLight.pdf",AM3SSCPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3SSCPDisPlotLight.svg",AM3SSCPDisPlotLight)

    AM3DIPSSCPDisPlotDark = Diplodocus.DiplodocusPlots.AM3_DIP_Combo_MomentumDistributionPlot("./AM3/ssc_test.jld2",sol,PhaseSpace,1e7,1e1,"l",plot_limits=((-15.0,7.0),(3.5,10.5)),theme=DiplodocusDark(),ele_err=false,pos_err=false,yticks=1:2:11)
    AM3DIPSSCPDisPlotLight = Diplodocus.DiplodocusPlots.AM3_DIP_Combo_MomentumDistributionPlot("./AM3/ssc_test.jld2",sol,PhaseSpace,1e7,1e1,"l",plot_limits=((-15.0,7.0),(3.5,10.5)),theme=DiplodocusLight(),ele_err=false,pos_err=false)
    Diplodocus.DiplodocusPlots.save("AM3DIPSSCPDisPlotDark.pdf",AM3DIPSSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3DIPSSCPDisPlotDark.svg",AM3DIPSSCPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3DIPSSCPDisPlotLight.pdf",AM3DIPSSCPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3DIPSSCPDisPlotLight.svg",AM3DIPSSCPDisPlotLight)


    DataDirectory = pwd()*"\\examples\\Data"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true,Bin_corrected=true);
    FluxM = BuildFluxMatrices(PhaseSpace);
    scheme = EulerStruct(sol.f[1],PhaseSpace,BigM,FluxM,false)

    SSCTimeScalePlotDark = Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol,[1,2,3,12,85,1002,10002],["Pho","Ele"];wide=false,plot_limits=((-14,8),(0.0,18.0)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),u_avg=true,p_timescale=true,plot_dt=true,logt=true,legend=true)
    SSCTimeScalePlotLight = Diplodocus.DiplodocusPlots.TimeScalePlot(scheme,sol,[1,2,3,12,85,1002,10002],["Pho","Ele"];wide=false,plot_limits=((-14,8),(0.0,18.0)),TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),u_avg=true,p_timescale=true,plot_dt=true,logt=true,legend=true)

    Diplodocus.DiplodocusPlots.save("SSCTimeScalePlotDark.pdf",SSCTimeScalePlotDark)
    Diplodocus.DiplodocusPlots.save("SSCTimeScalePlotDark.svg",SSCTimeScalePlotDark)
    Diplodocus.DiplodocusPlots.save("SSCTimeScalePlotLight.pdf",SSCTimeScalePlotLight)
    Diplodocus.DiplodocusPlots.save("SSCTimeScalePlotLight.svg",SSCTimeScalePlotLight)

    
    