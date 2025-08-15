using Diplodocus

# ==== Define domains of time and space ====== #

    #t_up::Float64 = 1e-11 # seconds * (σT*c)
    #t_low::Float64 = 0.0 # seconds * (σT*c)
    #t_num::Int64 = 100000
    #t_grid::String = "u"

    t_up::Float64 = log10(SIToCodeUnitsTime(1e8)) # -12e0 seconds * (σT*c)
    t_low::Float64 = log10(SIToCodeUnitsTime(1e0)) #-20e0 seconds * (σT*c)
    t_num::Int64 = 400
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
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[1e-4],Ani())];
    Forces::Vector{ForceType} = [CoordinateForce(),SyncRadReact(Ani(),1e-4),];

    # build phase space
    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# location of DataDirectory where Interaction Matrices are stored

    DataDirectory = pwd()*"\\examples\\Synchrotron\\Data"

# Load interaction matrices

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

# Set initial conditions

    Initial_Ele = Initial_PowerLaw(PhaseSpace,"Ele",1e3,1e6,0.0,1.0,0.0,2.0,2.0,1e6);
    Initial_Pho = Initial_Constant(PhaseSpace,"Pho",1,80,1,15,1,1,0.0);
    Initial = ArrayPartition(Initial_Ele,Initial_Pho);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "SyncAni_Test.jld2";
    fileLocation = pwd()*"\\examples\\Synchrotron\\Data";

    sol_Iso = Solve(Initial,scheme;save_steps=5,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol_Iso) = SolutionFileLoad(fileLocation,fileName);

    MomentumAndPolarAngleDistributionPlot(sol_Iso,"Ele",PhaseSpace,(1,400,800),order=1)
    MomentumAndPolarAngleDistributionPlot(sol_Iso,"Pho",PhaseSpace,(1,400,800),order=1)

    MomentumDistributionPlot(sol_Iso,["Pho","Ele"],PhaseSpace,Diplodocus.DiplodocusPlots.Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime)

    NumberDensityPlot(sol_Iso,PhaseSpace,species="Pho",theme=DiplodocusDark(),title=nothing)
    NumberDensityPlot(sol_Iso,PhaseSpace,species="Ele",theme=DiplodocusDark(),title=nothing)
    FracNumberDensityPlot(sol_Iso,PhaseSpace,species="Ele",theme=DiplodocusDark(),title=nothing)

    EnergyDensityPlot(sol_Iso,PhaseSpace,theme=DiplodocusDark(),title=nothing)
    FracEnergyDensityPlot(sol_Iso,PhaseSpace,theme=DiplodocusDark(),title=nothing,only_all=true)

    MomentumDistributionPlot(sol_Iso,["Ele","Pho"],PhaseSpace,Animated(),thermal=false,order=2,plot_limits=((-15,7),(2,10)),wide=true,TimeUnits=CodeToSIUnitsTime,filename="SyncElePhoAnimated.mp4")

# ====== Anisotropic Observations ====== # 

    Binary_list::Vector{BinaryStruct} = [];
    Emi_list::Vector{EmiStruct} = [EmiStruct("Ele","Ele","Pho","Sync",[1e-4],Ani())];
    Forces::Vector{ForceType} = [SyncRadReact(Ani(),1e-4),];

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

    DataDirectory = pwd()*"\\examples\\Synchrotron\\Data";   

    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=true);
    FluxM = BuildFluxMatrices(PhaseSpace);

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)
    fileName = "SyncAni_test.jld2";
    fileLocation = pwd()*"\\examples\\Synchrotron\\Data";

    sol_Ani = Solve(Initial,scheme;save_steps=5,progress=true,fileName=fileName,fileLocation=fileLocation);

# ====== Plot Observer angle dependence ======= #

    (PhaseSpace, sol_Ani) = SolutionFileLoad(fileLocation,fileName);

    CodeToSIUnitsTime(sol_Ani.t[52])
    ObserverFluxPlot(PhaseSpace,sol_Ani,52,[0.1,0.2,0.3,0.4,0.5],1.0,TimeUnits=CodeToSIUnitsTime,plot_limits=(-6.5,3.5,-0.5,5.5))

# ==== Saving plots for tutorial /paper ==== #

    
    SyncPDisPlotDark = MomentumDistributionPlot(sol_Iso,["Pho","Ele"],PhaseSpace,Diplodocus.DiplodocusPlots.Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusDark(),plot_limits=((-15,7),(2,10)))
    SyncPDisPlotLight = MomentumDistributionPlot(sol_Iso,["Pho","Ele"],PhaseSpace,Diplodocus.DiplodocusPlots.Static(),step=10,order=2,wide=true,TimeUnits=CodeToSIUnitsTime,theme=DiplodocusLight(),plot_limits=((-15,7),(2,10)))
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.pdf",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotDark.svg",SyncPDisPlotDark)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.pdf",SyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.svg",SyncPDisPlotLight)
    Diplodocus.DiplodocusPlots.save("SyncPDisPlotLight.png",SyncPDisPlotLight)

    # ==== AM3 comparison plots ==== # 

    AM3ElePDisPlotDark = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",1e8,1e0,"l",wide=true,plot_limits=((-15,7),(2,10)),theme=DiplodocusDark())
    AM3ElePDisPlotLight = Diplodocus.DiplodocusPlots.AM3_MomentumDistributionPlot("./AM3/syn_test_0-8.jld2",1e8,1e0,"l",wide=true,plot_limits=((-15,7),(2,10)),theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotDark.pdf",AM3ElePDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotDark.svg",AM3ElePDisPlotDark)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotLight.pdf",AM3ElePDisPlotLight)
    Diplodocus.DiplodocusPlots.save("AM3ElePDisPlotLight.svg",AM3ElePDisPlotLight)
    

    # ==== Observer flux plots ==== #
    
    ObsPlotDark = ObserverFluxPlot(PhaseSpace,sol_Ani,52,[0.1,0.2,0.3,0.4,0.5],1.0,TimeUnits=CodeToSIUnitsTime,plot_limits=(-6.5,3.5,-0.5,5.5),theme=DiplodocusDark())
    ObsPlotLight = ObserverFluxPlot(PhaseSpace,sol_Ani,52,[0.1,0.2,0.3,0.4,0.5],1.0,TimeUnits=CodeToSIUnitsTime,plot_limits=(-6.5,3.5,-0.5,5.5),theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("ObsPlotDark.pdf",ObsPlotDark)
    Diplodocus.DiplodocusPlots.save("ObsPlotDark.svg",ObsPlotDark)
    Diplodocus.DiplodocusPlots.save("ObsPlotLight.pdf",ObsPlotLight)
    Diplodocus.DiplodocusPlots.save("ObsPlotLight.svg",ObsPlotLight)
    Diplodocus.DiplodocusPlots.save("ObsPlotDark.png",ObsPlotDark)
    Diplodocus.DiplodocusPlots.save("ObsPlotLight.png",ObsPlotLight)
    

    # ===== Sync Timescale plot ==== #

    PhaseSpace.Grids.pzr_list[2] 
    PhaseSpace.Grids.dt[2]
    SIToCodeUnitsTime(1e0)

    I = Diplodocus.DiplodocusTransport.diag(FluxM.I_Flux)[1:px_num_list[1]*py_num_list[1]*pz_num_list[1]]
    I = reshape(I,(px_num_list[1],py_num_list[1],pz_num_list[1]))
    #I = dropdims(sum(I,dims=(2,3)),dims=(2,3)) # I flux 
    I = I[:,8,1] # I flux
    A = Diplodocus.DiplodocusTransport.diag(FluxM.Ap_Flux)[1:px_num_list[1]*py_num_list[1]*pz_num_list[1]]
    A = reshape(A,(px_num_list[1],py_num_list[1],pz_num_list[1]))
    A = A[:,1,1] # A flux

    test = I ./ A ./ PhaseSpace.Grids.dt[1]
    inv_test = ones(size(test)) ./ test

    SI_inv_test = CodeToSIUnitsTime.(inv_test)

    analytic = ones(length(SI_inv_test))
    @. analytic /= ((1e-4)^2 / (4pi * 1e-7) * 2/3 * PhaseSpace.Grids.dE_list[1] / (9.11e-31 * 2.99792458e8^2))
    analytic = CodeToSIUnitsTime.(analytic)

    PhaseSpace.Grids.pxr_list[1]

    x = zeros(Float64,length(PhaseSpace.Grids.mpx_list[1])) 
    for i = 1:80
        x[i] = log10(PhaseSpace.Grids.pxr_list[1][i+1]/PhaseSpace.Grids.pxr_list[1][i])
    end
    
    PhaseSpace.Grids.pxr_list[1]
    y = PhaseSpace.Grids.pxr_list[1][2]/PhaseSpace.Grids.pxr_list[1][1]
    hell = PhaseSpace.Grids.mpx_list[1] 
    dp = PhaseSpace.Grids.dpx_list[1]

    fig = Diplodocus.DiplodocusPlots.Figure()
    ax = Diplodocus.DiplodocusPlots.Axis(fig[1,1], xlabel="p", ylabel="(4pi)^2 times Sync timescale log10(s)", title="Sync timescale as a function of p")
    Diplodocus.DiplodocusPlots.scatter!(ax, log10.(PhaseSpace.Grids.mpx_list[1]), log10.(SI_inv_test))
    Diplodocus.DiplodocusPlots.scatter!(ax, log10.(PhaseSpace.Grids.mpx_list[1]), log10.(analytic))
    display(fig)

    
    
    PhaseSpace.Grids.pxr_list[1][40]/PhaseSpace.Grids.pxr_list[1][39]
    