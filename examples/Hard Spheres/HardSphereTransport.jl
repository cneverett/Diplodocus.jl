using Diplodocus

# ==== Define domains of time and space ====== #

    t_up::Float64 = 3.0 # seconds * (σT*c)
    t_low::Float64 = 0.0 # seconds * (σT*c)
    t_num::Int64 = 15000
    t_grid::String = "l"

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

# ===== Set Initial Conditions ================== #

    Initial = Initialise_Initial_Condition(PhaseSpace)
    Initial_Constant!(Initial,PhaseSpace,"Sph",pmin=10.0,pmax=13.0,umin=-0.25,umax=0.24,hmin=0.0,hmax=2.0,num_Init=1.0);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)

    fileName = "HardSphereTest.jld2";
    fileLocation = pwd()*"\\examples\\Hard Spheres\\Data";

    sol = Solve(Initial,scheme;save_steps=100,progress=true,fileName=fileName,fileLocation=fileLocation);

# ===== Load and Plot Results ================== # 

    (PhaseSpace, sol) = SolutionFileLoad(fileLocation,fileName);
    
    MomentumAndPolarAngleDistributionPlot(sol,"Sph",PhaseSpace,(0.0,10.0,1000.0),order=1)

    MomentumDistributionPlot(sol,"Sph",PhaseSpace,step=1,thermal=true,order=1,plot_limits=(-0.2,1.9,-2.1,0.8))
    
    IsThermalAndIsotropicPlot(sol,PhaseSpace)
    FracNumberDensityPlot(sol,PhaseSpace)
    FracEnergyDensityPlot(sol,PhaseSpace,species="Sph")

    NumberDensityPlot(sol,PhaseSpace)
    EnergyDensityPlot(sol,PhaseSpace)

    MomentumComboAnimation(sol,["Sph"],PhaseSpace;plot_limits_momentum=(-0.2,1.9,-2.1,0.8),filename="HardSphereMomentumComboAnimation.mp4",thermal=true)


    ####

    #Diplodocus.DiplodocusTransport.AllPlots_Ani(sol,PhaseSpace,numInit_list,engInit_list,tempInit_list,"test1.mp4";fps=12,istart=1,istop=1000,iframe=nothing,step=5)




# ==== Saving plots for tutorial /paper ==== #

    #=
    
    HardSphereFracEngDenPlotDark = FracEnergyDensityPlot(sol,PhaseSpace,species="Sph",theme=DiplodocusDark())
    Diplodocus.DiplodocusPlots.save("HardSphereFracEngDenPlotDark.svg",HardSphereFracEngDenPlotDark)
    Diplodocus.DiplodocusPlots.save("HardSphereFracEngDenPlotDark.pdf",HardSphereFracEngDenPlotDark)

    HardSphereFracEngDenPlotLight = FracEnergyDensityPlot(sol,PhaseSpace,species="Sph",theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("HardSphereFracEngDenPlotLight.svg",HardSphereFracEngDenPlotLight)
    Diplodocus.DiplodocusPlots.save("HardSphereFracEngDenPlotLight.pdf",HardSphereFracEngDenPlotLight)

    HardSphereFracNumDenPlotDark = FracNumberDensityPlot(sol,PhaseSpace,theme=DiplodocusDark())
    Diplodocus.DiplodocusPlots.save("HardSphereFracNumDenPlotDark.svg",HardSphereFracNumDenPlotDark)
    Diplodocus.DiplodocusPlots.save("HardSphereFracNumDenPlotDark.pdf",HardSphereFracNumDenPlotDark)

    HardSphereFracNumDenPlotLight = FracNumberDensityPlot(sol,PhaseSpace,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("HardSphereFracNumDenPlotLight.svg",HardSphereFracNumDenPlotLight)
    Diplodocus.DiplodocusPlots.save("HardSphereFracNumDenPlotLight.pdf",HardSphereFracNumDenPlotLight)

    HardSphereIsTAndIPlotDark = IsThermalAndIsotropicPlot(sol,PhaseSpace,theme=DiplodocusDark())
    Diplodocus.DiplodocusPlots.save("HardSphereIsTAndIPlotDark.svg",HardSphereIsTAndIPlotDark)
    Diplodocus.DiplodocusPlots.save("HardSphereIsTAndIPlotDark.pdf",HardSphereIsTAndIPlotDark)

    HardSphereIsTAndIPlotLight = IsThermalAndIsotropicPlot(sol,PhaseSpace,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("HardSphereIsTAndIPlotLight.svg",HardSphereIsTAndIPlotLight)
    Diplodocus.DiplodocusPlots.save("HardSphereIsTAndIPlotLight.pdf",HardSphereIsTAndIPlotLight)

    HardSpherePAndUDisPlotDark = MomentumAndPolarAngleDistributionPlot(sol,"Sph",PhaseSpace,Static(),(0.0,10.0,1000.0),order=1,theme=DiplodocusDark())
    Diplodocus.DiplodocusPlots.save("HardSpherePAndUDisPlotDark.svg",HardSpherePAndUDisPlotDark)
    Diplodocus.DiplodocusPlots.save("HardSpherePAndUDisPlotDark.pdf",HardSpherePAndUDisPlotDark)

    HardSpherePAndUDisPlotLight = MomentumAndPolarAngleDistributionPlot(sol,"Sph",PhaseSpace,Static(),(0.0,10.0,1000.0),order=1,theme=DiplodocusLight())
    Diplodocus.DiplodocusPlots.save("HardSpherePAndUDisPlotLight.svg",HardSpherePAndUDisPlotLight)
    Diplodocus.DiplodocusPlots.save("HardSpherePAndUDisPlotLight.pdf",HardSpherePAndUDisPlotLight)

    HardSpherePDisPlotDark = MomentumDistributionPlot(sol,["Sph"],PhaseSpace,Static(),step=10,thermal=true,legend=false,order=1,theme=DiplodocusDark(),plot_limits=(-0.2,1.9,-2.1,0.8))
    Diplodocus.DiplodocusPlots.save("HardSpherePDisPlotDark.svg",HardSpherePDisPlotDark)
    Diplodocus.DiplodocusPlots.save("HardSpherePDisPlotDark.pdf",HardSpherePDisPlotDark)

    HardSpherePDisPlotLight = MomentumDistributionPlot(sol,["Sph"],PhaseSpace,Static(),step=10,thermal=true,legend=false,order=1,theme=DiplodocusLight(),plot_limits=(-0.2,1.9,-2.1,0.8))
    Diplodocus.DiplodocusPlots.save("HardSpherePDisPlotLight.svg",HardSpherePDisPlotLight)
    Diplodocus.DiplodocusPlots.save("HardSpherePDisPlotLight.pdf",HardSpherePDisPlotLight)

    =#
    