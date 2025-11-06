using Diplodocus

# ==== Define domains of time and space ====== #

    t_up::Float64 = 3.0 # seconds * (σT*c)
    t_low::Float64 = 0.0 # seconds * (σT*c)
    t_num::Int64 = 15000
    t_grid::String = "l"

    time = TimeStruct(t_up,t_low,t_num,t_grid)

    space_coords = Cartesian() # x = x, y = y, z = z

    x_up::Float64 = 0.5
    x_low::Float64 = -0.5
    x_grid::String = "u"
    x_num::Int64 = 1

    y_up::Float64 = 0.5
    y_low::Float64 = -0.5
    y_grid::String = "u"
    y_num::Int64 = 1

    z_up::Float64 = 0.5
    z_low::Float64 = -0.5
    z_grid::String = "u"
    z_num::Int64 = 1

    space = SpaceStruct(space_coords,x_up,x_low,x_grid,x_num,y_up,y_low,y_grid,y_num,z_up,z_low,z_grid,z_num)

# ==== Define particles and their momentum space domains ====== #

    name_list::Vector{String} = ["name1","name2","name3","name4"];

    momentum_coords = Spherical() # px = p, py = u, pz = phi

    px_up_list::Vector{Float64} = [4.0,4.0,4.0,4.0,];
    px_low_list::Vector{Float64} = [-5.0,-5.0,-5.0,-5.0,];
    px_grid_list::Vector{String} = ["l","l","l","l",];
    px_num_list::Vector{Int64} = [72,72,72,72,];

    py_up_list::Vector{Float64} = [1.0,1.0,1.0,1.0,];
    py_low_list::Vector{Float64} = [-1.0,-1.0,-1.0,-1.0,];
    py_grid_list::Vector{String} = ["u","u","u","u",];
    py_num_list::Vector{Int64} = [8,8,8,8,];

    pz_up_list::Vector{Float64} = [2.0pi,2.0pi,2.0pi,2.0pi];
    pz_low_list::Vector{Float64} = [0.0,0.0,0.0,0.0];
    pz_grid_list::Vector{String} = ["u","u","u","u"];
    pz_num_list::Vector{Int64} = [1,1,1,1,];

    momentum = MomentumStruct(momentum_coords,px_up_list,px_low_list,px_grid_list,px_num_list,py_up_list,py_low_list,py_grid_list,py_num_list,pz_up_list,pz_low_list,pz_grid_list,pz_num_list,"upwind");

    Binary_list::Vector{BinaryStruct} = [BinaryStruct("name1","name2","name3","name4")];
    Emi_list::Vector{EmiStruct} = [EmiStruct("name1","name1","name2","type",[Ext],type)];
    Forces::Vector{ForceType} = [ForceStruct()];

    PhaseSpace = PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,Forces);

# ==== Build Interaction and Flux Matrices ====== #

    DataDirectory = pwd()*"\\Data"
    BigM = BuildBigMatrices(PhaseSpace,DataDirectory;loading_check=false);
    FluxM = BuildFluxMatrices(PhaseSpace);

# ===== Set Initial Conditions ================== #

    Initial = Initialise_Initial_Condition(PhaseSpace);
    Initial_Constant!(Initial,PhaseSpace,"name1",pmin=10.0,pmax=13.0,umin=-0.25,umax=0.24,hmin=0.0,hmax=2.0,num_Init=1.0);

# ===== Run the Solver ================== #

    scheme = EulerStruct(Initial,PhaseSpace,BigM,FluxM,false)

    fileName = "Tutorial0.jld2";
    fileLocation = pwd()*"\\Data";

    sol = Solve(Initial,scheme;save_steps=100,progress=true,fileName=fileName,fileLocation=fileLocation);

