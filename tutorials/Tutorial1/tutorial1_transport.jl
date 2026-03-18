using Diplodocus

Precision::DataType = Float32
Backend::BackendType = CPUBackend()

CHAR_number_density::Float64 = 1.0 # m^-3
CHAR_time = 1.0 / DiplodocusTransport.CONST_σT / DiplodocusTransport.CONST_c / CHAR_number_density # s

space_coords::CoordinateType = Cartesian()

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

name_list::Vector{String} = ["Sph",];

momentum_coords = Spherical() # px = p, py = u, pz = phi

px_up_list::Vector{Float64} = [4.0,]
px_low_list::Vector{Float64} = [-5.0,]
px_grid_list::Vector{String} = ["l",]
px_num_list::Vector{Int64} = [72,]

py_up_list::Vector{Float64} = [1.0,]
py_low_list::Vector{Float64} = [-1.0,]
py_grid_list::Vector{String} = ["u",]
py_num_list::Vector{Int64} = [9,]

pz_up_list::Vector{Float64} = [2.0*pi,]
pz_low_list::Vector{Float64} = [0.0,]
pz_grid_list::Vector{String} = ["u",]
pz_num_list::Vector{Int64} = [1,]

PhaseSpace = PhaseSpaceStruct()

Binary_list::Vector{BinaryStruct} = [BinaryStruct("Sph","Sph","Sph","Sph")]
Binary_Domain = CollisionDomain(PhaseSpace)

Emission_list::Vector{EmiStruct} = []
Forces_list::Vector{ForceType} = []

DataDirectory = joinpath(pwd(),"Data")

BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=false); # jl
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=false) # md

EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true); # jl
EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true) # md

FluxM = BuildFluxMatrices(PhaseSpace,Forces_list); # jl
FluxM = BuildFluxMatrices(PhaseSpace,Forces_list) # md

Initial = Initialise_Initial_Condition(PhaseSpace)
Injection = Initialise_Injection_Condition(PhaseSpace)

Initial_Constant!(Initial,PhaseSpace,"Sph",pmin=10.0,pmax=13.0,umin=-0.11,umax=0.11,hmin=0.0,hmax=2.0,num_Init=1.0)

dt_initial::Float64 = 1.0
t_save::Vector{Precision} = range(0.0,1e3,length=1001)

scheme = ForwardEulerStruct(PhaseSpace,Initial,Injection,BinM,EmiM,FluxM)

fileName = "tutorial1_output.jld2"
fileLocation = joinpath(pwd(),"Data")

sol = Solve(scheme,dt_initial,t_save;fileName=fileName,fileLocation=fileLocation,progress=true,Verbose=2);
