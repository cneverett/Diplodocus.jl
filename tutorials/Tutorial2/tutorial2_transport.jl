using Diplodocus

Precision::DataType = Float32
Backend::BackendType = CPUBackend()

CHAR_length::Float64 = 1e14 # m = 1e16 cm

space_coords::CoordinateType = Spherical(xp_BC=Escape())

x_up::Float64 = 1.0
x_low::Float64 = 0.0
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

name_list::Vector{String} = ["Ele","Pho"];

momentum_coords = Spherical() # px = p, py = u, pz = phi

px_up_list::Vector{Float64} = [9.0,9.0]
px_low_list::Vector{Float64} = [-3.0,-16.0]
px_grid_list::Vector{String} = ["l","l"]
px_num_list::Vector{Int64} = [288,600]

py_up_list::Vector{Float64} = [1.0,1.0]
py_low_list::Vector{Float64} = [-1.0,-1.0]
py_grid_list::Vector{String} = ["u","u"]
py_num_list::Vector{Int64} = [1,1]

pz_up_list::Vector{Float64} = [2.0*pi,2.0*pi]
pz_low_list::Vector{Float64} = [0.0,0.0]
pz_grid_list::Vector{String} = ["u","u"]
pz_num_list::Vector{Int64} = [1,1]

ElectroMagneticField::ElectroMagneticFieldStruct = ElectroMagneticField_Constant(parameters=[1e-3,1e0])

PhaseSpace = PhaseSpaceStruct()

Binary_list::Vector{BinaryStruct} = [BinaryStruct("Ele","Pho","Ele","Pho")]
Binary_Domain = CollisionDomain(PhaseSpace)

B_sampled::Vector{Float64} = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1];
Sync_Domain = CollisionDomain(PhaseSpace);
Emission_list::Vector{EmiStruct} = [EmiStruct(name1="Ele",name2="Ele",name3="Pho",EmiName="Sync",Ext_sampled=B_sampled,mode=Iso(),Force=true,Domain=Sync_Domain)]

Forces_list::Vector{ForceType} = []

DataDirectory = joinpath(pwd(),"Data")

BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=true); # jl
BinM = BuildBinaryMatrices(PhaseSpace,Binary_list,Binary_Domain,DataDirectory;Bin_corrected=true,Bin_sparse=true) # md

EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true); # jl
EmiM = BuildEmissionMatrices(PhaseSpace,Emission_list,DataDirectory;Emi_corrected=true,Emi_sparse=true) # md

FluxM = BuildFluxMatrices(PhaseSpace,Forces_list); # jl
FluxM = BuildFluxMatrices(PhaseSpace,Forces_list) # md

Initial = Initialise_Initial_Condition(PhaseSpace)
Injection = Initialise_Injection_Condition(PhaseSpace)

Injection_PowerLawExpDecay!(Injection,PhaseSpace,"Ele",pmin=10^1.9,pmax=10^2.1,index=2.0,num_Inj=1e3,rate_Inj=1.0)

dt_initial::Float64 = 0.1
t_save::Vector{Precision} = range(0.0,1e1,length=101)

scheme = ForwardEulerStruct(PhaseSpace,Initial,Injection,BinM,EmiM,FluxM)

fileName = "tutorial2_output.jld2"
fileLocation = joinpath(pwd(),"Data")

sol = Solve(scheme,dt_initial,t_save;fileName=fileName,fileLocation=fileLocation,progress=true,Verbose=2);
