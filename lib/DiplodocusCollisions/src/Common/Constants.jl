#= 
This defines all the values to be taken as constant in the evaluation of the Monte Carlo Integration of the S and T Matrices.
Including particle properties and some domain bounds
=#

# Physical constants taken from https://physics.nist.gov/cuu/Constants/index.html

const c::Float64 = Float64(299792458);     # Speed of light [m s-1]
const σT::Float64 = 6.6524587e-29;         # Thompson scattering cross section [m^2]
const ħ::Float64 = 1.0545718e-34;          # Reduced Planck constant [J s]
const μ0::Float64 = 4*π*1e-7;              # Permeability of free space [kg m s-2 A-2]
const RSph::Float64 = 4.601677272e-15;     # Sphere radius used for hard sphere collisions [m]
const q::Float64 = 1.60217662e-19;         # Elementary charge [C]

# Particle Masses
const mEle::Float64 = 9.10938356e-31;      # Mass of Electron [kg]
const mPos::Float64 = 9.10938356e-31;      # Mass of Positron [kg]
const mPho::Float64 = 0.0;                 # Mass of Photon [kg]
const mPro::Float64 = 1.6726219e-27;       # Mass of Proton [kg]
const mSph::Float64 = 1.6726219e-27;       # Mass of hard sphere [kg]

# Normalised Particle Masses (wrt electron mass)
const muEle::Float64 = 1.0;                # Reduced mass of Electron
const muPos::Float64 = 1.0;                # Reduced mass of Positron
const muPho::Float64 = 0.0;                # Reduced mass of Photon
const muPro::Float64 = 1836.1528;          # Reduced mass of Proton
const muSph::Float64 = 1836.1528;          # Reduced mass of hard sphere

# Normalised Particle Charges (wrt elementary charge)
const zEle::Float64 = -1.0;                # Charge of Electron 
const zPos::Float64 = 1.0;                 # Charge of Positron
const zPho::Float64 = 0.0;                 # Charge of Photon
const zPro::Float64 = 1.0;                 # Charge of Proton
const zSph::Float64 = 0.0;                 # Charge of hard sphere

# Domain bounds
const u_low::Float64 = -1.0;               # Lower bound for cos(theta)
const u_up::Float64 = 1.0;                 # Upper bound for cos(theta) 
const h_low::Float64 = 0.0;                # Lower bound for phi normalised by pi
const h_up::Float64 = 2.0;                 # Upper bound for phi normalised by pi

abstract type GridType end
struct UniformGrid <: GridType end
abstract type LogTenGridType end
struct LogTenGrid <: GridType end
abstract type BinaryGridType end
struct BinaryGrid <: GridType end
abstract type BoostGridType end
struct BoostGrid <: GridType end
abstract type ReBoostGridType end
struct ReBoostGrid <: GridType end