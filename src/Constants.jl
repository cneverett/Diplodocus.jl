#= 
This defines all the values to be taken as constant in the evaluation of the Monte Carlo Integration of the S and T Matrices.
Including particle properties and some domain bounds
=#

# Physical constants taken from https://physics.nist.gov/cuu/Constants/index.html

const CONST_c::Float64 = Float64(299792458);     # Speed of light [m s-1]
const CONST_σT::Float64 = 6.6524587e-29;         # Thompson scattering cross section [m^2]
const CONST_ħ::Float64 = 1.0545718e-34;          # Reduced Planck constant [J s]
const CONST_μ0::Float64 = 4*π*1e-7;              # Permeability of free space [kg m s-2 A-2]
const CONST_ε0::Float64 = 8.854187817e-12;       # Permittivity of free space [A^2 s^4 kg-1 m-3]
const CONST_RSph::Float64 = 4.601677272e-15;     # Sphere radius used for hard sphere collisions [m]
const CONST_q::Float64 = 1.60217662e-19;         # Elementary charge [C]
const CONST_kb::Float64 = 1.38064852e-23;        # Boltzmann constant [J K-1]

# Particle Masses
const CONST_mEle::Float64 = 9.10938356e-31;      # Mass of Electron [kg]
const CONST_mPos::Float64 = 9.10938356e-31;      # Mass of Positron [kg]
const CONST_mPho::Float64 = 0.0;                 # Mass of Photon [kg]
const CONST_mPro::Float64 = 1.6726219e-27;       # Mass of Proton [kg]
const CONST_mSph::Float64 = 1.6726219e-27;       # Mass of hard sphere [kg]

# Normalised Particle Masses (wrt electron mass)
const CONST_muEle::Float64 = 1.0;                # Reduced mass of Electron
const CONST_muPos::Float64 = 1.0;                # Reduced mass of Positron
const CONST_muPho::Float64 = 0.0;                # Reduced mass of Photon
const CONST_muPro::Float64 = 1836.1528;          # Reduced mass of Proton
const CONST_muSph::Float64 = 1836.1528;          # Reduced mass of hard sphere

# Normalised Particle Charges (wrt elementary charge)
const CONST_zEle::Float64 = -1.0;                # Charge of Electron 
const CONST_zPos::Float64 = 1.0;                 # Charge of Positron
const CONST_zPho::Float64 = 0.0;                 # Charge of Photon
const CONST_zPro::Float64 = 1.0;                 # Charge of Proton
const CONST_zSph::Float64 = 0.0;                 # Charge of hard sphere

# Common Length scales
const CONST_pcs::Float64 = 3.0857e16;           # Parsec in meters [m]
const CONST_rs_Sun::Float64 = 2.95325e3;        # Schwarzschild radius of the Sun [m]

# Domain bounds
const CONST_u0::Float64 = -1.0;                 # Lower bound for cos(theta)
const CONST_u1::Float64 = 1.0;                  # Upper bound for cos(theta) 
const CONST_h0::Float64 = 0.0;                  # Lower bound for phi normalised by pi
const CONST_h1::Float64 = 2.0;                  # Upper bound for phi normalised by pi