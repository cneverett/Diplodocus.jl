#= New approach to handling coordinates.
User needs only define:
    - spacetime type (e.g. Minkowski, Schwarzschild, Kerr, etc.)
    - coordinate type (e.g. Cartesian, Cylindrical, Spherical, etc.)
    - force-free electromagnetic field configuration (e.g. monopole, paraboloidal etc.)

From these, the following functions need be defined in the code:
    - metric function g_αβ(t,x,y,z) which returns the metric tensor at a given point in spacetime
    - Christoffel symbols function Γ^α_βγ(t,x,y,z) which returns the Christoffel symbols at a given point in spacetime
    - Local tetrad function e_a^α(t,x,y,z) which returns the local tetrad at a given point in spacetime in the form e_a^α = (T^α, X^α, Y^α, Z^α)
    - Inverse tetrad function e_α^a(t,x,y,z) which returns the inverse local tetrad at a given point in spacetime in the form e_α^a = transpose(-T_α, X_α, Y_α, Z_α).

From these, the Ricci rotation coefficients can be calculated using Automatic Differentiation. Then the Flux integrals can be evaluated for each point in phase space using numerical integration.
=#

# ========== Metrics ========== #
abstract type AbstractMetric end
struct Minkowski <: AbstractMetric end
@kwdef struct Schwarzschild{T} <: AbstractMetric 
    M::T = 1.0 # mass of black hole
end
@kwdef struct Kerr{T} <: AbstractMetric
    M::T = 1.0 # mass of black hole
    a::T = 0.0 # spin parameter
end

# ========== Boundary Conditions =====#
abstract type AbstractBoundaryCondition end

struct Periodic <: AbstractBoundaryCondition end # particles leaving at one coordinate boundary re-appear on the opposite coordinate boundary
struct Open <: AbstractBoundaryCondition end # particles may leave or enter at this flux boundary. Entering may be specified, could be constrained to only photons with other particles closed 
struct Closed <: AbstractBoundaryCondition end # no particles may leave or enter this boundary (no particle reflection when they hit this boundary)
struct Reflective <: AbstractBoundaryCondition end # particles reflect off this boundary (momentum component normal to boundary reverses sign) TODO: this reversal of sign may not be well defined for specific coordinate systems and momentum grids
struct Escape <: AbstractBoundaryCondition end # particles leave based on their momentum magnitude and not direction, useful for isotropic particle distributions in spherical regions.

#===== Modes ======#
abstract type AbstractMode end
struct Iso <: AbstractMode end
struct Axi <: AbstractMode end
struct Ani <: AbstractMode end

# ========== Coordinates ========== #
abstract type AbstractCoordinates end

@kwdef struct Cartesian <: AbstractCoordinates # x=x, y=y, z=z
    xp_BC::AbstractBoundaryCondition = Periodic()
    xm_BC::AbstractBoundaryCondition = Periodic()
    yp_BC::AbstractBoundaryCondition = Periodic()
    ym_BC::AbstractBoundaryCondition = Periodic()
    zp_BC::AbstractBoundaryCondition = Periodic()
    zm_BC::AbstractBoundaryCondition = Periodic() 
end
@kwdef struct Cylindrical <: AbstractCoordinates # x=ρ, y=ϑ, z=z
    xp_BC::AbstractBoundaryCondition = Open()
    xm_BC::AbstractBoundaryCondition = Closed() # ρ=0 is axis so should be closed
    yp_BC::AbstractBoundaryCondition = Periodic()
    ym_BC::AbstractBoundaryCondition = Periodic()
    zp_BC::AbstractBoundaryCondition = Periodic()
    zm_BC::AbstractBoundaryCondition = Periodic()
end
@kwdef struct Spherical <: AbstractCoordinates # x=r, y=θ, z=ψ
    xp_BC::AbstractBoundaryCondition = Open()
    xm_BC::AbstractBoundaryCondition = Closed() # r=0 is origin so should be closed
    yp_BC::AbstractBoundaryCondition = Closed()
    ym_BC::AbstractBoundaryCondition = Closed()
    zp_BC::AbstractBoundaryCondition = Periodic()
    zm_BC::AbstractBoundaryCondition = Periodic()
end
@kwdef struct ModifiedSpherical <: AbstractCoordinates # momentum space x=p, y=u, z=ϕ
    xp_BC::AbstractBoundaryCondition = Closed()
    xm_BC::AbstractBoundaryCondition = Closed()
    yp_BC::AbstractBoundaryCondition = Closed()
    ym_BC::AbstractBoundaryCondition = Closed()
    zp_BC::AbstractBoundaryCondition = Periodic()
    zm_BC::AbstractBoundaryCondition = Periodic()
end
@kwdef struct Paraboloidal <: AbstractCoordinates # x=ϕ, y=u, z=v (where u and v have dimensions of sqrt(length))
    xp_BC::AbstractBoundaryCondition = Periodic()
    xm_BC::AbstractBoundaryCondition = Periodic()
    yp_BC::AbstractBoundaryCondition = Open() 
    ym_BC::AbstractBoundaryCondition = Closed() # u = 0 is jet axis (z<0) so should be like cylindrical ρ axis
    zp_BC::AbstractBoundaryCondition = Open()
    zm_BC::AbstractBoundaryCondition = Closed() # v = 0 is jet axis (z>0) so should be like cylindrical ρ axis
end
