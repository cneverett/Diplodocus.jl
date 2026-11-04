"""
    ChristoffelComponents!(::AbstractMetric,::AbstractCoordinates,Γ::MArray{4,4,4},pos)

Returns the Christoffel symbol components ``Γ=Γ^α_{~βγ}`` at a point `pos` in spacetime defined by the coordinates `pos=(t,x,y,z)`. The Christoffel symbol components are stored in the static array Γ, which is modified in-place. The specific form of the Christoffel symbols depends on the type of metric and coordinate basis being used.
"""
ChristoffelComponents!(pos,Γ,metric::AbstractMetric,coordinates::AbstractCoordinates) = error("Christoffel function not defined for this metric $(typeof(metric)) and coordinate $(typeof(coordinates)).")

#===== Minkowski =======#
#=======================#
@inline function ChristoffelComponents!(pos::MVector{4,T},Γ::MArray{Tuple{4,4,4},Float64,3,64},::Minkowski,::Cartesian) where T

    return nothing

end

@inline function ChristoffelComponents!(pos::MVector{4,T},Γ::MArray{Tuple{4,4,4},Float64,3,64},::Minkowski,::Cylindrical) where T

    t = pos[1]
    ρ = pos[2]
    ϕ = pos[3]
    z = pos[4]

    Γ[2,3,3]::T = -ρ
    Γ[3,2,3]::T = Γ[3,3,2]::T = ρ == T(0) ? T(1) : T(1)/ρ

    return nothing

end

@inline function ChristoffelComponents!(pos::MVector{4,T},Γ::MArray{Tuple{4,4,4},Float64,3,64},::Minkowski,::Spherical) where T

    t = pos[1]
    r = pos[2]
    θ = pos[3]
    ϕ = pos[4]

    sθ,cθ = sincos(θ)

    Γ[2,3,3]::T = -r
    Γ[3,2,3]::T = Γ[3,3,2]::T = 1/r
    Γ[2,4,4]::T = -r*sθ^2
    Γ[4,2,4]::T = Γ[4,4,2]::T = 1/r
    Γ[3,4,4]::T = -sθ*cθ
    Γ[4,3,4]::T = Γ[4,4,3]::T = cθ/sθ

    return nothing

end

@inline function ChristoffelComponents!(pos::MVector{4,T},Γ::MArray{Tuple{4,4,4},Float64,3,64},::Minkowski,::Paraboloidal) where T

    t = pos[1]
    ϕ = pos[2]
    u = pos[3]
    v = pos[4]

    cond1 = (u^2+v^2) == zero(T)

    Γ[2,2,3]::T = u == zero(T) ? T(1) : T(1)/u
    Γ[2,3,2]::T = Γ[2,2,3]
    Γ[2,2,4]::T = v == zero(T) ? T(1) : T(1)/v
    Γ[2,4,2]::T = Γ[2,2,4]
    Γ[3,2,2]::T = cond1 ? T(-1) : -u*v^2/(u^2+v^2)
    Γ[3,3,3]::T = cond1 ? T(1) : u/(u^2+v^2)
    Γ[3,3,4]::T = cond1 ? T(1) : v/(u^2+v^2)
    Γ[3,4,3]::T = Γ[3,3,4]
    Γ[3,4,4]::T = cond1 ? T(-1) : -u/(u^2+v^2)
    Γ[4,2,2]::T = cond1 ? T(-1) : -u^2*v/(u^2+v^2)
    Γ[4,3,3]::T = cond1 ? T(-1) : -v/(u^2+v^2)
    Γ[4,3,4]::T = cond1 ? T(1) : u/(u^2+v^2)
    Γ[4,4,3]::T = Γ[4,3,4]
    Γ[4,4,4]::T = cond1 ? T(1) : v/(u^2+v^2)

    return nothing

end