"""
    CoordinateTransform(pos_from,pos_to,from::AbstractCoordinates,to::AbstractCoordinates)

Converts the coordinates `pos_from` in the `from` coordinate system to the `to` coordinate system, storing the result in `pos_to`.
"""
CoordinateTransform!(pos_from,pos_to,from::AbstractCoordinates,to::AbstractCoordinates) = error("Coordinate transform function not defined for coordinates $(typeof(from)) to coordinates $(typeof(to)).")


function CoordinateTransform!(pos_from::AbstractVector{T},pos_to::AbstractVector{T},::Paraboloidal,::Cartesian) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    x::T = u*v*cos(ϕ)
    y::T = u*v*sin(ϕ)
    z::T = 1/2 * (u^2 - v^2)

    pos_to[1] = t 
    pos_to[2] = x
    pos_to[3] = y
    pos_to[4] = z

    return nothing

end

function CoordinateTransform!(pos_from::AbstractVector{T},pos_to::AbstractVector{T},::Paraboloidal,::Cylindrical) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    ρ::T = u*v
    z::T = 1/2 * (u^2 - v^2)

    pos_to[1] = t
    pos_to[2] = ρ
    pos_to[3] = ϕ
    pos_to[4] = z

    return nothing

end

function CoordinateTransform!(pos_from::AbstractVector{T},pos_to::AbstractVector{T},::Paraboloidal,::Paraboloidal) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    pos_to[1] = t 
    pos_to[2] = ϕ
    pos_to[3] = u
    pos_to[4] = v

    return nothing

end