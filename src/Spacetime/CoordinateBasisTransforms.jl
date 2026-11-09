"""
    CoordinateBasisTransform(vec_from,vec_to,from::AbstractCoordinates,to::AbstractCoordinates)

Converts the 4-vector `vec_from` in the `from` coordinate basis at a position `pos_from` to the `to` coordinate basis, storing the result in `vec_to`.
"""
CoordinateBasisTransform!(vec_from,vec_to,pos_from,from::AbstractCoordinates,to::AbstractCoordinates) = error("Coordinate basis transform function not defined for coordinates $(typeof(from)) to coordinates $(typeof(to)).")

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::Cartesian,::Paraboloidal) where T

    t = pos_from[1]
    x = pos_from[2]
    y = pos_from[3]
    z = pos_from[4]

    r = sqrt(x^2+y^2+z^2)
    ρ = sqrt(x^2+y^2)

    #println("x: $x, y: $y, z: $z, r: $r, ρ: $ρ")

    vt = vec_from[1]
    vx = vec_from[2]
    vy = vec_from[3]
    vz = vec_from[4]

    vt = vt
    vϕ = -y/ρ^2 * vx + x/ρ^2 * vy
    vu = x*sqrt(r-z)/(2r*ρ) * vx + y*sqrt(r-z)/(2r*ρ) * vy + sqrt(r+z)/(2r) * vz
    vv = x*sqrt(r+z)/(2r*ρ) * vx + y*sqrt(r+z)/(2r*ρ) * vy - sqrt(r-z)/(2r) * vz

    #println("vt: $vt, vϕ: $vϕ, vu: $vu, vv: $vv")

    vec_to[1] = vt 
    vec_to[2] = vϕ
    vec_to[3] = vu
    vec_to[4] = vv

    return nothing

end

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::Paraboloidal,::Cartesian) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    sϕ,cϕ = sincos(ϕ)

    #println("x: $x, y: $y, z: $z, r: $r, ρ: $ρ")

    vt = vec_from[1]
    vϕ = vec_from[2]
    vu = vec_from[3]
    vv = vec_from[4]

    vt = vt
    vx = v*cϕ * vu + u*cϕ * vv - u*v*sϕ * vϕ
    vy = v*sϕ * vu + u*sϕ * vv + u*v*cϕ * vϕ
    vz = u * vu - v * vv

    #println("vt: $vt, vϕ: $vϕ, vu: $vu, vv: $vv")

    vec_to[1] = vt 
    vec_to[2] = vx
    vec_to[3] = vy
    vec_to[4] = vz

    return nothing

end

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::Cylindrical,::Paraboloidal) where T

    t = pos_from[1]
    ρ = pos_from[2]
    ϕ = pos_from[3]
    z = pos_from[4]

    r = sqrt(ρ^2+z^2)

    vt = vec_from[1]
    vρ = vec_from[2]
    vϕ = vec_from[3]
    vz = vec_from[4]

    vt = vt
    vu = sqrt(r-z)/(2r) * vρ + sqrt(r+z)/(2r) * vz
    vv = sqrt(r+z)/(2r) * vρ - sqrt(r-z)/(2r) * vz

    vec_to[1] = vt 
    vec_to[2] = vϕ
    vec_to[3] = vu
    vec_to[4] = vv

    return nothing

end

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::Paraboloidal,::Cylindrical) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    vt = vec_from[1]
    vϕ = vec_from[2]
    vu = vec_from[3]
    vv = vec_from[4]

    vt = vt
    vρ = v * vu + u * vv
    vz = u * vu - v * vv

    vec_to[1] = vt 
    vec_to[2] = vρ
    vec_to[3] = vϕ
    vec_to[4] = vz

    return nothing

end