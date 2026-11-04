"""
    RicciRotationComponents!(pos,ω,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad)

Returns the Ricci rotation coefficient components 
```math
ω=ω^a_{~bc}=-e_b^{~β}e_c^{~γ}∇_γe^a_{~β}=-e_b^{~β}e_c^{~γ}\\left(∂_γ - Γ^α_{~βγ}e^a_{~α}\\right)
```
at a point `pos` in spacetime defined by the coordinates `pos=(t,x,y,z)`. The Ricci rotation coefficient components are stored in the static array ω, which is modified in-place. The specific form of the Christoffel symbols depends on the type of metric and coordinate basis being used. Gradients of the tetrad components are calculated using automatic differentiation.
"""
function RicciRotationComponents!(pos::MVector{4,T},ω::MArray{Tuple{4,4,4},Float64,3,64},metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) where T

    fill!(ω, zero(T))

    e = MMatrix{4,4,T,16}(zeros(Float64,4,4))
    inve = MMatrix{4,4,T,16}(zeros(Float64,4,4))
    Γ = MArray{Tuple{4,4,4},T,3,64}(zeros(Float64,4,4,4))
    ∂inve = MArray{Tuple{4,4,4},T,3,64}(zeros(Float64,4,4,4))

    ChristoffelComponents!(pos,Γ,metric,coordinates)
    TetradComponents!(pos,e,metric,coordinates,tetrad)
    InverseTetradComponents!(pos,inve,metric,coordinates,tetrad)

    @inline f! = function (inve_local, xvec)
        InverseTetradComponents!(metric,coordinates,tetrad,inve_local,xvec)
        return nothing
    end

    cfg = ForwardDiff.JacobianConfig(f!, inve, pos)
    ForwardDiff.jacobian!(∂inve, f!, inve, pos, cfg)

    @inbounds for a in 1:4, b in 1:4, c in 1:4 
        tmp = zero(T)
        for β in 1:4, γ in 1:4
            tmp -= e[b,β]*e[c,γ]*∂inve[β,a,γ]
            for α in 1:4
                tmp += e[b,β]*e[c,γ]*Γ[α,β,γ]*inve[α,a]  
            end
        end
        ω[a,b,c] = tmp
    end

    return nothing

end


function RicciRotationComponents!(pos::MVector{4,T},ω::MArray{Tuple{4,4,4},Float64,3,64},metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad,e::MMatrix{4,4,T,16},inve::MMatrix{4,4,T,16},Γ::MArray{Tuple{4,4,4},T,3,64},∂inve::MArray{Tuple{4,4,4},T,3,64},cfg::ForwardDiff.JacobianConfig) where T

    # Non-allocating version of RicciRotationComponents! that takes pre-allocated arrays for the tetrad, inverse tetrad, Christoffel symbols, the jacobian of the inverse tetrad, and the ForwardDiff JacobianConfig.

    fill!(ω, zero(T))

    ChristoffelComponents!(pos,Γ,metric,coordinates)
    TetradComponents!(pos,e,metric,coordinates,tetrad)
    InverseTetradComponents!(pos,inve,metric,coordinates,tetrad)

    ForwardDiff.jacobian!(∂inve, func!, inve, pos, cfg)

    @inbounds for a in 1:4, b in 1:4, c in 1:4 
        tmp = zero(T)
        for β in 1:4, γ in 1:4
            tmp -= e[b,β]*e[c,γ]*∂inve[β,a,γ]
            for α in 1:4
                tmp += e[b,β]*e[c,γ]*Γ[α,β,γ]*inve[α,a]  
            end
        end
        ω[a,b,c] = tmp
    end

    return nothing

end

"""
    CoordinateForceSpaceIntegrand!(pos,CFSpaceArray,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad,e,inve,Γ,∂inve,cfg)

Returns the spacetime array associated with calculations of coordinate force fluxes `CFSpaceArray`. This is equal to the product of the Ricci rotation coefficient components `ω` and the spacetime volume element `χ` 
```math
\\text{CFSpaceArray}=ω^a_{~bc}χ=-e_b^{~β}e_c^{~γ}∇_γe^a_{~β}=-e_b^{~β}e_c^{~γ}\\left(∂_γ - Γ^α_{~βγ}e^a_{~α}\\right)χ
```
at a point `pos` in spacetime defined by the coordinates `pos=(t,x,y,z)`. The resulting components are stored in the static array `CFSpaceArray`, which is modified in-place. The specific form of the Christoffel symbols depends on the type of metric and coordinate basis being used. Gradients of the tetrad components are calculated using automatic differentiation.

This function is non-allocating so must be supplied with pre-allocated arrays for the tetrad, inverse tetrad, Christoffel symbols, the jacobian of the inverse tetrad, and the ForwardDiff JacobianConfig. This is to avoid allocations when this function is called repeatedly during numerical integration of the coordinate force fluxes.
"""
function CoordinateForceSpaceIntegrand!(pos::MVector{4,Float64},CFSpaceArray::MArray{Tuple{4,4,4},Float64,3,64},metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad,e::MMatrix{4,4,T,16},inve::MMatrix{4,4,T,16},Γ::MArray{Tuple{4,4,4},T,3,64},∂inve::MArray{Tuple{4,4,4},T,3,64},func!,cfg::ForwardDiff.JacobianConfig) where T

    fill!(CFSpaceArray, zero(T))

    ChristoffelComponents!(pos,Γ,metric,coordinates)
    TetradComponents!(pos,e,metric,coordinates,tetrad)
    InverseTetradComponents!(pos,inve,metric,coordinates,tetrad)

    fill!(∂inve, zero(T))

    ForwardDiff.jacobian!(∂inve, func!, inve, pos, cfg)

    χ::T = VolumeElement(pos,metric,coordinates)

    @inbounds for a in 1:4, b in 1:4, c in 1:4 
        tmp = zero(T)
        for β in 1:4, γ in 1:4
            tmp -= e[b,β]*e[c,γ]*∂inve[β,a,γ]
            for α in 1:4
                tmp += e[b,β]*e[c,γ]*Γ[α,β,γ]*inve[α,a]  
            end
        end
        CFSpaceArray[a,b,c] = iszero(χ) ? zero(T) : -tmp * χ # -Γ*χ to get the correct sign
    end

    return nothing

end