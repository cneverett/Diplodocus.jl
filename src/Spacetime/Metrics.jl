"""
    MetricComponents!(pos,g,::AbstractMetric,::AbstractCoordinates,)

Returns the metric components ``g=g_{αβ}`` at a given point `pos` in spacetime defined by the coordinates ``pos=(t,x,y,z)``. The metric components are stored in the static matrix g, which is modified in-place. The specific form of the metric components depends on the type of spacetime and coordinate basis being used.
"""
MetricComponents!(pos,g,metric::AbstractMetric,coordinates::AbstractCoordinates) = error("Metric function not defined for this metric $(typeof(metric)) and coordinate $(typeof(coordinates)).")

"""
    VolumeElement!(pos,::AbstractMetric,::AbstractCoordinates)

Returns the volume element `χ=sqrt(-det(g))`, where `g` is the metric, at a given point `pos` in spacetime defined by the coordinates ``pos=(t,x,y,z)``. The metric is not supplied as the volume element can be evaluated directly.
"""
VolumeElement!(pos,metric::AbstractMetric,coordinates::AbstractCoordinates) = error("Volume element function not defined for this metric $(typeof(metric)) and coordinate $(typeof(coordinates)).")

#===== Minkowski =======#
    #====== Cartesian ========#
        @inline function MetricComponents!(pos::MVector{4,T},g::MMatrix{4,4,T,16},::Minkowski,::Cartesian) where T

            g[1,1]::T = -1.0
            g[2,2]::T = 1.0
            g[3,3]::T = 1.0
            g[4,4]::T = 1.0

            return nothing

        end
        @inline function VolumeElement(pos::MVector{4,T},::Minkowski,::Cartesian) where T
            return T(1.0)
        end

    #====== Cylindrical ========#    
        @inline function MetricComponents!(pos::MVector{4,T},g::MMatrix{4,4,T,16},::Minkowski,::Cylindrical) where T

            t = pos[1]
            ρ = pos[2]
            ϕ = pos[3]
            z = pos[4]

            g[1,1]::T = -1.0
            g[2,2]::T = 1.0
            g[3,3]::T = ρ^2
            g[4,4]::T = 1.0

            return nothing 

        end
        @inline function VolumeElement(pos::MVector{4,T},::Minkowski,::Cylindrical) where T

            t = pos[1]
            ρ = pos[2]
            ϕ = pos[3]
            z = pos[4]

            return T(ρ)
        end
    
    #====== Spherical ========#
        @inline function MetricComponents!(pos::MVector{4,T},g::MMatrix{4,4,T,16},::Minkowski,::Spherical) where T

            t = pos[1]
            r = pos[2]
            θ = pos[3]
            ϕ = pos[4]

            g[1,1]::T = -1.0
            g[2,2]::T = 1.0
            g[3,3]::T = r^2
            g[4,4]::T = r^2*sin(θ)^2

            return nothing

        end
        @inline function VolumeElement(pos::MVector{4,T},::Minkowski,::Spherical) where T

            t = pos[1]
            r = pos[2]
            θ = pos[3]
            ϕ = pos[4]

            return T(r^2*sin(θ))
        end

    #====== Paraboloidal ========#
        @inline function MetricComponents!(pos::MVector{4,T},g::MMatrix{4,4,T,16},::Minkowski,::Paraboloidal) where T

            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]

            g[1,1]::T = -1.0
            g[2,2]::T = u^2*v^2
            g[3,3]::T = u^2+v^2
            g[4,4]::T = u^2+v^2

            return nothing

        end
        @inline function VolumeElement(pos::MVector{4,T},::Minkowski,::Paraboloidal) where T

            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]

            return T(u*v*(u^2+v^2))
        end