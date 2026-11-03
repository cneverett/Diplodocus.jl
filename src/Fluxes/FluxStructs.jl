"""
    FluxMatricesStruct()

A struct for storing the flux matrices associated with the simulation.
"""
struct FluxMatricesStruct{T<:Union{Float32,Float64}}

    Ap_Flux::Vector{T}                  # Forward Boundary time flux 
    Am_Flux::Vector{T}                  # Backward Boundary time flux
    
    X_Flux::SparseMatrixCSC{T,Int32}    # sum of space fluxes 
    P_Flux::SparseMatrixCSC{T,Int32}    # sum of momentum fluxes 

    Vol::Vector{T}                      # SpaceTime volume element

    B_Flux::SparseMatrixCSC{T,Int32}    # B Flux through x boundaries
    C_Flux::SparseMatrixCSC{T,Int32}    # C Flux through y boundaries
    D_Flux::SparseMatrixCSC{T,Int32}    # D Flux through z boundaries

    I_Flux::SparseMatrixCSC{T,Int32}    # I Flux through px boundaries
    J_Flux::SparseMatrixCSC{T,Int32}    # J Flux through py boundaries
    K_Flux::SparseMatrixCSC{T,Int32}    # K Flux through pz boundaries

end

#========= Force =========#
#=========================#

abstract type AbstractForce end
abstract type AnalyticForce <: AbstractForce end
abstract type NumericalForce <: AbstractForce end
abstract type Space3TensorForce <: NumericalForce end
abstract type SpaceMatrixForce <: NumericalForce end
abstract type SpaceVectorForce <: NumericalForce end
abstract type SpaceScalarForce <: NumericalForce end


@kwdef struct CoordinateForce <: AnalyticForce
    Domain::Union{Vector{Int64},Nothing} = nothing
end

@kwdef struct SyncRadReact <: AnalyticForce
    mode::AbstractMode = Ani()
    B::Union{Float64,Nothing} = nothing
    B_sampled::Union{Vector{Float64},Nothing} = nothing
    Domain::Union{Vector{Int64},Nothing} = nothing
end

"""
    LocalLorentzForceZerothOrder()

The Lorentz force acting on a particle as measured by a local observer who measures E parallel B may be expanded in powers of the gyro-displacement. This struct represents the zeroth-order contribution to that force which corresponds to a homogenous electromagnetic field. 
"""
@kwdef struct LocalLorentzForceZerothOrder <: SpaceMatrixForce
    Domain::Union{Vector{Int64},Nothing} = nothing
end

"""
    LocalLorentzForce()

The Lorentz force acting on a particle as measured by a local observer who measures E parallel B and B is in the e_3 direction. 
"""
@kwdef struct LocalLorentzForce <: SpaceMatrixForce
    Domain::Union{Vector{Int64},Nothing} = nothing
end

@kwdef struct FirstOrderGuidingCentre <: SpaceVectorForce 
    Domain::Union{Vector{Int64},Nothing} = nothing
end

struct ExB <: AbstractForce
    E0::Float64
    B0::Float64
end

struct GradBInvZDecay <: AbstractForce
end