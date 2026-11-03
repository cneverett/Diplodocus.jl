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
abstract type Space3TensorForce <: AbstractForce end
abstract type SpaceMatrixForce <: AbstractForce end
abstract type SpaceVectorForce <: AbstractForce end
abstract type SpaceScalarForce <: AbstractForce end
abstract type AnalyticForce <: AbstractForce end

@kwdef struct CoordinateForce <: AbstractForce
    Domain::Union{Vector{Int64},Nothing} = nothing
end

@kwdef struct SyncRadReact <: AnalyticForce
    mode::AbstractMode = Ani()
    B::Union{Float64,Nothing} = nothing
    B_sampled::Union{Vector{Float64},Nothing} = nothing
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