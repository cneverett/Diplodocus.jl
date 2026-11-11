abstract type AbstractInteraction end

"""
    BinaryInteraction(name1,name2,name3,name4)

A struct to contain information pertaining to a particular emissive interaction. 
# Fields
- `name1::String`: name of incoming particle 1
- `name2::String`: name of incoming particle 2
- `name3::String`: name of outgoing particle 3
- `name4::String`: name of outgoing particle 4
"""
struct BinaryInteraction <: AbstractInteraction
    name1::String
    name2::String
    name3::String
    name4::String
end

"""
    EmissiveInteraction()

A struct to contain information pertaining to a particular emissive interaction. 
# Fields
- `name1::String`: name of incoming particle (absorbed)
- `name2::String`: name of outgoing particle 1 (emitted)
- `name3::String`: name of outgoing particle 2 (emitted)
- `EmiName::String`: name of the emissive interaction e.g. "Sync" for synchrotron radiation
- `Ext_sampled::Vector{Float64}`: vector of sampled values of the external parameter for which the emissive interaction matrices are built e.g. sampled magnetic field values for synchrotron radiation
- `mode::AbstractMode`: mode of the emissive interaction e.g. "Iso", "Axi" or "Ani"
- `Force::Bool`: if name1==name2 the emissive interaction may be treated as reactive force e.g. radiation reaction for synchrotron. In which case this force can be applied directly, within the emission matrix or separately as a regular force.
- `Domain::Union{Vector{Int64},Nothing}`: vector of spatial indices (offset_space) referring to spatial sub-domains in which to apply the emissive interaction. If `nothing` the interaction will be applied to all spatial sub-domains.
"""
@kwdef struct EmissiveInteraction <: AbstractInteraction
    name1::String
    name2::String
    name3::String
    EmiName::String
    Ext_sampled::Vector{Float64}
    mode::AbstractMode
    Force::Bool = true
    Domain::Union{Vector{Int64},Nothing} = nothing
end

"""
    BinaryMatricesStruct()

A struct for storing the big matrices associated with binary interactions in the simulation.
"""
struct BinaryMatricesStruct{T<:Union{Float32,Float64}}
    
    M_Bin::AbstractMatrix{T}  # big matrix for binary interactions
    Binary_list::Vector{BinaryInteraction} # list of binary interactions
    Domain::Union{Vector{Int64},Nothing} # domain of the binary interaction in the state vector

end

struct BinaryMatricesStructPatankar{T<:Union{Float32,Float64}}
    
    Gijk::AbstractMatrix{T}  # gain matrix for binary interactions
    Liij::AbstractMatrix{T}  # loss matrix for binary interactions
    Binary_list::Vector{BinaryInteraction} # list of binary interactions
    Domain::Union{Vector{Int64},Nothing} # domain of the binary interaction in the state vector

end

struct BinaryMatricesStructPatankarSymmetric{T<:Union{Float32,Float64}}
    
    Gijk::AbstractMatrix{T}  # gain matrix for binary interactions
    Lijk::AbstractMatrix{T}  # loss matrix for binary interactions
    Liij::AbstractMatrix{T}  # loss matrix for binary interactions
    Liik::AbstractMatrix{T}  # loss matrix for binary interactions
    Binary_list::Vector{BinaryInteraction} # list of binary interactions
    Domain::Union{Vector{Int64},Nothing} # domain of the binary interaction in the state vector

end

"""
    EmissionMatricesStruct()

A struct for storing the big matrices associated with emission interactions and their corresponding reaction forces if applicable in the simulation.
"""
struct EmissionMatricesStruct{T<:Union{Float32,Float64}}
    
    M_Emi::Vector{Union{Matrix{T},SparseMatrixCSC{T,Int32}}}  # big matrix for emission interactions
    Emission_list::Vector{EmissiveInteraction} # list of emission interactions

end