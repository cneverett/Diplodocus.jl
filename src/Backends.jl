abstract type AbstractBackend end

struct CPUBackend <: AbstractBackend end
struct CUDABackend <: AbstractBackend end

# zeros functions
function Base.zeros(::CPUBackend, T::Type{<:AbstractFloat}, dims::Int...)
    return zeros(T, dims...)
end
function Base.zeros(::CUDABackend, T::Type{<:AbstractFloat}, dims::Int...)
    #@assert T == Float32 "CUDA backend only supports Float32 type."
    return CUDA.zeros(T, dims...)
end

# ones functions
function Base.ones(::CPUBackend, T::Type{<:AbstractFloat}, dims::Int...)
    return ones(T, dims...)
end
function Base.ones(::CUDABackend, T::Type{<:AbstractFloat}, dims::Int...)
    #@assert T == Float32 "CUDA backend only supports Float32 type."
    return CUDA.ones(T, dims...)
end

# sparse conversion 
function SparseArrays.sparse(::CUDABackend,Precision::Type{<:AbstractFloat},SA::SparseMatrixCSC{T,I}) where {T<:AbstractFloat,I<:Integer}
    return CUDA.CUSPARSE.CuSparseMatrixCSR(Precision.(SA))
end

# sparse conversion 
function SparseArrays.sparse(::CPUBackend,Precision::Type{<:AbstractFloat},SA::SparseMatrixCSC{T,I}) where {T<:AbstractFloat,I<:Integer}
    return SparseMatrixCSC{Precision}{Int32}(Precision.(SA))
end

# copy conversions
function Base.copy(::CUDABackend,Precision::Type{<:AbstractFloat}, A::AbstractArray{T}) where {T<:AbstractFloat}
    return CUDA.CuArray{Precision}(Precision.(copy(A)))
end

function Base.copy(::CPUBackend,Precision::Type{<:AbstractFloat}, A::AbstractArray{T}) where {T<:AbstractFloat}
    return Array{Precision}(Precision.(copy(A)))
end


function device_array(::CPUBackend,Precision::Type{<:AbstractFloat}, A::AbstractArray{T}) where {T<:AbstractFloat}
    # returns a new CPU array of type Precision with the same values as A
    return Array{Precision}(Precision.(A))
end

function device_array(::CUDABackend,Precision::Type{<:AbstractFloat}, A::AbstractArray{T}) where {T<:AbstractFloat}
    # returns a new GPU array of type Precision with the same values as A
    return CUDA.CuArray{Precision}(Precision.(A))
end