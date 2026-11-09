################### Modified ExponetialUtilities.jl functions for Arnoldi ####################

# old function (slightly modified) from ExponentialUtilities.jl, kept for reference
@inline function arnoldi_orthogonalize!(y,V,H::AbstractMatrix{U},j::Integer,lo::Integer,hi::Integer) where U
    @inbounds for i in lo:hi
        α = U(dot(@view(V[:, i]), y))
        H[i, j] = α
        axpy!(-α, @view(V[:, i]), y)
    end
    return nothing
end

#=@inline function arnoldi_orthogonalize!(y,V,H::AbstractMatrix{U},j::Integer,lo::Integer,hi::Integer) where U
    # this version is only correct if the columns of V are actually orthogonal, use with caution
    @inbounds begin 
        @views mul!(H[lo:hi,j], transpose(V[:,lo:hi]),y,one(U),zero(U)) # H = V' * y
        @views mul!(y,V[:,lo:hi],H[lo:hi,j],-one(U),one(U)) # y = y - V * H
    end
end=#

@inline function arnoldi_reorthogonalize!(y,V,H::AbstractMatrix{U},j::Integer,lo::Integer,hi::Integer) where U
    @inbounds for i in lo:hi
        α = U(dot(@view(V[:, i]), y))
        H[i, j] += α
        axpy!(-α, @view(V[:, i]), y)
    end
end

@inline function arnoldi_reorthogonalize_block!(y, V, H::AbstractMatrix{U},G,c, j::Integer, lo::Integer, hi::Integer) where U
    @views begin
        W = V[:, lo:hi]

        mul!(c[lo:hi], transpose(W), y)          # use transpose(W) if real-only
        mul!(G[lo:hi, lo:hi], transpose(W), W)

        ldiv!(UnitLowerTriangular(G[lo:hi, lo:hi]), c[lo:hi]) # forward substitution

        H[lo:hi, j] .+= c[lo:hi]
        mul!(y, W, c[lo:hi], -one(U), one(U))
    end
end


function firststep!(Ks::KrylovSubspace, V, H, b)
    return @inbounds begin
        Ks.beta = norm(b)
        if !iszero(Ks.beta)
            @. V[:, 1] = b / Ks.beta
        end
    end
end

function arnoldi_step!(j::Integer, iop::Integer, A::AT,V::AbstractMatrix{T}, H::AbstractMatrix{U};reorthogonalize::Bool = false,remove_drift::Bool = false) where {AT, T, U}

    x, y = @view(V[:, j]), @view(V[:, j + 1])
    mul!(y, A, x) # y = A * x

    lo = max(1, j - iop + 1)

    # First orthogonalization pass
    arnoldi_orthogonalize!(y, V, H, j, lo, j) 

    # Second pass (full reorthogonalization over the same window)
    if reorthogonalize
        arnoldi_reorthogonalize!(y, V, H, j, lo, j) # G = V' * V, c = V' * y
    end

    # remove energy drift from orthogonalisation
    #drift_before = abs(sum(y)) / (sqrt(length(y))*norm(y))
    if remove_drift
        α = sum(y) / length(y)
        y .-= α
    end

    #drift_after = abs(sum(y)) / (sqrt(length(y))*norm(y))
    #println("drift after = $drift_after, drift before = $drift_before, j = $j")

    β = H[j + 1, j] = norm(y)
    if !iszero(β)
        @. y /= β
    end
    return β
end

function arnoldi!(Ks::KrylovSubspace{T1, U}, A::AT, b; tol::Real = 1.0e-7,m::Int = min(Ks.maxiter, size(A, 1)),opnorm = nothing, iop::Int = 0,reorthogonalize::Bool = false, remove_drift::Bool = false, init::Int = 0, t::Number = NaN, mu::Number = NaN, l::Int = -1) where {T1 <: Number, U <: Number, AT}

    Ks.wasbreakdown = false

    m > Ks.maxiter ? resize!(Ks, m) : Ks.m = m

    @inbounds V, H = ExponentialUtilities.getV(Ks), ExponentialUtilities.getH(Ks)
    fill!(H, zero(eltype(H)))
    fill!(V, zero(eltype(V)))

    if iszero(init)
        ExponentialUtilities.firststep!(Ks::KrylovSubspace, V, H, b)
        init = 1
    end

    iszero(Ks.beta) && return Ks
    iszero(iop) && (iop = m)

    #G = similar(Ks.H, Ks.m, Ks.m)
    #c = similar(Ks.H, Ks.m)

    for j in init:m
        beta = arnoldi_step!(j, iop, A, V, H; reorthogonalize = reorthogonalize,remove_drift=remove_drift)
        if beta < tol
            Ks.m = j
            Ks.wasbreakdown = true
            break
        end
    end

    return Ks
end


################ Experimental block Arnoldi implementation ################



function arnoldi_scalar!(A::MT1, u::VT1, V::MT2, H::MT2, w::VT2; reorth::Bool = false) where {T1,T2,MT1<:AbstractMatrix{T1},VT1<:AbstractVector{T1},MT2<:AbstractMatrix{T2},VT2<:AbstractVector{T2}}

    n = length(u)
    m = size(H, 2) - 1

    @assert size(A, 1) == n && size(A, 2) == n
    @assert size(V, 1) == n && size(V, 2) >= m + 1
    @assert size(H, 1) >= m + 1
    @assert length(w) == n

    fill!(H, zero(eltype(H)))
    fill!(V, zero(eltype(V)))

    β0::T2 = norm(u)
    iszero(β0) && throw(ArgumentError("Initial vector u has zero norm"))

    copyto!(@view(V[:, 1]), u)
    rmul!(@view(V[:, 1]), inv(β0)) 

    #println("β0: $β0, sum(u) = $(sum(u)), norm(u) = $(norm(u))")

    hkp1k = zero(T2)

    @inbounds for k in 1:m

        Vk  = @view V[:, 1:k]
        vk  = @view V[:, k]
        vkp = @view V[:, k+1]
        hk  = @view H[1:k, k]

        Vkt = transpose(Vk)

        mul!(vkp, A, vk)
        mul!(hk, Vkt, vkp)
        mul!(vkp, Vk, hk, -one(T2), one(T2))

        #println("norm(w): ", norm(vkp))
        #println("hk: ", hk)

        if reorth # TODO: This currently doesn't work!!!!
            mul!(hk, Vkt, vkp)
            mul!(vkp, Vk, hk, -one(T2), one(T2))
        end
        
        #println("norm(w): ", norm(vkp))
        #println("hk: ", hk)

        hkp1k = norm(vkp)

        if hkp1k < 1e2 * eps(T2) # if the next basis vector is too small, stop the iteration
            println("happy")
            fill!(vkp, zero(T2))
            break
        end

        fill!(@view(H[k+1, k]), hkp1k)

        rmul!(vkp, inv(hkp1k))

    end

    # fill last entry for error estimate 
    fill!(@view(H[1, m+1]), one(T2))
    # remove last residual for error estimate form 
    fill!(@view(H[m+1, m]), zero(T2))

    return m, hkp1k, β0
end

function arnoldi_block4!(A::MT1, v::VT1, V::MT2, H::MT2, W::MT2, R::MT2; reorth::Bool = false) where {T1,T2,MT1<:AbstractMatrix{T1},VT1<:AbstractVector{T1},MT2<:AbstractMatrix{T2}}

    T = eltype(V)
    n = length(v)
    s = 4

    #@assert size(W, 1) == n && size(W, 2) == s
    #@assert size(R1, 1) == s && size(R1, 2) == s
    #@assert size(V, 1) == n
    #@assert size(V, 2) % s == 0
    #@assert size(H, 1) == size(H, 2) == size(V, 2)

    nb = (size(V, 2)-1) ÷ s   # total number of 4-column basis block

    # Initial block: [v, Av, A^2v, A^3v]
    @views begin
        copyto!(W[:, 1], v)
        for j in 2:s
            mul!(W[:, j], A, W[:, j-1])
        end
    end

    F0 = qr!(W)
    if W isa Matrix 
        Q = Matrix(F0.Q)
    elseif W isa CuArray
        Q = CuArray(F0.Q)
    else
        error("Unsupported array type for W")
    end
    @views V[:, 1:s] .= Q

    # last column for error estimate
    @views H[1:s, end] .= F0.R[:, 1]

    # first block of H
    @views mul!(W, A, V[:, 1:s])
    @views mul!(H[1:s, 1:s], adjoint(V[:, 1:s]), W)

    # Generate remaining blocks
    for j in 1:nb-1
        c1 = (j - 1) * s + 1
        c2 = j * s
        n1 = c2 + 1
        n2 = (j + 1) * s

        Vprev = @view V[:, 1:c2]
        Vj    = @view V[:, c1:c2]

        # W = A * current block (all 4 columns)
        mul!(W, A, Vj)

        # Orthogonalize against all previous basis vectors
        Hproj = @view H[1:c2, c1:c2]
        mul!(Hproj, adjoint(Vprev), W)
        mul!(W, Vprev, Hproj, -one(T), one(T))

        if reorth
            mul!(Hproj, adjoint(Vprev), W)
            mul!(W, Vprev, Hproj, -one(T), one(T))
        end

        Fs = qr!(W)
        if W isa Matrix 
            Q = Matrix(Fs.Q)
        elseif W isa CuArray
            Q = CuArray(Fs.Q)
        else
            error("Unsupported array type for W")
        end

        @views V[:, n1:n2] .= Q
        @views H[n1:n2, c1:c2] .= Fs.R

        if j == nb - 1
            # last column for error estimate
            R .= Fs.R
        end
    end

    return R
end