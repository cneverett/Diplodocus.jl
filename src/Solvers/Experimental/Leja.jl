function leja_interval(A::AbstractMatrix{T}) where {T<:Real}
    n, m = size(A)
    n == m || throw(ArgumentError("matrix must be square"))

    lo = Inf
    hi = -Inf

    @inbounds for i in 1:n
        r = zero(T)
        for j in 1:n
            j == i && continue
            r += abs(A[i, j])
        end
        c = A[i, i]
        lo = min(lo, c - r)
        hi = max(hi, c + r)
    end

    return lo, hi
end

# Sparse CSC matrices
function leja_interval(A::SparseMatrixCSC{T,Ti}) where {T<:Real,Ti<:Integer}
    n, m = size(A)
    n == m || throw(ArgumentError("matrix must be square"))

    rad = zeros(T, n)
    nz = nonzeros(A)
    rv = rowvals(A)

    @inbounds for col in 1:n
        for p in nzrange(A, col)
            row = rv[p]
            row == col && continue
            rad[row] += abs(nz[p])
        end
    end

    lo = Inf
    hi = -Inf

    @inbounds for i in 1:n
        c = A[i, i]
        lo = min(lo, c - rad[i])
        hi = max(hi, c + rad[i])
    end

    return lo, hi
end

# ----------------------------
# Main Leja routine
#
# Computes phi_1(dt*A)*v
# ----------------------------
function phi1_leja_dense(A,v,Leja_nodes;dt,tau_min = nothing,mu_est = nothing,gamma = 1.5,maxdeg = 80,tol = 1e-10,use_dense_mu2 = false,use_directional_mu = true)

    n = size(A, 1)
    @assert size(A, 2) == n
    @assert length(v) == n

    # 1. Estimate tau_min if not provided
    if tau_min === nothing
        tau_min = estimate_tau_min_from_diag(A)
    end
    #println(norm(A,Inf), " ", tau_min)
    #tau_min = 1/norm(A,Inf)

    # 2. Estimate positive growth endpoint
    if mu_est === nothing
        if use_dense_mu2
            # More reliable, but costs eigmax of symmetric dense matrix
            mu_est = estimate_mu2_dense(A)
        elseif use_directional_mu
            # Cheap and vector-specific
            mu_est = estimate_mu_from_v(A, v)
        else
            # No positive extension
            mu_est = 0.0
        end
    end

    a, b = leja_interval(A)
    a *= dt
    b *= 0.0 #dt
    #b = -10.0
    a -= 10.0
    #b = 0.0#10.0

    # 3. Build interval for z = dt*lambda
    a = -gamma * dt / tau_min
    b = max(0.0, dt * mu_est)

    if !(b > a)
        error("Invalid Leja interval: [$a, $b]")
    end

    # 4. Shift and scale interval [a,b] to [-2,2]
    #
    # z = q + theta*xi
    #
    q = (a + b) / 2
    theta = (b - a) / 4

    # Bhat = (dt*A - q*I)/theta
    # We do not need to form I explicitly for the recurrence,
    # but for dense A this is fine either way.
    I_n = Matrix{eltype(A)}(I, n, n)
    Bhat = (dt .* A .- q .* I_n) ./ theta

    # 5. Leja nodes xi on [-2,2]
    xi = Leja_nodes

    # 6. Interpolate g(xi) = phi_1(q + theta*xi)
    gvals = [phi1_scalar(q + theta * x) for x in xi]

    # Divided differences with respect to xi, not z
    coeffs = divided_differences(xi, gvals)

    # 7. Newton recurrence applied to Bhat
    #
    # p_m(Bhat)v =
    # c_1 v
    # + c_2 (Bhat - xi_1 I)v
    # + c_3 (Bhat - xi_2 I)(Bhat - xi_1 I)v
    # + ...
    #
    y = coeffs[1] .* v
    w = copy(v)

    last_inc_norm = Inf
    err_est = Inf
    degree_used = 0

    for m in 2:maxdeg
        # w <- (Bhat - xi[m-1] I) w
        w = Bhat * w .- xi[m-1] .* w

        inc = coeffs[m] .* w
        y_new = y .+ inc

        inc_norm = norm(inc)

        # Cheap heuristic Leja error indicator.
        # max of last two increments is slightly safer than one increment.
        err_est = (inc_norm + last_inc_norm) / norm(y_new)

        y = y_new
        last_inc_norm = inc_norm
        degree_used = m - 1

        if m  > 64 && err_est < tol
            return y, (
                err_est = err_est,
                degree_used = degree_used,
                interval = (a, b),
                tau_min = tau_min,
                mu_est = mu_est,
                q = q,
                theta = theta,
            )
        end
    end

    return y, (
        err_est = err_est,
        degree_used = degree_used,
        interval = (a, b),
        tau_min = tau_min,
        mu_est = mu_est,
        q = q,
        theta = theta,
        warning = "maximum degree reached before tolerance",
    )
end


# ----------------------------
# Scalar phi_1(z) = (exp(z)-1)/z
# ----------------------------
function phi1_scalar(z::Real)
    if abs(z) < 1e-8
        # Taylor expansion near zero
        return 1 + z/2 + z^2/6 + z^3/24 + z^4/120 + z^5/720
    else
        return expm1(z) / z
    end
end

# ----------------------------
# Approximate tau_min from A diagonal
#
# If A has decay rates on the diagonal, e.g. a_ii ≈ -1/tau_i,
# then the largest-magnitude negative diagonal entry gives
#
# tau_min ≈ 1 / max_i(-a_ii)
# ----------------------------
function estimate_tau_min_from_diag(A)
    neg_rates = [-A[i,i] for i in axes(A,1) if A[i,i] < 0]

    if isempty(neg_rates)
        error("No negative diagonal entries found; cannot estimate tau_min this way.")
    end

    fastest_rate = maximum(neg_rates)
    return 1 / fastest_rate
end

# ----------------------------
# Cheap growth estimate from current vector v
#
# mu_v = v'Av / v'v
#
# This is not an upper bound for mu_2, but is a cheap
# direction-dependent growth estimate.
# ----------------------------
function estimate_mu_from_v(A, v)
    Av = A * v
    return dot(v, Av) / dot(v, v)
end

# ----------------------------
# Safer but more expensive dense estimate:
#
# mu_2 = lambda_max((A + A')/2)
#
# Since A is dense, this may be acceptable for moderate size.
# ----------------------------
function estimate_mu2_dense(A)
    H = Symmetric((A + A') / 2)
    return eigmax(H)
end

# ----------------------------
# Simple approximate Leja nodes on [-2,2]
#
# For production use, you should use a proper precomputed Leja sequence.
# This greedy construction is acceptable for experimentation.
# ----------------------------
function leja_nodes_interval(n::Int; a=-2.0, b=2.0, grid_size=20_000)
    grid = collect(range(a, b, length=grid_size))
    nodes = Vector{Float64}(undef, n)

    # Start at right endpoint; common for real Leja sequences
    nodes[1] = b

    logprod = zeros(Float64, grid_size)

    for k in 2:n
        # Update log product: log prod_j |x - xi_j|
        @inbounds for i in eachindex(grid)
            logprod[i] += log(abs(grid[i] - nodes[k-1]) + eps())
        end

        idx = argmax(logprod)
        nodes[k] = grid[idx]
    end

    return nodes
end

# ----------------------------
# Newton divided differences
#
# Input:
#   z[j] = interpolation nodes in physical z-plane
#   f[j] = phi1(z[j])
#
# Output:
#   c[j] = Newton coefficients
# ----------------------------
function divided_differences(z, f)
    c = copy(f)
    n = length(z)

    for j in 2:n
        for i in n:-1:j
            c[i] = (c[i] - c[i-1]) / (z[i] - z[i-j+1])
        end
    end

    return c
end