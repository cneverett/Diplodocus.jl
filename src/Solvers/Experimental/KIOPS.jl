# -----------------------------------------------------------------------------
# Workspace specialised for exponential Rosenbrock-Euler:
#     exp(tau*A)u0 + phi_1(tau*A)u1
# Typical use: A = h*J, u1 = h*g, tau = 1.
# Assumptions: real Float32/Float64, p = 1, task1 = false, one output time.
# -----------------------------------------------------------------------------

#=mutable struct KIOPSRosenbrockWorkspace{T,AM<:AbstractMatrix,AV<:AbstractVector}
    n::Int
    mmin::Int
    mmax::Int
    orth_len::Int

    # Core KIOPS storage for the augmented dimension n+1.
    V::AM                 # (n+1) x (mmax+1)
    H::AM                 # (mmax+1) x (mmax+1)
    w::AV                 # output/current physical vector, length n
    u_phi::AV             # n x 1, stores scaled u1 forcing vector
    av::AV                # length n, A*v work vector

    # Dense active Hessenberg exponential storage, max size (mmax+1)^2.
    Hexp::AM
    Hscaled::AM
    F::AM
    A2::AM
    A4::AM
    A6::AM
    A8::AM
    Umat::AM
    Vmat::AM
    tmp1::AM
    tmp2::AM
    lhs::AM
    rhs::AM
    Id::AM
    absM::AM
    colsum::AM

    # Taylor action workspace
    y1::AV
    y2::AV
    e1::AV
    ej::AV
    tmpv::AV

    # norm tmp 
    normtmp::AM
end

function KIOPSRosenbrockWorkspace(A::AbstractMatrix, u::AbstractVector;
                                  mmin::Integer = 10,
                                  mmax::Integer = 128,
                                  orth_len::Integer = 2)
    T = eltype(A)
    n = size(A, 1)
    kmax = Int(mmax) + 1

    mat(r, c) = begin
        X = similar(A, T, r, c)
        fill!(X, zero(T))
        X
    end
    vec(r) = begin
        x = similar(A, T, r)
        fill!(x, zero(T))
        x
    end

    Id = mat(kmax, kmax)
    copyto!(Id, Matrix{T}(I, kmax, kmax))

    return KIOPSRosenbrockWorkspace{T,typeof(mat(1, 1)),typeof(vec(1))}(
        n, Int(mmin), Int(mmax), Int(orth_len),
        mat(n + 1, kmax),      # V
        mat(kmax, kmax),       # H
        vec(n),                # w
        vec(n),                # u_phi
        vec(n),                # av
        mat(kmax, kmax),       # Hexp
        mat(kmax, kmax),       # Hscaled
        mat(kmax, kmax),       # F
        mat(kmax, kmax),       # A2
        mat(kmax, kmax),       # A4
        mat(kmax, kmax),       # A6
        mat(kmax, kmax),       # A8
        mat(kmax, kmax),       # Umat
        mat(kmax, kmax),       # Vmat
        mat(kmax, kmax),       # tmp1
        mat(kmax, kmax),       # tmp2
        mat(kmax, kmax),       # lhs
        mat(kmax, kmax),       # rhs
        Id,                    # Id
        mat(kmax, kmax),       # absM
        mat(1, kmax),           # colsum
        vec(kmax),             # y1
        vec(kmax),             # y2
        vec(kmax),             # e1
        vec(kmax),             # ej
        vec(kmax),             # tmpv
        mat(1,1)               # normtmp
    )
end=#

# -----------------------------------------------------------------------------
# In-place Padé scaling-and-squaring exponential on the active k x k block only.
# This keeps the original full-matrix exp(H) path, not the e1-action variant.
# -----------------------------------------------------------------------------

function expm_kiops_pade!(F::M, A::M,
                          ws::KIOPSRosenbrockWorkspace{T},
                          k::Integer) where {T,M<:AbstractMatrix{T}}
    nrm = _matrix_one_norm!(ws, A, k)

    if T === Float32
        theta3 = 4.258730016922831f-1
        theta5 = 1.880152677804762f0
        theta7 = 3.925724783138660f0
        #if nrm <= theta3
        #    return _pade3!(F, A, ws, k)
        #elseif nrm <= theta5
        #    return _pade5!(F, A, ws, k)
        #else
            s = max(0, ceil(Int, log2(nrm / theta7)))
            @views ws.Hscaled[1:k, 1:k] .= A[1:k, 1:k] ./ T(2)^s
            _pade7!(F, ws.Hscaled, ws, k)
            T1 = @view ws.tmp1[1:k, 1:k]
            F1 = @view F[1:k, 1:k]
            for _ in 1:s
                mul!(T1, F1, F1)
                F1 .= T1
            end
            return F
        #end
    else
        theta3  = 1.495585217958292e-2
        theta5  = 2.539398330063230e-1
        theta7  = 9.504178996162932e-1
        theta9  = 2.097847961257068e0
        theta13 = 5.371920351148152e0
        if nrm <= theta3
            return _pade3!(F, A, ws, k)
        elseif nrm <= theta5
            return _pade5!(F, A, ws, k)
        elseif nrm <= theta7
            return _pade7!(F, A, ws, k)
        elseif nrm <= theta9
            return _pade9!(F, A, ws, k)
        else
            s = max(0, ceil(Int, log2(nrm / theta13)))
            @views ws.Hscaled[1:k, 1:k] .= A[1:k, 1:k] ./ T(2)^s
            _pade13!(F, ws.Hscaled, ws, k)
            for _ in 1:s
                @views mul!(ws.tmp1[1:k, 1:k], F[1:k, 1:k], F[1:k, 1:k])
                @views F[1:k, 1:k] .= ws.tmp1[1:k, 1:k]
            end
            return F
        end
    end
end

function _finish_pade!(F::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    @views ws.lhs[1:k, 1:k] .= ws.Vmat[1:k, 1:k] .- ws.Umat[1:k, 1:k]
    @views ws.rhs[1:k, 1:k] .= T(2) .* ws.Umat[1:k, 1:k]
    L = @view ws.lhs[1:k, 1:k]
    R = @view ws.rhs[1:k, 1:k]
    luL = lu(L)
    #@views ldiv!(ws.rhs[1:k, 1:k], lu!(ws.lhs[1:k, 1:k]), ws.rhs[1:k, 1:k])
    ldiv!(R,luL,R)
    @views F[1:k, 1:k] .= R .+ ws.Id[1:k, 1:k]
    return F
end

function _pade3!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    @views begin
        A2 = ws.A2[1:k, 1:k]; U = ws.Umat[1:k, 1:k]; V = ws.Vmat[1:k, 1:k]
        mul!(A2, A[1:k, 1:k], A[1:k, 1:k])
        U .= T(1) .* A2 .+ T(60) .* ws.Id[1:k, 1:k]
        mul!(ws.tmp1[1:k, 1:k], A[1:k, 1:k], U)
        U .= ws.tmp1[1:k, 1:k]
        V .= T(12) .* A2 .+ T(120) .* ws.Id[1:k, 1:k]
    end
    return _finish_pade!(F, ws, k)
end

function _pade5!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    @views begin
        A2 = ws.A2[1:k, 1:k]; A4 = ws.A4[1:k, 1:k]; U = ws.Umat[1:k, 1:k]; V = ws.Vmat[1:k, 1:k]
        mul!(A2, A[1:k, 1:k], A[1:k, 1:k])
        mul!(A4, A2, A2)
        U .= T(1) .* A4 .+ T(420) .* A2 .+ T(15120) .* ws.Id[1:k, 1:k]
        mul!(ws.tmp1[1:k, 1:k], A[1:k, 1:k], U)
        U .= ws.tmp1[1:k, 1:k]
        V .= T(30) .* A4 .+ T(3360) .* A2 .+ T(30240) .* ws.Id[1:k, 1:k]
    end
    return _finish_pade!(F, ws, k)
end

function _pade7!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    A11 = @view A[1:k,1:k];
    A2 = @view ws.A2[1:k, 1:k]
    A4 = @view ws.A4[1:k, 1:k]
    A6 = @view ws.A6[1:k, 1:k]
    U = @view ws.Umat[1:k, 1:k] 
    V = @view ws.Vmat[1:k, 1:k]
    T1 = @view ws.tmp1[1:k, 1:k]
    I1 = @view ws.Id[1:k, 1:k]

    mul!(A2, A11, A11)
    mul!(A4, A2, A2)
    mul!(A6, A2, A4)
    U .= T(1) .* A6 .+ T(1512) .* A4 .+ T(277200) .* A2 .+ T(8648640) .* I1
    mul!(T1, A11, U)
    U .= T1
    V .= T(56) .* A6 .+ T(25200) .* A4 .+ T(1995840) .* A2 .+ T(17297280) .* I1

    return _finish_pade!(F, ws, k)
end

function _pade9!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    @views begin
        A2 = ws.A2[1:k, 1:k]; A4 = ws.A4[1:k, 1:k]; A6 = ws.A6[1:k, 1:k]; A8 = ws.A8[1:k, 1:k]
        U = ws.Umat[1:k, 1:k]; V = ws.Vmat[1:k, 1:k]
        mul!(A2, A[1:k, 1:k], A[1:k, 1:k])
        mul!(A4, A2, A2)
        mul!(A6, A2, A4)
        mul!(A8, A4, A4)
        U .= T(1) .* A8 .+ T(3960) .* A6 .+ T(2162160) .* A4 .+ T(302702400) .* A2 .+ T(8821612800) .* ws.Id[1:k, 1:k]
        mul!(ws.tmp1[1:k, 1:k], A[1:k, 1:k], U)
        U .= ws.tmp1[1:k, 1:k]
        V .= T(90) .* A8 .+ T(110880) .* A6 .+ T(30270240) .* A4 .+ T(2075673600) .* A2 .+ T(17643225600) .* ws.Id[1:k, 1:k]
    end
    return _finish_pade!(F, ws, k)
end

function _pade13!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}, k::Integer) where {T,M<:AbstractMatrix{T}}
    @views begin
        A2 = ws.A2[1:k, 1:k]
        A4 = ws.A4[1:k, 1:k]
        A6 = ws.A6[1:k, 1:k]
        U = ws.Umat[1:k, 1:k]
        V = ws.Vmat[1:k, 1:k]
        tmp1 = ws.tmp1[1:k, 1:k]
        tmp2 = ws.tmp2[1:k, 1:k]
        Ak = A[1:k, 1:k]
        Ik = ws.Id[1:k, 1:k]
    end

        mul!(A2, Ak, Ak)
        mul!(A4, A2, A2)
        mul!(A6, A2, A4)

        tmp1 .= T(1) .* A6 .+ T(16380) .* A4 .+ T(40840800) .* A2
        mul!(tmp2, A6, tmp1)
        tmp1 .= tmp2 .+ T(33522128640) .* A6 .+ T(10559470521600) .* A4 .+ T(1187353796428800) .* A2 .+ T(32382376266240000) .* Ik
        mul!(U, Ak, tmp1)

        tmp1 .= T(182) .* A6 .+ T(960960) .* A4 .+ T(1323241920) .* A2
        mul!(tmp2, A6, tmp1)
        V .= tmp2 .+ T(670442572800) .* A6 .+ T(129060195264000) .* A4 .+ T(7771770303897600) .* A2 .+ T(64764752532480000) .* Ik
    return _finish_pade!(F, ws, k)
end

# -----------------------------------------------------------------------------
# Specialised in-place KIOPS for Rosenbrock-Euler.
# Computes ws.w = exp(tau*A)u0 + phi_1(tau*A)u1.
# For usual Exprb-Euler: A = h*J, u1 = h*g, tau = 1.
# -----------------------------------------------------------------------------

function kiops_roe!(ws::KIOPSRosenbrockWorkspace{T}, tau_out::Real,
                    A::AbstractMatrix, u0::AbstractVector, u1::AbstractVector;
                    tol::Real = 1.0e-7,
                    m_init::Integer = ws.mmax) where {T}
    n = ws.n
    mmin = ws.mmin
    mmax = ws.mmax
    orth_len = ws.orth_len
    m = Int(max(mmin, min(m_init, mmax)))

    V = ws.V
    H = ws.H
    w = ws.w
    u_phi = ws.u_phi
    av = ws.av

    #fill!(V, zero(T))
    #fill!(H, zero(T))
    #fill!(w, zero(T))
    #fill!(u_phi, zero(T))
    #fill!(av, zero(T))

    # Scale the phi_1 input as in KIOPS, but p=1 only.
    @. av = abs(u1)
    normU = Float64(_host_scalar(sum(av)))
    if normU > 0.0
        ex = ceil(Int, log2(normU))
        nu = T(exp2(-ex))
        mu = T(exp2(ex))
    else
        nu = one(T)
        mu = one(T)
    end
    @. u_phi = nu * u1

    sgn = sign(Float64(tau_out))
    tau_now = 0.0
    tau_end = abs(Float64(tau_out))
    tau = tau_end
    j = 0

    @. w = u0

    if tau_end > 1.0
        gamma = 0.2
        gamma_mmax = 0.1
    else
        gamma = 0.9
        gamma_mmax = 0.6
    end
    delta = 1.4

    step = 0
    krystep = 0
    ireject = 0
    reject = 0
    exps = 0

    oldm = -1
    oldtau = NaN
    omega = NaN
    orderold = true
    kestold = true
    order = 0.0
    kest = 2.0
    happy = false
    last_nrm = 0.0
    beta = 0.0
    m_new = m

    while tau_now < tau_end
        if j == 0
            fill!(H, zero(T))

            wnorm2 = Float64(_host_scalar(dot(w, w)))
            beta = sqrt(wnorm2 + Float64(mu * mu))

            @views V[1:n, 1] .= w ./ T(beta)
            _set_one_entry!(V, n + 1, 1, mu / T(beta))
        end

        while j < m
            j += 1

            @views mul!(av, A, V[1:n, j])
            #@views V[1:n, j + 1] .= av
            #display(@views V[n + 1:n + 1, j])
            #mul!(@view(V[1:n, j + 1]), u_phi, @view(V[n + 1:n + 1, j]), one(T), one(T))
            @views @. V[1:n, j + 1] = V[n + 1:n + 1, j] * u_phi + av
            #@views LinearAlgebra.axpby!(V[n + 1:n + 1, j],u_phi,one(T),V[1:n, j + 1])
            fill!(@view(V[n + 1:n + 1, j + 1:j + 1]), zero(T))

            for i in max(1, j - orth_len + 1):j
                hij = _host_scalar(dot(@view(V[:, i]), @view(V[:, j + 1])))
                _set_one_entry!(H, i, j, hij)
                @views axpy!(-hij, V[:, i], V[:, j + 1])
            end

            last_nrm = Float64(_host_scalar(norm(@view(V[:, j + 1]))))
            if last_nrm < tol
                happy = true
                break
            end

            _set_one_entry!(H, j + 1, j, T(last_nrm))
            @views V[:, j + 1] ./= T(last_nrm)
            krystep += 1
        end

        # KIOPS error estimator: add the temporary H[1,j+1] coupling, compute exp(tau*H),
        # then restore H[j+1,j].  Only the active (j+1)x(j+1) block is exponentiated.
        _set_one_entry!(H, 1, j + 1, one(T))
        nrm = happy ? 0.0 : last_nrm
        _set_one_entry!(H, j + 1, j, zero(T))

        kactive = j + 1
        fill!(ws.Hexp, zero(T))
        @views ws.Hexp[1:kactive, 1:kactive] .= H[1:kactive, 1:kactive]
        @views ws.Hexp[1:kactive, 1:kactive] .*= T(sgn * tau)

        # action for output vector
        #_set_e1!(ws.e1, 1)
        #expm_kiops_taylor_action!(ws.y1, ws.Hexp, ws.e1, ws, kactive)
        # action for KIOPS error estimator
        #_set_ej!(ws.ej, j+1)
        #expm_kiops_taylor_action!(ws.y2, ws.Hexp, ws.ej, ws, kactive)

        expm_kiops_pade!(ws.F, ws.Hexp, ws, kactive)
        exps += 1
        _set_one_entry!(H, j + 1, j, T(nrm))

        if happy
            omega = 0.0
            happy = false
            m_new = m
            tau_new = min(max(tau_end - (tau_now + tau), 0.0), tau)
        else
            ferr = _host_entry(ws.F, j, j + 1)
            #ferr = _host_entry(ws.y2,j)
            err = abs(beta * nrm * Float64(ferr))
            oldomega = omega
            omega = tau_end * err / (tau * Float64(tol))

            if m == oldm && tau != oldtau && ireject >= 1
                order = max(1.0, log(omega / oldomega) / log(tau / oldtau))
                orderold = false
            elseif orderold || ireject == 0
                orderold = true
                order = j / 4
            else
                orderold = true
            end

            if m != oldm && tau == oldtau && ireject >= 1
                kest = max(1.1, (omega / oldomega)^(1 / (oldm - m)))
                kestold = false
            elseif kestold || ireject == 0
                kestold = true
                kest = 2.0
            else
                kestold = true
            end

            remaining_time = omega > delta ? tau_end - tau_now : tau_end - (tau_now + tau)
            same_tau = min(remaining_time, tau)

            if omega == 0
                tau_opt = remaining_time
                m_opt = m
            else
                tau_opt = tau * (gamma / omega)^(1 / order)
                tau_opt = min(remaining_time, max(tau / 5, min(5 * tau, tau_opt)))
                m_opt = ceil(Int, j + log(omega / gamma) / log(kest))
                #m_opt = max(mmin,min(mmax,max(floor(Int, 3m / 4),min(m_opt, ceil(Int, 4m / 3)))))
                m_opt = max(mmin,min(mmax,m_opt))
            end

            if j == mmax
                if omega > delta
                    m_new = j
                    tau_new = tau * (gamma_mmax / omega)^(1 / order)
                    tau_new = min(tau_end - tau_now, max(tau / 5, tau_new))
                else
                    tau_new = tau_opt
                    m_new = m
                end
            else
                m_new = m_opt
                tau_new = same_tau
            end
        end

        if omega <= delta
            reject += ireject
            step += 1

            @views mul!(w, V[1:n, 1:j], ws.F[1:j, 1])
            #@views mul!(w, V[1:n, 1:j], ws.y1[1:j])
            @views w .*= T(beta)

            tau_now += tau
            j = 0
            ireject = 0
        else
            ireject += 1
            _set_one_entry!(H, 1, j + 1, zero(T))
        end

        oldtau = tau
        tau = tau_new
        oldm = m
        m = 64#m_new
        
        if ireject >= 10
            @warn "KIOPS Rosenbrock-Euler failed to converge after 10 rejections."
            m_new = 65 # triggers timestep rejection and reduction in dt 
            break
        end
    end

    stats = (steps = step,
             rejected = reject,
             krylov_steps = krystep,
             exponentials = exps,
             m_final = m_new#=m=#)
    return w, stats
end

# Compatibility convenience: accept an n x 2 U with columns [u0, u1].
function kiops_roe!(ws::KIOPSRosenbrockWorkspace{T}, tau_out::Real,
                    A::AbstractMatrix, U::AbstractMatrix;
                    tol::Real = 1.0e-7,
                    m_init::Integer = ws.mmin) where {T}
    @views return kiops_roe!(ws, tau_out, A, U[:, 1], U[:, 2];
                             tol = tol,
                             m_init = m_init)
end

# Allocating wrapper for quick tests only. Prefer kiops_roe! in timesteppers.
function kiops_roe(tau_out::Real, A::AbstractMatrix, u0::AbstractVector, u1::AbstractVector;
                   tol::Real = 1.0e-7,
                   mmin::Integer = 10,
                   mmax::Integer = 128,
                   m_init::Integer = mmin,
                   orth_len::Integer = 2)
    ws = KIOPSRosenbrockWorkspace(A, u0;
                                  mmin = mmin,
                                  mmax = mmax,
                                  orth_len = orth_len)
    w, stats = kiops_roe!(ws, tau_out, A, u0, u1;
                          tol = tol,
                          m_init = m_init)
    return w, stats.m_final, stats
end


##### TAYLOR 

function _matrix_one_norm!(ws, A::AbstractMatrix, k::Integer)
    @views ws.absM[1:k, 1:k] .= abs.(A[1:k, 1:k])
    @views sum!(ws.colsum[:, 1:k], ws.absM[1:k, 1:k])
    return Float64(_host_scalar(maximum(@view(ws.colsum[:, 1:k]))))
end

# Precomputed reciprocal factorial coefficients up to degree 12.
# Kept as tuples for type stability.
@inline function _taylor_coeffs(::Type{Float32})
    return (Float32(1),
            Float32(1),
            Float32(1)/Float32(2),
            Float32(1)/Float32(6),
            Float32(1)/Float32(24),
            Float32(1)/Float32(120),
            Float32(1)/Float32(720),
            Float32(1)/Float32(5040),
            Float32(1)/Float32(40320),
            Float32(1)/Float32(362880),
            Float32(1)/Float32(3628800),
            Float32(1)/Float32(39916800),
            Float32(1)/Float32(479001600))
end

@inline function _taylor_coeffs(::Type{Float64})
    return (Float64(1),
            Float64(1),
            Float64(1)/Float64(2),
            Float64(1)/Float64(6),
            Float64(1)/Float64(24),
            Float64(1)/Float64(120),
            Float64(1)/Float64(720),
            Float64(1)/Float64(5040),
            Float64(1)/Float64(40320),
            Float64(1)/Float64(362880),
            Float64(1)/Float64(3628800),
            Float64(1)/Float64(39916800),
            Float64(1)/Float64(479001600))
end

@inline _taylor_degree(::Type{Float32}) = 8
@inline _taylor_degree(::Type{Float64}) = 12

"""
    expm_kiops_taylor!(F, A, ws, k)

Experimental drop-in replacement for the Padé small-matrix exponential used by
KIOPS. It computes the active k×k block of `exp(A)` by scaling-and-squaring with
an unrolled Horner-form truncated Taylor polynomial.

This is intended mainly for benchmarking against the Padé version.
"""
function expm_kiops_taylor!(F::AbstractMatrix, A::AbstractMatrix,
                            ws, k::Integer)
    T = eltype(A)
    nrm = _matrix_one_norm!(ws, A, k)

    # Conservative scaling so the truncated Taylor polynomial sees a small norm.
    # This is intentionally simple for benchmarking.
    θ = T === Float32 ? Float32(0.5) : Float64(0.5)
    s = max(0, ceil(Int, log2(nrm / Float64(θ))))

    @views begin
        Ak  = A[1:k, 1:k]
        Hs  = ws.Hscaled[1:k, 1:k]
        Fk  = F[1:k, 1:k]
        Tmp = ws.tmp1[1:k, 1:k]
        Ik  = ws.Id[1:k, 1:k]
    end

    scale = ldexp(one(T), -s)
    Hs .= Ak .* scale

    coeffs = _taylor_coeffs(T)
    m = _taylor_degree(T)

    # Start from c_m * I.
    Fk .= coeffs[m + 1] .* Ik

    # Horner loop for the truncated Taylor polynomial.
    for j = m:-1:1
        mul!(Tmp, Hs, Fk)
        Fk .= Tmp
        Fk .+= coeffs[j] .* Ik
    end

    # Squaring phase.
    for _ in 1:s
        mul!(Tmp, Fk, Fk)
        Fk .= Tmp
    end

    return F
end

@inline function _set_e1!(e1::AbstractVector, k::Integer)
    fill!(e1, zero(eltype(e1)))
    @views fill!(e1[1:1], one(eltype(e1)))
    return e1
end
@inline function _set_ej!(ej::AbstractVector, k::Integer)
    fill!(ej, zero(eltype(ej)))
    @views fill!(ej[k:k], one(eltype(ej)))
    return ej
end

# Compute dest = p(A) * src on the active k×k block using Horner form,
# where p(z) = sum_{i=0}^m z^i / i!.
function _taylor_apply_poly!(dest::AbstractVector, A::AbstractMatrix,
                             src::AbstractVector, ws::KIOPSRosenbrockWorkspace,
                             k::Integer)
    T = eltype(A)
    coeffs = _taylor_coeffs(T)
    m = _taylor_degree(T)

    @views begin
        Ak = A[1:k, 1:k]
        dk = dest[1:k]
        sk = src[1:k]
        tk = ws.tmpv[1:k]

        dk .= coeffs[m + 1] .* sk
        for j = m:-1:1
            mul!(tk, Ak, dk)
            dk .= tk
            dk .+= coeffs[j] .* sk
        end
    end

    return dest
end

# Drop-in replacement in spirit of the Padé small-matrix exponential, but
# returns only exp(A)*e1 on the active block.
function expm_kiops_taylor_action!(y::AbstractVector, A::AbstractMatrix,e::AbstractVector,
                                   ws::KIOPSRosenbrockWorkspace, k::Integer)
    T = eltype(A)
    nrm = _matrix_one_norm!(ws, A, k)

    # Conservative scaling; keeps the truncated Taylor kernel stable.
    θ = T === Float32 ? Float32(0.5) : Float64(0.5)
    s = max(0, ceil(Int, log2(nrm / Float64(θ))))
    scale = ldexp(one(T), -s)

    @views begin
        Ak = A[1:k, 1:k]
        Hs = ws.Hscaled[1:k, 1:k]
        yk = y[1:k]

        Hs .= scale .* Ak

        # First evaluate p(H/2^s)e1.
        _taylor_apply_poly!(y, Hs, e, ws, k)

        # Repeated squaring, but keep the computation as actions only.
        for _ in 1:s
            _taylor_apply_poly!(ws.tmpv, Hs, y, ws, k)
            yk .= ws.tmpv[1:k]
        end
    end

    return y
end
