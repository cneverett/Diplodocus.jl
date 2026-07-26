# -----------------------------------------------------------------------------
# Workspace specialised for exponential Rosenbrock-Euler:
#     exp(tau*A)u0 + phi_1(tau*A)u1
# Typical use: A = h*J, u1 = h*g, tau = 1.
# Assumptions: real Float32/Float64, p = 1, task1 = false, one output time.
# -----------------------------------------------------------------------------

mutable struct KIOPSRosenbrockWorkspace{T,AM<:AbstractMatrix,AV<:AbstractVector}
    n::Int
    mmin::Int
    mmax::Int
    orth_len::Int

    # Core KIOPS storage for the augmented dimension n+1.
    V::AM                 # (n+1) x (mmax+1)
    H::AM                 # (mmax+1) x (mmax+1)
    w::AV                 # output/current physical vector, length n
    u_phi::AM             # n x 1, stores scaled u1 forcing vector
    av::AV                # length n, A*v work vector
    Atilde::AM            # (n+1) x (n+1), augmented matrix for Rosenbrock-Euler

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

    panelW  :: AM   # (n+1) x panel_max   raw candidates / QR buffer
    panelG  :: AM   # orth_len x panel_max projection against existing basis
    panelR  :: AM   # panel_max x panel_max QR triangular factor
    panelTmp:: AM   # optional scratch, same size as panelG or panelW if needed
end

function KIOPSRosenbrockWorkspace(A::AbstractMatrix, u::AbstractVector; mmin::Integer = 10, mmax::Integer = 128,
orth_len::Integer = 4,panel_max::Integer = 4)

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
        mat(n, 1),             # u_phi
        vec(n),                # av
        mat(n + 1, n + 1),     # Atilde
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
        mat(1,1),               # normtmp
        mat(n + 1, panel_max), # panelW
        mat(orth_len, panel_max), # panelG
        mat(panel_max, panel_max), # panelR
        mat(kmax, panel_max) # panelTmp
    )
end

# -----------------------------------------------------------------------------
# Backend-generic helpers.
# -----------------------------------------------------------------------------

@inline function _matrix_one_norm!(ws::KIOPSRosenbrockWorkspace, A::AbstractMatrix)
    ws.absM .= abs.(A)
    sum!(ws.colsum, ws.absM)
    return Float64(_host_scalar(maximum(ws.colsum)))
end

@inline function _norm!(out::AbstractMatrix, in::AbstractVector)
    out .= sqrt.(sum(abs2, in, dims=1))    # elementwise sqrt on-device
    return nothing
end

function _set_one_entry!(A::AbstractMatrix, i::Integer, j::Integer, value)
    fill!(@view(A[i, j]), value)
    return nothing
end

function _set_one_entry!(A::AbstractVector, i::Integer, value)
    fill!(@view(A[i]), value)
    return nothing
end

function _host_entry(A::AbstractMatrix, i::Integer, j::Integer)
    return Array(@view(A[i, j]))[1]
end

function _host_entry(A::AbstractVector, i::Integer)
    return Vector(@view(A[i]))[1]
end

_host_scalar(x::Number) = x
_host_scalar(x) = Array(x)[1]


# -----------------------------------------------------------------------------
# In-place Padé scaling-and-squaring exponential on the active k x k block only.
# This keeps the original full-matrix exp(H) path, not the e1-action variant.
# -----------------------------------------------------------------------------
function expm_kiops_pade!(F::M, A::M,ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    
    #nrm = _matrix_one_norm!(ws, A)
    nrm = norm(A,1)
    #println("nrm = $nrm, 1-norm: $(norm(A,1))")

    if T === Float32
        theta3 = 4.258730016922831f-1
        theta5 = 1.880152677804762f0
        theta7 = 3.925724783138660f0
        theta9 = 6.42276f0 # more digits?
        if nrm <= theta3
            return _pade3!(F, A, ws)
        elseif nrm <= theta5
            return _pade5!(F, A, ws)
        else
            s = max(0, ceil(Int, log2(nrm / theta7)))
            println("nrm = $nrm, theta7 = $theta7, s = $s")
            if isodd(s) # force s to be even to allow mul step below to always give F without need for mem copies
                s += 1
            end
            ws.Hscaled .= A ./ T(2)^s
            #println("nrmH: $(norm(ws.Hscaled)), maxH: $(maximum(ws.Hscaled)), minH: $(minimum(ws.Hscaled))")
            _pade7!(F, ws.Hscaled, ws)
            T1 = ws.tmp1
            for sv in 1:s
                if isodd(sv)
                    mul!(T1, F, F)
                else
                    mul!(F, T1, T1)
                end
                #println("nrmF: $(norm(F)), maxF: $(maximum(F)), minF: $(minimum(F))")
                #mul!(T1, F, F)
                #copyto!(F,T1)
            end
            return F
        end
    else
        theta3  = 1.495585217958292e-2
        theta5  = 2.539398330063230e-1
        theta7  = 9.504178996162932e-1
        theta9  = 2.097847961257068e0
        theta13 = 5.371920351148152e0
        if nrm <= theta3
            return _pade3!(F, A, ws)
        elseif nrm <= theta5
            return _pade5!(F, A, ws)
        elseif nrm <= theta7
            return _pade7!(F, A, ws)
        elseif nrm <= theta9
            return _pade9!(F, A, ws)
        else
            s = max(0, ceil(Int, log2(nrm / theta13))) + 16
            println("nrm = $nrm, theta13 = $theta13, s = $s")
            ws.Hscaled .= A ./ T(2)^s
            _pade13!(F, ws.Hscaled, ws)
            for _ in 1:s
                mul!(ws.tmp1, F, F)
                F .= ws.tmp1
            end
            return F
        end
    end
end

function _finish_pade!(F::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    ws.lhs .= ws.Vmat .- ws.Umat
    #ws.rhs .= T(2) .* ws.Umat
    ws.rhs .= ws.Vmat .+ ws.Umat
    L = ws.lhs
    R = ws.rhs
    #ldiv!(R, lu!(L),R)
    #luL = lu(L)
    #ldiv!(R,luL,R)
    #F .= R .+ ws.Id

    LAPACK.gesv!(L, R)
    F .= R #.+ ws.Id
    return F
end

function _pade3!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    U .= T(1) .* A2 .+ T(60) .* ws.Id
    mul!(ws.tmp1, A, U)
    U .= ws.tmp1
    V .= T(12) .* A2 .+ T(120) .* ws.Id
    return _finish_pade!(F, ws)
end

function _pade5!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    mul!(A4, A2, A2)
    U .= T(1) .* A4 .+ T(420) .* A2 .+ T(15120) .* ws.Id
    mul!(ws.tmp1, A, U)
    U .= ws.tmp1
    V .= T(30) .* A4 .+ T(3360) .* A2 .+ T(30240) .* ws.Id
    return _finish_pade!(F, ws)
end

function _pade7!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    A6 = ws.A6
    U = ws.Umat 
    V = ws.Vmat
    T1 = ws.tmp1
    I1 = ws.Id

    mul!(A2, A, A)
    mul!(A4, A2, A2)
    mul!(A6, A2, A4)
    U .= T(1) .* A6 .+ T(1512) .* A4 .+ T(277200) .* A2 .+ T(8648640) .* I1
    mul!(T1, A, U)
    U .= T1
    V .= T(56) .* A6 .+ T(25200) .* A4 .+ T(1995840) .* A2 .+ T(17297280) .* I1

    return _finish_pade!(F, ws)
end

function _pade9!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    A6 = ws.A6
    A8 = ws.A8
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    mul!(A4, A2, A2)
    mul!(A6, A2, A4)
    mul!(A8, A4, A4)
    U .= T(1) .* A8 .+ T(3960) .* A6 .+ T(2162160) .* A4 .+ T(302702400) .* A2 .+ T(8821612800) .* ws.Id
    mul!(ws.tmp1, A, U)
    U .= ws.tmp1
    V .= T(90) .* A8 .+ T(110880) .* A6 .+ T(30270240) .* A4 .+ T(2075673600) .* A2 .+ T(17643225600) .* ws.Id
    return _finish_pade!(F, ws)
end

function _pade13!(F::M, A::M, ws::KIOPSRosenbrockWorkspace{T}) where {T,M<:AbstractMatrix{T}}

        A2 = ws.A2
        A4 = ws.A4
        A6 = ws.A6
        U = ws.Umat
        V = ws.Vmat
        tmp1 = ws.tmp1
        tmp2 = ws.tmp2
        Ik = ws.Id

        mul!(A2, A, A)
        mul!(A4, A2, A2)
        mul!(A6, A2, A4)

        tmp1 .= T(1) .* A6 .+ T(16380) .* A4 .+ T(40840800) .* A2
        mul!(tmp2, A6, tmp1)
        tmp1 .= tmp2 .+ T(33522128640) .* A6 .+ T(10559470521600) .* A4 .+ T(1187353796428800) .* A2 .+ T(32382376266240000) .* Ik
        mul!(U, A, tmp1)

        tmp1 .= T(182) .* A6 .+ T(960960) .* A4 .+ T(1323241920) .* A2
        mul!(tmp2, A6, tmp1)
        V .= tmp2 .+ T(670442572800) .* A6 .+ T(129060195264000) .* A4 .+ T(7771770303897600) .* A2 .+ T(64764752532480000) .* Ik
    return _finish_pade!(F, ws)
end

# -----------------------------------------------------------------------------
# Specialised in-place KIOPS for Rosenbrock-Euler.
# Computes ws.w = exp(tau*A)u0 + phi_1(tau*A)u1.
# For usual Exprb-Euler: A = h*J, u1 = h*g, tau = 1.
# -----------------------------------------------------------------------------

function kiops_roe_noview!(ws::KIOPSRosenbrockWorkspace{T}, tau_out::Real,
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
    Atilde = ws.Atilde
    normtmp = ws.normtmp

    #fill!(V, zero(T))
    #fill!(H, zero(T))
    #fill!(w, zero(T))
    #fill!(u_phi, zero(T))
    #fill!(av, zero(T))

    # Scale the phi_1 input as in KIOPS, but p=1 only.

    #normU = sum(abs,u1)
    #if normU > 0.0
    #    ex = ceil(Int, log2(normU))
    #    nu = T(exp2(-ex))
    #    mu = T(exp2(ex))
    #else
    #    nu = one(T)
    #    mu = one(T)
    #end

    # form augmented matrix Atilde = [A u_phi; 0 0] for Rosenbrock-Euler
    @views copyto!(Atilde[1:n, 1:n], A)
    @views copyto!(Atilde[1:n, n+1:n+1], #=nu *=# u1)

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

        fill!(H, zero(T))
        fill!(V, zero(T))

        if j == 0

            @views V[1:n, 1] .= w 
            _set_one_entry!(V, n + 1, 1, T(1))

            beta = norm(@view(V[:,1])) #_norm!(normtmp, @view(V[:,1]))
            @view(V[:,1]) ./= beta #normtmp

        end

        while j < m
            j += 1

            @inbounds begin

            # version with augmented matrix
            mul!(@view(V[:,j+1]), Atilde, @view(V[:,j]))

            i0 = max(1,j-orth_len+1)
            Vblk = @view V[:,i0:j]
            vnew = @view V[:,j+1]
            hblk = @view H[i0:j,j]

            Vblkt = transpose(Vblk)

            mul!(hblk, Vblkt, vnew) 
            mul!(vnew, Vblk, hblk, -one(T), one(T))  # vnew .= vnew - Vblk*hblk
            
            _norm!(normtmp, @view(V[:,j+1]))

            
            last_nrm = Float64(_host_scalar(normtmp))
            if last_nrm < tol
                happy = true
                println("happy")
                break
            end

            copyto!(@view(H[j+1,j]), normtmp)
            @view(V[:,j+1]) ./= normtmp

            end # inbounds
            krystep += 1
        end

        # KIOPS error estimator: add the temporary H[1,j+1] coupling, compute exp(tau*H),
        # then restore H[j+1,j].  Only the active (j+1)x(j+1) block is exponentiated.
        _set_one_entry!(H, 1, j + 1, one(T))
        #nrm = happy ? 0.0 : last_nrm
        _set_one_entry!(H, j + 1, j, zero(T)) # resets residual norm in H

        #kactive = j + 1
        #fill!(ws.Hexp, zero(T))
        ws.Hexp .= H
        ws.Hexp .*= T(sgn * tau)

        #display(H)
        #error("stop")

        expm_kiops_pade!(ws.F, ws.Hexp, ws)
        exps += 1
        #_set_one_entry!(H, j + 1, j, T(nrm))

        #display(ws.F)

        nrm = Float64(_host_scalar(normtmp))

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

            println("err = $err, beta = $beta, nrm = $nrm, ferr = $ferr")

            if isnan(err)
                @warn "KIOPS Rosenbrock-Euler encountered NaN error estimate."
                tau = ldexp(tau, -1)  # halve the timestep
                continue
            end

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
        m = 16#32 # 64#m_new
        
        if ireject >= 10
            @warn "KIOPS Rosenbrock-Euler failed to converge after 10 rejections."
            m_new = 17#33#65 # triggers timestep rejection and reduction in dt 
            break
        end
    end

    stats = (steps = step,rejected = reject,krylov_steps = krystep,exponentials = exps,m_final = m_new#=m=#)
    return w, stats
end

# -----------------------------------------------------------------------------
# Outer Rosenbrock/KIOPS loop: unchanged structure, but the inner Arnoldi loop
# is now panel-based.
# -----------------------------------------------------------------------------
function kiops_roe_panel!(
    ws::KIOPSRosenbrockWorkspace{T},
    tau_out::Real,
    A::AbstractMatrix,
    u0::AbstractVector,
    u1::AbstractVector;
    tol::Real = 1.0e-7,
    m_init::Integer = ws.mmax,
    panel_size::Integer = 4,     # use 4 or 8
    reorth::Bool = false,
) where {T}

    n = ws.n
    mmin = ws.mmin
    mmax = ws.mmax
    m = Int(max(mmin, min(m_init, mmax)))

    V = ws.V
    H = ws.H
    w = ws.w
    u_phi = ws.u_phi
    av = ws.av
    Atilde = ws.Atilde
    normtmp = ws.normtmp

    fill!(V, zero(T))
    fill!(H, zero(T))
    fill!(Atilde, zero(T))

    # Same scaling logic you already have
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

    @views copyto!(Atilde[1:n, 1:n], A)
    @views copyto!(Atilde[1:n, n+1:n+1], nu * u1)

    sgn = sign(Float64(tau_out))
    tau_now = 0.0
    tau_end = abs(Float64(tau_out))
    tau = tau_end

    j = 0

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
    last_nrm = 0.0
    beta = 1.0
    m_new = m

    # form first basis vector v0 = [u0; 1] with some scaling
    v0 = vcat(u0,mu)
    beta = norm(v0)
    v0 ./= beta

    while tau_now < tau_end

        fill!(H, zero(T))

        nrm = arnoldi_block!(ws, Atilde,v0,2;orth_len=2)

        #c0 = @views transpose(ws.V[:,1:4]) * v0
        #c0 .*= T(beta)
        #c = CUDA.zeros(T,m)
        #@views c[1:4] .= c0
        # ---------------------------------------------------------------------
        # Same reduced exponential / error-estimation path as before
        # ---------------------------------------------------------------------
        #@views copyto!(H[1:m,m+1],c)
        _set_one_entry!(H, 1, m+1, one(T))
        _set_one_entry!(H, m + 1, m, zero(T))

        ws.Hexp .= H
        ws.Hexp .*= T(sgn * tau)

        Hm = @view H[1:m, 1:m]

        @show any(!isfinite, Hm)
        @show maximum(abs, Hm)
        @show minimum(abs, Hm)

        
        println("norm H: $(norm(H)), $sgn, $tau")
        display(H)
        #error("stop")

        expm_kiops_pade!(ws.F, ws.Hexp, ws)
        exps += 1

        display(ws.F)

        #nrm = norm(@view(V[:, m + 1]))
        ferr = _host_entry(ws.F, m, m + 1)
        err = abs(beta * nrm * Float64(ferr))

        println("err = $err, beta = $beta, nrm = $nrm, ferr = $ferr")

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
            m_opt = max(mmin, min(mmax, m_opt))
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

        if omega <= delta
            reject += ireject
            step += 1

            #wtmp = ws.F[1:m,1:m] * c 
            #@views mul!(w, V[1:n, 1:m], wtmp)

            @views mul!(w, V[1:n, 1:m], ws.F[1:m,1])
            w .*= T(beta)

            tau_now += tau
            j = 0
            ireject = 0
        else
            ireject += 1
            _set_one_entry!(H, 1, m + 1, zero(T))
        end

        oldtau = tau
        tau = tau_new
        oldm = m
        m = 64  # keep your current policy if you want

        if ireject >= 10
            @warn "KIOPS Rosenbrock-Euler failed to converge after 10 rejections."
            m_new = 65
            break
        end
    end

    stats = (
        steps = step,
        rejected = reject,
        krylov_steps = krystep,
        exponentials = exps,
        m_final = m_new,
    )

    return w, stats
end

function kiops_roe_panel_new!(
    ws::KIOPSRosenbrockWorkspace{T},
    tau_out::Real,
    A::AbstractMatrix,
    u0::AbstractVector,
    u1::AbstractVector;
    tol::Real = 1.0e-7,
    m_init::Integer = ws.mmax,
    panel_size::Integer = 4,     # use 4 or 8
    reorth::Bool = false,
) where {T}

    n = ws.n
    mmin = ws.mmin
    mmax = ws.mmax
    m = 16 # Int(max(mmin, min(m_init, mmax)))

    w = ws.w
    av = ws.av
    Atilde = ws.Atilde

    fill!(Atilde, zero(T))


    @views copyto!(Atilde[1:n, 1:n], A)
    @views copyto!(Atilde[1:n, n+1:n+1], u1)

    sgn = sign(Float64(tau_out))
    tau_now = 0.0
    tau_end = abs(Float64(tau_out))
    tau = tau_end

    j = 0

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
    last_nrm = 0.0
    beta = 1.0
    m_new = m

    # form first basis vector v0 = [u0; 1] with some scaling
    v0 = vcat(u0,T(1))

    while tau_now < tau_end

        H = ws.H
        V = ws.V

        β, rho = build_panels!(ws,Atilde, v0, 4; s = 4)

        display(H)
        println("tau = $tau, beta = $β, rho = $rho, maxH = $(maximum(H)), minH = $(minimum(H))")

        ws.Hexp .= ws.H
        ws.Hexp .*= T(sgn * tau)

        #Hexp = zeros(Float64,size(ws.Hexp))
        #copyto!(Hexp, ws.Hexp)

        #copyto!(ws.F,exp(Hexp))

        expm_kiops_pade!(ws.F, ws.Hexp, ws)
        exps += 1

        #display(ws.F)

        #nrm = norm(@view(V[:, m + 1]))
        ferr = _host_entry(ws.F, m, m + 1)
        err = abs(β * rho * Float64(ferr))

        println("err = $err, beta = $β , nrm = $rho, ferr = $ferr")

        if isnan(err)
            
            @warn "KIOPS Rosenbrock-Euler encountered NaN error estimate."
            #tau = ldexp(tau, -1)  # halve the timestep
            #continue
            m_new = 17 
            break
        end

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
            m_opt = max(mmin, min(mmax, m_opt))
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

        if omega <= delta
            reject += ireject
            step += 1

            #wtmp = ws.F[1:m,1:m] * c 
            #@views mul!(w, V[1:n, 1:m], wtmp)

            @views mul!(w, V[1:n, 1:m], ws.F[1:m,1])
            w .*= T(β)

            #w .+= u0

            tau_now += tau
            j = 0
            ireject = 0
        else
            ireject += 1
            _set_one_entry!(H, 1, m + 1, zero(T))
        end

        oldtau = tau
        tau = tau_new
        oldm = m
        m = 16#32  # keep your current policy if you want

        if ireject >= 10
            @warn "KIOPS Rosenbrock-Euler failed to converge after 10 rejections."
            m_new = 17 # 65
            break
        end
    end

    stats = (
        steps = step,
        rejected = reject,
        krylov_steps = krystep,
        exponentials = exps,
        m_final = m_new,
    )

    return w, stats
end

# -----------------------------------------------------------------------------
# Panel Arnoldi with QR orthonormal panels.
#
# Idea:
#   1) Build a small candidate panel W on-device
#   2) Project the panel seed against only the previous 2 basis vectors
#      (your IOC(2) rule) to keep the local recurrence cheap
#   3) QR the panel to get orthonormal columns
#   4) Update the reduced matrix block from projection + R
#
# Important: this is NOT the old scalar orthogonalization loop.
# -----------------------------------------------------------------------------
 function arnoldi_block!(ws::KIOPSRosenbrockWorkspace{T},Atilde::AbstractMatrix,v0::AbstractVector,s::Int;npanels::Int = round(Int,(size(ws.V,2)-1)/s),orth_len::Int = 4,tol::Real = 1.0e-12) where {T}

    n = length(v0)-1

    V = ws.V
    H = ws.H
    W = @view ws.panelW[:, 1:s]
    R = @view ws.panelR[1:s, 1:s]
    G = @view ws.panelG[1:orth_len, 1:s]

    # first block 

    @views V[:,1] .= v0
    #=for k = 2:s
        @views mul!(W[:,k], Atilde, W[:,k-1])
    end

    #Vprev = @view V[:, 1]
    #Gs = @view G[1:1, 1:3]
    #Ws = @view W[:, 1:3]
    #Gs .= transpose(Vprev) * Ws
    #Ws .-= Vprev * Gs

    # orthonormalise first block 
    F = qr!(W)     # thin QR
    if W isa Matrix 
        Q = Matrix(F.Q)
    elseif W isa CuArray
        Q = CuArray(F.Q)
    else
        error("Unsupported array type for W")
    end
    @views V[:, 1:s] .= Q=#

    # block recurance

    old_start = 1
    old_end = 1
    new_start = old_start + 1
    new_end = old_end + s - 1 # for first block

    for j in 1:npanels

        if j == 1
            b = 1
            sb = s-1
            old_start = 1
            old_end = 1
            new_start = old_end + 1
            new_end = new_start + sb - 1
        else
            b = s
            sb = s
            old_start = (j-2)*s + 1
            old_end = (j-1)*s
            new_start = old_end + 1
            new_end = new_start + s - 1
        end

        Vprev = @view V[:, old_start:old_end]

        # Form the next raw power panel:
        # [A^(4(j-1))b, A^(4(j-1)+1)b, A^(4(j-1)+2)b, A^(4(j-1)+3)b]

        @views mul!(W[:,1],Atilde,Vprev[:,end])
        for k = 2:sb
            @views mul!(W[:,k], Atilde, W[:,k-1])
        end

        # IOP(q=4): orthogonalize against previous block (last 4 vectors)
        mul!(@view(G[1:b,1:sb]), transpose(Vprev), @view(W[:,1:sb]))
        @view(W[:,1:sb]) .-= Vprev * @view(G[1:b,1:sb])

        # QR to get the new orthonormal block
        F = qr!(@view(W[:,1:sb]))     # thin QR
        if W isa Matrix 
            Q = Matrix(F.Q)
        elseif W isa CuArray
            Q = CuArray(F.Q)
        else
            error("Unsupported array type for W")
        end


        # Store explicit block Hessenberg coefficients
        # diagonal block coupling from previous block
        if j == 1 
            @views V[:, new_start:new_end] .= Q
            Vseed = @view V[:, 1:s]
            @view(H[1:s,1:s]) .= Vseed' * (Atilde * Vseed)
            #display(@views H[1:s,1:s])
        elseif j == nblocks+1
            @views H[old_start:old_end, old_start:old_end] .= G[1:b,:]
            # subdiagonal block from QR of the residual
            @views w = V[:, end-1]
            @views aw = Atilde * w
            @views h = V[:, 1:end-2]' * aw
            @views r = aw - V[:, 1:end-2] * h
            display(size(r))
            display(size(@view(H[new_start:new_start, old_end:old_end])))
            nrm = norm(r)
            println("nrm = ", nrm)
            _set_one_entry!(H, new_start, old_end, nrm)
            @views V[:,end] .= r ./ nrm

            return nrm
        else
            @views V[:, new_start:new_end] .= Q
            @views H[old_start:old_end, old_start:old_end] .= G[1:b,:]
            # subdiagonal block from QR of the residual
            @views H[new_start:new_end, old_start:old_end] .= F.R
        end
    end

    return nrm 

end

# Panel-4 Krylov basis with QR panels.
# Basis layout:
#   V[:, 1:4]   = first panel
#   V[:, 5:8]   = second panel
#   ...
#
# H layout:
#   H[1:k, 1:k] = reduced operator for the built basis
#   H[k+1, k]   = scalar residual used in the KIOPS-style error estimate
#   H has size (k+1) x (k+1), so it can also be reused as the augmented
#   small matrix in the timestep function.

function build_panels!(ws,A, u0, npanels::Int; s::Int = 4)
    
    T = promote_type(eltype(A), eltype(u0))
    n = length(u0)
    k = s * npanels

    V = ws.V
    H = ws.H
    W = ws.panelW
    g = ws.panelG

    β = norm(u0)
    q1 = u0 / β
    @views V[:, 1] .= q1

    # ------------------------------------------------------------
    # First panel: [q1, q2, q3, q4]
    # q2:q4 come from powers of A*q1, then are orthonormalized
    # against q1 with QR.
    # ------------------------------------------------------------
    @views W[:, 1] .= A * q1
    for j = 2:s
        @views W[:, j] .= A * W[:, j - 1]
    end

    g1 = @view g[1:1, 1:s-1]

    # Remove the q1 component from the three candidate vectors
    @views mul!(g1, transpose(q1), W[:, 1:s-1])          # 1 x 3
    @views mul!(W[:, 1:s-1], q1, g1, -1, 1)

    F = qr!(view(W, :, 1:s-1))
    @views V[:, 2:s] .= Matrix(F.Q)              # n x 3

    V1 = @view V[:, 1:s]
    @views H[1:s, 1:s] .= transpose(V1) * (A * V1)

    # ------------------------------------------------------------
    # Later panels:
    # Wraw = [A w, A^2 w, A^3 w, A^4 w]
    # where w is the last vector of the previous panel
    # Orthogonalize only against the previous panel, then QR.
    # Store the panel recurrence Wraw = Vprev*C + Vnew*R.
    # ------------------------------------------------------------
    for p = 2:npanels
        prev_start = s * (p - 2) + 1
        prev_end   = s * (p - 1)
        new_start  = s * (p - 1) + 1
        new_end    = s * p

        Vprev = @view V[:, prev_start:prev_end]

        # Candidate panel from powers of the last vector of the previous panel
        @views mul!(W[:, 1], A, Vprev[:, end])
        for j = 2:s
            @views mul!(W[:,j], A, W[:,j-1])
        end

        # Orthogonalize only against the previous panel
        mul!(g, transpose(Vprev) , W)
        mul!(W, Vprev, g, -1, 1)

        # QR gives the next basis panel
        F = qr!(W)
        #Q .= CuMatrix(Fp.Q)
        #@views V[:, new_start:new_end] .= Q
        @views V[:, new_start:new_end] .= Matrix(F.Q)

        # Panel recurrence:
        # Wraw = Vprev * C + Vnew * R
        H[prev_start:prev_end, new_start:new_end] .= g
        H[new_start:new_end, new_start:new_end] .= F.R
    end

    # ------------------------------------------------------------
    # Scalar residual for the KIOPS-style error estimate
    # (panel analogue: use the last basis vector of the last panel)
    # ------------------------------------------------------------
    vlast = @view V[:, k]
    av    = A * vlast
    @views hlast = V[:, 1:k]' * av
    @views r     = av - V[:, 1:k] * hlast
    rho   = norm(r)

    #H[k + 1, k] = rho

    # Augmented column for the φ1 trick
    fill!(@view(H[1, k + 1]), one(T))

    return β, rho
end

# -----------------------------------------------------------------------------
# Helper: initialize the first Krylov basis vector safely.
# This removes the j=0 indexing hazard.
# -----------------------------------------------------------------------------
function init_krylov_basis!(
    ws::KIOPSRosenbrockWorkspace{T},
    u0::AbstractVector,
    u1::AbstractVector,
    mu::T,
) where {T}

    V = ws.V
    H = ws.H
    w = ws.w
    n = ws.n

    fill!(H, zero(T))
    @. w = u0

    wnorm2 = Float64(_host_scalar(dot(w, w)))
    beta = sqrt(wnorm2 + Float64(mu * mu))

    @views V[1:n, 1] .= w ./ T(beta)
    _set_one_entry!(V, n + 1, 1, mu / T(beta))

    return beta
end
