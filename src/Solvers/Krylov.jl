    mutable struct ERBEKrylovStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

            PhaseSpace::PhaseSpaceStruct
            Precision::Type{T}

            Binary_Interactions::Bool
            Emission_Interactions::Bool

            DistributionDomainMask::Union{Vector{Int64},Nothing}
            DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}
            ActiveDomain::Vector{Int64}

            M_Bin::MBT
            Bin_Domain::BD

            M_Emi::Vector{Union{MT,SMT}}
            A_Flux::SMT
            invA_Flux::SMT                  # inv Ap flux for time stepping
            X_Flux::SMT
            P_Flux::SMT

            invImMP::SMT                      # (I-dt*A^{-1}(M_Emi-P_Flux))^{-1} for momentum update

            Vol::Vector{T}
            invA::Vector{T}                 # vector of diagonal entries of invA_Flux for each spatial point (used for scaling)

            Adaptive::Bool
            Implicit::Bool
            dt0::T
            Cr::T
            n_cut::T

            step::Int64

            M_Bin_Mul_Step::MT             # temporary array for matrix multiplication of binary terms
            M_Bin_Mul_Step_reshape::VT     # temporary array for reshaped matrix multiplication of binary terms
            f_init::VT                     # initial distribution function (used by solver to define output struct)
            f::VT                          # current distribution function
            fstep::VT                      # distribution function after a step
            df::VT                         # change in distribution function
            df_Bin::VT                     # change in distribution function due to binary interactions
            df_Emi::VT                     # change in distribution function due to emission interactions
            df_Flux::VT                    # change in distribution function due to fluxes
            df_Inj::VT                     # change in distribution function due to injection of particles
            df_tmp::VT                     # temporary array the size of f for CFL calculations
            f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
            df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

            F::VT                          # vector for implicit solve residuals
            J::MT                          # Jacobian matrix for implicit solve 
            Jsparse::SMT                   # sparse Jacobian matrix for implicit solve
            D::VT             # Diagonal scaling matrix
            Dinv::VT          # Inverse of diagonal scaling matrix
            ϕ::MT                          # matrix of ϕ functions for exponential Rosenbrock method

            fold::VT                        # local distribution function from previous step 
            fout::VT                        # local distribution function from current step
            fscale::VT                      # scaling vector for exponential Rosenbrock method
            δ::VT                          # temporary vector for exponential Rosenbrock method

            E::VT                           # energy vector for correcting step 

            # Kyrlov Subspace specific caches these will be Float64 precision
            m::Int64                        # dimension of Krylov subspace
            V::AbstractMatrix{Float64}                           # Krylov subspace basis vectors
            Hs::AbstractMatrix{Float64}                          # Hessenberg matrix for scalar Arnoldi
            Hsexp::AbstractMatrix{Float64}                       # exponential of Hessenberg matrix for scalar Arnoldi
            w::AbstractVector{Float64}                          # temporary vector for scalar Arnoldi
            Hb::AbstractMatrix{Float64}                          # Hessenberg matrix for block Arnoldi
            Hbexp::AbstractMatrix{Float64}                       # exponential of Hessenberg matrix for block Arnoldi
            W::AbstractMatrix{Float64}                           # temporary matrix for block Arnoldi
            R0::AbstractMatrix{Float64}                          # matrix for first block residual for block Arnoldi
            Hmmp1::AbstractMatrix{Float64}                       # Hessenberg matrix for block Arnoldi with last residual block

            # Pade Exponential workspace 
            PadeExpWorkspace_s::PadeExpWorkspaceStruct{Float64} # for scalar Arnoldi
            PadeExpWorkspace_b::PadeExpWorkspaceStruct{Float64} # for block Arnoldi


            dt_guess::Vector{T}             # vector of dt guesses for adaptive time stepping


            function ERBEKrylovStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,m::Int64=128)

                Backend = getfield(Main,Symbol("Backend"))
                Precision = getfield(Main,Symbol("Precision"))

                Momentum = PhaseSpace.Momentum
                Spacetime = PhaseSpace.Spacetime
                x_num = Spacetime.x_num
                y_num = Spacetime.y_num
                z_num = Spacetime.z_num
                px_num_list = Momentum.px_num_list
                py_num_list = Momentum.py_num_list
                pz_num_list = Momentum.pz_num_list
                dE_list = PhaseSpace.Grids.dE_list

                n_space = x_num*y_num*z_num
                n_momentum =PhaseSpace.Grids.n_momentum
                momentum_offset_species = PhaseSpace.Grids.momentum_species_offset

                E = zeros(Backend,Precision,n_momentum)
                Etmp = zeros(Precision,n_momentum)
                for species in eachindex(PhaseSpace.name_list)
                    px_num = px_num_list[species]
                    py_num = py_num_list[species]
                    pz_num = pz_num_list[species]
                    dE = dE_list[species]
                    for px in 1:px_num
                        for py in 1:py_num
                            for pz in 1:pz_num
                                idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                                Etmp[idx] = Precision(dE[px])
                            end
                        end
                    end
                end

                E = Backend === CuArray ? CuArray(Etmp) : Etmp

                @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

                Binary_Interactions = !isempty(BinM.Binary_list) && !isnothing(BinM.Domain)
                Emission_Interactions = !isempty(EmiM.Emission_list)

                Bin_Domain = BinM.Domain

                if Binary_Interactions
                    M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
                    M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
                else
                    M_Bin_Mul_Step = zeros(Backend,Precision,0,0)
                    M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,0)
                end
                df = zeros(Backend,Precision,length(Initial))
                df_Bin = zeros(Backend,Precision,length(Initial))
                df_Emi = zeros(Backend,Precision,length(Initial))
                df_Flux = zeros(Backend,Precision,length(Initial))
                df_tmp = zeros(Backend,Precision,length(Initial))

                fstep = zeros(Backend,Precision,length(Initial))
                F = zeros(Backend,Precision,n_momentum)
                J = zeros(Backend,Precision,n_momentum,n_momentum)
                Jsparse = sparse(zeros(Backend,Precision,n_momentum,n_momentum))
                fold = zeros(Backend,Precision,n_momentum)
                fout = zeros(Backend,Precision,n_momentum)
                fscale = zeros(Backend,Precision,n_momentum)
                δ = zeros(Backend,Precision,n_momentum)

                D = one(Precision) ./ copy(E)
                Dinv = copy(E)

                # Krylov Subspace specific caches 
                # Common 
                V = zeros(Backend,Float64,n_momentum,m+1)
                # Scalar Arnoldi 
                Hs = zeros(Backend,Float64,m+1,m+1)
                Hsexp = zeros(Backend,Float64,m+1,m+1)
                w = zeros(Backend,Float64,n_momentum)
                # Block Arnoldi
                Hb = zeros(Backend,Float64,m+1,m+1)
                Hbexp = zeros(Backend,Float64,m+1,m+1)
                W = zeros(Backend,Float64,n_momentum,4)
                R0 = zeros(Backend,Float64,4,4)
                Hmmp1 = zeros(Backend,Float64,4,4)

                PadeExpWorkspace_s = PadeExpWorkspaceStruct{Float64}(Hsexp)
                PadeExpWorkspace_b = PadeExpWorkspaceStruct{Float64}(Hbexp)
                

                Vol = FluxM.Vol

                f_init = convert(Vector{Precision},Initial)
                M_Bin = Precision.(BinM.M_Bin)
                X_Flux = Precision.(FluxM.X_Flux)
                P_Flux = Precision.(FluxM.P_Flux)
                A_Flux = Precision.(spdiagm(FluxM.Ap_Flux)) # diagonal matrix of Ap flux for Modified Patankar Euler method
                invA_Flux = Precision.(spdiagm(1 ./ FluxM.Ap_Flux)) # invert Ap Flux for time stepping
                f = convert(Vector{Precision},copy(Initial))
                df_Inj = convert(Vector{Precision},copy(Injection))

                # Making invA = vector of diagonal entries of invA_Flux for each spatial point (used for scaling)
                invA = zeros(Precision,n_space)
                for off_space in 0:n_space-1

                    start_idx = n_momentum*off_space + 1

                    invA[off_space+1] = invA_Flux[start_idx,start_idx]

                end


                #= Making invMP = (I-dt*A^{-1}(M_Emi-P_Flux))^{-1}
                  This assuming only emissive interactions coming from one set of specie to another (e.g. electron to photon but not photon to electron) (I-A^{-1}(M_Emi-P_Flux)) has a block triangular structure:
                   I-MP = [A 0]
                          [B C]
                   where A corresponds to the non-emissive species and C corresponds to the emissive species. This means we can invert (I-MP) as:
                   inv(I-MP) = [A^{-1}          0    ]
                               [-C^{-1}BA^{-1} C^{-1}]    

                =#
                ImMP = I - (dt_initial/2)*invA_Flux*(#=M_Emi=# - P_Flux)
                #invImMP = spzeros(Precision,size(P_Flux))
                invImMP_rows = Int64[]
                invImMP_cols = Int64[]
                invImMP_vals = Precision[]
                momentum_offset = [momentum_offset_species ; n_momentum]
                for space in 0:n_space-1
                    off_space = space*n_momentum
                    # diagonal blocks
                    #= for block diagonals Bii = ii component of the inverse matrix 
                          Bii = Aii^{-1} where Aii is the ii component of (I-MP) = I-dt*A^{-1}(M_Emi-P_Flux)
                    =#
                    for speciesi in eachindex(PhaseSpace.name_list)
                        pi_low = momentum_offset[speciesi] + off_space + 1
                        pi_up = momentum_offset[speciesi+1] + off_space

                        #invImMP_view = @view(invImMP[pi_low:pi_up,pi_low:pi_up])
                        ImMP_view = @view(ImMP[pi_low:pi_up,pi_low:pi_up])

                        #invImMP_view .= sparse(inv(ImMP_view))
                        rows, cols, vals = findnz(sparse(inv(ImMP_view)))
                        append!(invImMP_rows, rows .+ (pi_low - 1))
                        append!(invImMP_cols, cols .+ (pi_low - 1))
                        append!(invImMP_vals, vals)
                    end
                end
                invImMP = sparse(invImMP_rows, invImMP_cols, invImMP_vals, size(P_Flux,1), size(P_Flux,2))

                if Backend isa CUDABackend
                    f_init = CuArray(f_init)
                    if M_Bin isa AbstractSparseArray
                        M_Bin = CuSparseMatrixCSC(M_Bin)
                    else
                        M_Bin = CuArray(M_Bin)
                    end
                    X_Flux = CuSparseMatrixCSC(X_Flux)
                    P_Flux = CuSparseMatrixCSC(P_Flux)
                    A_Flux = CuSparseMatrixCSC(A_Flux)
                    invA_Flux = CuSparseMatrixCSC(invA_Flux)
                    f = CuArray(f)
                    df_Inj = CuArray(df_Inj)
                    invImMP = CuSparseMatrixCSC(invImMP)
                    D = cu(D)
                    Dinv = cu(Dinv)
                end

                # Build new MEmi
                if Backend isa CPUBackend
                    M_Emi = Vector{Union{Matrix{Precision},SparseMatrixCSC{Precision,Int64}}}(undef,n_space)
                    for off_space in 1:n_space
                        if isassigned(EmiM.M_Emi,off_space)
                            M_Emi[off_space] = Precision.(EmiM.M_Emi[off_space])
                        end
                    end
                elseif Backend isa CUDABackend
                    M_Emi = Vector{Union{CuMatrix{Precision},CuSparseMatrixCSC{Precision,Int32}}}(undef,n_space)
                    for off_space in 1:n_space
                        if isassigned(EmiM.M_Emi,off_space)
                            M_Emi[off_space] = cu(EmiM.M_Emi[off_space])
                        end
                    end
                end

                if !isnothing(DistributionDomainMask)
                    f_mask = ones(Precision,length(Initial))
                    for off_space_idx in DistributionDomainMask
                        for species_idx in eachindex(PhaseSpace.name_list)
                        LocationSpeciesToStateVector(f_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= Precision(0.0)
                        end
                    end
                    if Backend isa CUDABackend
                        f_mask = CuArray(f_mask)
                    end
                else
                    f_mask = nothing
                end

                if !isnothing(DeltaDistributionDomainMask)
                    df_mask = ones(Precision,length(Initial))
                    for off_space_idx in DeltaDistributionDomainMask
                        for species_idx in eachindex(PhaseSpace.name_list)
                        LocationSpeciesToStateVector(df_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= Precision(0.0)
                        end
                    end
                    if Backend isa CUDABackend
                        df_mask = CuArray(df_mask)
                    end
                else
                    df_mask = nothing
                end

                # cut initial values that are smaller than n_cut 
                    @. f_init = ifelse(f_init<=n_cut,zero(eltype(f_init)),f_init)
                    @. f = ifelse(f<=n_cut,zero(eltype(f)),f)

                dt_guess = zeros(Precision, n_space)
                fill!(dt_guess,ldexp(dt_initial, -4)) # initial guess for dt is 1/64 of the initial dt

                ###### Actually Build the Struct with Concrete Types ######

                self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(X_Flux),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

                self.PhaseSpace = PhaseSpace
                self.Implicit = true
                self.Precision = Precision
                self.Adaptive = Adaptive
                self.dt0 = PhaseSpace.Spacetime.dt0
                self.Cr = zero(Precision)
                self.n_cut = n_cut
                self.step = 0
                self.Binary_Interactions = Binary_Interactions
                self.Emission_Interactions = Emission_Interactions
                self.DistributionDomainMask = DistributionDomainMask
                self.DeltaDistributionDomainMask = DeltaDistributionDomainMask
                if isnothing(DistributionDomainMask)
                    self.ActiveDomain = InclusiveDomainMask(PhaseSpace)
                else
                    self.ActiveDomain = setdiff(InclusiveDomainMask(PhaseSpace),DistributionDomainMask)
                end
                self.PhaseSpace = PhaseSpace
                self.Bin_Domain = BinM.Domain
                self.f_init = f_init
                self.M_Bin = M_Bin
                self.M_Emi = M_Emi
                self.X_Flux = X_Flux
                self.P_Flux = P_Flux
                self.A_Flux = A_Flux
                self.invA_Flux = invA_Flux
                self.invA = invA
                self.Vol = Vol
                self.f = f
                self.fstep = fstep
                self.df_Inj = df_Inj
                self.M_Bin_Mul_Step = M_Bin_Mul_Step
                self.M_Bin_Mul_Step_reshape = M_Bin_Mul_Step_reshape
                self.df = df
                self.df_Bin = df_Bin
                self.df_Emi = df_Emi
                self.df_Flux = df_Flux
                self.df_tmp = df_tmp
                self.f_mask = f_mask
                self.df_mask = df_mask

                self.F = F
                self.J = J
                self.Jsparse = Jsparse
                self.D = D 
                self.Dinv = Dinv 
                self.fold = fold
                self.fout = fout
                self.fscale = fscale
                self.δ = δ
                self.m = m

                self.E = E

                self.invImMP = invImMP

                self.V = V
                self.Hs = Hs
                self.Hsexp = Hsexp
                self.w = w
                self.Hb = Hb
                self.Hbexp = Hbexp
                self.W = W
                self.R0 = R0
                self.Hmmp1 = Hmmp1

                self.PadeExpWorkspace_s = PadeExpWorkspace_s
                self.PadeExpWorkspace_b = PadeExpWorkspace_b

                self.dt_guess = dt_guess

                return self
            end

    end


"""
    ERBEKrylov(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. 

"""
function (method::ERBEKrylovStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

    # will we reached the next t_save?

    t_next = t_start + dt
    if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
        save = true
    elseif t_next >= t_stop
        adaptive_factor *= (t_stop - t_start) / dt # adjust adaptive factor for final time step to ensure we end exactly at t_stop 
        dt = t_stop - t_start
        save = true
    else
        save = false
    end

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = similar(method.f)

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # Half injection 

        #@. method.fstep += method.df_Inj * dt_scale / 2

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end 

    # binary update

        update_momentum!(method,dt_scale)

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end

    # half injection 

        #@. method.fstep += method.df_Inj * dt_scale / 2

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end


function update_momentum!(method::ERBEKrylovStruct,dt::T) where T

    Precision = method.Precision
    
    n_momentum = method.PhaseSpace.Grids.n_momentum
    #n_space = method.PhaseSpace.Grids.n_space
    momentum_species_offset = method.PhaseSpace.Grids.momentum_species_offset
    name_list = method.PhaseSpace.name_list
    num_species = length(name_list)
    
    fold = method.fold
    fout = method.fout
    J = method.J
    Jsparse = sparse(copy(J))
    #Jtmp = copy(J)
    F = method.F
    fscale = method.fscale #::Vector{Precision} = zeros(Precision,n_momentum)

    D = method.D #Diagonal(ones(Precision,n_momentum)) # 1/E
    Dinv = method.Dinv #Diagonal(ones(Precision,n_momentum)) # E
    E = method.E

    δ = method.δ #::Vector{Precision} = zeros(Precision,n_momentum)

    m = method.m # dimension of Krylov subspace for exponential action approximation, can adjust based on problem size and desired accuracy


    EmiTrue::Bool = true

    for off_space in method.ActiveDomain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        df_Inj = @view(method.df_Inj[start_idx:end_idx])

        has_injection = sum(df_Inj) > zero(Precision)

        if !has_injection && sum(fstep) == zero(Precision) 
            continue
        end

        @inbounds vol = method.Vol[off_space+1]
        @inbounds invA = method.invA[off_space+1]
        @inbounds dt_guess = method.dt_guess[off_space+1]

        if method.Binary_Interactions && off_space in method.Bin_Domain

            if isassigned(method.M_Emi, off_space+1)
                EmiTrue = true
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0 # t is ratio of current t to dt. i.e. t=1.0 is end of timestep

            dt_local = dt_guess # initial guess for local time step, based on previous time steps
            dt_old = dt_local
            dt_next = dt_local

            # k's scale dt_local
            kE = 1.0
            kH = 1.0
            k = 1.0

            while t < 1.0

                dtscale = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision))
                    # Form J
                    @. J = Precision(2) * method.M_Bin_Mul_Step
                    #fill!(J,zero(Precision))
                    #fill!(method.M_Bin_Mul_Step,zero(Precision)) 
                    if EmiTrue
                        @. J += M_Emi
                        @. method.M_Bin_Mul_Step += M_Emi
                    end

                    @. J *= dtscale * invA

                    #println("max M_Emi: ", maximum(M_Emi), " min M_Emi : ", minimum(M_Emi))
                    #println("max J: ", maximum(J), " min J : ", minimum(J))
                    #println("max fold: ", maximum(fold), " min fold : ", minimum(fold))
                    #println("max df_Inj: ", maximum(df_Inj), " min df_Inj : ", minimum(df_Inj))

                    # Form F
                    mul!(F,method.M_Bin_Mul_Step,fold)
                    @. F *= invA

                    # Energy error estimate using F 
                    #ηEest = abs(dot(E,F) * dtscale) / dot(E,fold)
                    #println("Energy error estimate: ", ηEest)

                    @. F += #=A *=# df_Inj
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end
                    lmul!(dtscale,F)

                    # Diagonal scaling of J and F to improve conditioning for Krylov subspace approximation
                    J .*= (D)' # right mul
                    J .*= (Dinv) # left mul
                    F .*= (Dinv) # left mul

                    # TODO: define J and F as views of Aaug to avoid copying data
                    #Aaug = method.Aaug
                    #fill!(Aaug,zero(Precision))
                    #@view(Aaug[1:n_momentum,1:n_momentum]) .= J
                    #@view(Aaug[1:n_momentum,n_momentum+1]) .= F

                    #println("sum(F): $(sum(F)), max(F): $(maximum(F)), min(F): $(minimum(F))")
                    #println("sum(J): $(sum(J)), max(J): $(maximum(J)), min(J): $(minimum(J))")

                    # scalar
                        _, hmp1m, β = arnoldi_scalar!(J,F,method.V,method.Hs,method.w;reorth=false) # TODO: reoth doesn't work
                    # block
                        #Rlast = arnoldi_block4!(J,F,method.V,method.Hb,method.W,method.R0;reorth=false)

                        #display(method.Hb)
                        #display(method.V)
                        #error("Stopping for debugging new, dtscale = $dtscale, $hmp1m")

                    order = 1.0 # energy error is order 1
                    ηHtarget = 1e-16
                    ηEtarget = Precision == Float32 ? 1e-4 : 1e-6
                    k = min(1.0,(1.0-t)*dt/dt_local) # initial guess for k, limited to not overshoot final time 
                    k_old = 1.0
                    k_next = 1.0
                    q_init = m/4
                    q = q_init
                    ηH_old = 1.0
                    γ = 0.6
                    substeps = 0
                    ηE = Inf 
                    ηH = Inf
                    non_finite = false

                    while (!isfinite(ηE) || !isfinite(ηH)) || (ηE > ηEtarget || ηH > 1.4ηHtarget) 

                        # scalar
                            expm_pade!(k, method.Hsexp, method.Hs, method.PadeExpWorkspace_s)
                            #display(method.Hs)
                            #display(method.Hsexp)
                            #display(exp(method.Hs))
                            ferr = CUDA.@allowscalar method.Hsexp[m,m+1]
                            #println("dtscale: ", dtscale, " k: ", k, " ferr: ", ferr)
                            ηH = abs(k * dtscale * β * hmp1m * ferr) # error estimate from Krylov subspace approximation of exponential action
                            @views mul!(δ, method.V[1:n_momentum, 1:m], method.Hsexp[1:m, m+1],β*k,zero(Precision)) # δ is energy change
                            δ .*= D # scale δ back to original space
                        

                        # block
                        #=    expm_pade!(k, method.Hbexp, method.Hb, method.PadeExpWorkspace_s)
                            ferr = @view(method.Hbexp[m-3:m,m+1])
                            ηH = k * norm(Rlast * ferr)
                            @views mul!(δ, method.V[1:n_momentum, 1:m], method.Hbexp[1:m, m+1],one(Precision),zero(Precision)) # δ is energy change
                            δ .*= D # scale δ back to original space

                            #display(method.Hb)
                            #display(method.Hbexp)
                        =#

                        @. fout = fold + δ
                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))

                        #println("max delta: $(maximum(δ)), min delta: $(minimum(δ)), sum delta: $(sum(δ))")
                        #println("max fold: $(maximum(fold)), min fold: $(minimum(fold)), sum fold: $(sum(fold))")
                        #println("max fout: $(maximum(fout)), min fout: $(minimum(fout)), sum fout: $(sum(fout))")

                        #error("end for debugging, dtscale = $dtscale, errest = $ηH, eng err: $ηE")

                        # KIOPS adaptive time 
                        #=if substeps == 0 || substeps == 1
                            q = q_init
                            q_old = q 
                        else
                            q_old = q
                            q = log(k/k_old) / log(ηH/ηH_old)
                        end=#

                        kE = k*(ηEtarget/(ηE+eps(ηEtarget)))^(1.0/(order+1)) # k from energy error estimate 
                        kH = k*(k*γ/((ηH+eps(ηHtarget))/ηHtarget))^(1/q) #(0.77 * 0.5(tanh(log10(ηH/ηHtarget)+1)+1)) # TODO: update this with a more specific method

                        k_old = k

                        if !isfinite(ηE) || !isfinite(ηH)
                            @warn "Rejecting step due to non-finite error estimate: ηE: $ηE, ηH: $ηH reducing time step by 0.5"
                            k *= 0.5
                            non_finite = true
                            break
                        elseif ηE < ηEtarget && ηH < 1.4ηHtarget # good step
                            # define how large to guess the next step can be
                            k_next = min(sqrt(kE),sqrt(kH),1.5) # max k is 2.0 or 1.0-t to avoid overshooting final time
                            #println("Accepting step with k: ", k, " k_next: ", k_next, " t: ", t, " dt_local: ", dt_local," ηE: ", ηE, " ηH: ", ηH, " dtscale: ", dtscale)
                            break
                        else # 
                            k = min(0.9*kE,kH)
                            k = max(k,0.5*k_old)
                            #println("Rejecting step with k: ", k, " t: ", t, " dt_local: ", dt_local," ηE: ", ηE, " ηH: ", ηH)
                        end

                        substeps += 1

                    end

                    if non_finite
                        dt_local = 0.5*dt_local
                        continue # skip to next spatial point if non-finite error estimate
                    end

                    t += k*dt_local / dt
                    println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", k*dt_local)
                    dt_old = k*dt_local # used for updating dt_guess
                    dt_local = k_next*dt_local

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

            @inbounds method.dt_guess[off_space+1] = dt_old # update dt_guess with last plus one dt (avoid last as this could be limited by (1.0-t))

        else #

            if isassigned(method.M_Emi, off_space+1)
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = Precision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    if EmiTrue
                        copyto!(Jsparse,M_Emi)
                        Jsparse *= dtscale * invA
                    else
                        fill!(Jsparse,zero(Precision))
                        dropzeros!(Jsparse)
                    end

                    # Form F
                    mul!(F,Jsparse,fold)
                    @. F += #=A *=# df_Inj * dtscale
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end

                    Jsparse .*= D' # right mul
                    Jsparse .*= Dinv # left mul
                    F .*= Dinv # left mul

                    #@. J64 = Float64(J)

                    if has_injection 
                        #@. F64 = Float64(F)
                        #arnoldi!(Ks64,J64,F64;m=m,reorthogonalize=true)
                        arnoldi!(Ks,Jsparse,F;m=m,reorthogonalize=true)
                    else # no injection, just linear Jacobian so use exp over phi
                        #f64 = Float64.(Dinv * fold)
                        #arnoldi!(Ks64,J64,f64;m=m,reorthogonalize=true)
                        mul!(fscale,Dinv,fold)
                        arnoldi!(Ks,Jsparse,fscale;m=m,reorthogonalize=true)
                    end

                    V = ExponentialUtilities.getV(Ks)[:,1:end-1]
                    H = ExponentialUtilities.getH(Ks)[1:end-1,1:end]
                    #V = ExponentialUtilities.getV(Ks64)[:,1:end-1]
                    #H = ExponentialUtilities.getH(Ks64)[1:end-1,1:end]

                    norm_unscaled = norm(V' * V - I)
                    if isnan(norm_unscaled) # no J or F means no Krylov subspace generated, so just accept step as is and move on,
                        break
                    end

                    if #=cond((I - dt_local*H)) > 1f3 ||=#  norm(V' * V - I) > 1f-5 
                        dt_local = dt_local * 0.5
                        @warn "Krylov subspace has poor orthogonalisation $(norm(V' * V - I)), reducing time step to $dt_local"
                        continue
                    end

                    # Compute φ functions of H
                    if has_injection 
                        phiv!(ϕ,k,Ks,1;cache=ϕcache,correct=true,errest=false) # TODO: This allocates
                        @. δ = D * @view(ϕ[:,2]) * k
                        @. fout = fold + δ
                    else # no injection, just linear Jacobian so use exp over phi
                        expv!(δ,k#=*invA=#,Ks)
                        @. fout = D * δ
                    end


                    # adaptive time stepping

                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        if has_injection
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * dtscale)
                            Eold = dot(E, fold .+ df_Inj * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        #println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        if isinf(ηE)
                            @warn "Energy error is Inf or NaN, may be unstable, consider reducing time step or adjusting ηtarget"
                            dt_local *= convert(typeof(dt_local), 1/64)
                            continue
                        end

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        ηtarget = 1e-16

                        kE = (1e-4/(ηE+eps(1e-4)))^(1.0/(order+1)) # k from energy error estimate
                        if has_injection
                            kϕmax = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                            kϕ = min(kE,kϕmax) # max k is 2.0
                            errest = Inf
                            while errest > ηtarget # k from ϕv errestimate, can be more strict than energy error estimate
                                _, errest = phiv!(ϕ,kϕ,Ks,1;cache=ϕcache,correct=true,errest=true) # TODO: This allocates
                                println("ϕv error estimate during: ", errest)
                                if errest > ηtarget
                                    kϕ *= 0.77 * 0.5(tanh(log10(errest)-log10(ηtarget)+1)+1)
                                end
                            end
                            println("ϕv error estimate after: ", errest)
                            k = min(kE,kϕ,2.0)
                            println("kE: ", kE, " kϕ: ", kϕ, " errest: ", errest)
                        else
                            k = min(kE,2.0)
                        end
                        #k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        dt_local = k*dt_old
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                    if k != 1.0
                        if has_injection 
                            phiv!(ϕ,k,Ks,1;cache=ϕcache,correct=true)
                            mul!(δ,D,@view(ϕ[:,2]),k#=*invA=#,zero(Precision))
                            @. fout = fold + δ
                        else # no injection, just linear Jacobian so use exp over phi
                            expv!(δ,k#=*invA=#,Ks)
                            @. fout = D.diag * δ
                        end

                        # energy after correction  
                        if has_injection 
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * k * dtscale)
                            Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error after adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        #=if ηE > 2e-4
                            @warn "Energy error in EXPRB1 step is large $ηE, may be unstable, consider reducing time step or adjusting ηtarget"
                        end=#
                        
                    end

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        end

    end

    return nothing

end