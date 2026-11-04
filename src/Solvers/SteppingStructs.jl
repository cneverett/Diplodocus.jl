abstract type AbstractSteppingMethod <: Function end
abstract type ImplicitSteppingMethod <: AbstractSteppingMethod end
abstract type ExplicitSteppingMethod <: AbstractSteppingMethod end

##### FORWARD EULER ######
    mutable struct ForwardEulerStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},MET<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ExplicitSteppingMethod

        PhaseSpace::PhaseSpaceStruct

        Binary_Interactions::Bool
        Emission_Interactions::Bool

        M_Bin::MBT
        Bin_Domain::BD

        M_Emi::MET

        F_Flux::SMT
        invAp_Flux::VT                  # inv Ap flux for time stepping
        Vol::Vector{T}

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
        df::VT                         # change in distribution function
        df_Bin::VT                     # change in distribution function due to binary interactions
        df_Emi::VT                     # change in distribution function due to emission interactions
        df_Flux::VT                    # change in distribution function due to fluxes
        df_Inj::VT                     # change in distribution function due to injection of particles
        df_tmp::VT                     # temporary array the size of f for CFL calculations
        f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
        df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

        function ForwardEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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

            n_space = x_num*y_num*z_num
            n_momentum = 0
            for i in eachindex(px_num_list)
                n_momentum += px_num_list[i]*py_num_list[i]*pz_num_list[i]
            end

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

            Vol = FluxM.Vol

            f_init = convert(Vector{Precision},Initial)
            M_Bin = BinM.M_Bin
            M_Emi = EmiM.M_Emi
            F_Flux = FluxM.X_Flux + FluxM.P_Flux # sum of space and momentum fluxes
            invAp_Flux = 1 ./ FluxM.Ap_Flux # invert Ap Flux
            f = convert(Vector{Precision},copy(Initial))
            df_Inj = convert(Vector{Precision},copy(Injection))
            if Backend isa CUDABackend
                f_init = CuArray(f_init)
                if M_Bin isa AbstractSparseArray
                    M_Bin = CuSparseMatrixCSC(M_Bin)
                else
                    M_Bin = CuArray(M_Bin)
                end
                if M_Emi isa AbstractSparseArray
                    M_Emi = CuSparseMatrixCSC(M_Emi)
                else
                    M_Emi = CuArray(M_Emi)
                end
                F_Flux = CuSparseMatrixCSC(F_Flux) # sum of space and momentum fluxes
                invAp_Flux = CuArray(invAp_Flux)
                f = CuArray(f)
                df_Inj = CuArray(df_Inj)
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

            ###### Actually Build the Struct with Concrete Types ######

            self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(M_Emi),typeof(F_Flux),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

            self.Implicit = false
            self.Adaptive = Adaptive
            self.dt0 = PhaseSpace.Spacetime.dt0
            self.Cr = zero(Precision)
            self.n_cut = n_cut
            self.step = 0
            self.Binary_Interactions = Binary_Interactions
            self.Emission_Interactions = Emission_Interactions
            self.PhaseSpace = PhaseSpace
            self.Bin_Domain = BinM.Domain
            self.f_init = f_init
            self.M_Bin = M_Bin
            self.M_Emi = M_Emi
            self.F_Flux = F_Flux
            self.invAp_Flux = invAp_Flux
            self.Vol = Vol
            self.f = f
            self.df_Inj = df_Inj
            self.M_Bin_Mul_Step = M_Bin_Mul_Step
            self.M_Bin_Mul_Step_reshape = M_Bin_Mul_Step_reshape # Thanks to Emma Godden for fixing a bug here
            self.df = df
            self.df_Bin = df_Bin
            self.df_Emi = df_Emi
            self.df_Flux = df_Flux
            self.df_tmp = df_tmp
            self.f_mask = f_mask
            self.df_mask = df_mask

            return self
        end

    end

##### FORWARD SYMPLECTIC EULER ######
    mutable struct ForwardSymplecticEulerStruct{T<:AbstractFloat} <: ExplicitSteppingMethod

        PhaseSpace::PhaseSpaceStruct

        Binary_Interactions::Bool
        Emission_Interactions::Bool

        M_Bin::AbstractMatrix{T}
        Bin_Domain::Union{Vector{Int64},Nothing}

        M_Emi::AbstractMatrix{T}

        X_Flux::AbstractSparseArray{T,<:Integer,2}
        P_Flux::AbstractSparseArray{T,<:Integer,2}
        invAp_Flux::AbstractVector{T}                   # inv Ap flux for time stepping
        Vol::AbstractVector{T}

        Adaptive::Bool
        Implicit::Bool
        dt0::T
        p_cut::T

        M_Bin_Mul_Step::AbstractMatrix{T}               # temporary array for matrix multiplication of binary terms
        M_Bin_Mul_Step_reshape::AbstractVector{T}       # temporary array for reshaped matrix multiplication of binary terms
        f_init::AbstractVector{T}                       # initial distribution function (used by solver to define output struct)
        f::AbstractVector{T}                            # current distribution function
        df::AbstractVector{T}                           # change in distribution function
        df_Momentum::AbstractVector{T}                  # change in distribution function due to momentum fluxes 
        df_Space::AbstractVector{T}                     # change in distribution function due to spatial fluxes
        df_Bin::AbstractVector{T}                       # change in distribution function due to binary interactions
        df_Emi::AbstractVector{T}                       # change in distribution function due to emission interactions
        df_XFlux::AbstractVector{T}                     # change in distribution function due to spatial fluxes
        f_Space::AbstractVector{T}                      # temporary array for f after spatial fluxes
        df_PFlux::AbstractVector{T}                     # change in distribution function due to momentum fluxes
        f_Momentum::AbstractVector{T}                      # temporary array for f after momentum fluxes
        df_Inj::AbstractVector{T}                       # change in distribution function due to injection of particles
        df_tmp::AbstractVector{T}                       # temporary array the size of f for CFL calculations 
        f_tmp::AbstractVector{T}                        # temporary array the size of f for symplectic Euler calculations  

        function ForwardSymplecticEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,p_cut::Float64=1e-45)

            Backend = getfield(Main,Symbol("Backend"))
            Precision = getfield(Main,Symbol("Precision"))

            @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

            self = new{Precision}()

            self.Adaptive = Adaptive
            self.Implicit = false
            self.dt0 = convert(Precision,PhaseSpace.Spacetime.dt0)
            self.p_cut = convert(Precision,p_cut)

            self.Binary_Interactions = !isempty(BinM.Binary_list)
            self.Emission_Interactions = !isempty(EmiM.Emission_list)

            self.PhaseSpace = PhaseSpace

            self.Bin_Domain = BinM.Domain

            invAp_Flux = 1 ./ FluxM.Ap_Flux # invert Ap Flux

            if Backend isa CPUBackend
                self.f_init = convert(Vector{Precision},Initial)
                self.M_Bin = BinM.M_Bin
                self.M_Emi = EmiM.M_Emi
                self.X_Flux = FluxM.X_Flux
                self.P_Flux = FluxM.P_Flux
                self.invAp_Flux = invAp_Flux
                self.Vol = FluxM.Vol
                self.f = convert(Vector{Precision},Initial)
                self.df_Inj = convert(Vector{Precision},Injection)
            elseif Backend isa CUDABackend
                self.f_init = CuArray(Initial)
                if BinM.M_Bin isa AbstractSparseMatrix
                    self.M_Bin = CuSparseMatrixCSC(BinM.M_Bin)
                else
                    self.M_Bin = CuArray(BinM.M_Bin)
                end
                if EmiM.M_Emi isa AbstractSparseMatrix
                    self.M_Emi = CuSparseMatrixCSC(EmiM.M_Emi)
                else
                    self.M_Emi = CuArray(EmiM.M_Emi)
                end
                self.X_Flux = CuSparseMatrixCSC(FluxM.X_Flux)
                self.P_Flux = CuSparseMatrixCSC(FluxM.P_Flux)
                self.invAp_Flux = CuArray(invAp_Flux)
                self.Vol = CuArray(FluxM.Vol)
                self.f = CuArray(Initial)
                self.df_Inj = CuArray(Injection)
            else
                error("Backend type not recognized.")
            end  

            Momentum = PhaseSpace.Momentum
            Spacetime = PhaseSpace.Spacetime
            x_num = Spacetime.x_num
            y_num = Spacetime.y_num
            z_num = Spacetime.z_num
            px_num_list = Momentum.px_num_list
            py_num_list = Momentum.py_num_list
            pz_num_list = Momentum.pz_num_list

            n_space = x_num*y_num*z_num
            n_momentum = sum(sum(px_num_list.*py_num_list.*pz_num_list))

            if !isempty(BinM.Binary_list)
                self.M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
                self.M_Bin_Mul_Step_reshape = reshape(self.M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
            end
            self.df = zeros(Backend,Precision,length(Initial))
            self.df_Momentum = zeros(Backend,Precision,length(Initial))
            self.df_Space = zeros(Backend,Precision,length(Initial))
            self.df_Bin = zeros(Backend,Precision,length(Initial))
            self.df_Emi = zeros(Backend,Precision,length(Initial))
            self.df_XFlux = zeros(Backend,Precision,length(Initial))
            self.df_PFlux = zeros(Backend,Precision,length(Initial))
            self.f_Space = zeros(Backend,Precision,length(Initial))
            self.f_Momentum = zeros(Backend,Precision,length(Initial))
            self.df_tmp = zeros(Backend,Precision,length(Initial))
            self.f_tmp = zeros(Backend,Precision,length(Initial))

            return self
        end

    end



##### Exponential Rosenbrock Euler ######

    mutable struct ExponentialRosenbrockEulerKrylovStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,Int32,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

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
            MEmiPFlux::Vector{Union{MT,SMT}} # M_Emi - P_Flux for each spatial point  

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
            Jtmp::MT                       # temporary Jacobian matrix for implicit solve
            Jsparse::SMT                   # sparse Jacobian matrix for implicit solve
            D::VT             # Diagonal scaling matrix
            Dinv::VT          # Inverse of diagonal scaling matrix
            ϕ::MT                          # matrix of ϕ functions for exponential Rosenbrock method

            fold::VT                        # local distribution function from previous step 
            fout::VT                        # local distribution function from current step
            fscale::VT                      # scaling vector for exponential Rosenbrock method
            δ::VT                          # temporary vector for exponential Rosenbrock method

            KsB#::KrylovSubspace{Float64,Float64,AbstractMatrix{Float64}}      # Krylov subspace for exponential Rosenbrock method (higher precision for more accuracy)
            KsL
            mB::Int64                        # dimension of Krylov subspace for binary
            mL::Int64                        # dimension of Krylov subspace for linear
            ϕcacheB::ExponentialUtilities.PhivCache{useview,T} where useview # cache for ϕ functions for binary 
            ϕcacheL::ExponentialUtilities.PhivCache{useview,T} where useview # cache for ϕ functions for linear

            E::VT                           # energy vector for correcting step (length of n_momentum) 
            E_long::VT                      # energy vector for correcting step (length of n_momentum*n_space)

            dt_guess::Vector{T}             # guess for dt based on CFL condition

            function ExponentialRosenbrockEulerKrylovStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,mB::Int64=192,mL::Int64=64)

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

                E_long_tmp = zeros(Backend,Precision,n_momentum*n_space)
                for space in 1:n_space
                    @view(E_long_tmp[(space-1)*n_momentum+1:space*n_momentum]) .= Etmp
                end

                E = Backend === CuArray ? CuArray(Etmp) : Etmp
                E_long = Backend === CuArray ? CuArray(E_long_tmp) : E_long_tmp

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
                ϕ = zeros(Backend,Precision,n_momentum,2)
                ϕcacheB = ExponentialUtilities.PhivCache(ϕ,mB,1)
                ϕcacheL = ExponentialUtilities.PhivCache(ϕ,mL,1)

                D = one(Precision) ./ copy(E)
                Dinv = copy(E)

                # use higher precision for Krylov subspace to avoid numerical issues
                if Backend isa CUDABackend
                    KsB = KrylovSubspace{Float64,Float64,CuArray{Float64,2}}(n_momentum,mB)
                    KsL = KrylovSubspace{Float64,Float64,CuArray{Float64,2}}(n_momentum,mL)
                else
                    KsB = KrylovSubspace{Float64,Float64,Array{Float64,2}}(n_momentum,mB)
                    KsL = KrylovSubspace{Float64,Float64,Array{Float64,2}}(n_momentum,mL)
                end

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
                invImMP_rows = Int32[]
                invImMP_cols = Int32[]
                invImMP_vals = Precision[]
                momentum_offset = [momentum_offset_species ; n_momentum]
                #=for space in 0:n_space-1
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
                        append!(invImMP_rows, Int32.(rows .+ (pi_low - 1)))
                        append!(invImMP_cols, Int32.(cols .+ (pi_low - 1)))
                        append!(invImMP_vals, Precision.(vals))
                    end
                    # off-diagonal blocks
                    #=
                        for i>j and Bij = the ij components of the inverse matrix

                            Bij = -Bii \sum_{k=j}^{i-1} Aik Bkj

                    =#
                    #=for speciesi in eachindex(PhaseSpace.name_list)
                        pi_low = momentum_offset[speciesi] + off_space + 1
                        pi_up = momentum_offset[speciesi+1] + off_space


                        for speciesj in 1:speciesi-1
                            pj_low = momentum_offset[speciesj] + off_space + 1
                            pj_up = momentum_offset[speciesj+1] + off_space

                            
                            #Bij = @view(invImMP[pi_low:pi_up,pj_low:pj_up])

                            Bij = zeros(Precision, pi_up - pi_low + 1, pj_up - pj_low + 1)

                            for speciesk in speciesj:speciesi-1
                                pk_low = momentum_offset[speciesk] + off_space + 1
                                pk_up = momentum_offset[speciesk+1] + off_space

                                Aik = @view(ImMP[pi_low:pi_up,pk_low:pk_up])
                                Bkj = @view(invImMP[pk_low:pk_up,pj_low:pj_up])

                                #Bij .-= sparse(Aik * Bkj) 
                                Bij .-= Aik * Bkj
                            end

                            # TODO: check order of rows and cols here, might need to swap them 
                            rows, cols, vals = findnz(sparse(inv(Bij)))
                            append!(invImMP_rows, rows .+ pi_low - 1)
                            append!(invImMP_cols, cols .+ pj_low - 1)
                            append!(invImMP_vals, vals)
                        end

 
                    end=#
                end=#
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

                # Build new MEmi and MEmiPFlux
                if Backend isa CPUBackend
                    M_Emi = Vector{Union{Matrix{Precision},SparseMatrixCSC{Precision,Int32}}}(undef,n_space)
                    MEmiPFlux = Vector{Union{Matrix{Precision},SparseMatrixCSC{Precision,Int32}}}(undef,n_space)
                    for off_space in 0:n_space-1
                        start_idx = n_momentum*off_space+1
                        end_idx = n_momentum*(off_space+1)
                        in_Binary = Binary_Interactions && in(off_space,BinM.Domain)
                        if in_Binary 
                            MEmiPFlux[off_space+1] = zeros(Precision,n_momentum,n_momentum)
                            MEmiPFlux[off_space+1] .-= @view(P_Flux[start_idx:end_idx,start_idx:end_idx])
                        else
                            MEmiPFlux[off_space+1] = spzeros(Precision,Int32,n_momentum,n_momentum)
                            MEmiPFlux[off_space+1] .-= @view(P_Flux[start_idx:end_idx,start_idx:end_idx])
                        end

                        if isassigned(EmiM.M_Emi,off_space+1) # already defined as sparse or dense 
                            M_Emi[off_space+1] = Precision.(EmiM.M_Emi[off_space+1])
                            MEmiPFlux[off_space+1] .+= Precision.(EmiM.M_Emi[off_space+1])
                        end

                        if !in_Binary 
                            dropzeros!(MEmiPFlux[off_space+1])
                        end
                    end
                elseif Backend isa CUDABackend
                    M_Emi = Vector{Union{CuMatrix{Precision},CuSparseMatrixCSC{Precision,Int32}}}(undef,n_space)
                    MEmiPFlux = Vector{Union{CuMatrix{Precision},CuSparseMatrixCSC{Precision,Int32}}}(undef,n_space)
                    for off_space in 0:n_space-1
                        start_idx = n_momentum*off_space+1
                        end_idx = n_momentum*(off_space+1)
                        in_Binary = Binary_Interactions && in(off_space,BinM.Domain)
                        if in_Binary 
                            MEmiPFlux[off_space+1] = zeros(CUDABackend(),Precision,n_momentum,n_momentum)
                            MEmiPFlux[off_space+1] .-= @view(P_Flux[start_idx:end_idx,start_idx:end_idx])
                        else
                            MEmiPFlux[off_space+1] = CuSparseMatrixCSC(spzeros(Precision,Int32,n_momentum,n_momentum))
                            MEmiPFlux[off_space+1] .-= @view(P_Flux[start_idx:end_idx,start_idx:end_idx])
                        end
                        if isassigned(EmiM.M_Emi,off_space+1)
                            M_Emi[off_space+1] = cu(EmiM.M_Emi[off_space+1])
                            MEmiPFlux[off_space+1] .+= Precision.(EmiM.M_Emi[off_space+1])
                        end

                        if !in_Binary 
                            dropzeros!(MEmiPFlux[off_space+1])
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

                # cut initial values that are smaller than n_cut in both number and energy 
                @. f_init = ifelse(f_init<=n_cut && f_init * E_long < n_cut,zero(eltype(f_init)),f_init)
                @. f = ifelse(f<=n_cut && f * E_long < n_cut,zero(eltype(f)),f)

                dt_guess = zeros(Precision, n_space)
                fill!(dt_guess,ldexp(dt_initial, -5)) # initial guess for dt is 1/32 of the initial dt

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
                self.Jtmp = copy(J)
                self.Jsparse = Jsparse
                self.D = D 
                self.Dinv = Dinv 
                self.fold = fold
                self.fout = fout
                self.fscale = fscale
                self.δ = δ
                self.KsB = KsB
                self.KsL = KsL
                self.ϕ = ϕ
                self.ϕcacheB = ϕcacheB
                self.ϕcacheL = ϕcacheL
                self.mB = mB
                self.mL = mL

                self.E = E
                self.E_long = E_long

                self.invImMP = invImMP
                self.MEmiPFlux = MEmiPFlux

                self.dt_guess = dt_guess

                return self
            end

    end

##### Exponential Rosenbrock Euler Mixed Backend ######

    struct WorkerPoolStruct{T} 
        jobs::Channel{Tuple{Int,T,Channel{Nothing}}}
        tasks::Vector{Task}
    end

    mutable struct ExponentialRosenbrockEulerKrylovMixedStruct{HT<:AbstractFloat,DT<:AbstractFloat,HVT<:Vector{HT},DVT<:Union{Vector{DT},CuArray{DT, 1, CUDACore.DeviceMemory}},HMT<:Matrix{HT},DMT<:Union{Matrix{DT},CuArray{DT, 2, CUDACore.DeviceMemory}},HSMT<:SparseMatrixCSC{HT, Int32},DSMT<:Union{CuSparseMatrixCSR{DT, Int32},SparseMatrixCSC{DT,Int32}},BD<:Union{Vector{Int64},Nothing},FD<:Union{DVT,Nothing},DFD<:Union{DVT,Nothing}} <: ImplicitSteppingMethod

        # mixed backend so only MBin multiplication is done on GPU and then transfered back to CPU for the rest of the calculations
        # i.e. Jacobian and F are generated on GPU then transfered back if in a BinaryDomain, if not they are taken from the stored sparse arrays on CPU.

        nworkers::Int64
        WorkerPool::WorkerPoolStruct{HT}

        PhaseSpace::PhaseSpaceStruct
        HostPrecision::Type{HT}
        DevicePrecision::Type{DT}
        LocalMomentumUpdate::Bool

        BinaryInteractions::Bool
        EmissionInteractions::Bool

        DistributionDomainMask::Union{Vector{Int64},Nothing}
        DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}
        ActiveDomain::Vector{Int64}

        MBin::DSMT                    # GPU
        BinaryDomain::BD

        MEmi::Vector{Union{DMT,HSMT}}  # on GPU if dense and CPU is sparse
        #AFlux_d::DVT                     # on CPU
        invAFlux_d::DVT                  # inv Ap flux for time stepping
        XFlux_d::DSMT                    # on GPU
        PFlux_d::DSMT                    # on GPU

        MEmiPFlux::Vector{Union{DMT,HSMT}}# M_Emi - P_Flux for each spatial point, on GPU if dense and CPU if sparse

        Vol::Vector{HT}
        invA::Vector{HT}                 # vector of diagonal entries of invA_Flux for each spatial point (used for scaling)

        Adaptive::Bool
        Implicit::Bool
        dt0::HT
        dt_min::HT                        # suggestion for smallest timestep in form 2^n
        n_cut::HT
        n_cut_d::DT

        step::Int64
        momentum_step::Vector{Int64}    # vector indicating if momentum should be updated for a given step    
        HasInjection::Vector{Bool}     # vector of bools to indicate if a spatial cell has injection       

        MBinMulStep::Vector{DMT}            # GPU temporary array for matrix multiplication of binary terms
        MBinMulStepReshape::Vector{DVT}    # GPU temporary array for reshaped matrix multiplication of binary terms
        f_init::HVT                     # initial distribution function (used by solver to define output struct)
        f::DVT           # current distribution function on GPU
        fstep_d::DVT     # distribution function after a step on GPU
        df::HVT          # change in distribution function
        df_d::DVT        # change in distribution function on GPU
        df_Inj::HVT      # change in distribution function due to injection of particles
        df_Inj_d::DVT    # change in distribution function due to injection of particles on GPU
        df_tmp::HVT      # temporary array the size of f for CFL calculations
        df_tmp_d::DVT    # temporary array the size of f for CFL calculations on GPU
        f_mask_d::FD    # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
        df_mask_d::DFD   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

        F::Vector{HVT}          # vector for implicit solve residuals
        J::Vector{HMT}          # Jacobian matrix for implicit solve 
        JF::Vector{HVT}         # vector for J*F used for scaling
        F_d::Vector{DVT}        # vector for implicit solve residuals on GPU
        J_d::Vector{DMT}        # Jacobian matrix for implicit solve on GPU
        Jsparse::Vector{HSMT}   # sparse Jacobian matrix for implicit solve on CPU
        E::HVT                  # energy vector for correcting step (length of n_momentum)) on CPU
        E_d::DVT                # energy vector for correcting step (length of n_momentum)) on GPU
        E_long_d::DVT           # energy vector for correcting step (length of n_momentum*n_space) on GPU
        D::HVT                  # Diagonal scaling matrix
        D_d::DVT                # Diagonal scaling matrix on GPU
        Dinv::HVT               # Inverse of diagonal scaling matrix
        Dinv_d::DVT             # Inverse of diagonal scaling matrix on GPU
        S::Vector{HVT}                  # Second scaling vector (mutable)
        Sinv::Vector{HVT}               # Inverse of second scaling vector (mutable)
        ϕ::Vector{HMT}          # matrix of ϕ functions for exponential Rosenbrock method

        fold::Vector{HVT}       # local distribution function from previous step 
        fold_d::Vector{DVT}     # local distribution function from previous step on GPU
        fout::Vector{HVT}       # local distribution function from current step
        fout_d::Vector{DVT}     # local distribution function from current step on GPU
        fscale::Vector{HVT}     # scaling vector for exponential Rosenbrock method
        fscale_d::Vector{DVT}   # scaling vector for exponential Rosenbrock method on GPU
        δ::Vector{HVT}          # temporary vector for exponential Rosenbrock method
        δ_d::Vector{DVT}        # temporary vector for exponential Rosenbrock method on GPU

        KsB::Vector{KrylovSubspace{HT,HT,HT,HMT,HMT}}      # Krylov subspace for exponential Rosenbrock method (higher precision for more accuracy)
        KsL::Vector{KrylovSubspace{HT,HT,HT,HMT,HMT}}      # Krylov subspace for exponential Rosenbrock method (higher precision for more accuracy)
        mB::Int64                        # dimension of Binary Krylov subspace
        mL::Int64                        # dimension of Linear Krylov subspace
        ϕcacheB::Vector{ExponentialUtilities.PhivCache{true,HT}} # cache for ϕ functions for Binary terms
        ϕcacheL::Vector{ExponentialUtilities.PhivCache{true,HT}} # cache for ϕ functions for Linear terms

        dt_guess::Vector{HT}             # guess for dt based on CFL condition

        function ExponentialRosenbrockEulerKrylovMixedStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,mB::Int64=128,mL::Int64=32,nworkers::Int64=8,LocalMomentumUpdate::Bool=false)

            Backend = getfield(Main,Symbol("Backend"))
            DevicePrecision = getfield(Main,Symbol("Precision"))
            HostPrecision = Float64

            @assert DevicePrecision == Float32 || DevicePrecision == Float64 "Precision for device arrays must be either Float32 or Float64"

            if Backend isa CPUBackend 
                device_array_type = Matrix{DevicePrecision}
                device_vector_type = Vector{DevicePrecision}
            elseif Backend isa CUDABackend
                device_array_type = CuArray{DevicePrecision,2}
                device_vector_type = CuArray{DevicePrecision,1}
            end

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

            E = zeros(HostPrecision,n_momentum)
            for species in eachindex(PhaseSpace.name_list)
                px_num = px_num_list[species]
                py_num = py_num_list[species]
                pz_num = pz_num_list[species]
                dE = dE_list[species]
                for px in 1:px_num
                    for py in 1:py_num
                        for pz in 1:pz_num
                            idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                            E[idx] = HostPrecision(dE[px])
                        end
                    end
                end
            end
            E_d = zeros(Backend,DevicePrecision,n_momentum)
            copyto!(E_d,E)

            E_long_d = zeros(Backend,DevicePrecision,n_momentum*n_space)
            for space in 1:n_space
                copyto!(@view(E_long_d[(space-1)*n_momentum+1:space*n_momentum]), E)
            end

            BinaryInteractions = !isempty(BinM.Binary_list) && !isnothing(BinM.Domain)
            EmissionInteractions = !isempty(EmiM.Emission_list)

            BinaryDomain = BinM.Domain

            Vol = FluxM.Vol

            MBin = sparse(Backend,DevicePrecision,BinM.M_Bin)

            if Backend isa CUDABackend
                CUDA.pool_status() 
            end

            XFlux_d = sparse(Backend,DevicePrecision,FluxM.X_Flux)
            #AFlux = device_array(Backend,DevicePrecision,FluxM.Ap_Flux) # diagonal matrix of Ap flux 
            invAFlux_d = device_array(Backend,DevicePrecision,1 ./ FluxM.Ap_Flux) # invert Ap Flux for time stepping
            df_Inj = convert(Vector{HostPrecision},copy(Injection))
            df_Inj_d = device_array(Backend,DevicePrecision,df_Inj)

            f_init = convert(Vector{HostPrecision},copy(Initial))
            # cut initial values that are smaller than n_cut in both number and energy
            E_long_tmp = Vector{HostPrecision}(E_long_d) # host copy
            @. f_init = ifelse(f_init < n_cut && f_init * E_long_tmp < n_cut,zero(HostPrecision),f_init)

            f = zeros(Backend,DevicePrecision,length(Initial)) # initial f set in solver 

            if BinaryInteractions
                MBinMulStep = zeros(Backend,DevicePrecision,n_momentum,n_momentum)
                MBinMulStepReshape = reshape(MBinMulStep,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
                MBinMulStep_vec = Vector{device_array_type}(undef,nworkers)
                MBinMulStepReshape_vec = Vector{device_vector_type}(undef,nworkers)
                for i in 1:nworkers
                    MBinMulStep_vec[i] = zeros(Backend,DevicePrecision,n_momentum,n_momentum)
                    MBinMulStepReshape_vec[i] = reshape(MBinMulStep_vec[i],n_momentum^2)
                end
            else
                MBinMulStep = zeros(Backend,DevicePrecision,0,0)
                MBinMulStepReshape = reshape(MBinMulStep,0)
                MBinMulStep_vec = Vector{device_array_type}(undef,nworkers)
                MBinMulStepReshape_vec = Vector{device_vector_type}(undef,nworkers)
                for i in 1:nworkers
                    MBinMulStep_vec[i] = zeros(Backend,DevicePrecision,0,0)
                    MBinMulStepReshape_vec[i] = reshape(MBinMulStep_vec[i],0)
                end
            end
            df = zeros(HostPrecision,length(Initial))
            df_d = device_array(Backend,DevicePrecision,df)
            df_tmp = zeros(HostPrecision,length(Initial))
            df_tmp_d = device_array(Backend,DevicePrecision,df_tmp)

            fstep_d = zeros(Backend,DevicePrecision,length(Initial))
            F = zeros(HostPrecision,n_momentum)
            F_d = device_array(Backend,DevicePrecision,F)
            J = zeros(HostPrecision,n_momentum,n_momentum)
            J_d = device_array(Backend,DevicePrecision,J)
            Jsparse = sparse(zeros(HostPrecision,n_momentum,n_momentum))
            JF = zeros(HostPrecision,n_momentum)
            fold = zeros(HostPrecision,n_momentum)
            fold_d = device_array(Backend,DevicePrecision,fold)
            fout = zeros(HostPrecision,n_momentum)
            fout_d = device_array(Backend,DevicePrecision,fout)
            fscale = zeros(HostPrecision,n_momentum)
            fscale_d = device_array(Backend,DevicePrecision,fscale)
            δ = zeros(HostPrecision,n_momentum)
            δ_d = device_array(Backend,DevicePrecision,δ)
            ϕ = zeros(HostPrecision,n_momentum,3)

            D = one(HostPrecision) ./ copy(E)
            Dinv = copy(E)
            D_d = device_array(Backend,DevicePrecision,D)
            Dinv_d = device_array(Backend,DevicePrecision,Dinv)
            S = similar(D)
            Sinv = similar(Dinv)

            # Making invA = vector of diagonal entries of invA_Flux for each spatial point (used for scaling)
            invA = zeros(HostPrecision,n_space)
            for off_space in 0:n_space-1
                start_idx = n_momentum*off_space + 1
                invA[off_space+1] = 1 / FluxM.Ap_Flux[start_idx]
            end

            # Build new MEmi
            MEmi = Vector{Union{device_array_type,SparseMatrixCSC{HostPrecision,Int32}}}(undef,n_space)
            MEmiPFlux = Vector{Union{device_array_type,SparseMatrixCSC{HostPrecision,Int32}}}(undef,n_space)
            tmpmatrix = zeros(HostPrecision,n_momentum,n_momentum)
            tmpsparsematrix = spzeros(HostPrecision,Int32,n_momentum,n_momentum)
            for off_space in 0:n_space-1
                start_idx = n_momentum*off_space+1
                end_idx = n_momentum*(off_space+1)
                in_Binary = BinaryInteractions && in(off_space,BinaryDomain)
                if in_Binary # dense matrix for binary cells as no default CUDA for dense + sparse matrix addition
                    if Backend isa CUDABackend   
                    fill!(tmpmatrix,zero(HostPrecision))
                    tmpmatrix .-= @view(FluxM.P_Flux[start_idx:end_idx,start_idx:end_idx])
                    if isassigned(EmiM.M_Emi,off_space+1)
                        @assert EmiM.M_Emi[off_space+1] isa Matrix "Emission matrix must be dense for binary interactions"
                        tmpmatrix .+= EmiM.M_Emi[off_space+1]
                    end
                    MEmiPFlux[off_space+1] = copy(Backend,DevicePrecision,tmpmatrix)
                    elseif Backend isa CPUBackend
                        tmpsparsematrix.nzval .= zero(HostPrecision)
                        dropzeros!(tmpsparsematrix)
                        tmpsparsematrix .-= @view(FluxM.P_Flux[start_idx:end_idx,start_idx:end_idx])
                        if isassigned(EmiM.M_Emi,off_space+1)
                            @assert EmiM.M_Emi[off_space+1] isa Matrix "Emission matrix must be dense for binary interactions"
                            tmpsparsematrix .+= EmiM.M_Emi[off_space+1]
                        end
                        dropzeros!(tmpsparsematrix)
                        MEmiPFlux[off_space+1] = copy(tmpsparsematrix)
                    end
                else # sparse matrix for non-binary cells
                    tmpsparsematrix.nzval .= zero(HostPrecision)
                    dropzeros!(tmpsparsematrix)
                    tmpsparsematrix .-= @view(FluxM.P_Flux[start_idx:end_idx,start_idx:end_idx])
                    if isassigned(EmiM.M_Emi,off_space+1)
                        @assert EmiM.M_Emi[off_space+1] isa SparseMatrixCSC "Emission matrix must be sparse for non-binary interactions"
                        tmpsparsematrix .+= EmiM.M_Emi[off_space+1]
                    end
                    dropzeros!(tmpsparsematrix)
                    MEmiPFlux[off_space+1] = copy(tmpsparsematrix)
                end
                #=if isassigned(EmiM.M_Emi,off_space+1)
                    if EmiM.M_Emi[off_space+1] isa SparseMatrixCSC # no binary interactions so stay on GPU for momentum_update
                        M_Emi[off_space+1] = Precision.(EmiM.M_Emi[off_space+1])
                    else
                        M_Emi[off_space+1] = CuArray(Precision.(EmiM.M_Emi[off_space+1]))
                    end
                end=#
            end

            if isnothing(DistributionDomainMask)
                ActiveDomain = InclusiveDomainMask(PhaseSpace)
            else
                ActiveDomain = setdiff(InclusiveDomainMask(PhaseSpace),DistributionDomainMask)
            end

            # work out suggestion for smallest timestep
            dt_min = HostPrecision(Inf)
            for off_space in ActiveDomain
                start_idx = n_momentum*off_space+1
                end_idx = n_momentum*(off_space+1)
                X_Flux_diag_local = @view(diag(FluxM.X_Flux)[start_idx:end_idx])
                invA_local = invA[off_space+1]
                XdivA_local = X_Flux_diag_local .* invA_local
                dt_min_local = nextpow2_float(minimum(PhaseSpace.Spacetime.dt0 ./ XdivA_local))
                dt_min = min(dt_min,dt_min_local)
            end
            # work out momentum_step 
            momentum_step::Vector{Int64} = ones(Int64,n_space)
            if LocalMomentumUpdate
                fill!(momentum_step,256) # random large number so 2^n overflows back to zero without affecting minimum calculations below 
                for off_space in ActiveDomain
                    start_idx = n_momentum*off_space+1
                    end_idx = n_momentum*(off_space+1)
                    X_Flux_diag_local = @view(diag(FluxM.X_Flux)[start_idx:end_idx])
                    invA_local = invA[off_space+1]
                    XdivA_local = X_Flux_diag_local .* invA_local
                    # Δf/f = -X/A * Δt/Δt0 (Δt = Δt0 when evaluating X and A)
                    # Therefore Δtmin = min(-A * Δt0 / X) is the smallest timestep required for Δf/f < 1
                    # we want this in the form 2^n where 2^n > Δtmin is the smallest power of 2 greater than Δtmin (as space update is perfomred using two sets in the split operator scheme)
                    # here n = steps
                    steps = Int64(log2(nextpow2_float(minimum(PhaseSpace.Spacetime.dt0 ./ XdivA_local))))
                    momentum_step[off_space+1] = steps
                end
                # now work out minimum of steps 
                min_steps = minimum(momentum_step)
                println("minsteps: $min_steps")
                # we now want to scale n -> n - nmin so that each n is a multiple of the smallest timestep 
                momentum_step .-= min_steps 
                @. momentum_step = ifelse(momentum_step > 0, momentum_step-1, momentum_step) # do twice the number of required steps if > 2^0
                # now generate Δt intervals for each momentum step to be evaluated
                momentum_step = 2 .^ momentum_step
                # now when mod(method.step,momentum_steps) = 0 momentum update should be evaluated
            end

            # has injection? 
            HasInjection = zeros(Bool,n_space)
            for off_space in 0:n_space-1 
                idx_range = n_momentum*off_space+1:n_momentum*(off_space+1)
                df_Inj_local = @view(df_Inj[idx_range])
                HasInjection[off_space+1] = !all(iszero,df_Inj_local)
            end

            if !isnothing(DistributionDomainMask)
                f_mask = ones(HostPrecision,length(Initial))
                for off_space_idx in DistributionDomainMask
                    for species_idx in eachindex(PhaseSpace.name_list)
                    LocationSpeciesToStateVector(f_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= HostPrecision(0)
                    end
                end
                f_mask_d = device_array(Backend,DevicePrecision,f_mask)
            else
                f_mask_d = nothing
            end

            if !isnothing(DeltaDistributionDomainMask)
                df_mask = ones(HostPrecision,length(Initial))
                for off_space_idx in DeltaDistributionDomainMask
                    for species_idx in eachindex(PhaseSpace.name_list)
                    LocationSpeciesToStateVector(df_mask,PhaseSpace,off_space_idx=off_space_idx,species_index=species_idx) .= HostPrecision(0)
                    end
                end
                df_mask_d = device_array(Backend,DevicePrecision,df_mask)
            else
                df_mask_d = nothing
            end

            dt_guess = zeros(HostPrecision, n_space)
            fill!(dt_guess,ldexp(dt_initial, -5)) # initial guess for dt is 1/64 of the initial dt

            ###### Actually Build the Struct with Concrete Types ######

            self = new{HostPrecision,DevicePrecision,Vector{HostPrecision},typeof(f),Matrix{HostPrecision},typeof(J_d),SparseMatrixCSC{HostPrecision, Int32},typeof(MBin),typeof(BinaryDomain),typeof(f_mask_d),typeof(df_mask_d)}()

            self.nworkers = nworkers

            self.PhaseSpace = PhaseSpace
            self.HostPrecision = HostPrecision
            self.DevicePrecision = DevicePrecision
            self.LocalMomentumUpdate = LocalMomentumUpdate
 
            self.BinaryInteractions = BinaryInteractions
            self.EmissionInteractions = EmissionInteractions

            self.DistributionDomainMask = DistributionDomainMask
            self.DeltaDistributionDomainMask = DeltaDistributionDomainMask
            self.ActiveDomain = ActiveDomain

            self.MBin = MBin
            self.BinaryDomain = BinaryDomain

            self.MEmi = MEmi
            #self.AFlux = AFlux
            self.invAFlux_d = invAFlux_d
            self.XFlux_d = XFlux_d

            self.MEmiPFlux = MEmiPFlux

            self.Vol = Vol
            self.invA = invA

            self.Adaptive = Adaptive
            self.Implicit = true
            self.dt0 = PhaseSpace.Spacetime.dt0
            self.dt_min = dt_min
            self.n_cut = HostPrecision(n_cut)
            self.n_cut_d = DevicePrecision(n_cut)

            self.step = 0
            self.momentum_step = momentum_step
            self.HasInjection = HasInjection

            self.MBinMulStep = MBinMulStep_vec
            self.MBinMulStepReshape = MBinMulStepReshape_vec
            self.f_init = f_init
            self.f = f

            self.fstep_d = fstep_d

            self.df = df
            self.df_d = df_d

            self.df_Inj = df_Inj
            self.df_Inj_d = df_Inj_d

            self.df_tmp = df_tmp
            self.df_tmp_d = df_tmp_d

            self.f_mask_d = f_mask_d
            self.df_mask_d = df_mask_d

            self.F = [similar(F) for _ in 1:nworkers]
            self.J = [similar(J) for _ in 1:nworkers]
            self.JF = [similar(JF) for _ in 1:nworkers]
            self.F_d = [similar(F_d) for _ in 1:nworkers]
            self.J_d = [similar(J_d) for _ in 1:nworkers]
            self.Jsparse = [similar(Jsparse) for _ in 1:nworkers]
            self.E = E
            self.E_d = E_d
            self.E_long_d = E_long_d
            self.D = D 
            self.D_d = D_d
            self.Dinv = Dinv
            self.Dinv_d = Dinv_d
            self.S = [similar(S) for _ in 1:nworkers]
            self.Sinv = [similar(Sinv) for _ in 1:nworkers]

            self.ϕ = [similar(ϕ) for _ in 1:nworkers]

            self.fold = [similar(fold) for _ in 1:nworkers]
            self.fold_d = [similar(fold_d) for _ in 1:nworkers]
            self.fout = [similar(fout) for _ in 1:nworkers]
            self.fout_d = [similar(fout_d) for _ in 1:nworkers]
            self.fscale = [similar(fscale) for _ in 1:nworkers]
            self.fscale_d = [similar(fscale_d) for _ in 1:nworkers]
            self.δ = [similar(δ) for _ in 1:nworkers]
            self.δ_d = [similar(δ_d) for _ in 1:nworkers]

            self.KsB = [KrylovSubspace{Float64,Float64,Array{Float64,2}}(n_momentum,mB) for _ in 1:nworkers]
            self.KsL = [KrylovSubspace{Float64,Float64,Array{Float64,2}}(n_momentum,mL) for _ in 1:nworkers]
            self.mB = mB
            self.mL = mL
            self.ϕcacheB = [ExponentialUtilities.PhivCache(ϕ,mB,2) for _ in 1:nworkers]
            self.ϕcacheL = [ExponentialUtilities.PhivCache(ϕ,mL,2) for _ in 1:nworkers]

            self.dt_guess = dt_guess

            self.WorkerPool = WorkerPoolStruct(self,nworkers)

            return self
        end

    end
    
    function WorkerPoolStruct(method::ExponentialRosenbrockEulerKrylovMixedStruct,nwokers::Int)

        jobs = Channel{Tuple{Int,method.HostPrecision,Channel{Nothing}}}(length(method.ActiveDomain))
        tasks = [Threads.@spawn worker!(method,wid,jobs) for wid in 1:nwokers]

        return WorkerPoolStruct{method.HostPrecision}(jobs, tasks)
    end

    function close_pool!(pool::WorkerPoolStruct)
        close(pool.jobs)
        wait.(pool.tasks)
    end
