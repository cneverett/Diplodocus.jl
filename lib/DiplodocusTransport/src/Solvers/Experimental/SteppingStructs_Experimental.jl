##### Heun's Method ######
    mutable struct HeunStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},MET<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ExplicitSteppingMethod

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
    f_step::VT                     # distribution function after first step of Heun's method
    df::VT                         # change in distribution function
    df_step::VT                    # change in distribution function after first step of Heun's method
    df_Bin::VT                     # change in distribution function due to binary interactions
    df_Emi::VT                     # change in distribution function due to emission interactions
    df_Flux::VT                    # change in distribution function due to fluxes
    df_Inj::VT                     # change in distribution function due to injection of particles
    df_tmp::VT                     # temporary array the size of f for CFL calculations
    f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
    df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

    function HeunStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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
        df_step = zeros(Backend,Precision,length(Initial))
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
        f_step = convert(Vector{Precision},copy(Initial))
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
            f_step = CuArray(f_step)
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
        self.f_step = f_step
        self.df_Inj = df_Inj
        self.M_Bin_Mul_Step = M_Bin_Mul_Step
        self.M_Bin_Mul_Step_reshape = M_Bin_Mul_Step_reshape
        self.df = df
        self.df_step = df_step
        self.df_Bin = df_Bin
        self.df_Emi = df_Emi
        self.df_Flux = df_Flux
        self.df_tmp = df_tmp
        self.f_mask = f_mask
        self.df_mask = df_mask

        return self
    end

    end



##### Modified Patankar Euler Method ######
    mutable struct MPEStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

    PhaseSpace::PhaseSpaceStruct

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    M_Bin::MBT
    Bin_Domain::BD

    FMEmi::SMT                     # sum of space and momentum fluxes plus any emissive interactions 
    MEmi::SMT
    A_Flux::SMT
    Q::SMT                         # sparse matrix to be inverted 
    Qtemplate::SMT                 # sparse matrix to be inverted template (same sparsity pattern as Q to allow easy addition of sparse and dense matrices to build Q each time step)  
    Vol::Vector{T}
    invA_Flux::SMT                  # inv Ap flux for time stepping

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
    f_tmp::VT                      # temp distribution function 
    df::VT                         # change in distribution function
    df_step::VT                    # change in distribution function after first step of Heun's method
    df_Bin::VT                     # change in distribution function due to binary interactions
    df_Emi::VT                     # change in distribution function due to emission interactions
    df_Flux::VT                    # change in distribution function due to fluxes
    df_Inj::VT                     # change in distribution function due to injection of particles
    df_tmp::VT                     # temporary array the size of f for CFL calculations
    f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
    df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

    function MPEStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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

        Vol = FluxM.Vol

        f_init = convert(Vector{Precision},Initial)
        M_Bin = BinM.M_Bin
        M_Emi = EmiM.M_Emi
        # build FMEmi term
        nzX_Flux = findnz(FluxM.X_Flux)
        nzP_Flux = findnz(FluxM.P_Flux)
        nzM_Emi = findnz(M_Emi)
        FMEmi_I = vcat(nzX_Flux[1],nzP_Flux[1],nzM_Emi[1])
        FMEmi_J = vcat(nzX_Flux[2],nzP_Flux[2],nzM_Emi[2])
        FMEmi_V = vcat(nzX_Flux[3] .* Precision(-1),nzP_Flux[3] .* Precision(-1),nzM_Emi[3] .* Precision(1)) # minus on X and P flux values as these terms are on LHS of transport equation
        FMEmi = sparse(FMEmi_I,FMEmi_J,FMEmi_V,size(FluxM.X_Flux,1),size(FluxM.X_Flux,2))
        MEmi = M_Emi
        EmiM = nothing # for GC 
        # build A_Flux term
        A_Flux = spdiagm(FluxM.Ap_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
        invA_Flux = spdiagm(1 ./ FluxM.Ap_Flux) # invert Ap Flux for time stepping
        # build Q and Qtemplate terms
        block = sparse(ones(Precision,n_momentum,n_momentum))
        BlockDiag = blockdiag([block for _ in 1:n_space]...)
        nzBlockDiag = findnz(BlockDiag)
        nzFMEmi = findnz(FMEmi)
        nzAFlux = findnz(A_Flux)
        Qtemplate_I = vcat(nzBlockDiag[1],nzFMEmi[1])
        Qtemplate_J = vcat(nzBlockDiag[2],nzFMEmi[2])
        Qtemplate_V = vcat(nzBlockDiag[3] * zero(Precision),nzFMEmi[3]) # sets up  (F+M_Emi) with zeros in locations where MBin terms will go
        Qtemplate = sparse(Qtemplate_I,Qtemplate_J,Qtemplate_V,size(FluxM.X_Flux,1),size(FluxM.X_Flux,2))
        # TODO: simplify Qtemplate for cases where there are no emissive interactions or Binary interactions or just binary interactions at certain spatial locations
        Q = copy(Qtemplate)

        if Binary_Interactions
            M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
            f_ones = ones(Backend,Precision,n_momentum) # temporary array of ones for matrix multiplication to build Qtemplate
            for off_space in Bin_Domain
                start_idx = n_momentum*off_space+1
                end_idx = n_momentum*(off_space+1)
                @inbounds vol = Vol[off_space+1]
                mul!(M_Bin_Mul_Step_reshape,M_Bin,f_ones,vol,zero(eltype(f_ones))) 
                @. @view(Qtemplate[start_idx:end_idx,start_idx:end_idx]) += M_Bin_Mul_Step
            end
            fill!(M_Bin_Mul_Step,zero(Precision)) # zero out M_Bin_Mul_Step after using it to build Qtemplate
            fill!(M_Bin_Mul_Step_reshape,zero(Precision)) # zero out M_Bin_Mul_Step_reshape after using it to build Qtemplate
        else
            M_Bin_Mul_Step = zeros(Backend,Precision,0,0)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,0)
        end
        Qtemplate = invA_Flux * Qtemplate 
        display(Qtemplate)

        f_tmp = zeros(Backend,Precision,length(Initial))
        df = zeros(Backend,Precision,length(Initial))
        df_Bin = zeros(Backend,Precision,length(Initial))
        df_Emi = zeros(Backend,Precision,length(Initial))
        df_Flux = zeros(Backend,Precision,length(Initial))
        df_tmp = zeros(Backend,Precision,length(Initial))

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
            FMEmi = CuSparseMatrixCSC(FMEmi) # sum of space and momentum fluxes
            A_Flux = CuSparseMatrixCSC(A_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
            f = CuArray(f)
            df_Inj = CuArray(df_Inj)
            Qtemplate = CuSparseMatrixCSC(Qtemplate)
            Q = CuSparseMatrixCSC(Q)
            invA_Flux = CuArray(invA_Flux)
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

        GC.gc() # force garbage collection before building the struct to free up memory from temporary arrays

        ###### Actually Build the Struct with Concrete Types ######

        self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(FMEmi),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

        self.Implicit = true
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
        self.FMEmi = FMEmi
        self.A_Flux = A_Flux
        self.Vol = Vol
        self.invA_Flux = invA_Flux
        self.f = f
        self.f_tmp = f_tmp
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
        self.Q = Q
        self.Qtemplate = Qtemplate
        self.MEmi = MEmi

        return self
    end

    end

##### Symplectic Modified Patankar Euler Method ######
    mutable struct SymplecticMPEStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MET<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},LUT<:Union{LinearAlgebra.LU{T, <:AbstractMatrix{T}, <:AbstractVector{<:Integer}},Nothing},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

    PhaseSpace::PhaseSpaceStruct
    Precision::Type{T}

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    #M_Bin::MBT
    Gijk::MBT
    Lij::MBT
    Bin_Domain::BD

    M_Emi::MET
    A_Flux::SMT
    X_Flux::SMT
    P_Flux::SMT
    Q::MT                           # matrix to be inverted
    QLU::LUT                        # LU factorization of the matrix for implicit solving
    E::VT                           # energy vector for correcting Patankar energy error
    GMRESWorkspace::GmresWorkspace
    Vol::Vector{T}
    invA_Flux::SMT                  # inv Ap flux for time stepping

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
    f_step::VT                     # change distribution function after first step
    f_tmp::VT                      # temp distribution function 
    df::VT                         # change in distribution function
    df_step::VT                    # change in distribution function after first step
    df_Bin::VT                     # change in distribution function due to binary interactions
    df_Emi::VT                     # change in distribution function due to emission interactions
    df_Flux::VT                    # change in distribution function due to fluxes
    df_Inj::VT                     # change in distribution function due to injection of particles
    df_tmp::VT                     # temporary array the size of f for CFL calculations
    f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
    df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

    function SymplecticMPEStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStructPatankar,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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

        E = zeros(Backend,Precision,length(Initial))
        for species in eachindex(PhaseSpace.name_list)
            px_num = px_num_list[species]
            py_num = py_num_list[species]
            pz_num = pz_num_list[species]
            dE = dE_list[species]
            for px in 1:px_num
                for py in 1:py_num
                    for pz in 1:pz_num
                        idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                        E[idx] = dE[px]
                    end
                end
            end
        end

        n_space = x_num*y_num*z_num
        n_momentum = 0
        for i in eachindex(px_num_list)
            n_momentum += px_num_list[i]*py_num_list[i]*pz_num_list[i]
        end

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        Binary_Interactions = !isempty(BinM.Binary_list) && !isnothing(BinM.Domain)
        Emission_Interactions = !isempty(EmiM.Emission_list)

        Bin_Domain = BinM.Domain

        Vol = FluxM.Vol

        f_init = convert(Vector{Precision},Initial)
        #M_Bin = BinM.M_Bin
        Gijk = BinM.Gijk
        Lij = BinM.Liij
        M_Emi = EmiM.M_Emi
        M_Emi.nzval[abs.(M_Emi.nzval) .< 1e4eps(Float64)] .= zero(Precision)
        # build FMEmi term
        X_Flux = FluxM.X_Flux
        P_Flux = copy(FluxM.P_Flux)
        P_Flux[diagind(P_Flux)] .= zero(eltype(P_Flux)) # remove diagonal of P_Flux for Patankar method
        dropzeros!(P_Flux) # drop zeros from P_Flux after removing diagonal
        # build A_Flux term
        A_Flux = spdiagm(FluxM.Ap_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
        invA_Flux = spdiagm(1 ./ FluxM.Ap_Flux) # invert Ap Flux for time stepping

        if Binary_Interactions
            M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
            fill!(M_Bin_Mul_Step,zero(Precision)) # zero out M_Bin_Mul_Step after using it to build Qtemplate
            fill!(M_Bin_Mul_Step_reshape,zero(Precision)) # zero out M_Bin_Mul_Step_reshape after using it to build Qtemplate
        else
            M_Bin_Mul_Step = zeros(Backend,Precision,0,0)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,0)
        end

        df = zeros(Backend,Precision,length(Initial))
        df_Bin = zeros(Backend,Precision,length(Initial))
        df_Emi = zeros(Backend,Precision,length(Initial))
        df_Flux = zeros(Backend,Precision,length(Initial))
        df_tmp = zeros(Backend,Precision,length(Initial))
        df_Inj = convert(Vector{Precision},copy(Injection))

        f = convert(Vector{Precision},copy(Initial))
        @. f = ifelse(f<=n_cut,zero(eltype(f)),f)
        f_tmp = zeros(Backend,Precision,length(Initial))
        f_step = zeros(Backend,Precision,length(Initial))

        
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
            A_Flux = CuSparseMatrixCSC(A_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
            f = CuArray(f)
            df_Inj = CuArray(df_Inj)
            invA_Flux = CuArray(invA_Flux)
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

        # setup for building Q and its LU factorization
        Q = zeros(Backend,Precision,n_momentum,n_momentum)
        QLU = lu(Q+I)
        GMRESWorkspace = GmresWorkspace(n_momentum,n_momentum, typeof(f))

        GC.gc() # force garbage collection before building the struct to free up memory from temporary arrays

        ###### Actually Build the Struct with Concrete Types ######

        self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Emi),typeof(Gijk),typeof(X_Flux),typeof(QLU),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

        self.PhaseSpace = PhaseSpace
        self.Precision = Precision
        self.Implicit = true
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
        #self.M_Bin = M_Bin
        self.Gijk = Gijk
        self.Lij = Lij
        self.M_Emi = M_Emi
        self.A_Flux = A_Flux
        self.invA_Flux = invA_Flux
        self.X_Flux = X_Flux
        self.P_Flux = P_Flux
        self.Q = Q
        self.QLU = QLU
        self.E = E
        self.GMRESWorkspace = GMRESWorkspace
        self.Vol = Vol
        self.f = f
        self.f_tmp = f_tmp
        self.f_step = f_step
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

        return self
    end

    end

##### Symplectic Symmetric Modified Patankar Euler Method ######
    mutable struct SymplecticSymmetricMPEStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MET<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},LUT<:Union{LinearAlgebra.LU{T, <:AbstractMatrix{T}, <:AbstractVector{<:Integer}},Nothing},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

    PhaseSpace::PhaseSpaceStruct
    Precision::Type{T}

    Binary_Interactions::Bool
    Emission_Interactions::Bool

    #M_Bin::MBT
    Gijk::MBT
    Lijk::MBT
    Liij::MBT
    Liik::MBT
    Aijk::MBT
    Bin_Domain::BD

    M_Emi::MET
    A_Flux::SMT
    X_Flux::SMT
    P_Flux::SMT
    Q::MT                           # matrix to be inverted
    QLU::LUT                        # LU factorization of the matrix for implicit solving
    E::VT                           # energy vector for correcting Patankar energy error
    GMRESWorkspace::GmresWorkspace
    Vol::Vector{T}
    invA_Flux::SMT                  # inv Ap flux for time stepping

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
    f_step::VT                     # change distribution function after first step
    f_tmp::VT                      # temp distribution function 
    df::VT                         # change in distribution function
    df_step::VT                    # change in distribution function after first step
    df_Bin::VT                     # change in distribution function due to binary interactions
    df_Emi::VT                     # change in distribution function due to emission interactions
    df_Flux::VT                    # change in distribution function due to fluxes
    df_Inj::VT                     # change in distribution function due to injection of particles
    df_tmp::VT                     # temporary array the size of f for CFL calculations
    f_mask::FD                     # mask for spatial domain f (1 for points in domain, 0 for points outside domain) 
    df_mask::DFD                   # mask for spatial domain df (1 for points in domain, 0 for points outside domain)  

    function SymplecticSymmetricMPEStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStructPatankarSymmetric,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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

        E = zeros(Backend,Precision,length(Initial))
        for species in eachindex(PhaseSpace.name_list)
            px_num = px_num_list[species]
            py_num = py_num_list[species]
            pz_num = pz_num_list[species]
            dE = dE_list[species]
            for px in 1:px_num
                for py in 1:py_num
                    for pz in 1:pz_num
                        idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                        E[idx] = dE[px]
                    end
                end
            end
        end

        n_space = x_num*y_num*z_num
        n_momentum = 0
        for i in eachindex(px_num_list)
            n_momentum += px_num_list[i]*py_num_list[i]*pz_num_list[i]
        end

        @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

        Binary_Interactions = !isempty(BinM.Binary_list) && !isnothing(BinM.Domain)
        Emission_Interactions = !isempty(EmiM.Emission_list)

        Bin_Domain = BinM.Domain

        Vol = FluxM.Vol

        f_init = convert(Vector{Precision},Initial)
        #M_Bin = BinM.M_Bin
        Gijk = BinM.Gijk
        Lijk = BinM.Lijk
        Liij = BinM.Liij
        Liik = BinM.Liik
        Aijk = zeros(Precision,size(Gijk)) # preallocate Aijk as the same size as Gijk and Lijk
        @. Aijk = Gijk + Lijk
        #for i in 1:n_momentum, j in 1:n_momentum, k in 1:n_momentum
        #    Aijk[(j-1)*n_momentum+(i-1)+1,k] *= E[i] / E[j] / E[k]
            #println("i: $i, j: $j, k: $k")
        #end
        M_Emi = EmiM.M_Emi
        M_Emi.nzval[abs.(M_Emi.nzval) .< 1e4eps(Float64)] .= zero(Precision)
        dropzeros!(M_Emi) # drop zeros from M_Emi after setting small values to zero
        # build FMEmi term
        X_Flux = FluxM.X_Flux
        P_Flux = copy(FluxM.P_Flux)
        P_Flux[diagind(P_Flux)] .= zero(eltype(P_Flux)) # remove diagonal of P_Flux for Patankar method
        dropzeros!(P_Flux) # drop zeros from P_Flux after removing diagonal
        # build A_Flux term
        A_Flux = spdiagm(FluxM.Ap_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
        invA_Flux = spdiagm(1 ./ FluxM.Ap_Flux) # invert Ap Flux for time stepping

        if Binary_Interactions
            M_Bin_Mul_Step = zeros(Backend,Precision,n_momentum,n_momentum)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2) # Thanks to Emma Godden for fixing a bug here
            fill!(M_Bin_Mul_Step,zero(Precision)) # zero out M_Bin_Mul_Step after using it to build Qtemplate
            fill!(M_Bin_Mul_Step_reshape,zero(Precision)) # zero out M_Bin_Mul_Step_reshape after using it to build Qtemplate
        else
            M_Bin_Mul_Step = zeros(Backend,Precision,0,0)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,0)
        end

        df = zeros(Backend,Precision,length(Initial))
        df_Bin = zeros(Backend,Precision,length(Initial))
        df_Emi = zeros(Backend,Precision,length(Initial))
        df_Flux = zeros(Backend,Precision,length(Initial))
        df_tmp = zeros(Backend,Precision,length(Initial))
        df_Inj = convert(Vector{Precision},copy(Injection))

        f = convert(Vector{Precision},copy(Initial))
        @. f = ifelse(f<=n_cut,zero(eltype(f)),f)
        f_tmp = zeros(Backend,Precision,length(Initial))
        f_step = zeros(Backend,Precision,length(Initial))

        
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
            A_Flux = CuSparseMatrixCSC(A_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
            f = CuArray(f)
            df_Inj = CuArray(df_Inj)
            invA_Flux = CuArray(invA_Flux)
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

        # setup for building Q and its LU factorization
        Q = zeros(Backend,Precision,n_momentum,n_momentum)
        QLU = lu(Q+I)
        GMRESWorkspace = GmresWorkspace(n_momentum,n_momentum, typeof(f))

        GC.gc() # force garbage collection before building the struct to free up memory from temporary arrays

        ###### Actually Build the Struct with Concrete Types ######

        self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Emi),typeof(Gijk),typeof(X_Flux),typeof(QLU),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

        self.PhaseSpace = PhaseSpace
        self.Precision = Precision
        self.Implicit = true
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
        #self.M_Bin = M_Bin
        self.Gijk = Gijk
        self.Lijk = Lijk
        self.Liij = Liij
        self.Liik = Liik
        self.Aijk = Aijk
        self.M_Emi = M_Emi
        self.A_Flux = A_Flux
        self.invA_Flux = invA_Flux
        self.X_Flux = X_Flux
        self.P_Flux = P_Flux
        self.Q = Q
        self.QLU = QLU
        self.E = E
        self.GMRESWorkspace = GMRESWorkspace
        self.Vol = Vol
        self.f = f
        self.f_tmp = f_tmp
        self.f_step = f_step
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

        return self
    end

    end



##### Backwards Euler Method ######

    mutable struct BackwardEulerStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},MET<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

            PhaseSpace::PhaseSpaceStruct
            Precision::Type{T}

            Binary_Interactions::Bool
            Emission_Interactions::Bool

            M_Bin::MBT
            Bin_Domain::BD

            M_Emi::MET
            A_Flux::SMT
            invA_Flux::SMT                  # inv Ap flux for time stepping
            X_Flux::SMT
            P_Flux::SMT

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

            function BackwardEulerStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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

                fstep = zeros(Backend,Precision,length(Initial))
                F = zeros(Backend,Precision,n_momentum)
                J = zeros(Backend,Precision,n_momentum,n_momentum)

                Vol = FluxM.Vol

                f_init = convert(Vector{Precision},Initial)
                M_Bin = BinM.M_Bin
                M_Emi = EmiM.M_Emi
                X_Flux = FluxM.X_Flux
                P_Flux = FluxM.P_Flux
                A_Flux = spdiagm(FluxM.Ap_Flux) # diagonal matrix of Ap flux for Modified Patankar Euler method
                invA_Flux = spdiagm(1 ./ FluxM.Ap_Flux) # invert Ap Flux for time stepping
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
                    X_Flux = CuSparseMatrixCSC(X_Flux)
                    P_Flux = CuSparseMatrixCSC(P_Flux)
                    invA_Flux = CuArray(invA_Flux)
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

                self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(M_Emi),typeof(X_Flux),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

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
                self.PhaseSpace = PhaseSpace
                self.Bin_Domain = BinM.Domain
                self.f_init = f_init
                self.M_Bin = M_Bin
                self.M_Emi = M_Emi
                self.X_Flux = X_Flux
                self.P_Flux = P_Flux
                self.A_Flux = A_Flux
                self.invA_Flux = invA_Flux
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

                return self
            end

    end

##### ES2 With Correction Method ######

    mutable struct ES2Struct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},MET<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

            PhaseSpace::PhaseSpaceStruct
            Precision::Type{T}

            Binary_Interactions::Bool
            Emission_Interactions::Bool

            M_Bin::MBT
            Bin_Domain::BD

            M_Emi::MET
            A_Flux::SMT
            invA_Flux::SMT                  # inv Ap flux for time stepping
            X_Flux::SMT
            P_Flux::SMT

            invImMP::SMT                      # (I-dt*A^{-1}(M_Emi-P_Flux))^{-1} for momentum update

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

            E::VT                           # energy vector for correcting step 
            N::MT                           # matrix of number densities for each species for correcting step [N1, N2, ...]

            function ES2Struct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing)

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
                N = zeros(Backend,Precision,length(PhaseSpace.name_list),n_momentum)
                for species in eachindex(PhaseSpace.name_list)
                    px_num = px_num_list[species]
                    py_num = py_num_list[species]
                    pz_num = pz_num_list[species]
                    dE = dE_list[species]
                    for px in 1:px_num
                        for py in 1:py_num
                            for pz in 1:pz_num
                                idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                                E[idx] = Precision(dE[px])
                                N[species,idx] = Precision(1.0)
                            end
                        end
                    end
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

                fstep = zeros(Backend,Precision,length(Initial))
                F = zeros(Backend,Precision,n_momentum)
                J = zeros(Backend,Precision,n_momentum,n_momentum)

                Vol = FluxM.Vol

                f_init = convert(Vector{Precision},Initial)
                M_Bin = Precision.(BinM.M_Bin)
                M_Emi = Precision.(EmiM.M_Emi)
                X_Flux = Precision.(FluxM.X_Flux)
                P_Flux = Precision.(FluxM.P_Flux)
                A_Flux = Precision.(spdiagm(FluxM.Ap_Flux)) # diagonal matrix of Ap flux for Modified Patankar Euler method
                invA_Flux = Precision.(spdiagm(1 ./ FluxM.Ap_Flux)) # invert Ap Flux for time stepping
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
                    X_Flux = CuSparseMatrixCSC(X_Flux)
                    P_Flux = CuSparseMatrixCSC(P_Flux)
                    invA_Flux = CuArray(invA_Flux)
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

                #= Making invMP = (I-dt*A^{-1}(M_Emi-P_Flux))^{-1}
                  This assuming only emissive interactions coming from one set of specie to another (e.g. electron to photon but not photon to electron) (I-A^{-1}(M_Emi-P_Flux)) has a block triangular structure:
                   I-MP = [A 0]
                          [B C]
                   where A corresponds to the non-emissive species and C corresponds to the emissive species. This means we can invert (I-MP) as:
                   inv(I-MP) = [A^{-1}          0    ]
                               [-C^{-1}BA^{-1} C^{-1}]    

                =#

                ImMP = I - (dt_initial/2)*invA_Flux*(M_Emi - P_Flux)
                invImMP = spzeros(Precision,size(P_Flux))
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

                        invImMP_view = @view(invImMP[pi_low:pi_up,pi_low:pi_up])
                        ImMP_view = @view(ImMP[pi_low:pi_up,pi_low:pi_up])

                        invImMP_view .= sparse(inv(ImMP_view))
                    end
                    # off-diagonal blocks
                    #=
                        for i>j and Bij = the ij components of the inverse matrix

                            Bij = -Bii \sum_{k=j}^{i-1} Aik Bkj

                    =#
                    for speciesi in eachindex(PhaseSpace.name_list)
                        pi_low = momentum_offset[speciesi] + off_space + 1
                        pi_up = momentum_offset[speciesi+1] + off_space


                        for speciesj in 1:speciesi-1
                            pj_low = momentum_offset[speciesj] + off_space + 1
                            pj_up = momentum_offset[speciesj+1] + off_space

                            
                            Bij = @view(invImMP[pi_low:pi_up,pj_low:pj_up])

                            for speciesk in speciesj:speciesi-1
                                pk_low = momentum_offset[speciesk] + off_space + 1
                                pk_up = momentum_offset[speciesk+1] + off_space

                                Aik = @view(ImMP[pi_low:pi_up,pk_low:pk_up])
                                Bkj = @view(invImMP[pk_low:pk_up,pj_low:pj_up])

                                Bij .-= sparse(Aik * Bkj) 
                            end
                        end

 
                    end
                end

                # cut initial values that are smaller than n_cut 
                    @. f_init = ifelse(f_init<=n_cut,zero(eltype(f_init)),f_init)
                    @. f = ifelse(f<=n_cut,zero(eltype(f)),f)

                ###### Actually Build the Struct with Concrete Types ######

                self = new{Precision,typeof(f),typeof(M_Bin_Mul_Step),typeof(M_Bin),typeof(M_Emi),typeof(X_Flux),typeof(Bin_Domain),typeof(f_mask),typeof(df_mask)}()

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
                self.PhaseSpace = PhaseSpace
                self.Bin_Domain = BinM.Domain
                self.f_init = f_init
                self.M_Bin = M_Bin
                self.M_Emi = M_Emi
                self.X_Flux = X_Flux
                self.P_Flux = P_Flux
                self.A_Flux = A_Flux
                self.invA_Flux = invA_Flux
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

                self.E = E
                self.N = N

                self.invImMP = invImMP

                return self
            end

    end

##### Experimental Krylov method #####
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

#######################

    mutable struct ExpRBKIOPSStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

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
            g::VT                          # vector for implicit solve residuals
            U::MT                          # matrix for implicit solve residuals
            Jsparse::SMT                   # sparse Jacobian matrix for implicit solve
            D::VT             # Diagonal scaling matrix
            Dinv::VT          # Inverse of diagonal scaling matrix
            ϕ::MT                          # matrix of ϕ functions for exponential Rosenbrock method

            fold::VT                        # local distribution function from previous step 
            fout::VT                        # local distribution function from current step
            fscale::VT                      # scaling vector for exponential Rosenbrock method
            δ::VT                          # temporary vector for exponential Rosenbrock method

            Ks#::KrylovSubspace{T,T,T,MT,AbstractMatrix{T}}      # Krylov subspace for exponential Rosenbrock method
            m::Int64                        # dimension of Krylov subspace
            ϕcache::ExponentialUtilities.PhivCache{useview,T} where useview # cache for ϕ functions

            E::VT                           # energy vector for correcting step
            
            KIOPS_workspace::KIOPSRosenbrockWorkspace{T,MT,VT} # workspace for KIOPS method

            dt_guess::Vector{T}             # vector of dt guesses for adaptive time stepping

            function ExpRBKIOPSStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,m::Int64=128)

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
                g = zeros(Backend,Precision,n_momentum)
                Jsparse = sparse(zeros(Backend,Precision,n_momentum,n_momentum))
                fold = zeros(Backend,Precision,n_momentum)
                fout = zeros(Backend,Precision,n_momentum)
                fscale = zeros(Backend,Precision,n_momentum)
                δ = zeros(Backend,Precision,n_momentum)
                ϕ = zeros(Backend,Precision,n_momentum,2)
                ϕcache = ExponentialUtilities.PhivCache(ϕ,m,1)

                D = one(Precision) ./ copy(E)
                Dinv = copy(E)

                if Backend isa CUDABackend
                    Ks = KrylovSubspace{Precision,Precision,CuArray{Precision,2}}(n_momentum,m)
                else
                    Ks = KrylovSubspace{Precision,Precision,Array{Precision,2}}(n_momentum,m)
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

                            
                            Bij = @view(invImMP[pi_low:pi_up,pj_low:pj_up])

                            for speciesk in speciesj:speciesi-1
                                pk_low = momentum_offset[speciesk] + off_space + 1
                                pk_up = momentum_offset[speciesk+1] + off_space

                                Aik = @view(ImMP[pi_low:pi_up,pk_low:pk_up])
                                Bkj = @view(invImMP[pk_low:pk_up,pj_low:pj_up])

                                Bij .-= sparse(Aik * Bkj) 
                            end
                        end


                    end=#
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

                KIOPS_workspace = KIOPSRosenbrockWorkspace(J, F;mmin = 10,mmax = m,orth_len = 4)

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
                self.g = g
                self.Jsparse = Jsparse
                self.D = D 
                self.Dinv = Dinv 
                self.fold = fold
                self.fout = fout
                self.fscale = fscale
                self.δ = δ
                self.Ks = Ks
                self.ϕ = ϕ
                self.ϕcache = ϕcache
                self.m = m

                self.E = E

                self.invImMP = invImMP

                self.KIOPS_workspace = KIOPS_workspace

                self.dt_guess = dt_guess

                return self
            end

    end

    mutable struct ExponentialRosenbrockEulerLejaStruct{T<:AbstractFloat,VT<:AbstractVector{T},MT<:AbstractMatrix{T},MBT<:AbstractMatrix{T},SMT<:AbstractSparseArray{T,<:Integer,2},BD<:Union{Vector{Int64},Nothing},FD<:Union{VT,Nothing},DFD<:Union{VT,Nothing}} <: ImplicitSteppingMethod

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
            D::Diagonal{T,VT}              # Diagonal scaling matrix
            Dinv::Diagonal{T,VT}           # Inverse of diagonal scaling matrix
            ϕ::MT                          # matrix of ϕ functions for exponential Rosenbrock method

            fold::VT                        # local distribution function from previous step 
            fout::VT                        # local distribution function from current step
            fscale::VT                      # scaling vector for exponential Rosenbrock method
            δ::VT                          # temporary vector for exponential Rosenbrock method

            Ks::KrylovSubspace{T,T,T,MT,MT}      # Krylov subspace for exponential Rosenbrock method
            m::Int64                        # dimension of Krylov subspace
            ϕcache::ExponentialUtilities.PhivCache{useview,T} where useview # cache for ϕ functions

            Leja_nodes::Vector{T}                 # Leja nodes for exponential Rosenbrock method

            E::VT                           # energy vector for correcting step 

            function ExponentialRosenbrockEulerLejaStruct(PhaseSpace::PhaseSpaceStruct,Initial::Vector{Float64},Injection::Vector{Float64},BinM::BinaryMatricesStruct,EmiM::EmissionMatricesStruct,FluxM::FluxMatricesStruct;Adaptive::Bool=false,dt_initial::Float64=1.0,n_cut::Float64=1e-45,DistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,DeltaDistributionDomainMask::Union{Vector{Int64},Nothing}=nothing,m::Int64=128)

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
                for species in eachindex(PhaseSpace.name_list)
                    px_num = px_num_list[species]
                    py_num = py_num_list[species]
                    pz_num = pz_num_list[species]
                    dE = dE_list[species]
                    for px in 1:px_num
                        for py in 1:py_num
                            for pz in 1:pz_num
                                idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                                E[idx] = Precision(dE[px])
                            end
                        end
                    end
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

                fstep = zeros(Backend,Precision,length(Initial))
                F = zeros(Backend,Precision,n_momentum)
                J = zeros(Backend,Precision,n_momentum,n_momentum)
                fold = zeros(Backend,Precision,n_momentum)
                fout = zeros(Backend,Precision,n_momentum)
                fscale = zeros(Backend,Precision,n_momentum)
                δ = zeros(Backend,Precision,n_momentum)
                ϕ = zeros(Backend,Precision,n_momentum,2)
                ϕcache = ExponentialUtilities.PhivCache(ϕ,m,1)

                D = Diagonal(one(Precision) ./ copy(E))
                Dinv = Diagonal(copy(E))

                if Backend isa CUDABackend
                    Ks = KrylovSubspace{Precision,Precision,CuArray{Precision,2}}(n_momentum,m)
                else
                    Ks = KrylovSubspace{Precision,Precision,Array{Precision,2}}(n_momentum,m)
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
                    X_Flux = CuSparseMatrixCSC(X_Flux)
                    P_Flux = CuSparseMatrixCSC(P_Flux)
                    invA_Flux = CuArray(invA_Flux)
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

                #= Making invMP = (I-dt*A^{-1}(M_Emi-P_Flux))^{-1}
                  This assuming only emissive interactions coming from one set of specie to another (e.g. electron to photon but not photon to electron) (I-A^{-1}(M_Emi-P_Flux)) has a block triangular structure:
                   I-MP = [A 0]
                          [B C]
                   where A corresponds to the non-emissive species and C corresponds to the emissive species. This means we can invert (I-MP) as:
                   inv(I-MP) = [A^{-1}          0    ]
                               [-C^{-1}BA^{-1} C^{-1}]    

                =#

                ImMP = I - (dt_initial/2)*invA_Flux*(#=M_Emi=# - P_Flux)
                invImMP = spzeros(Precision,size(P_Flux))
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

                        invImMP_view = @view(invImMP[pi_low:pi_up,pi_low:pi_up])
                        ImMP_view = @view(ImMP[pi_low:pi_up,pi_low:pi_up])

                        invImMP_view .= sparse(inv(ImMP_view))
                    end
                    # off-diagonal blocks
                    #=
                        for i>j and Bij = the ij components of the inverse matrix

                            Bij = -Bii \sum_{k=j}^{i-1} Aik Bkj

                    =#
                    for speciesi in eachindex(PhaseSpace.name_list)
                        pi_low = momentum_offset[speciesi] + off_space + 1
                        pi_up = momentum_offset[speciesi+1] + off_space


                        for speciesj in 1:speciesi-1
                            pj_low = momentum_offset[speciesj] + off_space + 1
                            pj_up = momentum_offset[speciesj+1] + off_space

                            
                            Bij = @view(invImMP[pi_low:pi_up,pj_low:pj_up])

                            for speciesk in speciesj:speciesi-1
                                pk_low = momentum_offset[speciesk] + off_space + 1
                                pk_up = momentum_offset[speciesk+1] + off_space

                                Aik = @view(ImMP[pi_low:pi_up,pk_low:pk_up])
                                Bkj = @view(invImMP[pk_low:pk_up,pj_low:pj_up])

                                Bij .-= sparse(Aik * Bkj) 
                            end
                        end

 
                    end
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

                Leja_nodes = leja_nodes_interval(m; a=-2.0, b=2.0)

                # cut initial values that are smaller than n_cut 
                    @. f_init = ifelse(f_init<=n_cut,zero(eltype(f_init)),f_init)
                    @. f = ifelse(f<=n_cut,zero(eltype(f)),f)

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
                self.D = D 
                self.Dinv = Dinv 
                self.fold = fold
                self.fout = fout
                self.fscale = fscale
                self.δ = δ
                self.Ks = Ks
                self.ϕ = ϕ
                self.ϕcache = ϕcache
                self.m = m

                self.Leja_nodes = Leja_nodes

                self.E = E

                self.invImMP = invImMP

                return self
            end

    end