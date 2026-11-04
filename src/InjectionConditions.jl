"""
    InitialiseInjectionCondition(PhaseSpace)

Returns a zero vector with elements for the injection conditions of all particles at all positions in phase space. 
"""
function InitialiseInjectionCondition(PhaseSpace::PhaseSpaceStruct)

    px_num_list = PhaseSpace.Momentum.px_num_list
    py_num_list = PhaseSpace.Momentum.py_num_list
    pz_num_list = PhaseSpace.Momentum.pz_num_list
    x_num = PhaseSpace.Spacetime.x_num
    y_num = PhaseSpace.Spacetime.y_num
    z_num = PhaseSpace.Spacetime.z_num

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*x_num*y_num*z_num

    return zeros(Float64,m)

end


"""
    InjectionPowerLaw!(Injection,PhaseSpace,species,index,pmin,pmax;kwargs...)

Modifies the injection state vector `Injection` with a power law distribution generated using `InitialPowerLaw!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionPowerLaw!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,index::Float64,pmin::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::AbstractFloat=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end
    
    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialPowerLaw!(tmp,PhaseSpace,species,index,pmin,pmax;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method="hcubature",samples=32)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end


"""
    InjectionPowerLawExpDecay!(Injection,PhaseSpace,species,index,pmin,pmax;kwargs...)

Modifies the injection state vector `Injection` with a power law distribution with exponential cut-off generated using `InitialPowerLawExpDecay!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionPowerLawExpDecay!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,index::Float64,pmin::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::AbstractFloat=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end
    
    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialPowerLawExpDecay!(tmp,PhaseSpace,species,index,pmin,pmax;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method="hcubature",samples=32)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end

"""
    InjectionBoostedPowerLaw!(Injection,PhaseSpace,species,Gamma,index,pmin,pmax;kwargs...)

Modifies the injection state vector `Injection` with a power law distribution generated using `InitialBoostedPowerLaw!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionBoostedPowerLaw!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,Gamma::Float64,index::Float64,pmin::Float64,pmax::Float64;num_Inj::Float64=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end

     # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialBoostedPowerLaw!(tmp,PhaseSpace,species,Gamma,index,pmin,pmax;num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method="hcubature",samples=32)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing

end

"""
    InjectionConstant!(Injection,PhaseSpace,species,pmin,pmax,umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a constant distribution generated using `InitialConstant!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionConstant!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::T,pmax::T,umin::T,umax::T,hmin::T,hmax::T,num_Inj::AbstractFloat=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing) where T <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialConstant!(tmp,PhaseSpace,species;pmin=pmin,pmax=pmax,umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end

"""
    InjectionMaxwellJuttner!(Injection,PhaseSpace,species,T;umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a Maxwell-Juttner distribution generated using `InitialMaxwellJuttner!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionMaxwellJuttner!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::Float64=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialMaxwellJuttner!(tmp,PhaseSpace,species,T;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method=method,samples=samples)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end

"""
    InjectionBlackBody!(Injection,PhaseSpace,species,T;umin,umax,hmin,hmax,num_Inj,rate_Inj)

Modifies the injection state vector `Injection` with a Black-Body distribution generated using `InitialBlackBody!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionBlackBody!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::Union{Nothing,Float64}=nothing,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialBlackBody!(tmp,PhaseSpace,species,T;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method=method,samples=samples)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    println("sum injection: ",sum(rate_Inj .* tmp))

    return nothing
end

"""
    InjectionKappaExpDecay!(Injection,PhaseSpace,species,κ,T,pmax;kwargs...)

Modifies the injection state vector `Injection` with a kappa distribution generated using `InitialKappaExpDecay!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionKappa!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,κ::Float64,T::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::AbstractFloat=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end
    
    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialKappa!(tmp,PhaseSpace,species,κ,T,pmax;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method="hcubature",samples=32)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end


"""
    InjectionKappaExpDecay!(Injection,PhaseSpace,species,κ,T,pmax;kwargs...)

Modifies the injection state vector `Injection` with a kappa distribution with exponential cut-off generated using `InitialKappaExpDecay!` scaled by the rate of injection `rate_Inj`. Such that particles with that distribution and number density `num_Init` are injected at a rate of `rate_Inj`.
"""
function InjectionKappaExpDecay!(Injection::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,κ::Float64,T::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Inj::AbstractFloat=1.0,rate_Inj::Union{Float64,Vector{Float64}}=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32,p_cut::Float64=1e2*pmax) where F<:AbstractFloat
    
    if rate_Inj isa Vector{Float64}
        @assert length(rate_Inj) == length(Injection) "rate_Inj and Injection must have the same length if rate_Inj is a vector"
    end

    # Create temporary vector to hold initial conditions
    tmp = zeros(eltype(Injection), length(Injection))

    # Add initial conditions to temporary vector
    InitialKappaExpDecay!(tmp,PhaseSpace,species,κ,T,pmax;umin=umin,umax=umax,hmin=hmin,hmax=hmax,num_Init=num_Inj,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx,off_space_idx=off_space_idx,method="hcubature",samples=32,p_cut=p_cut)

    #tr = PhaseSpace.Grids.tr
    #dt0 = tr[2] - tr[1]
    # dt0 now taken to be 1.0 by default
    # Scale by rate (rate scaled by dt0 to convert to per time step)
    @. Injection += rate_Inj * tmp

    return nothing
end