"""
    InitialiseInitialCondition(PhaseSpace)

Returns a zero vector with elements for the initial conditions of all particles at all positions in phase space. 
"""
function InitialiseInitialCondition(PhaseSpace::PhaseSpaceStruct)
    
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
    InitialPowerLaw!(Initial,PhaseSpace,species,index,pmin,pmax;kwargs...)

Modifies the initial state vector `Initial` with a power law distribution with `index` for `species`. A power-law distribution is typically defined by N(E) ∝ E^(-index). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index) * p/E = pE^(-index-1). Averaging this over a cell gives f(p)_avg = [E^(1-index)/(1-index)]/[p] where [] denote evaluation at the cell bounds.
"""
function InitialPowerLaw!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,index::Float64,pmin::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::Float64=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32)  where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_PowerLaw,index,pmin,pmax,mass,umin=umin,umax=umax,hmin=hmin,hmax=hmax)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    InitialPowerLawExpDecay!(Initial,PhaseSpace,species,index,pmin,pmax;kwargs...)

Modifies the initial state vector `Initial` with a power law distribution with `index` for `species`, with an exponential cut-off. A power-law distribution is typically defined by N(E) ∝ E^(-index)exp(-E/Emax). N(E) = f(E) therefore f(p) for a power-law distribution is given by f(p) = f(E)*dE/dp = E^(-index)exp(-E/Emax) * p/E = pE^(-index-1)exp(-E/Emax).
"""
function InitialPowerLawExpDecay!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,index::Float64,pmin::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::F=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32)  where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_PowerLawExpDecay,index,pmin,pmax,mass,umin=umin,umax=umax,hmin=hmin,hmax=hmax)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    InitialBoostedPowerLaw!(Initial,PhaseSpace,species,Gamma,index,pmin,pmax;kwargs...)

Takes an isotropic power-law distribution, with minimum momentum `pmin`, maximum momentum `pmax` and `index` in some frame propagating with Lorentz factor `Gamma` in the local z-direction and modifies the initial state vector (distribution), with a number density of `num_Init`.
"""
function InitialBoostedPowerLaw!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,Gamma::Float64,index::Float64,pmin::Float64,pmax::Float64;num_Init::Float64=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    Grids = PhaseSpace.Grids
    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    mass = Grids.mass_list[species_index]
    pu = Float64(p_up_list[species_index])
    pl = Float64(p_low_list[species_index])
    p_num = p_num_list[species_index]
    p_grid = p_grid_list[species_index]
    u_num = u_num_list[species_index]
    u_grid = u_grid_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]

    f0_3D_species = zeros(Float64,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    # Set Rapidity 
    w::Float64 = acosh(Gamma)
    cw = cosh(w)
    sw = sinh(w)
    # de-boost pmin and pmax 
    Emin::Float64 = sqrt(mass^2+pmin^2)
    Emax::Float64 = sqrt(mass^2+pmax^2)
    p1 = sqrt((cw*Emin-sw*pmin)^2-mass^2)
    p2 = sqrt((cw*Emin+sw*pmin)^2-mass^2)
    p3 = sqrt((cw*Emax-sw*pmax)^2-mass^2)
    p4 = sqrt((cw*Emax+sw*pmax)^2-mass^2)
    pmin_UB::Float64 = min(p1,p2,p3,p4)
    pmax_UB::Float64 = max(p1,p2,p3,p4)
    Emin_UB::Float64 = sqrt(mass^2+pmin_UB^2)
    Emax_UB::Float64 = sqrt(mass^2+pmax_UB^2)

    pmin_index = location(pl,pu,p_num,pmin_UB,p_grid)
    pmax_index = location(pl,pu,p_num,pmax_UB,p_grid)

    if pmax_index > p_num
        println("Boost Gamma: ",Gamma," w:",w)
        println("p1,p2,p3,p4: ",p1,", ",p2,", ",p3,", ",p4)
        println("Emin: ",Emin," Emax: ",Emax)
        println("Unboosted pmin: ",pmin_UB," Unboosted pmax: ",pmax_UB)
        println("Unboosted Emin: ",Emin_UB," Unboosted Emax: ",Emax_UB)
        error("Unboosted pmax index exceeds momentum grid size. p_max: ",pmax_UB)
    end

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_Boosted,Gamma,mass,Distribution_PowerLaw,index,pmin,pmax,mass)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    InitialConstant!(Initial,PhaseSpace,species,pmin,pmax,umin,umax,hmin,hmax,num_Init)

Divides the initial number density `num_Init` equally among momentum-space bins in the range of `pmin` to `pmax`, `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values. This is then applies to the initial state vector `Initial`.
"""
function InitialConstant!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String;pmin::T,pmax::T,umin::T=-1.0,umax::T=1.0,hmin::T=0.0,hmax::T=2.0,num_Init::AbstractFloat=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing) where T <: Union{Float32,Float64,Int64} where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum

    name_list = PhaseSpace.name_list
    p_up_list = Momentum.px_up_list
    p_low_list = Momentum.px_low_list
    p_grid_list = Momentum.px_grid_list
    p_num_list = Momentum.px_num_list
    u_grid_list = Momentum.py_grid_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    h_grid_list = Momentum.pz_grid_list

    #Grids = PhaseSpace.Grids
    #dp_list = Grids.dpx_list
    #du_list = Grids.dpy_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    #dp = dp_list[species_index]
    #du = du_list[species_index]
    f0_3D_species = zeros(Float32,p_num_list[species_index],u_num_list[species_index],h_num_list[species_index])

    pu = p_up_list[species_index]
    pl = p_low_list[species_index]
    p_grid = p_grid_list[species_index]
    p_num = p_num_list[species_index]
    u_grid = u_grid_list[species_index]
    u_num = u_num_list[species_index]
    h_num = h_num_list[species_index]
    h_grid = h_grid_list[species_index]

    type = zero(T)
    if (typeof(type) == Float32) || (typeof(type) == Float64) 
        pmin_index = location(pl,pu,p_num,Float64(pmin),p_grid)
        pmax_index = location(pl,pu,p_num,Float64(pmax),p_grid)
        umin_index = location(CONST_u0,CONST_u1,u_num,Float64(umin),u_grid)
        umax_index = location(CONST_u0,CONST_u1,u_num,Float64(umax),u_grid)
        hmin_index = location(CONST_h0,CONST_h1,h_num,Float64(hmin),h_grid)
        hmax_index = location(CONST_h0,CONST_h1,h_num,Float64(hmax),h_grid)
    elseif typeof(type)==Int64
        pmin_index = pmin
        pmax_index = pmax
        umin_index = umin
        umax_index = umax
        hmin_index = hmin
        hmax_index = hmax
    end

    # set values and normalise to initial number density (in m^{-3})
    for px in pmin_index:pmax_index, py in umin_index:umax_index, pz in hmin_index:hmax_index 
        f0_3D_species[px,py,pz] = 1e0
    end

    num = sum(f0_3D_species)
    f0_3D_species *= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    InitialMaxwellJuttner(PhaseSpace,species,T,umin,umax,hmin,hmax,num_Init)

Modeifies the initial state vector `Initial` with a Maxwell-Juttner distribution for `species` of temperature `T` in Kelvin with a number density of `num_Init` and angular range `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function InitialMaxwellJuttner!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::AbstractFloat=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_MaxwellJuttner,T,mass,umin=umin,umax=umax,hmin=hmin,hmax=hmax)
    
    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing

end

"""
    InitialBlackBody(PhaseSpace,species,T,umin,umax,hmin,hmax,num_Init)

Modeifies the initial state vector `Initial` with a Black-Body distribution for `species` of temperature `T` in Kelvin with a number density of `num_Init` and angular range `umin` to `umax` and `hmin to hmax`. These ranges may be defined as either grid indices or physical values.
"""
function InitialBlackBody!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::Union{Nothing,Float64}=nothing,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32) where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    CHAR_number_density = Characteristic.CHAR_number_density

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_BlackBody,T,umin=umin,umax=umax,hmin=hmin,hmax=hmax)
    
    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    if isnothing(num_Init) # standard Blackbody distribution which is optically thin
        num_Init = 16pi * 1.202 * (2pi*CONST_kb*T/(CONST_ħ *CONST_c))^3 / CHAR_number_density
    end
    f0_3D_species .*= num_Init/num
    println("Initial number density: ",num_Init," code units")

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing

end

"""
    InitialKappa!(Initial,PhaseSpace,species,κ,T,pmax;kwargs...)

Modifies the initial state vector `Initial` with a kappa distribution with `κ` and temperature `T` for `species`. 
"""
function InitialKappa!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,κ::Float64,T::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::F=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32)  where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_Kappa,κ,T,pmax,mass,umin=umin,umax=umax,hmin=hmin,hmax=hmax)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end

"""
    InitialKappaExpDecay!(Initial,PhaseSpace,species,κ,T,pmax;kwargs...)

Modifies the initial state vector `Initial` with a kappa distribution with `κ` and temperature `T` for `species`, with an exponential cut-off. 
"""
function InitialKappaExpDecay!(Initial::Vector{F},PhaseSpace::PhaseSpaceStruct,species::String,κ::Float64,T::Float64,pmax::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,num_Init::F=1.0,x_idx=nothing,y_idx=nothing,z_idx=nothing,off_space_idx=nothing,method="hcubature",samples=32,p_cut::Float64=1e2*pmax)  where F<:AbstractFloat

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    name_list = PhaseSpace.name_list

    if species == "Pos" && isnothing(findfirst(==("Pos"),name_list))
        species = "Ele"
    end
    species_index = findfirst(==(species),name_list)
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]
    p_num = Momentum.px_num_list[species_index]
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]

    f0_3D_species = zeros(Float64,p_num,u_num,h_num)

    DistributionToDIPIntegration!(f0_3D_species,pr,ur,hr,method,samples,Distribution_KappaExpDecay,κ,T,pmax,mass,umin=umin,umax=umax,hmin=hmin,hmax=hmax,p_cut=p_cut)

    # set values and normalise to initial number density (in m^{-3})
    num = sum(f0_3D_species)
    f0_3D_species .*= num_Init/num

    f0_species = reshape(f0_3D_species,p_num*u_num*h_num)

    if isnothing(off_space_idx) && (isnothing(x_idx) || isnothing(y_idx) || isnothing(z_idx))
        off_space_idx = 0
    end

    if !isnothing(off_space_idx)
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,off_space_idx=off_space_idx)
    else
        Initial_local = LocationSpeciesToStateVector(Initial,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
    end

    Initial_local .+= convert(typeof(Initial),f0_species)

    return nothing
end


