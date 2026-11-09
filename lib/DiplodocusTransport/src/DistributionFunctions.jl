"""
    MaxwellJuttner_Distribution(PhaseSpace,species,T;n=1.0)

Generates a Maxwell-Juttner distribution in momentum space: 
```math
f(\\vec{p})\\mathrm{d}^3\\vec{p} = f(p,u,ϕ)\\mathrm{d}p\\mathrm{d}u\\mathrm{d}ϕ = \\frac{p^2}{4π m^3c^3θK_2(1/θ)}e^{-γ(p)/θ}\\mathrm{d}p\\mathrm{d}u\\mathrm{d}ϕ
```
where ``γ(p) = \\sqrt{1+(p/mc)^2}`` and ``1/θ = m c^2/(k_B T)``. 

As f(p,u,ϕ)ΔpΔuΔϕ is the normalised distribution function in a single momentum-space bin. The value of the distribution function in each bin is given by the mean value of the Maxwell-Juttner distribution at the mean momentum of the bin. 

"""
function MaxwellJuttner_Distribution(PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;n::Float64=1e0)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    species_index = findfirst(x->x==species,name_list)

    dp = Grids.dpx_list[species_index]
    #du = Grids.dpy_list[species_index] # angle averaged therefore not used and dudh = 4π
    #dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    mass = Grids.mass_list[species_index]

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    MJPoints = zeros(Float64,size(meanp))
    # besselk doesn't do well with small T besselkx = besselk * e^x does better.
    invtheta = mEle*c^2/(kb*T)
    @. MJPoints = (meanp^2)*(invtheta/(2*mass^2)) * (1/besselkx(2,mass*invtheta)) * exp((mass-sqrt(mass^2+meanp^2))*invtheta) * dp

    num = sum(MJPoints)
    MJPoints = MJPoints ./ num .*n # scale to number density

    return MJPoints

end

function BlackBody_Distribution(PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;n::Float64=1e0)
    # Generates the height of the BB distribution at positions of the mean momentum per bin

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    species_index = findfirst(x->x==species,name_list)

    dp = Grids.dpx_list[species_index]
    #du = Grids.dpy_list[species_index] # angle averaged therefore not used and dudh = 4π
    #dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    mass = Grids.mass_list[species_index]

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23
    h = 6.62607015e-34

    BBPoints = zeros(Float64,size(meanp))
    invtheta = mEle*c^2/(kb*T)
    @. BBPoints = (meanp^2) / (exp(meanp*invtheta)-1) * dp

    num = sum(BBPoints)
    BBPoints = BBPoints ./ num .*n # scale to number density

    return BBPoints

end





### ============================================================ ###
# inline Distribution functions for different distributions in the form f(vec{p}) for integration with trapezium rule

"""
    Distribution_MaxwellJuttner(px,py,pz,T,m)

Generates a Maxwell-Juttner distribution f(\vec{p}): 
```math
f(\\vec{p})} = \\frac{1}{4π m^3c^3θK_2(1/θ)}e^{-γ(p)/θ}
```
where ``γ(p) = \\sqrt{1+(p/mc)^2}`` and ``1/θ = m c^2/(k_B T)``, and f is normalised by (mEle c)^3
"""
@inline function Distribution_MaxwellJuttner(px::Float64,py::Float64,pz::Float64,T::Float64,m::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0)

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    if (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        # calculate only if within angular bounds
        # besselk doesn't do well with small T besselkx = besselk * e^x does better.
        invtheta = m*mEle*c^2/(kb*T)
        return (invtheta/(4*pi*m^3)) * (1/besselkx(2,invtheta)) * exp((1-sqrt(1+(px/m)^2))*invtheta)
    end

    return 0.0

end

"""
    Distribution_KappaExpDecay(px,py,pz,κ,p_max,m)

Generates a Kappa distribution with exponential decay f(\vec{p}): 
```math
f(\\vec{p})} = \\left(1+\\frac{\\gamma(p)-1}{κθ}\\right)^{-(κ+1)}e^{-E/Emax}
```
where ``γ(p) = \\sqrt{1+(p/mc)^2}`` and ``1/θ = m c^2/(k_B T)``. 
"""
@inline function Distribution_KappaExpDecay(px::Float64,py::Float64,pz::Float64,κ::Float64,T::Float64,p_max::Float64,m::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0,p_cut::Float64=1e2*p_max)

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    if px <= p_cut && (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        # added 1e2 above is cut for exponential tail
        E = sqrt(m^2 + px^2)
        θ = (kb*T)/(m*mEle*c^2)
        return (1+(E/m-1)/(κ*θ))^(-κ-1) * exp(-E/sqrt(m^2+p_max^2)) 
    end

    return 0.0

end

"""
    Distribution_Kappa(px,py,pz,κ,p_max,m)

Generates a Kappa distribution f(\vec{p}): 
```math
f(\\vec{p})} = \\left(1+\\frac{\\gamma(p)-1}{κθ}\\right)^{-(κ+1)}e^{-E/Emax}
```
where ``γ(p) = \\sqrt{1+(p/mc)^2}`` and ``1/θ = m c^2/(k_B T)``. 
"""
@inline function Distribution_Kappa(px::Float64,py::Float64,pz::Float64,κ::Float64,T::Float64,p_max::Float64,m::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0)

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    if px <= p_max && (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        # added 1e2 above is cut for exponential tail
        E = sqrt(m^2 + px^2)
        θ = (kb*T)/(m*mEle*c^2)
        return (1+(E-1)/(κ*θ))^(-κ-1) 
    end

    return 0.0

end

@inline function Distribution_BlackBody(px::Float64,py::Float64,pz::Float64,T::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0)
    # Generates the height of the BB distribution at positions of the mean momentum per bin

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23
    
    if (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        # calculate only if within angular bounds
        # exp(x)-1 = expm1(x) does better for small x
        invtheta = mEle*c^2/(kb*T)
        return 1.0 / (expm1(px*invtheta))
    end

    return 0.0

end

"""
    Distribution_PowerLaw(px,py,pz,index,p_min,p_max,m)

Generates an isotropic Power-law distribution dN/dE ∝ E^(-index) ∴ f(\vec{p}): 
```math
f(\\vec{p})} \\propto \\frac{E^{(-index-1)}}{4π p}
```
where ``E(p) = \\sqrt{m^2+p^2}``. This power law is taken to start at momentum `p_min` and end at `p_max`.
"""
@inline function Distribution_PowerLaw(px::Float64,py::Float64,pz::Float64,index::Float64,p_min::Float64,p_max::Float64,m::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0)

    if (px >= p_min && px <= p_max) && (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        return (sqrt(m^2 + px^2))^(-index - 1) / (4 * pi * px)
    end

    return 0.0

end

"""
    Distribution_PowerLawExpDecay(px,py,pz,index,p_min,p_max,m)

Generates an isotropic Power-law distribution with exponential decay dN/dE ∝ E^(-index)exp(-E/Emax) ∴ f(\vec{p}): 
```math
f(\\vec{p})} \\propto \\frac{E^{(-index-1)}}{4π p}\\exp(-E/Emax)
```
where ``E(p) = \\sqrt{m^2+p^2}``. This power law is taken to start at momentum `p_min` and end at `p_max`.
"""
@inline function Distribution_PowerLawExpDecay(px::Float64,py::Float64,pz::Float64,index::Float64,p_min::Float64,p_max::Float64,m::Float64;umin::Float64=-1.0,umax::Float64=1.0,hmin::Float64=0.0,hmax::Float64=2.0)

    if (px >= p_min && px <= p_max*1e2) && (py >= umin && py <= umax) && (pz >= hmin*pi && pz <= hmax*pi)
        # added 1e2 above is cut for exponential tail
        E = sqrt(m^2 + px^2)
        return E^(-index - 1) / (4 * pi * px) * exp(-E/sqrt(m^2+p_max^2)) # here we take Emax = pmax for simplicity, as for highly relativistic particles E ~ p
    end

    return 0.0

end

"""
    Distribution_Boosted(px,py,pz,Gamma,DistributionFunction,Parameters...)

Generates an isotropic Power-law distribution dN/dE ∝ E^(-index) ∴ f(\vec{p}): 
```math
f(\\vec{p})} \\propto \\frac{E^{(-index-1)}}{4π p}
```
where ``E(p) = \\sqrt{m^2+p^2}``, and the power law is taken to start at momentum `p_min` and end at `p_max`. In one frame in which it is at rest, and boosts it by `Gamma` into the observer frame where it is moving along the local z-direction in momentum space.
"""
@inline function Distribution_Boosted(px::Float64,py::Float64,pz::Float64,Gamma::Float64,m::Float64,DistributionFunction,Parameters...;kwargs...)

    cw::Float64 = Gamma
    sw::Float64 = sqrt(Gamma^2 - 1.0)

    E::Float64 = sqrt(px^2 + m^2)

    E_prime::Float64 = cw*E - sw*px*py
    pxpy_prime::Float64 = cw*px*py - sw*E  # p*u in boosted frame
    px_prime::Float64 = sqrt(E_prime^2 - m^2)
    py_prime::Float64 = pxpy_prime / px_prime # u in boosted frame
    # pz unaffected by boost

    return DistributionFunction(px_prime,py_prime,pz,Parameters...;kwargs...)

end


#### OLD METHODS

#=
"""
    DistributionToDIP_TrapeziumIntegration(f::Array{Float64,3},DistributionFunction::Function,pxr::Vector{Float64},pyr::Vector{Float64},pzr::Vector{Float64},Parameters...)

Integrates a distribution function over momentum space using trapezium integration to generate the discrete form of the distribution function used in DiplodocusTransport.jl. By default this uses 10 samples for each direction in momentum space per bin, but this can be increased using the `samples` keyword argument. `Parameters...` are any additional parameters required by the `DistributionFunction`.

"""
function DistributionToDIP_TrapeziumIntegration!(f::Array{Float64,3},pxr::Vector{Float64},pyr::Vector{Float64},pzr::Vector{Float64},DistributionFunction::Function,Parameters...;kwargs...)

    # Integrates a distribution function over momentum space using trapezium integration
    # f is the 3D array to be filled with the distribution function values
    # DistributionFunction is a function that takes Parameters... and returns the distribution function f(vec{p})

    samples::Int64 = 16 # number of samples per direction in momentum space per bin

    for px in 1:length(pxr)-1, py in 1:length(pyr)-1, pz in 1:length(pzr)-1

        xr = range(pxr[px],pxr[px+1],length=samples)
        yr = range(pyr[py],pyr[py+1],length=samples)    
        zr = range(pzr[pz],pzr[pz+1],length=samples)

        pxp = pxr[px+1]
        pxm = pxr[px]
        pyp = pyr[py+1]
        pym = pyr[py]
        pzp = pzr[pz+1]
        pzm = pzr[pz]

        for (ix,x) in enumerate(xr), (iy,y) in enumerate(yr), (iz,z) in enumerate(zr)
    
            if  (ix == samples || ix == 1) && (iy == samples || iy == 1) && (iz == samples || iz == 1) # corners 
                weight = 0.125 # 1/8
                #println("$ix , $iy , $iz, corner")
            elseif   (!(ix == samples || ix == 1) && (iy == samples || iy == 1) && (iz == samples || iz == 1)) || ((ix == samples || ix == 1) && !(iy == samples || iy == 1) && (iz == samples || iz == 1)) || ((ix == samples || ix == 1) && (iy == samples || iy == 1) && !(iz == samples || iz == 1)) # bounding edge
                weight = 0.25 # 1/4
                #println("$ix , $iy , $iz, edge")
            elseif ((ix == samples || ix == 1) && !(iy == samples || iy == 1) && !(iz == samples || iz == 1)) || (!(ix == samples || ix == 1) && (iy == samples || iy == 1) && !(iz == samples || iz == 1)) || (!(ix == samples || ix == 1) && !(iy == samples || iy == 1) && (iz == samples || iz == 1)) # bounding surface 
                weight = 0.5 # 1/2
                #println("$ix , $iy , $iz, surface")
            elseif (!(ix == samples || ix == 1) && !(iy == samples || iy == 1) && !(iz == samples || iz == 1)) # interior points
                weight = 1.0
                #println("$ix , $iy , $iz, interior")
            else
                error("Error in trapezium integration weight calculation")
            end

            dis = DistributionFunction(x,y,z,Parameters...;kwargs...)

            f[px,py,pz] += weight*dis * x^2 # x^2 from spherical coords

        end

        f[px,py,pz] *= (pxp - pxm)*(pyp - pym)*(pzp - pzm) / ((samples-1)^3)

    end

    if sum(f) == 0.0
        error("Distribution function integration revealed no non-zero values. Check parameters.")
    end

    return nothing

end
=#