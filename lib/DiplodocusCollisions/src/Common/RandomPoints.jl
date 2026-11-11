#=
Contains functions for generating random points
=#

"""
    RPointLogMomentum!(pu,pl,pv,num)

Edits the first element of `pv` with a random real-space momentum value between ``10^{pl}`` and ``10^{pu}``. This sample is chosen by first randomly picking a momentum bin in the range `1:num` and then uniformly sampling a momentum point in real-space (rather than log10 space) between pl and pu which are the momentum values at start and end of that bin. Sampling is done such there will be a constant number of points per momentum-space volume. As the momentum space between ``10^{pl}`` and ``10^{pu}`` it is a spherical shell hence the correct sampling is ``p = (U*(10^{pu})^3+(1-U)*(10^{pl})^3)^{1/3}`` with uniform ``U ∈ [0~~1]``.

Assumes ``f(x,p,μ)=f(x,\\vec{p})*(2πp^2)=const`` in bin, therefore momentum space volume element is ``\\mathrm{d}p`` and as such uniform sampling corresponds to ``U*10^{u}+(1-U)*10^{l}`` where ``U`` is a uniform random number between 0 and 1.

If instead ``f(x,\\vec{p})=const`` in bin, momentum space volume element is ``p^2 \\mathrm{d}p`` and uniform sampling corresponds to ``(10^pu)*\\sqrt[3]{U+(1-U)*10^{3pl-3pu}}`` where ``U`` is a uniform random number between 0 and 1.
"""
@inline function RPointLogMomentum!(pv::Vector{Float64},pu::Float64,pl::Float64,num::Int64) 
    # Inputs a momentum vector and momentum bounds and mutates first of said vector
    bin = rand(1:num)
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float64)
    #pv[1] = (10^u)*cbrt(U+(1-U)*1f3^(l-u)) 
    pv[1] = U*10^(u)+(1-U)*10^(l)  # if instead want to sample space uniformly.

    return nothing
    
end

@inline function RPointLogMomentum!(pv::Vector{Float64},pu::Float64,pl::Float64,num::Int64,bin::Int64) 
    # Inputs a momentum vector and momentum bounds and mutates first of said vector
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float64)
    #pv[1] = (10^u)*cbrt(U+(1-U)*1f3^(l-u)) 
    pv[1] = U*10^(u)+(1-U)*10^(l)  # if instead want to sample space uniformly.

    return nothing
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (cos(theta)) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi). 
"""
@inline function RPointSphereCosThetaPhi!(a::Vector{Float64}) 
    # Inputs a 4 element vector [p, cos(theta), phi/pi,theta/pi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi points are normalised by pi

    u::Float64 = rand(Float64)
    v::Float64 = rand(Float64)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]
    a[3] = 2*u       # phi bound by [0,2) 
    a[4] = acos(a[2])/pi # theta bound by [0,1]

    return nothing
    
end

"""
    RPointSphereThetaPhiBounds!()

Assigns the second (cos(theta)) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi), within the bounds of ``u_low ≤ u ≤ u_up`` and ``h_low ≤ h ≤ h_up``. 
"""
@inline function RPointSphereCosThetaPhiBounds!(a::Vector{Float64},u_low,u_up,h_low,h_up) 
    # Inputs a 4 element vector [p, cos(theta), phi/pi,theta/pi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi points are normalised by pi

    u::Float64 = rand(Float64)
    v::Float64 = rand(Float64)

    a[2] = u_low + (u_up-u_low)*v     # cos(theta) bound by u_low ≤ u ≤ u_up, where u = cos(theta)
    a[3] = h_low + (h_up-h_low)*u       # phi bound by h_low ≤ h ≤ h_up, normalised by pi 
    a[4] = acos(a[2])/pi # theta bound by [0,1]

    return nothing
    
end

"""
    RPointSphereTheta!()

Assigns the second (cos(theta)) element of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta). 
"""
@inline function RPointSphereCosTheta!(a::Vector{Float64}) 
    # Inputs a 3 element vector [p, cos(theta)] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta)
    # phi points are normalised by pi

    v::Float64 = rand(Float64)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]

    return nothing
    
end

"""
    betaVec!(βv,p1v,p2v)

Mutates the components of the centre of momentum velocity vector `βv` with components `[β,u,phi,γ,γβ] in terms of the incident state vectors `p1v` and `p1v`
"""
@inline function betaVec!(βv::Vector{Float64},p1v,p2v,m1,m2)

    p1 = p1v[1]
    p2 = p2v[1]

    ct1 = p1v[2]
    ct2 = p2v[2]
    st1 = sqrt(1-ct1^2)
    st2 = sqrt(1-ct2^2)

    (sh1,ch1) = sincospi(p1v[3])
    (sh2,ch2) = sincospi(p2v[3])
    ch1h2 = cospi(p1v[3]-p2v[3])

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    βv[1] = (1/(E1+E2)) * sqrt(p1^2 + p2^2 + 2*p1*p2*(ct1*ct2 + ch1h2*st1*st2))

    βv[2] = (1/(E1+E2)) * (p1*ct1 + p2*ct2) / βv[1] # cos(theta) bounded by [-1,1]
    x = p1*st1*ch1 + p2*st2*ch2
    y = p1*st1*sh1 + p2*st2*sh2
    βv[3] = mod2pi(atan(y,x))/pi # atan(y,x), note 1/(E1+E2) is common factor in y,x so ignored. phi bounded by [0,2]

    #gamma
    βv[4] = (E1+E2)
    βv[4] /= sqrt((m1+m2)^2 + 2*E1s*E2s+2*E1s*m2+2*E2s*m1-2*p1*p2*(ct1*ct2+ch1h2*st1*st2))

    # beta gamma
    βv[5] = sqrt(p1^2+p2^2+2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
    βv[5] /= sqrt((m1+m2)^2 + 2*E1s*E2s+2*E1s*m2+2*E2s*m1-2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
    
end

"""
    RPointSphereBoost!(βv,pvCOM)

Takes the random points on the sphere generated in the centre of momentum frame `pCv=[ct_deboosted,ct,h]` and boosts them back to the lab frame using `βv` to modify the lab frame vector `pLv`.
"""
@inline function RPointSphereBoost!(pLv::Vector{Float64},βv::Vector{Float64},pCv::Vector{Float64})

    ctC = pCv[2]
    stC = sqrt(1-ctC^2)
    (shC,chC) = sincospi(pCv[3])

    ctβ = βv[2]
    stβ = sqrt(1-ctβ^2)
    (shβ,chβ) = sincospi(βv[3])

    pLv[2] = ctC*ctβ - stC*stβ*chC
    x = -stC*shC*shβ + chβ*(stC*chC*ctβ + ctC*stβ)
    y = stC*shC*chβ + shβ*(stC*chC*ctβ + ctC*stβ)

    pLv[3] = mod2pi(atan(y,x))/pi

    # deboost ctC to lab frame for probability calculation
    γ = βv[4]
    βγ = βv[5]
    #pCv[1] = (ctC+βv[1]) / (1+βv[1]*ctC)
    pCv[1] = (γ*ctC+βγ) / (γ+βγ*ctC)

end

"""
    pdfBoost(βv,pCv)

returns the probability of sampling a point given by `pCv` dependent on the boost `βv`.
"""
@inline function pdfBoost(βv::Vector{Float64},ctC::Float64)

    # assumes pv in lab frame is valid and un-mirrored from the boosted frame
    β=βv[1]
    γ=βv[4]
    βγ=βv[5]

    #prob = (1-β^2)
    #prob /= 2*(1-β*ctC)^2
    
    prob = (γ^2-βγ^2)
    prob /= 2*(γ-βγ*ctC)^2
    
    return prob

end


"""
    WeightedFactors(p1v,p2v,m1,m2,scale)

Returns the weighting rapidities `w3` and `w4` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. `t` and `h` are the angles of the COM velocity direction, and the weights `w3` and `w4` are rapidities based on the expected angular spread of the particles 3 and 4 due to the incoming states of particles 1 and 2. 
"""
@inline function WeightedFactors(p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64,sBig::Float64,sSmol::Float64,scale::Float64) 

    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    (st1::Float64,ct1::Float64) = sincospi(p1v[4])
    (st2::Float64,ct2::Float64) = sincospi(p2v[4])
    (sh1::Float64,ch1::Float64) = sincospi(p1v[3])
    (sh2::Float64,ch2::Float64) = sincospi(p2v[3])
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    s::Float64 = sBig + sSmol
    E1::Float64 = sqrt(p1v[1]^2+m1^2)
    E2::Float64 = sqrt(p2v[1]^2+m2^2)

    # COM frame vector and angles
    γC::Float64 = (E1+E2)/sqrt(s)
    if γC < 1.0
        γC = 1.0 # avoid numerical issues
    end
    wC::Float64 = acosh(γC)
    if (pSmol::Float64 = p2/p1) < sqrt(eps(Float64))
        #z = sqrt(1e0+pSmol^2+2*pSmol*(ct1*ct2+ch1h2*st1*st2))
        a = (ct1*ct2+ch1h2*st1*st2)
        b = ct1
        c = ct2 
        # 1/sqrt(1+smol^2) approx 1-smol^2/2
        val =  b + (c-a*b)*pSmol + ((-1+3*a^2)*b/2-a*c)*pSmol^2 #+ (3*a*b/2-5*a^3*b/2-c+3*a^2*c/2)*pSmol^3 
        x = st1*ch1 + pSmol*st2*ch2
        y = st1*sh1 + pSmol*st2*sh2
    elseif (pSmol = p1/p2) < sqrt(eps(Float64))
        #z = sqrt(1e0+pSmol^2+2*pSmol*(ct1*ct2+ch1h2*st1*st2))
        a = (ct1*ct2+ch1h2*st1*st2)
        b = ct1
        c = ct2 
        # 1/sqrt(1+smol^2) approx 1-smol^2/2
        val =  c + (b-a*c)*pSmol + ((-1+3*a^2)*c/2-a*b)*pSmol^2 #+ (-b+3*a*c/2-5*a^3*c/2+3*a^2*b/2)*pSmol^3 
        x = pSmol*st1*ch1 + st2*ch2
        y = pSmol*st1*sh1 + st2*sh2
    else
        z = sqrt(p1^2+p2^2+2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
        x = p1*st1*ch1 + p2*st2*ch2
        y = p1*st1*sh1 + p2*st2*sh2
        val  = (p1*ct1+p2*ct2)/z
    end 
    t::Float64 = acos(val)/pi
    h::Float64 = mod2pi(atan(y,x))/pi   

    # outgoing COM frame momentum
    pC::Float64 = InvariantFluxSmall(sSmol,sBig,m3,m4)/sqrt(s)

    # pre allocate types
    w3Limit::Float64 = 0e0
    w4Limit::Float64 = 0e0

    # Lab Frame Angle limits
    # Due to conservation laws there is a limit on the lab frame scattering angle with respect to the COM frame direction. This angle is then mapped to a rapidity that "boosts" the lab frame angle samples such that 50% of samples lie within this angle limit.
    if m3 > m4
        w4Limit = 0e0
        tmp =  pC/(m3*sinh(wC))
        if tmp < 1e0 
            if tmp > sqrt(eps(Float64)) 
                w3Limit = atanh(sqrt(1-tmp^2))
            else # small tmp
                w3Limit = log(2e0)-log(tmp)-tmp^2/4
            end
            #w3Limit = min(w3Limit,18.7e0) # for tmp < 1e-8 sqrt=0 due to float precision, 18.7e0 is maximum value of w3 to this precision
        else
            w3Limit = 0e0
        end
    elseif m3 < m4
        w3Limit = 0e0
        tmp = pC/(m4*sinh(wC))
        if tmp < 1e0
            if tmp > sqrt(eps(Float64)) 
                w4Limit = atanh(sqrt(1-tmp^2))
            else # small tmp
                w4Limit = log(2e0)-log(tmp)-tmp^2/4
            end
            #w4Limit = min(w4Limit,18.7e0)
        else
            w4Limit = 0e0
        end
    else
        w3Limit = 0e0
        w4Limit = 0e0
    end

    sCOMrest = max((m1+m2)^2,(m3+m4)^2)
    wScale::Float64 = asinh(sqrt(sSmol/(sCOMrest))) # rough measure of how energetic the interaction is over the COM frame energy

    #w3::Float64 = min(w3Limit+wC+scale*wScale,18e0)
    #w4::Float64 = min(w4Limit+wC+scale*wScale,18e0)
    if w3Limit != 0e0
        w3 = (scale)*w3Limit #+ scale
    else
        w3 = scale*wC #+scale*wScale
    end
    if w4Limit != 0e0
        w4 = (scale)*w4Limit #+ scale
    else
        w4 = scale*wC #+scale*wScale
    end
    w3 = min(w3,18.0)
    w4 = min(w4,18.0)
    #if w3 == 18.0 || w4 == 18.0
    #    println("Warning: Weighting rapidity capped at 18.0, consider increasing scale or checking kinematics.")
    #end
    #w3::Float64 = w3Limit#+scale*wC #+scale*(wScale)
    #w4::Float64 = w4Limit#+scale*wC #+scale*(wScale) 

    return (w3,w4,t,h)
    
end


"""
    WeightedFactorsEmission!(p1v,m1,scale)

Returns the weighting rapidity `w` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. Weight is dependant on the energy of the emitting particle `p1v`. With angles `t` and `h` being the angles of the emitting particle. 
"""
@inline function WeightedFactorsEmission(p1v::Vector{Float64},m1::Float64,scale::Float64) 

    E1::Float64 = sqrt(p1v[1]^2+m1^2)
    gamma::Float64 = E1/m1

    w = scale*acosh(gamma)
    t = p1v[4]
    h = p1v[3]

    #println("E1 = $E1, gamma = $gamma, w = $w")

    return (w,t,h)
    
end

"""
    RPointSphereWeighted!()

Assigns randomly sampled angles on the sphere (cos(theta) and phi) weighed by a doppler boosting by rapidity `w`, returning the probability `prob` for such a sample and mutating the vector `a` with elements `[p, cos(theta), phi, theta]` with angles normalised by pi. 
"""
@inline function RPointSphereWeighted!(a::Vector{Float64},w::Float64) 
    # Inputs a 5 element vector [p, cos(theta), phi, theta, prob] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html using the doppler boost formula as a weighting with a rapidity `w`. 
    # phi points are normalised by pi
    # prob is then the probability of sampling the random points on the sphere given the doppler boosting.

    #=
        Angle transformations, where u is in lab frame and uB is in the boosted frame (and uniformly sampled with uB = 2v-1):
            u = (cosh(w)*uB+sinh(w))/(cosh(w)+sinh(w)*uB)
            phi = phiB
        This can be re-written as: 
            t = 2*atan(exp(-w)*tan(tB/2))
            where t = acos(u)

        Probability of sampling
            P(u) = 1/(cosh(w)-sinh(w)u)^2 = (cosh(w)+sinh(w)*uB)^2 = Doppler^2 i.e. solid angle element transformation dOmega = Doppler^2 dOmegaB
        Putting uB = 2v-1 for uniform sampling 
            P(u) = (e^w*v+(1-v)/e^w)^2

    =#

    v::Float64 = rand(Float64)
    #tB::Float64 = acos(2*v-1)/pi # "boosted" theta

    #costBdiv2sqr::Float64 = v
    #sintBdiv2sqr::Float64 = 1-v
    tantBdiv2::Float64 = sqrt((1e0-v)/v)


    h::Float64 = 2*rand(Float64) # phi points are normalised by pi
    #t::Float64 = 2*atan(exp(-w)*tanpi(tB/2))/pi  # "unboosted" theta divided by pi
    x::Float64 = exp(-w)*tantBdiv2
    t::Float64 = 2*atan(exp(-w)*tantBdiv2)/pi  # "unboosted" theta divided by pi
    
    #sintdiv2sqr::Float64 = x^2/(1+x^2)
    #costdiv2sqr::Float64 = 1/(1+x^2)

    ew::Float64 = exp(w)
    #prob::Float64 = (ew*sinpi(t/2)^2+cospi(t/2)^2/ew)^-2
    #prob::Float64 = (ew*sintdiv2sqr+costdiv2sqr/ew)^-2
    #prob::Float64 = (ew*cospi(tB/2)^2+sinpi(tB/2)^2/ew)^2
    #prob::Float64 = (ew*costBdiv2sqr+sintBdiv2sqr/ew)^2
    prob::Float64 = (ew*v+(1-v)/ew)^2

    
    a[2] = cospi(t)
    a[3] = h
    a[4] = t

    #=if w >= 7e0
        println("")
        println("w = $w")
        println("t = $t")
        println("tB = $tB")
        println("h = $h")
        println("prob = $prob")
        println("")
    end  =#

    return prob #/ (2*pi)
    
end

"""
    WeightedFactorsSync!()

Returns the weighting rapidity `w` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. 
"""
@inline function WeightedFactorsSync(pv::Vector{Float64},m::Float64,scale::Float64) 
    E1::Float64 = sqrt(pv[1]^2+m^2) # non dimensional units energy is gamma * m

    w = scale*acosh(E1/m) #*rand(Float64)
    t = pv[4]
    h = pv[3]

    return (w,t,h)
    
end