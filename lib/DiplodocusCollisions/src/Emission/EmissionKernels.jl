"""
    SyncKernel(p3v,p1v,m1,z1,B)

Returns the emission rate for a single photon ``p3v`` state emitted by a charged particle in state ``p1v`` with charge ``z1`` relative to the fundamental charge and mass ``m1`` relative to the mass of the electron, in a uniform magnetic field ``B``.
"""
function SyncKernel(p3v::Vector{Float64},p1v::Vector{Float64},m1::Float64,z1::Float64,Ext::Float64)

    # p3 is Photon
    # p1 is Charged Particle
    
    B::Float64 = Ext # B field in Tesla
    n_int::Int128 = 1
    tol::Float64 = 5e-4

    p3::Float64 = p3v[1]
    p1::Float64 = p1v[1]
    t1::Float64 = p1v[4]
    t3::Float64 = p3v[4]
    st1::Float64,ct1::Float64 = sincospi(t1)
    st3::Float64,ct3::Float64 = sincospi(t3)

    E1::Float64 = sqrt(p1^2 + m1^2)

    if st3 == 0.0
        return 0.0 # to avoid division by zero
    end

    Jfactor1 = (E1*ct3-p1*ct1)/(st3) # code breaks if st3 = 0 FIX
    Jfactor2 = p1*st1

    n::Float64 = abs((mEle^2*c^2)/(z1*ħ*q*B)) * p3 * (E1-p1*ct3*ct1)
    if n < 0.5/tol
        n_int = round(Int128,n)
    end

    y = p1 * st1 *st3 / (E1-p1*ct1*ct3) # y=x/n
    #println(y)

    # characteristic frequency
    #ω0 = abs((z1*q*B))/(E1*mEle)
    #println("critical photon momentum: "*string(ħ*ω0/(mEle*c^2)*E1^3))   
        
    if n > 1e3 # || n > 1e6 && 1-y < 1e-3 # large argument approximation
        # approximation for J's to second order 
        if abs(y) < sqrt(eps(Float64)) # y too close to 1 for numerical precision so calculate z=1-y and assume z is small
            if eps(p1^2) > 0.01m1^2 # small E1-p1
                if abs(t1-t3) < 100sqrt(eps(Float64)) # t1 ≈ t3 
                    e = (m1^2/p1^2)/st3^2 - (m1^2/p1^2)*(t1-t3)*ct3/st3^2+(1/st3^2+(m1^2/p1^2)*(ct3^2-2)/(2st3^4))*((t1-t3)^2)
                else
                    e = 1.0 - (st1^2*st3^2)/(1-ct1*ct3)^2 - (m1^2/p1^2)*st1^2*st3^2/(1.0 - ct1*ct3)^3
                end
            else # t1 ≈ t3 
                e = 1.0 - (p1^2*st3^4)/(E1-p1*ct3^2)^2 - 2*(t1-t3)*(E1-p1)*ct3*st3^3/(E1-p1*ct3^2)^3
            end
        else
            e = 1-y^2
        end
        # Approximation taken from Sokolov et al. 1969
        if e <= 0.0
            #println("p1=$p1, E1=$E1, ct1=$ct1, ct3=$ct3")
            #println("e=$e, y=$y, n=$n")
            #@warn "y value is greater than 1, invalid for synchrotron emission, skipping this point."
            return 0.0
        end
        K13 = besselk(1/3,n*e^(3/2)/3)
        K23 = besselk(2/3,n*e^(3/2)/3)
        J1 = ((sqrt(e))/(pi*sqrt(3)))*(K13 +(e/10)*(K13-2*n*e^(3/2)*K23))
        J2 = (e/(pi*sqrt(3)))*(K23 + (e/5)*(2*K23-(1/(e^(3/2)*n)+n*e^(3/2))*K13))    
    elseif n_int >= 1 && abs(n-n_int)/n_int < tol # for tol = 5e-4 this always true for n > 1e3 
        #J1 = (n*y/2)^n/Bessels.Gamma(n+1)
        #J2 = (1/2)*(n*y/2)^(n-1)/Bessels.Gamma(n)
        J1 = besselj(n_int,n_int*y)
        J2 = 1/2 * (besselj(n_int-1,n_int*y) - besselj(n_int+1,n_int*y))
    else
        J1 = 0e0
        J2 = 0e0
    end

    val = (p3/E1)*((Jfactor1*J1)^2+(Jfactor2*J2)^2)
    #println(val)

    factor = (abs(z1/B))*(3*c^4*mEle^5)/(4*pi*ħ^3*μ0*q^3) # synchrotron emission rate divided by c*σT

    #println(factor)

    return val*factor
    
end

#SyncKernel([1e-14,0.6],[1e1,0.5],1,1,1e-4)
#n = 7.36e10
#e = 1-(0.0045)^2
#besselk(1/3,n*e^(3/2)/3)