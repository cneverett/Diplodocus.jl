#= Functions for the S and T integration functions =#

"""
    CosThetaValue(p1v,p2v)

Returns the cosine of the angle between two momentum vectors `p1v` and `p2v` of format [p,cos(theta),phi/pi,theta/pi]. To aid in floating point precision, 1 is subtracted from the returned cosine value. i.e. CosTheta12m1 = cos(Theta12) - 1.0. This is done as it is often the case that two momentum states are very close together.
"""
function ThetaValue(p1v::Vector{Float64},p2v::Vector{Float64})

    t1::Float64 = p1v[4]
    h1::Float64 = p1v[3]
    t2::Float64 = p2v[4]
    h2::Float64 = p2v[3]
    
    if (td::Float64=pi*abs(t1-t2)/2) < 1e-6 && (hd::Float64=pi*abs(h1-h2)/2) < 1e-6
        ta::Float64 = (t1 + t2)/2
        CosTheta12m1 = 2*td^2*(hd^2-1)-2*hd^2*sinpi(ta)^2
    else
        CosTheta12m1 = cospi(t1)*cospi(t2) + sinpi(t1)*sinpi(t2)*cospi(h1-h2) - 1e0
    end

    return CosTheta12m1

end

"""
    LossValue(p1v,p2v,sigma,mu1,mu2)

returns `LossVal` based on initial momentum states `p1v` and `p2v` and cross section `sigma` based on particle selection.
```math
Loss_\\text{val} = \\frac{1}{p^0_1p^0_2}\\sigma(s)F_12(s)
```
If initial state fails `sCheck`, i.e. cannot generate a physical output state, LossVal is set to 0e0. 
Assumes f(x,p,u,ϕ)=f(x,vec{p})p^2=constant over bin
"""
function LossValue(p1v::Vector{Float64},p2v::Vector{Float64},sigma::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2

    m12::Float64 = m1^2
    m22::Float64 = m2^2

    LossVal::Float64 = 0e0
    
    #= This method leads to errors when p/m < 10^-6 due to Floating point precision, causing s<(m+m)^2 and InvFlux to return a complex number error 
    E1 = sqrt(p1^2+m12)
    E2 = sqrt(p2^2+m22) 
    s = m12+m22 + 2*E1*E2 - 2*p1*p2*(ct1*ct2+ch1h2*st1*st2)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...)) =#

    Es1::Float64 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es1s::Float64 = Es1/p1
    Es2::Float64 = (p2^2)/(sqrt(m22+p2^2)+m2)
    Es2s::Float64 = Es2/p2

    sBig::Float64 = max((m1+m2)^2,(m3+m4)^2)
    #sSmol::Float64 = ((m1+m2)^2-sBig) + 2*p1*p2*(-ctheta12 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)
    sSmol::Float64 = 0.0
    if m1 == 0.0 && m2 == 0.0
        sSmol += 2*p1*p2*(1.0-ctheta12) - sBig
    elseif m1 == 0.0 && m2 != 0.0
        sSmol += m2^2 + 2*p1*m22/(sqrt(m22+p2^2)+p2) + 2*p1*p2*(1-ctheta12) - sBig
    elseif m1 != 0.0 && m2 == 0.0
        sSmol += m1^2 + 2*p2*m12/(sqrt(m12+p1^2)+p1) + 2*p2*p1*(1-ctheta12) - sBig
    else 
        sSmol += m1^2 + m2^2 + 2*(p1*m22/(sqrt(m22+p2^2)+p2) + p2*m12/(sqrt(m12+p1^2)+p1) +m22/(sqrt(m22+p2^2)+p2)*m12/(sqrt(m12+p1^2)+p1)) + 2*p1*p2*(1-ctheta12) - sBig
    end
   

    if sCheck(sSmol,sBig,m1,m2,m3,m4) # check if s value is valid for interaction

        E1::Float64 = Es1 + m1
        E2::Float64 = Es2 + m2
        
        LossVal = (1/E1)*(1/E2)*(InvariantFluxSmall(sSmol,sBig,m1,m2))*sigma(sSmol,sBig)
        if LossVal==Inf || LossVal < 0e0
            println("")
            println("p1v = $p1v")
            println("p2v = $p2v")
            println("LossVal = $LossVal")
            println("sSmol = $sSmol")
            println("sBig = $sBig")
            error("ST1 Inf or -ve#")
        end

    else # if not valid set T value to zero

        LossVal = 0e0
        
    end

    return LossVal, sBig, sSmol

end


"""
    GainValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)

Returns `GainVal` based on initial momentum states `p1v` and `p2v` and final state `p3v` and differential cross section `dsigmadt` based on particle selection 12->34. 
```math
Gain_\\text{val}=\\frac{\\mathrm{d}\\sigma_{12|34}}{\\mathrm{d}t}\\frac{\\mathcal{F}_{12}^2}{\\pi\\left|p^0_3(p_1\\cos\\Theta_{31}+p_2\\cos\\Theta_{32})-p_3(p^0_1+p^0_2)\\right|}\frac{p_3^2}{p^0_1p^0_2}.
``` 
Assumes f(x,p,u,ϕ)=f(x,vec{p})p^2=constant over bin
"""
function GainValue3(p3v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    # cos(A-B) = 1.0 - b, so "a" part is 1.0, "b" part is the cos(A-B)-a
    b::Float64 = pi^2*(p1v[3]-p2v[3])^2/2
    ch1h2a::Float64, ch1h2b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p1v[3]-p2v[3])-1.0)
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    p3::Float64 = p3v[1]
    st3::Float64,ct3::Float64 = sincospi(p3v[4])
    b = pi^2*(p3v[3]-p1v[3])^2/2
    ch3h1a::Float64, ch3h1b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p3v[3]-p1v[3])-1.0)
    b = pi^2*(p3v[3]-p2v[3])^2/2
    ch3h2a::Float64, ch3h2b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p3v[3]-p2v[3])-1.0)
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2
    
    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float64 = Es3/p3
    E3::Float64 = Es3 + m3

    b = pi^2*(p1v[4]-p2v[4])^2/2
    ctheta12a::Float64, ctheta12b::Float64, ctheta12c::Float64 =  b < eps(Float64) ? (1.0, -b, ch1h2b*st1*st2) : (1.0,cospi(p1v[4]-p2v[4])-1.0, ch1h2b*st1*st2)
    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2

    b = pi^2*(p3v[4]-p1v[4])^2/2
    ctheta13a::Float64, ctheta13b::Float64, ctheta13c::Float64 = b < eps(Float64) ? (1.0, -b, ch3h1b*st3*st1) : (1.0,cospi(p3v[4]-p1v[4])-1.0, ch3h1b*st3*st1)
    ctheta13::Float64 = ct3*ct1+ch3h1*st3*st1

    b = pi^2*(p3v[4]-p2v[4])^2/2
    ctheta23a::Float64, ctheta23b::Float64, ctheta23c::Float64 = b < eps(Float64) ? (1.0, -b, ch3h2b*st3*st2) : (1.0,cospi(p3v[4]-p2v[4])-1.0, ch3h2b*st3*st2)
    ctheta23::Float64 = ct3*ct2+ch3h2*st3*st2
    
    # TODO: mak this more accurate for large separations of p1 and p3 consistent with the new definitions of ctheta a,b,c above
    deltacorrect::Float64 = Es1*p3 - Es3*p1*ctheta13
    deltacorrect += m1*p3 - m3*p1*ctheta13
    deltacorrect += Es2*p3 - Es3*p2*ctheta23    
    deltacorrect += m2*p3 - m3*p2*ctheta23

    # t = tBig + tSmol
    tBig::Float64 = (m3-m1)^2
    #tSmol::Float64 = -2*(m1*Es3 + m3*Es1 + Es3*Es1 - p3*p1*(ct3*ct1+ch3h1*st3*st1))
    #tSmol::Float64 = -2*m1*m3 - 2*E1*E3 + 2*p1*p3*ctheta13
    #tSmol::Float64 = 2*p3*p1*(ctheta13 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3)
    tSmol::Float64 = 0.0
    if m1 == 0.0 && m3 == 0.0
        tSmol += 2*p1*p3*(ctheta13b + ctheta13c)
    elseif m1 == 0.0 && m3 != 0.0
        #E = sqrt(p^2+m^2) = p + (m^2)/(sqrt(m^2+p^2)+p)
        #E-p = (m^2)/(sqrt(m^2+p^2)+p) 
        tSmol += -2*(m32/(sqrt(m32+p3^2)+p3))*p1 + 2*p3*p1*(ctheta13b + ctheta13c)
    elseif m1 != 0.0 && m3 == 0.0
        tSmol += -2*(m12/(sqrt(m12+p1^2)+p1))*p3 + 2*p3*p1*(ctheta13b + ctheta13c)
    else 
        tSmol += -2*(E1*E3) + 2*(p1*p3+m1*m3) + 2*p1*p3*(ctheta13b + ctheta13c)
    end


    # u = uBig + uSmol
    uBig::Float64 = (m2-m3)^2
    #uSmol::Float64 = m12+m22+m32+m42 - sBig - tBig - uBig - sSmol - tSmol  # this leads to Float64 issues better to calculate directly
    #uSmol::Float64 = -2*(m3*Es2 + m2*Es3 + Es2*Es3 - p2*p3*(ct2*ct3+ch3h2*st2*st3))
    #uSmol::Float64 = 2*p2*p3*(ctheta23 - Es2s*Es3s - m3*Es2s/p3 - m2*Es3s/p2)
    uSmol::Float64 = 0.0
    if m2 == 0.0 && m3 == 0.0
        uSmol += 2*p2*p3*(ctheta23b + ctheta23c)
    elseif m2 == 0.0 && m3 != 0.0
        #E = sqrt(p^2+m^2) = p + (m^2)/(sqrt(m^2+p^2)+p)
        #E-p = (m^2)/(sqrt(m^2+p^2)+p) 
        uSmol += -2*(m32/(sqrt(m32+p3^2)+p3))*p2 + 2*p3*p2*(ctheta23b + ctheta23c)
    elseif m2 != 0.0 && m3 == 0.0
        uSmol += -2*(m22/(sqrt(m22+p2^2)+p2))*p3 + 2*p3*p2*(ctheta23b + ctheta23c)
    else 
        uSmol += -2*(E2*E3) + 2*(p2*p3+m2*m3) + 2*p3*p2*(ctheta23b + ctheta23c)
    end

    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,sBig,m1,m2))/pi

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("stu check")
    end=#

    GainVal = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p3^2/(deltacorrect*sign(deltacorrect)))

    if (GainVal==Inf || GainVal == -Inf || GainVal < 0e0)
        println("")
        println("GainVal = $GainVal")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol, tSmol = $tSmol, uSmol = $uSmol, sSmol+tSmol+uSmol = $(sSmol+tSmol+uSmol), sBig+tBig+uBig = $(sBig+tBig+uBig)")
        println("ch1h2 = $ch1h2, ch3h1 = $ch3h1, ch3h2 = $ch3h2")
        println("ctheta12 = $ctheta12, ctheta13 = $ctheta13, ctheta23 = $ctheta23")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("GainVal Inf")  
    end

    #=if abs(sSmol+tSmol+uSmol) > 1e-15
        println("sSmol+tSmol+uSmol = $(sSmol+tSmol+uSmol)")
    end=#

    return GainVal

end


"""
    GainValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)

Returns `GainVal` based on initial momentum states `p1v` and `p2v` and final state `p4v` and differential cross section `dsigmadt` based on particle selection 12->34.  
Assumes f(x,p,μ)=constant over bin
"""
function GainValue4(p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    # cos(A-B) = 1.0 - b, so "a" part is 1.0, "b" part is the cos(A-B)-a
    b = pi^2*(p1v[3]-p2v[3])^2/2
    ch1h2a::Float64, ch1h2b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p1v[3]-p2v[3])-1.0)
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    p4::Float64 = p4v[1]
    st4::Float64,ct4::Float64 = sincospi(p4v[4])
    b = pi^2*(p4v[3]-p1v[3])^2/2
    ch4h1a::Float64, ch4h1b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p4v[3]-p1v[3])-1.0)
    b = pi^2*(p4v[3]-p2v[3])^2/2
    ch4h2a::Float64, ch4h2b::Float64 = b < eps(Float64) ? (1.0, -b) : (1.0,cospi(p4v[3]-p2v[3])-1.0)
    ch4h1::Float64 = cospi(p4v[3]-p1v[3])
    ch4h2::Float64 = cospi(p4v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    Es4::Float64 = m4 != 0e0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    Es4s::Float64 = Es4/p4
    E4::Float64 = Es4 + m4

    b = pi^2*(p1v[4]-p2v[4])^2/2
    ctheta12a::Float64, ctheta12b::Float64, ctheta12c::Float64 =  b < eps(Float64) ? (1.0, -b, ch1h2b*st1*st2) : (1.0,cospi(p1v[4]-p2v[4])-1.0, ch1h2b*st1*st2)
    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2 

    b = pi^2*(p4v[4]-p1v[4])^2/2
    ctheta14a::Float64, ctheta14b::Float64, ctheta14c::Float64 = b < eps(Float64) ? (1.0, -b, ch4h1b*st4*st1) : (1.0,cospi(p4v[4]-p1v[4])-1.0, ch4h1b*st4*st1)
    ctheta14::Float64 = ct4*ct1+ch4h1*st4*st1

    b = pi^2*(p4v[4]-p2v[4])^2/2
    ctheta24a::Float64, ctheta24b::Float64, ctheta24c::Float64 = b < eps(Float64) ? (1.0, -b, ch4h2b*st4*st2) : (1.0,cospi(p4v[4]-p2v[4])-1.0, ch4h2b*st4*st2)
    ctheta24::Float64 = ct4*ct2+ch4h2*st4*st2

    # TODO: make this more accurate for large separation of p1 p4 consistent with the new definitions of ctheta a,b,c above
    deltacorrect::Float64 = Es1*p4 - Es4*p1*ctheta14
    deltacorrect += m1*p4 - m4*p1*ctheta14
    deltacorrect += Es2*p4 - Es4*p2*ctheta24
    deltacorrect += m2*p4 - m4*p2*ctheta24

    # u = uBig + uSmol
    uBig::Float64 = (m4-m1)^2
    #uSmol::Float64 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    #uSmol::Float64 = 2*p4*p1*(ctheta14 - Es4s*Es1s - m1*Es4s/p1 - m4*Es1s/p4)
    #uSmol::Float64 = -2*E1*E4 + 2*(p1*p4+m1*m4) + 2*p4*p1*(ctheta14b + ctheta14c)
    uSmol::Float64 = 0.0
    if m1 == 0.0 && m4 == 0.0
        uSmol += 2*p4*p1*(ctheta14b + ctheta14c)
    elseif m1 == 0.0 && m4 != 0.0
        #E = sqrt(p^2+m^2) = p + (m^2)/(sqrt(m^2+p^2)+p)
        #E-p = (m^2)/(sqrt(m^2+p^2)+p) 
        uSmol += -2*(m42/(sqrt(m42+p4^2)+p4))*p1 + 2*p4*p1*(ctheta14b + ctheta14c)
    elseif m1 != 0.0 && m4 == 0.0
        uSmol += -2*(m12/(sqrt(m12+p1^2)+p1))*p4 + 2*p4*p1*(ctheta14b + ctheta14c)
    else 
        uSmol += -2*(E1*E4) + 2*(p1*p4+m1*m4) + 2*p4*p1*(ctheta14b + ctheta14c)
    end

    # t = tBig + tSmol
    tBig::Float64 = (m2-m4)^2
    #tSmol::Float64 = m12+m22+m32+m42 - sBig - uBig - tBig - sSmol - uSmol # Leads to Float64 issues, better to calculate directly 
    #tSmol::Float64 = -2*(m4*Es2 + m2*Es4 + Es2*Es4 - p2*p4*(ct2*ct4+ch4h2*st2*st4))
    #tSmol::Float64 = 2*p2*p4*(ctheta24 - Es2s*Es4s - m4*Es2s/p4 - m2*Es4s/p2)
    #tSmol::Float64 = -2*E2*E4 + 2*(p2*p4+m2*m4) + 2*p2*p4*(ctheta24b + ctheta24c)
    tSmol::Float64 = 0.0
    if m2 == 0.0 && m4 == 0.0
        tSmol += 2*p2*p4*(ctheta24b + ctheta24c)
    elseif m2 == 0.0 && m4 != 0.0
        #E = sqrt(p^2+m^2) = p + (m^2)/(sqrt(m^2+p^2)+p)
        #E-p = (m^2)/(sqrt(m^2+p^2)+p) 
        tSmol += -2*(m42/(sqrt(m42+p4^2)+p4))*p2 + 2*p4*p2*(ctheta24b + ctheta24c)
    elseif m2 != 0.0 && m4 == 0.0
        tSmol += -2*(m22/(sqrt(m22+p2^2)+p2))*p4 + 2*p4*p2*(ctheta24b + ctheta24c)
    else 
        tSmol += -2*(E2*E4) + 2*(p2*p4+m2*m4) + 2*p2*p4*(ctheta24b + ctheta24c)
    end

    #p2*p4 -E2*E4 +m2*m4
    
    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,sBig,m1,m2))/pi

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        error("stu check")
    end=#

    GainVal = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p4^2/(deltacorrect*sign(deltacorrect)))

    if (GainVal==Inf || GainVal == -Inf || GainVal < 0e0)
        println("")
        println("GainVal = $GainVal")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol, tSmol = $tSmol, uSmol = $uSmol, sSmol+tSmol+uSmol = $(sSmol+tSmol+uSmol), sBig+tBig+uBig = $(sBig+tBig+uBig)")
        println("ch1h2 = $ch1h2, ch4h1 = $ch4h1, ch4h2 = $ch4h2")
        println("ctheta12 = $ctheta12, ctheta14 = $ctheta14, ctheta24 = $ctheta24")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p4v = $p4v")
        error("GainVal Inf")  
    end

    #=if abs(sSmol+tSmol+uSmol) > 1e-15
        println("sSmol+tSmol+uSmol = $(sSmol+tSmol+uSmol)")
    end=#

    return GainVal

end

"""
    InvariantFlux(s,mu12,mu22)

returns the value of the invariant flux with 's' Mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvariantFlux(s::Float64,mu12::Float64,mu22::Float64)

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/2

end

"""
    InvariantFluxSmall(sSmol,sBig,m1,m2)

returns the value of the invariant flux with smalled 's' Mandelstram variable (sSmol = s - max((m1+m2)^2,(m3+m4)^2))
"""
function InvariantFluxSmall(sSmol::Float64,sBig::Float64,m1::Float64,m2::Float64)
    # Better accuracy for small s

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|

    if sBig == (m1+m2)^2
        # lambda(s,m1^2,m2^2) = sSmol*(sSmol+4*m1*m2)
        return sqrt(sSmol*(sSmol+4*m1*m2))/2
    else
        return sqrt(sSmol*(sSmol+2(sBig-(m1^2+m2^2)))+(sBig^2-2sBig*(m1^2+m2^2)+(m1^2-m2^2)^2))/2
    end

end

"""
    InvariantFlux2(s,mass12,mass22)

returns the value of the squared invariant flux with 's' Mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvariantFlux2(s::Float64,mu12::Float64,mu22::Float64)

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/4

end

"""
    InvariantFluxSmall(sSmol,mass1,mass2)

returns the value of the squared invariant flux with smalled 's' Mandelstram variable (sSmol = s - (m1+m2)^2)
"""
function InvariantFlux2Small(sSmol::Float64,sBig::Float64,m1::Float64,m2::Float64)
    # Better accuracy for small s

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    if sBig == (m1+m2)^2
        # lambda(s,m1^2,m2^2) = sSmol*(sSmol+4*m1*m2)
        return sSmol*(sSmol+4*m1*m2)/4
    else
        return (sSmol*(sSmol+2(sBig-(m1^2+m2^2)))+(sBig^2-2sBig*(m1^2+m2^2)+(m1^2-m2^2)^2))/4
    end

end


  

