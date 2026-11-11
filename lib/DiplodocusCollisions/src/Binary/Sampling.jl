"""
    UniformSampling!()

Uniform MC sampling of outgoing momentum states in the un-boosted Lab frame.
"""
function UniformSampling!(pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},#=Sval::Float64,Svalp::Float64,=#dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,3},SAtallyView4::AbstractArray{UInt32,3},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64})

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # generate random p direction in Lab frame for use in both p3 and p4 calculations
    RPointSphereCosThetaPhi!(pv)

    # === p3 === #
    #set random p3 direction 
    @view(p3v[1:3]) .= pv
    @view(p3pv[1:3]) .= pv

    # Calculate p3 value
    (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

    # S Array Tallies
    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states.
    #if NumStates != 0
        u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
        u3locMirror = location(u_low,u_up,u3_num,-p3v[2],u3_grid)
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        h3locMirror = location(h_low,h_up,h3_num,mod(p3v[3]+1e0,2e0),h3_grid)
        SAtallyView3[end,u3loc,h3loc] += UInt32(1)
        SAtallyView3[end,u3locMirror,h3locMirror] += UInt32(1)
    #end

    # Calculate S Array totals
    if NumStates == 1
        if p3_physical
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3loc,u3loc,h3loc] += Sval
            SAtallyView3[p3loc,u3loc,h3loc] += UInt32(1)
        end
    end

    if NumStates == 2
        if p3_physical
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3loc,u3loc,h3loc] += Sval
            SAtallyView3[p3loc,u3loc,h3loc] += UInt32(1)
        end
        if p3p_physical
            u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
            p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
            Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp
            SAtallyView3[p3ploc,u3ploc,h3ploc] += UInt32(1)
        end
    end

    # === p4 === #
    #set random p4 direction 
    @view(p4v[1:3]) .= pv
    @view(p4pv[1:3]) .= pv

    # Calculate p4 value
    (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p2v,p1v,mu2,mu1,mu4,mu3)

    # S Array Tallies
    # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states.
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
    h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
    SAtallyView4[end,u4loc,h4loc] += UInt32(1)
    SAtallyView4[end,u4locMirror,h4locMirror] += UInt32(1)


    # Calculate S Array totals
    if NumStates == 1
        if p4_physical
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView4[p4loc,u4loc,h4loc] += Sval
            SAtallyView4[p4loc,u4loc,h4loc] += UInt32(1)
        end
    end

    if NumStates == 2
        if p4_physical
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView4[p4loc,u4loc,h4loc] += Sval
            SAtallyView4[p4loc,u4loc,h4loc] += UInt32(1)
        end
        if p4p_physical
            u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
            h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
            p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
            Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView4[p4ploc,u4ploc,h4ploc] += Svalp
            SAtallyView4[p4ploc,u4ploc,h4ploc] += UInt32(1)
        end
    end

end

"""
    ImportanceSampling4!()

Importance MC sampling of outgoing momentum states weighted towards the direction of the incoming particle with highest momentum.
"""
function ImportanceSampling4!(p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,LocalGainTotal3::Array{Float64,3},LocalGainTotal4::Array{Float64,3},LocalGainTally3::Array{UInt32,3},LocalGainTally4::Array{UInt32,3},SmallParameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},WeightFactors::Tuple{Float64,Float64,Float64,Float64})


    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,m1,m2,m3,m4) = SmallParameters

    # weighting parameters
    (w3::Float64,w4::Float64,t::Float64,h::Float64) = WeightFactors
   
    prob3::Float64 = RPointSphereWeighted!(p3v,w3)    
    prob4::Float64 = RPointSphereWeighted!(p4v,w4)
    RotateToLab!(p3v,p4v,t,h)
    @. p3pv = p3v
    @. p4pv = p4v

    
    # === p3 === #

    # Calculate p3 value
    (p_physical,pp_physical,NumStates) = Momentum3Value2!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4)


    # S Array Tallies
    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states.
    #if NumStates != 0
        u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        LocalGainTally3[end,u3loc,h3loc] += UInt32(NumStates)# UInt32(1)
    #end

    # Calculate S Array totals
    if NumStates == 1
        if p_physical
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            Sval = SValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob3)
            LocalGainTotal3[p3loc,u3loc,h3loc] += Sval/prob3
            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
        end
    end

    if NumStates == 2
        if p_physical
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            Sval = SValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob3)
            LocalGainTotal3[p3loc,u3loc,h3loc] += Sval/prob3
            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
        end
        if pp_physical
            u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
            p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
            Svalp = SValue3(p3pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob3)
            LocalGainTotal3[p3ploc,u3ploc,h3ploc] += Svalp/prob3
            LocalGainTally3[p3ploc,u3ploc,h3ploc] += UInt32(1)
        end
    end

    # === p4 === #

    # Calculate p4 value
    (p_physical,pp_physical,NumStates) = Momentum3Value2!(p4v,p4pv,p2v,p1v,m2,m1,m4,m3)

    # S Array Tallies
    # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states.
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    #u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
    #h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
    LocalGainTally4[end,u4loc,h4loc] += UInt32(NumStates) #UInt32(1)
    #LocalGainTally4[end,u4locMirror,h4locMirror] += UInt32(1)

    # Calculate S Array totals
    if NumStates == 1
        if p_physical
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            Sval = SValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
            LocalGainTotal4[p4loc,u4loc,h4loc] += Sval/prob4
            LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
        end
    end

    if NumStates == 2
        if p_physical
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            Sval = SValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
            LocalGainTotal4[p4loc,u4loc,h4loc] += Sval/prob4
            LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
        end
        if pp_physical
            u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
            h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
            p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
            Svalp = SValue4(p4pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
            LocalGainTotal4[p4ploc,u4ploc,h4ploc] += Svalp/prob4
            LocalGainTally4[p4ploc,u4ploc,h4ploc] += UInt32(1)
        end
    end



    return nothing
end

"""
    ImportanceSampling()

Importance MC sampling of outgoing momentum states in the boosted centre-of-momentum frame. i.e. importance sampling weighted by the de-boosting of angles from the COM frame to the Lab frame.
"""
function ImportanceSampling(pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},βv::Vector{Float64},pCv::Vector{Float64}#=,Sval::Float64,Svalp::Float64=#,dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

        # Unpack parameters
        (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

        # Unpack locations NOT NEEDED
        #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate COM frame velocity vector 
    betaVec!(βv,p1v,p2v,mu1,mu2)
    # generate random p direction in COM frame and boost back to Lab frame for use in both p3 and p4 calculations
    RPointSphereCosThetaPhi!(pCv)
    RPointSphereBoost!(pv,βv,pCv)

    # === p3 === #
    #set random p3 direction 
    p3v .= pv
    p3pv .= pv

    # Calculate p3 value
    (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

    # S Array Tallies
    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). As this sampling is conducted in the COM frame, it would give the incorrect statistics to include the mirrored states. Therefore only the un-mirrored states are included.
    # un-mirrored states have sign(p3v[2]) == sign(pv[2]) 

    #if NumStates != 0
        u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        SAtallyView3[u3loc,h3loc] += UInt32(1)
    #end

    # Calculate S Array totals
    if NumStates == 1
        if p3_physical
            if sign(p3v[2]) == sign(pv[2])
                p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3loc,u3loc,h3loc] += Sval/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                    u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                    u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                end
            end
        end
    end

    if NumStates == 2
        if p3_physical
            if sign(p3v[2]) == sign(pv[2])
                p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3loc,u3loc,h3loc] += Sval/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                    u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                    u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                end
            end
        end
        if p3p_physical
            if sign(p3pv[2]) == sign(pv[2])
                u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3ploc,h3ploc] = max(p3MaxView[u3ploc,h3ploc],p3pv[1])
                    u3MinView[p3ploc,h3ploc] = min(u3MinView[p3ploc,h3ploc],p3pv[2])
                    u3MaxView[p3ploc,h3ploc] = max(u3MaxView[p3ploc,h3ploc],p3pv[2])
                end
            end
        end
    end

    # === p4 === #
    #set random p4 direction 
    p4v .= pv
    p4pv .= pv

    # Calculate p4 value
    (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p2v,p1v,mu2,mu1,mu4,mu3)

    # S Array Tallies
    # For each u4,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u4 and a rotation of h4 by pi i.e. mod(h4+1,2). As this sampling is conducted in the COM frame, it would give the incorrect statistics to include the mirrored states. Therefore only the un-mirrored states are included.
    # un-mirrored states have sign(p4v[2]) == sign(pv[2]) 
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    # Calculate S Array totals
    if NumStates == 1
        if p4_physical
            if sign(p4v[2]) == sign(pv[2])
                p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4loc,u4loc,h4loc] += Sval/prob
                if MinMax
                    p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                    u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                    u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                end
            end
        end
    end

    if NumStates == 2
        if p4_physical
            if sign(p4v[2]) == sign(pv[2])
                p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4loc,u4loc,h4loc] += Sval/prob
                if MinMax
                    p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                    u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                    u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                end
            end
        end
        if p4p_physical
            if sign(p4pv[2]) == sign(pv[2])
                u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4ploc,u4ploc,h4ploc] += Svalp/prob
                if MinMax
                    p4MaxView[u4ploc,h4ploc] = max(p4MaxView[u4ploc,h4ploc],p4pv[1])
                    u4MinView[p4ploc,h4ploc] = min(u4MinView[p4ploc,h4ploc],p4pv[2])
                    u4MaxView[p4ploc,h4ploc] = max(u4MaxView[p4ploc,h4ploc],p4pv[2])
                end
            end
        end
    end

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end

function ImportanceSampling2(p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},Tval::Float64,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate importance sampled points 
    importance!(p3v,p4v,p1v,p2v,p1cv,p1cBv,mu1,mu2,mu3,mu4)
    prob3 = p3v[5]
    prob4 = p4v[5]

    p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
    p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView3[u3loc,h3loc] += UInt32(1)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    Sval3,Sval4 = SValueTest(p3v,p4v,Tval,p1v,p2v,mu1,mu2,mu3,mu4)
    SAtotalView3[p3loc,u3loc,h3loc] += Sval3*prob3
    SAtotalView4[p4loc,u4loc,h4loc] += Sval4*prob4

    #=if p4v[1] <= 1e-8
        println("p4 = $(p4v[1])")
        println("sp = $(Sval4*prob4)")
        println("p4loc = $p4loc")

    end=#

    #println("prob = $prob3")
    #println("sval3 = $Sval3")

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end

function ImportanceSampling3(p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate importance sampled points 
    tSmol,sSmol = importance!(p3v,p4v,p1v,p2v,p1cv,p1cBv,mu1,mu2,mu3,mu4)
    prob3 = p3v[5]
    prob4 = p4v[5]

    p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
    p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView3[u3loc,h3loc] += UInt32(1)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    Sval3,Sval4 = SValueTest2(p3v,p4v,dsigmadt,tSmol,sSmol,p1v,p2v,mu1,mu2,mu3,mu4)
    SAtotalView3[p3loc,u3loc,h3loc] += Sval3*prob3
    SAtotalView4[p4loc,u4loc,h4loc] += Sval4*prob4

    #=if p4v[1] <= 1e-8
        println("p4 = $(p4v[1])")
        println("sp = $(Sval4*prob4)")
        println("p4loc = $p4loc")

    end=#

    #println("prob = $prob3")
    #println("sval3 = $Sval3")

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end


function Inv_cdfElePhoElePho(sS::Float64)

    # generates a random number x between [0,1] and returns the inverse cdf of the Compton cross section i.e. tSmol value.

    #cdf(tSmol) is
    #(((tS + sS*(sS + tS)) * (8*tS + sS*(16 + sS^3 + 2*sS^2*(9 + tS) + tS*(16 + tS) + sS*(32 + tS*(11 + tS)))))/(2*(1 + sS)^2 *(sS + tS)) + sS*(-8 + (-4 + sS)*sS)*log(1 + sS + tS + tS/sS))/(sS^3*(8/sS + (sS*(2 + sS))/(2*(1 + sS)^2) + ((-8 + (-4 + sS)*sS)*log(1 + sS))/sS^2))

    # tSmol can be normalised by sS^2/(sS+1) to be tSmolN, then cos(ts) = 1-2tSmolN.
    # Therefore prob(cos(ts)) = |-2*prob(tSmolN)|

    x = rand(Float64)

    #logsS=log10(sS)
    logsS=log10(sS)

    # Inverse cdf(tSmol) can be approximate by to regions
    if logsS < 1e0 #0.8

        tSmolN = (x - 1)
        #tS *= sS^2/(sS + 1)

        # more floating point safe than acos()
        ts = 2*atan(sqrt(-tSmolN),sqrt(tSmolN+1))/pi
        
        prob_cts = 1e0 #2e0 # factor of 2 absorbed into cos(theta) probability distribution

    else #logsS >= log10

        tSmolN = x^(logsS)-1

        # more floating point safe than acos()
        ts = 2*atan(sqrt(-tSmolN),sqrt(tSmolN+1))/pi

        prob = (1+tSmolN)^(1/logsS-1) 
        prob /= logsS 

        prob_cts = prob #2e0 # factor of 2 absorbed into cos(theta) probability distribution

    end

    tSmol = tSmolN * sS^2/(sS + 1)

    if ts < -1.0 #tS < -sS^2/(sS+1)
        println("tS = $tS")
        println("sS = $sS")
        println("tS min = $(-sS^2/(sS+1))")
        error("tS < sS^2/(sS+1)")
    end

    #println("tS = $(tS*sS^2/(sS+1))")

    return (tSmol,ts,prob_cts)

end

function RPointPhi()

    phi = 2*rand(Float64)

    return phi

end

function p1cv!(p1cv::Vector{Float64},p1cBv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,sS::Float64)

    p1 = p1v[1]
    ct1 = p1v[2]
    st1 = sqrt(1e0 - ct1^2)

    p2 = p2v[1]
    ct2 = p2v[2]
    st2 = sqrt(1e0 - ct2^2)

    sh1h2, ch1h2 = sincospi(p1v[3]-p2v[3])

    m12 = m2^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    val1 = (ct1*ct2+ch1h2*st1*st2)
    val2 = 2*p1*p2*val1
    val3 = val2 + p1^2 + p2^2

    w = asinh(sqrt(val3)/sqrt(sS+(m1+m2)^2))

    # calculate p1cv value
    if (y=p2/p1) <= 1e-6
        t1c = 2*atan((1-val1^2)*y^2/8)/pi
    elseif y >= 1e6
        t1c = 2*atan((1-val1)/(2*(1+val1)) + (val1-1)/(1+val1)/y)/pi
    else
        y = p1+p2*val1
        y /= sqrt(val3)
        t1c = acos(y)/pi
    end

    p1cv[2] = cospi(t1c)
    p1cv[4] = t1c

    y = sh1h2*st1*st2*sqrt(val3)
    x = ct2*(p1+p2*val1)-ct1*(p2+p1*val1)
    p1cv[3] = mod(atan(y,x)/pi,2)

    # Boost p1c to CM frame
    #LorentzBoost!(p1cBv,p1cv,m1,w)

    return nothing

end

function tstar(tS::Float64,sS::Float64,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # calculate tstar value

    #cts = (sS+(m1+m2)^2)*(2*tS+sS+(2*m12+m32-m42-4*m1*m3+2*m1*m2)) + (m12-m22)*(m32-m42)
    #cts /= 4*F12*F34

    #cts= ((m1 - m2)*(m1 + m2)*(m3 - m4)*(m3 + m4) + ((m1 + m2)^2 + sS)*(2*m1*(m1 + m2) - 4*m1*m3 + m3^2 - m4^2 + sS + 2*tS))/(sqrt(sS*(4*m1* m2 + sS))*sqrt((m3^2 - m4^2)^2 - 2*(m3^2 + m4^2)*((m1 + m2)^2 + sS) + ((m1 + m2)^2 + sS)^2))

    # cts with m values put in for compton
        #cts = 2*tS/sS^2 + 2*tS/sS + 1 
    # this is potentially very small negative +1 so subject to float issues
    # use 1-cts = 2*sin^2(ts/2) = -2*tS*(sS+1)/sS^2
    # use 1+cts = 2*cos^2(ts/2) = 2(tS*(sS+1)+sS^2)/sS^2
        #onemcts = -2*tS*(sS+1)/sS^2
        #ts = 2*asin(sqrt(onemcts/2))   
        #ts = 2*atan(-tS*(sS+1),tS*(sS+1)+sS^2)/pi

        # tS normalised by sS^2/(sS+1)
        ts = 2*atan(sqrt(-tS),sqrt(tS+1))/pi

    #println(tS)
    #println(sS)
    #println(sS^2/(sS+1))
    #=if abs(cts) > 1e0
        println("tS=$tS")
        println("sS=$sS")
        println("tS max = $(sS^2/(sS+1))")
        println("cts = $cts")
        error("cts > 1e0")
    end=#
    #if cts > 1e0
    #    cts = 1e0
    #elseif cts < -1e0
    #    cts = -1e0
    #end

    return ts 

end

function p3p4v!(p3v::Vector{Float64},p4v::Vector{Float64},p1cBv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},ts::Float64,hs::Float64,sS::Float64,m1::Float64,m2::Float64,m3::Float64,m4::Float64,p3B::Float64)

    t1cB = p1cBv[4]
    (st1cB,ct1cB) = sincospi(t1cB)

    (sts,cts) = sincospi(ts)

    #p1 = p1v[1]
    #p2 = p2v[1]

    #ct1 = p1v[2]
    #ct2 = p2v[2]
    #st1 = sqrt(1-ct1^2)
    #st2 = sqrt(1-ct2^2)

    #(sh1,ch1) = sincospi(p1v[3])
    #(sh2,ch2) = sincospi(p2v[3])
    #ch1h2 = cospi(p1v[3]-p2v[3])

    #m12 = m1^2
    #m22 = m2^2

    #Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    #E1::Float64 = Es1 + m1
    #Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    #E2::Float64 = Es2 + m2

    #val1 = (ct1*ct2+ch1h2*st1*st2)
    #val2 = 2*p1*p2*val1
    #val3 = val2 + p1^2 + p2^2
    #val3 = 2*p1*p2*(val1-1) + (p1+p2)^2

    #w = asinh(sqrt(val3)/sqrt(sS+(m1+m2)^2))

    # calculate ct3cs and ct4cs
    #t3cs = acos(cts*ct1cs+cospi(hs)*sts*st1cs)/pi
    if ts <= 1e-5 && t1cB <= 1e-5
        t3cB = sqrt(-2*cospi(hs)*ts*t1cB+ts^2+t1cB^2) # /pi not needed as accounted for in approximation
    else
        t3cB = acos(cts*ct1cB+cospi(hs)*sts*st1cB)/pi
    end

    # calculate h3,h4 from ts, hs
    A = atan(sts*sinpi(hs),cts*st1cB-cospi(hs)*sts*ct1cB)/pi
    #A1 = atan((sts*sinpi(hs))/(cts*st1cB-cospi(hs)*sts*ct1cB))/pi
    #A2 = acos((cts*st1cB-cospi(hs)*sts*ct1cB)/sinpi(t3cB))/pi
    h3cB = mod(p1cBv[3] + A,2)

    # Boosted frame angles and momenta
    p3v[4] = t3cB
    p3v[2] = cospi(p3v[4])
    p4v[4] = 1-t3cB
    p4v[2] = cospi(p4v[4])
    p3v[3] = h3cB
    p4v[3] = mod(h3cB+1e0,2)
    #p3B = InvariantFluxSmall(sS,m3,m4)/sqrt(sS+(m1+m2)^2)
    #p4B = p3B
    p3v[1] = p3B
    p4v[1] = p3B # in COM frame p3B = p4B

    # Doppler Factors
    #DopplerFactor!(p3v,m3,-w)
    #DopplerFactor!(p4v,m4,-w)

    # De-boost p3 and p4 (modifies p3v and p4v)
    #LorentzBoost!(p3v,p3v,m3,-w)
    #LorentzBoost!(p4v,p4v,m4,-w)

    # p3 value
    #m32 = m3^2
    #m42 = m4^2
    #p3s = InvariantFluxSmall(sS,m3,m4)/sqrt(sS+(m1+m2)^2)
    #E3s = sqrt(p3s^2+m32)
    #testE3s = 1/2*(sS+(m12+m32+2*m1*m2+m22-m42))/sqrt(sS+(m1+m2)^2)
    #p3s = sqrt(E3s^2-m32)
    
    #p3v[1] = γ*p3s*ct3cs+γβ*E3s
    #p3c = p3s*(exp(w)*cospi(t3cB/2)^2+exp(-w)*sinpi(t3cB/2)^2)+γβ*m32/(E3s+p3s)
    #p3c = p3c/ct3c*sign(ct3c)
    #p3v[1] = p3c 

    # p4 value
    #p4s = p3s
    #E4s = sqrt(p4s^2+m42)

    #p4v[1] = -γ*p4s*ct3cs+γβ*E4s
    #p4c = p4s*(exp(w)*sinpi(t3cs/2)^2+exp(-w)*cospi(t3cs/2)^2)+γβ*m42/(E4s+p4s)
    #p4c = p4c/ct4c#*sign(ct4c)
    #p4v[1] = p4c

    # rotate to lab frame
    #=βc = sqrt(val3)#/(E1+E2) 
    tβ = acos((p1*ct1 + p2*ct2)/βc)/pi #/ (E1+E2) # add an approximation here for p1/p2 << 1 etc
    if (y=p1/p2) <= 1e-6
        tβ = acos(ct2 + st2*(-ch1h2*ct2*st1 + ct1*st2)*y)/pi
    elseif y >= 1e6
        tβ = acos(ct1 + st1*(ct2*st1 - ch1h2*ct1*st2)/y)/pi
    else
        tβ = acos((p1*ct1 + p2*ct2)/βc)/pi
    end
    (stβ,ctβ) = sincospi(tβ)
    x = p1*st1*ch1 + p2*st2*ch2
    y = p1*st1*sh1 + p2*st2*sh2
    hβ = mod(atan(y,x)/pi,2)
    (shβ,chβ) = sincospi(hβ)
    (st3c,ct3c) = sincospi(p3v[4])
    (st4c,ct4c) = sincospi(p4v[4])
    (sh3c,ch3c) = sincospi(p3v[3])
    (sh4c,ch4c) = sincospi(p4v[3])
    (sh3hc,ch3hc) = sincospi(p3v[3]-hβ) 
    (sh4hc,ch4hc) = sincospi(p4v[3]-hβ) 

    p3v[4] = acos(ct3c*ctβ + st3c*stβ*ch3hc)/pi
    p3v[2] = cospi(p3v[4])
    #x = -st3c*sh3c*shβ + chβ*(st3c*ch3c*ctβ + ct3c*stβ)
    x = -sh3c*shβ*st3c+chβ*(ch3c*ctβ*st3c + ct3c*stβ)
    #y = st3c*sh3c*chβ + shβ*(st3c*ch3c*ctβ + ct3c*stβ)
    y = chβ*sh3c*st3c+shβ*(ch3c*ctβ*st3c+ct3c*stβ)
    p3v[3] = mod(atan(y,x)/pi,2)

    p4v[4] = acos(ct4c*ctβ + st4c*stβ*ch4hc)/pi
    p4v[2] = cospi(p4v[4])
    #x = -st4c*sh4c*shβ + chβ*(st4c*ch4c*ctβ + ct4c*stβ)
    x = -sh4c*shβ*st4c+chβ*(ch4c*ctβ*st4c + ct4c*stβ)
    #y = st4c*sh4c*chβ + shβ*(st4c*ch4c*ctβ + ct4c*stβ)
    y = chβ*sh4c*st4c+shβ*(ch4c*ctβ*st4c+ct4c*stβ)
    p4v[3] = mod(atan(y,x)/pi,2)
    =#



    #=if p1 >= 1e4 && p2 <= 1e0 && log10(sS) > 4.0
        println("ss = $sS")
        println("t1 = $(acos(p1v[2])/pi)")
        println("t2 = $(acos(p2v[2])/pi)")
        println("t1c = $(acos(p1cv[2])/pi)")
        println("t1cs = $t1cs")
        println("ts = $ts")
        println("t3cs = $t3cs")
        println("t3 = $(acos(p3v[2])/pi)")
        println("t4 = $(acos(p4v[2])/pi)")
        println("")
    end=#


end

function importance!(p3v::Vector{Float64},p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # calculate random tS and hs and then ts

        # pre-defining terms for efficiency 
        p1::Float64 = p1v[1]
        p2::Float64 = p2v[1]
    
        ct1::Float64 = p1v[2] #cospi(p1v[2])
        ct2::Float64 = p2v[2] #cospi(p2v[2]) 
    
        st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
        st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])
        (sh1,ch1) = sincospi(p1v[3])
        (sh2,ch2) = sincospi(p2v[3])
    
        ch1h2::Float64 = cospi(p1v[3]-p2v[3])
    
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

        val1 = (ct1*ct2+ch1h2*st1*st2)
        val2 = 2*p1*p2*val1
        val3 = val2 + p1^2 + p2^2
        #println("val3 = $val3")
        #val3 = 2*p1*p2*(val1-1e0) + (p1+p2)^2
        #println("val3 = $val3")

    sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol::Float64 = 2*p1*p2*(-val1 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    w = asinh(sqrt(val3)/sqrt(sSmol+(m1+m2)^2))

    (tSmol,ts,prob) = Inv_cdfElePhoElePho(sSmol)
    hs = RPointPhi()
    #ts = tstar(tSmol,sSmol,m1,m2,m3,m4)

    if (z=p1/p2) <= 1e-8 # p1 small
        tβ = acos(ct2 + st2*(-ch1h2*ct2*st1 + ct1*st2)*z)/pi
        x = z*st1*ch1 + st2*ch2
        y = z*st1*sh1 + st2*sh2
    elseif z >= 1e8 # p1 large
        tβ = acos(ct1 + st1*(ct2*st1 - ch1h2*ct1*st2)/z)/pi
        x = st1*ch1 + st2*ch2/z
        y = st1*sh1 + st2*sh2/z
    else
        β = sqrt(val3)
        tβ = acos((p1*ct1 + p2*ct2)/β)/pi
        x = p1*st1*ch1 + p2*st2*ch2
        y = p1*st1*sh1 + p2*st2*sh2
    end
    #(stβ,ctβ) = sincospi(tβ)
    hβ = mod(atan(y,x)/pi,2)

    # calculate p1cv
    @view(p1cv[1:3]) .= p1v
    p1cv[4] = acos(p1v[2])/pi

    p2cv::Vector{Float64} = zeros(Float64,4)
    p2cBv::Vector{Float64} = zeros(Float64,4)
    @view(p2cv[1:3]) .= p2v
    p2cv[4] = acos(p2v[2])/pi

    RotateToCentre!(p1cv,tβ,hβ)

    RotateToCentre!(p2cv,tβ,hβ)
    #println("p1cv = $(p1cv)")
    #println("p2cv = $(p2cv)")
    #println("p1v = $(p1v)")
    #println("p2v= $(p2v)")
    #println("tβ = $(tβ)")
    #println("hβ = $(hβ)")
    p1B = InvariantFluxSmall(sSmol,m1,m2)/sqrt(sSmol+sBig)
    p3B = InvariantFluxSmall(sSmol,m3,m4)/sqrt(sSmol+sBig)
    #println("p1B = $(p1B)")
    ForwardLorentzBoost!(p2cBv,p2cv,m2,w,p1B)
    #println("p2cBv = $(p2cv)")
    ForwardLorentzBoost!(p1cBv,p1cv,m1,w,p1B)
    #println("p1cBv = $(p1cBv)")
    #p1cv!(p1cv,p1cBv,p1v,p2v,m1,m2,sSmol)
    #println("p1cv = $(p1cv)")
    #println("")

    # calculate p3v and p4v in COM frame
    p3p4v!(p3v,p4v,p1cBv,p1v,p2v,ts,hs,sSmol,m1,m2,m3,m4,p3B)

    # energy check in COM
        #=
        eng_error = 1-(sqrt(m3^2+p3v[1]^2)+sqrt(m4^2+p4v[1]^2))/(sqrt(m1^2+p1cBv[1]^2)+sqrt(m2^2+p2cBv[1]^2))

        mom_z_error = (p1cBv[1]*p1cBv[2]+p2cBv[1]*p2cBv[2])-(p3v[1]*p3v[2]+p4v[1]*p4v[2])/(p1cBv[1]*p1cBv[2]+p2cBv[1]*p2cBv[2])

        mom_x_error = (p1cBv[1]*sin(acos(p1cBv[2]))*cospi(p1cBv[3])+p2cBv[1]*sin(acos(p2cBv[2]))*cospi(p2cBv[3]))-(p3v[1]*sin(acos(p3v[2]))*cospi(p3v[3])+p4v[1]*sin(acos(p4v[2]))*cospi(p4v[3]))/(p1cBv[1]*sin(acos(p1cBv[2]))*cospi(p1cBv[3])+p2cBv[1]*sin(acos(p2cBv[2]))*cospi(p2cBv[3]))
        #mom_x_error = 1-(p3v[1]*sinpi(p3v[4])*cospi(p3v[3])+p4v[1]*sinpi(p4v[4])*cospi(p4v[3]))/(p1v[1]*sin(acos(p1v[2]))*cospi(p1v[3])+p2v[1]*sin(acos(p2v[2]))*cospi(p2v[3]))
        if eng_error > 1e-6 || (mom_z_error > 1e-6 && (1-abs(p2cBv[1]*p2cBv[2]/(p1cBv[1]*p1cBv[2]))) > 1e-6) #|| mom_x_error > 1e-4
            println("COM")
            println("eng_error = $eng_error")
            println("mom_z_error = $mom_z_error")
            println("mom_x_error = $mom_x_error")
            println("w=$w")
            println("p1 = $(p1cBv[1])")
            println("ct1 = $(p1cBv[2])")
            println("h1 = $(p1cBv[3])")
            println("p2 = $(p2cBv[1])")
            println("ct2 = $(p2cBv[2])")
            println("h2 = $(p2cBv[3])")
            println("p3 = $(p3v[1])")
            println("ct3 = $(p3v[2])")
            println("h3 = $(p3v[3])")
            println("p4 = $(p4v[1])")
            println("ct4 = $(p4v[2])")
            println("h4 = $(p4v[3])")
            println("$(p3v[1]*p3v[2]+p4v[1]*p4v[2])")
            println("unboosed")
            println("p1/m = $(p1/m1)")
            println("t1 = $(p1cv[4])")
            println("p2 = $(p2)")
            println("t2 = $(p2cv[4])")
            println("")
        end
        
        =#

    # Doppler Factors
    DopplerFactor!(p3v,m3,w)
    DopplerFactor!(p4v,m4,w)

    # De-boost p3 and p4 (modifies p3v and p4v)
    BackwardsLorentzBoost!(p3v,p3v,m3,w)
    BackwardsLorentzBoost!(p4v,p4v,m4,w)

    # Rotate p3 p4 back to lab (modifies p3v and p4v)
    RotateToLab!(p3v,tβ,hβ)
    RotateToLab!(p4v,tβ,hβ)

        # energy check in lab
        #=
        eng_error = 1-(sqrt(m3^2+p3v[1]^2)+sqrt(m4^2+p4v[1]^2))/(E1+E2)
        mom_z_error = ((p1v[1]*p1v[2]+p2v[1]*p2v[2])-(p3v[1]*p3v[2]+p4v[1]*p4v[2]))/(p1v[1]*p1v[2]+p2v[1]*p2v[2])
        mom_x_error = 1-(p3v[1]*sin(acos(p3v[2]))*cospi(p3v[3])+p4v[1]*sin(acos(p4v[2]))*cospi(p4v[3]))/(p1v[1]*sin(acos(p1v[2]))*cospi(p1v[3])+p2v[1]*sin(acos(p2v[2]))*cospi(p2v[3]))
        #mom_x_error = 1-(p3v[1]*sinpi(p3v[4])*cospi(p3v[3])+p4v[1]*sinpi(p4v[4])*cospi(p4v[3]))/(p1v[1]*sin(acos(p1v[2]))*cospi(p1v[3])+p2v[1]*sin(acos(p2v[2]))*cospi(p2v[3]))
        if eng_error > 1e-4 || mom_z_error > 1e-4 #|| mom_x_error > 1e-2
            println("lab")
            println("eng_error = $eng_error")
            println("mom_z_error = $mom_z_error")
            println("mom_x_error = $mom_x_error")
            println("w=$w")
            println("p1 = $(p1v[1])")
            println("ct1 = $(p1v[2])")
            println("h1 = $(p1v[3])")
            println("p2 = $(p2v[1])")
            println("ct2 = $(p2v[2])")
            println("h2 = $(p2v[3])")
            println("p3 = $(p3v[1])")
            println("ct3 = $(p3v[2])")
            println("h3 = $(p3v[3])")
            println("p4 = $(p4v[1])")
            println("ct4 = $(p4v[2])")
            println("h4 = $(p4v[3])")
            println("")
        end
        =#
        

    p3v[5] /= prob
    p4v[5] /= prob

    return (tSmol,sSmol)

end


"""
    ForwardLorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64,pB::Float64)

    Calculates the forward Lorentz boost (in z direction) from the lab frame to the COM frame for the momentum vector `pv` of a particle with (normalised) mass `m` by the rapidity `w` (must be positive!), placing the results in the boosted vector `pBv``. Assumes boosted magnitude of momentum `pB` in the COM frame is known, therfore only computes the angle transformations.

    Momentum vectors have coponents [p,cos(theta),phi,theta]
"""
function ForwardLorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64,pB::Float64)

    if w < 0e0
        println("w=$w")
        error("w must be positive")
    end

    p::Float64 = pv[1]
    t::Float64 = pv[4]
    h::Float64 = pv[3]

    hB::Float64 = h # phi unaffected by boost

    if m == 0e0

        ew = exp(w)
        tB = 2*atan(ew*tanpi(t/2))/pi

    #=elseif t <= 1e-8 # small 
        # forth order is 1e-20 so beyond numerical precision

        # tmp is p*cosh(w)-sqrt(m^2+p^2)*sinh(w) # this changes sign
        tmp = (exp(w)*(p-sqrt(p^2+m^2))+exp(-w)*(p+sqrt(p^2+m^2)))/2
        # tmp2 is sqrt(m^2+p^2)*cosh(w)-p*sinh(w) # this is +ve
        #tmp2 = (exp(w)*(sqrt(p^2+m^2)-p)+exp(-w)*(sqrt(p^2+m^2)+p))/2
        #pBcheck = tmp*sign(tmp)
        #pBcheck += p*sinh(w)*tmp2/(2*tmp*sign(tmp)) * (pi^2*t^2)
        #pBcheck -= p*sinh(w)*(tmp2^3-m^2*tmp2+3*m^2*p*sinh(w)) / (24*tmp^3*sign(tmp)) * (pi^4*t^4)

        #=if isapprox(pBcheck,pB) == false
            println("small p")
            println("p = $p")
            println("w = $w")
            println("t = $t")
            println("pBcheck = $pBcheck")
            println("pB = $pB")
        end=#

        tB = sign(tmp)
        tB -= p^2*sign(tmp) / (2*tmp^2) * (pi*t)^2 
        #tB += p^2*(3*p^3*cosh(w)-2*p^3*cosh(3w)-sqrt(p^2+m^2)*sinh(w)*(-2*m^2+p^2+(2*m^2-4*p^2)*cosh(2w))) / (24*tmp^5*sign(tmp)) * (pi^4*t^4)

        if abs(tB) > 1e0
            println("tB = $tB")
            println("p = $p")
            println("w = $w")
            println("t = $t")
            println("pB = $pB")
            println("tmp= $tmp")
            error("tB > 1e0")
        end
        
        tB = acos(tB)/pi
    =#
    elseif (z=p/m) <= 1e-4 # small p approximation 
        # valid to at fifth order in p/m

        # pB check
        #pBcheck = m*sinh(w) 
        #pBcheck -= m*cospi(t)*cosh(w)*z 
        #pBcheck += m*(-cospi(2*t) + cosh(2*w))*csch(w)/4 * z^2 
        #pBcheck += m*(cospi(t)*coth(w)*csch(w)*sinpi(t)^2)/2 * z^3 
        #pBcheck -= m*((3*cospi(4*t) + (-2 - 4*cospi(2*t) + 2*cospi(4*t))*cosh(2*w) + cosh(4*w))*csch(w)^3)/64 * z^4

        #=if isapprox(pBcheck,pB) == false
            println("small p")
            println("z = $z")
            println("w = $w")
            println("t = $t")
            println("pBcheck = $pBcheck")
            println("pB = $pB")
        end=#

        tB = 2*sinh(w)/sinpi(t) / z
        tB -= 2*cosh(w)*cospi(t)/sinpi(t)
        tB += (sinpi(t)/sinh(w)+2*sinh(w)/sinpi(t))/2 * z

        tB = 2*atan(tB)/pi

        #tB = 2*atan(p*sinpi(t),2(sign(w)*sinh(w)))/pi

    elseif z >= 1e4 # large p approximation
        # valid to at fifth order in m/p

        # pB check 
        # tmp is cosh(w)-cospi(t)*sinh(w)
        tmp = exp(w)*sinpi(t/2)^2 + exp(-w)*cospi(t/2)^2
        #pBcheck = m*tmp * z 
        #pBcheck -= m*(1/tmp-cosh(w))/2 / z
        #tmp2 = exp(2*w)*sinpi(t/2)^2 + exp(-2*w)*cospi(t/2)^2
        #pBcheck += m*(tmp2/tmp^3-cosh(w))/8 /z^3

        #=if isapprox(pBcheck,pB) == false
            println("large p")
            println("z = $z")
            println("w = $w")
            println("t = $t")
            println("pBcheck = $pBcheck")
            println("pB = $pB")
        end=#

        ew = exp(w)
        ttd2 = tanpi(t/2)

        val = ew*ttd2
        val += ew*ttd2 * sinh(w)/(2*tmp) / z^2
        val += ew*ttd2 * sinh(w) * (-1 + 5*cosh(2*w) + 2*(2*cospi(t) + cospi(2*t))*sinh(w)^2 - (2 + 6*cospi(t))*sinh(2*w)) / (32*tmp)^3 / z^4

        #tB = 2*atan(ew*ttd2*sqrt(1+(m^2)/(p^2)*(sinh(w))/(cosh(w)-cospi(t)*sinh(w))))/pi

        tB = 2*atan(val)/pi

    else # no approximation

        if w >= 6e0 # not quite accurate
            # approximate acos(), keep to fifth order in exp(w)
            ew = exp(-w)
            val = 1e0
            val -= 2*p*sinpi(t)/(sqrt(p^2+m^2)-p*cospi(t))/pi * ew
            val += p*sinpi(t)*(-6*m^2+2*p^2*sinpi(t)^2)/(sqrt(p^2+m^2)-p*cospi(t))^3/(3*pi) * ew^3
            val -= 2*p*sinpi(t)*(p^2*sinpi(t)^2-5m^2)^2/(sqrt(p^2+m^2)-p*cospi(t))^5/(5*pi) * ew^5

            tB = val

            #tB = 2*atan(z)/pi

        #=elseif t <= 1e-8 # small 
            # forth order is 1e-20 so beyond numerical precision
    
            # tmp is p*cosh(w)-sqrt(m^2+p^2)*sinh(w) # this changes sign
            tmp = (exp(w)*(p-sqrt(p^2+m^2))+exp(-w)*(p+sqrt(p^2+m^2)))/2
            # tmp2 is sqrt(m^2+p^2)*cosh(w)-p*sinh(w) # this is +ve
            #tmp2 = (exp(w)*(sqrt(p^2+m^2)-p)+exp(-w)*(sqrt(p^2+m^2)+p))/2
            #pBcheck = tmp*sign(tmp)
            #pBcheck += p*sinh(w)*tmp2/(2*tmp*sign(tmp)) * (pi^2*t^2)
            #pBcheck -= p*sinh(w)*(tmp2^3-m^2*tmp2+3*m^2*p*sinh(w)) / (24*tmp^3*sign(tmp)) * (pi^4*t^4)
    
            #=if isapprox(pBcheck,pB) == false
                println("small p")
                println("p = $p")
                println("w = $w")
                println("t = $t")
                println("pBcheck = $pBcheck")
                println("pB = $pB")
            end=#
    
            tB = sign(tmp)
            println("sign(tmp) = $(sign(tmp))")
            tB -= p^2*sign(tmp) / (2*tmp^2) * (pi*t)^2 
            #tB += p^2*(3*p^3*cosh(w)-2*p^3*cosh(3w)-sqrt(p^2+m^2)*sinh(w)*(-2*m^2+p^2+(2*m^2-4*p^2)*cosh(2w))) / (24*tmp^5*sign(tmp)) * (pi^4*t^4)
    
            if abs(tB) > 1e0
                println("tB = $tB")
                println("p = $p")
                println("w = $w")
                println("t = $t")
                println("pB = $pB")
                println("tmp= $tmp")
                error("tB > 1e0")
            end
            
            tB = acos(tB)/pi
        =#
        else # no approximation

            y = pB+(exp(w)*(sqrt(p^2+m^2)-p*cospi(t))/2-exp(-w)*(p*cospi(t)+sqrt(p^2+m^2))/2)
            x = pB-(exp(w)*(sqrt(p^2+m^2)-p*cospi(t))/2-exp(-w)*(p*cospi(t)+sqrt(p^2+m^2))/2)

            if y<0e0
                #println("")
                #println("forwards boost")
                #println("warning: y<0,=$y, setting y=0e0")
                #println("x=$x")
                #println("y=$y")
                #println("pB=$pB")
                #println("w=$w")
                #println("p=$p")
                #println("t=$t")
                #println("")
                y = 0e0
                #tB = acos((-exp(w)*(sqrt(p^2+m^2)-p*cospi(t))/2+exp(-w)*(p*cospi(t)+sqrt(p^2+m^2))/2)/pB)/pi
                #println("tB = $tB")
            end

            if x < 0e0
                x=0e0
            end

            # acos(x) can lead to a domain error with x > 1 due to floating point precision, atan approach better 
            #tB = acos((cosh(w)*p*cospi(t)-sinh(w)*sqrt(p^2+m^2))/pB)/pi

            #tBtest = 2*atan(sqrt(-p*cosh(w)*cospi(t)+sqrt(p^2+m^2)*sinh(w)+pB),sqrt(+p*cosh(w)*cospi(t)-sqrt(p^2+m^2)*sinh(w)+pB))/pi

            tB = 2*atan(sqrt(y),sqrt(x))/pi

        end

    end

    pBv[1] = pB
    pBv[2] = cospi(tB)
    pBv[3] = hB
    pBv[4] = tB

    if isnan(pBv[1])
        println("pB=$pB")
        println("w=$w")
        println("p=$p")
        println("t=$t")
    end

end

"""
    BackwardsLorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64)

    Calculates the backwards Lorentz boost (in z direction) from the COM frame to the lab frame for the momentum vector `pv` of a particle with (normalised) mass `m` by the rapidity `w` (must be positive!), placing the results in the boosted vector `pBv``. 

    Momentum vectors have coponents [p,cos(theta),phi,theta]
"""
function BackwardsLorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64)

    if w < 0e0
        println("w=$w")
        error("w must be positive")
    end

    p::Float64 = pv[1]
    t::Float64 = pv[4]
    h::Float64 = pv[3]

    hB::Float64 = h # phi unaffected by boost

    if m == 0e0

        ew = exp(w)

        pB = p*(ew*cospi(t/2)^2+sinpi(t/2)^2/ew)
        tB = 2*atan(tanpi(t/2)/ew)/pi

    elseif (z=p/m) <= 1e-4 # small p approximation 
        # valid to at fifth order in p/m

        # pB
        pB = m*sinh(w) 
        pB += m*cospi(t)*cosh(w)*z 
        pB += m*(-cospi(2*t) + cosh(2*w))*csch(w)/4 * z^2 
        pB -= m*(cospi(t)*coth(w)*csch(w)*sinpi(t)^2)/2 * z^3 
        pB -= m*((3*cospi(4*t) + (-2 - 4*cospi(2*t) + 2*cospi(4*t))*cosh(2*w) + cosh(4*w))*csch(w)^3)/64 * z^4

        val = sinpi(t)*csch(w)/2 * z
        val -= sinpi(2*t)*coth(w)*csch(w)/4 * z^2
        val += sinpi(t)*csch(w)^3*(3+cospi(2*t)*(3+2*cosh(2*w)))/16 * z^3

        tB = 2*atan(val)/pi
        
        #tB = 2*atan(p*sinpi(t),2(sign(w)*sinh(w)))/pi

    elseif z >= 1e4 # large p approximation
        # valid to at fifth order in m/p

        # pB 
        tmp = exp(w)*cospi(t/2)^2 + exp(-w)*sinpi(t/2)^2
        pB = m*tmp * z 
        pB -= m*(1/tmp-cosh(w))/2 / z
        tmp2 = exp(2*w)*cospi(t/2)^2 + exp(-2*w)*sinpi(t/2)^2
        pB += m*(tmp2/tmp^3-cosh(w))/8 /z^3      

        ew = exp(-w)
        ttd2 = tanpi(t/2)

        val = ew*ttd2
        val -= ew*sinpi(t)*sinh(w) / ((1+cospi(t))*tmp) /2 / z^2
        val += ew*ttd2 * sinh(w) * (-1 + 5*cosh(2*w) + 2*(2*cospi(t) + cospi(2*t))*sinh(w)^2 - (2 + 6*cospi(t))*sinh(2*w)) / (-2*sinh(w)*cospi(t)+2*cosh(w))^3 / z^4

        #tB = 2*atan(ew*ttd2*sqrt(1+(m^2)/(p^2)*(sinh(w))/(cosh(w)-cospi(t)*sinh(w))))/pi

        tB = 2*atan(val)/pi

    else # no approximation

        test = (cosh(w)*sqrt(p^2+m^2)-sinh(w)*p*cospi(t))^2-m^2
        if test < 0e0
            println("")
            println("test =$test")
            println("w=$w")
            println("p=$p")
            println("t=$t")
            println("")
        end

        pB = sqrt((cosh(w)*sqrt(p^2+m^2)+sinh(w)*p*cospi(t))^2-m^2)

        if pB < 0e0
            println("pB < 0e0, $pB, setting pB=0e0")
            pB = 0e0
        end
        
        if w >= 6e0
            # approximate atan() to fifth order in exp(w)
            ew = exp(w)
            val = p*sinpi(t)/(sqrt(p^2+m^2)+p*cospi(t)) / ew
            val += m^2*p*sinpi(t)/(sqrt(p^2+m^2)+p*cospi(t))^3 / ew^3
            val += m^2*p*sinpi(t)*(m^2-p^2*sinpi(t)^2)/(sqrt(p^2+m^2)+p*cospi(t))^5 / ew^5
            #val += m^2*p*sinpi(t)*()/(sqrt(p^2+m^2)+p*cospi(t))^7 / ew^7

            tB = 2*atan(val)/pi

        else # no approximation

            y = pB+(-exp(w)*(sqrt(p^2+m^2)+p*cospi(t))/2+exp(-w)*(sqrt(p^2+m^2)-p*cospi(t))/2)
            x = pB-(-exp(w)*(sqrt(p^2+m^2)+p*cospi(t))/2+exp(-w)*(sqrt(p^2+m^2)-p*cospi(t))/2)

            #y = pB + sqrt(p^2+m^2)*sinh(w) - p*cosh(w)*cospi(t) 
            #x = pB + p*cosh(w)*cospi(t) - sqrt(p^2+m^2)*sinh(w)

            if y<0e0
                #println("")
                #println("backwards boost")
                #println("warning: y<0,=$y, setting y=0e0")
                #println("x=$x")
                #println("y=$y")
                #println("pB=$pB")
                #println("w=$w")
                #println("p=$p")
                #println("t=$t")
                #println("")
                y = 0e0
                #tB = acos((exp(w)*(sqrt(p^2+m^2)p*cospi(t))/2-exp(-w)*(-p*cospi(t)+sqrt(p^2+m^2))/2)/pB)/pi
                #println("tB = $tB")
            end

            if x < 0e0
                x=0e0
            end

            # acos(x) can lead to a domain error with x > 1 due to floating point precision, atan approach better 
            #tB = acos((cosh(w)*p*cospi(t)+sinh(w)*sqrt(p^2+m^2))/pB)/pi

            #tBtest = 2*atan(sqrt(-p*cosh(w)*cospi(t)+sqrt(p^2+m^2)*sinh(w)+pB),sqrt(+p*cosh(w)*cospi(t)-sqrt(p^2+m^2)*sinh(w)+pB))/pi

            tB = 2*atan(sqrt(y),sqrt(x))/pi

        end

    end

    pBv[1] = pB
    pBv[2] = cospi(tB)
    pBv[3] = hB
    pBv[4] = tB

    if isnan(pBv[1])
        println("pB=$pB")
        println("w=$w")
        println("p=$p")
        println("t=$t")
    end

end

"""
    Weight!(pv::Vector{Float64},w::Float64)

Computes the doppler factor aka cosine angle jacobian for a particle of mass `m` between the CM frame and lab frame, for a boost with rapididty `w` in the z direction (w must be +ve). The result is stored in the fifth entry of the vector `pv`. Which is expected to have components [p,cos(theta),phi,theta,doppler].
"""
function Weight!(pv::Vector{Float64},m::Float64,w::Float64)

    p::Float64 = pv[1]
    t::Float64 = pv[4]

    if m == 0e0

        ew = exp(w)

        DF = 1/(ew*cospi(t/2)^2+sinpi(t/2)^2/ew)^2

    elseif (z=p/m) <= 1e-4 # small p approximation 
        # valid to at least second order in p/m

        # COULD be more accurate
        DF = cospi(t)/sinh(w)^2 * z^2
        DF -= (4-3*sinpi(t)^2)*cosh(w)/(2*sinh(w)^3) * z^3

    elseif z >= 1e4 # large p approximation
        # valid to at least second order in m/p

        # COULD be more accurate
        ew = exp(w)
        # tmp is cosh(w)+cospi(t)*sinh(w)
        tmp = ew*cospi(t/2)^2+sinpi(t/2)^2/ew
        # tmp2 is cosh(2w)-cospi(t)*sinh(2w)
        tmp2 = exp(2*w)*sinpi(t/2)^2+cospi(t/2)^2*exp(-2*w)
        DF = 1/tmp^2
        DF += (-1+2*tmp2-tmp^2)/2/tmp^4 / z^2


    else # no approximation

        #=if w >= 8e0

            ew = exp(w)
            tmp = (sqrt(p^2+m^2)+p*cospi(t))

            DF = 4*p^2*(p+sqrt(p^2+m^2)*cospi(t))/tmp^3 / ew^2
            DF += (8*p^2*(2*m^2*p - p^3 + (m^2 - p^2)*sqrt(m^2 + p^2)*cospi(t) + (-m^2*p + p^3)*cospi(t)^2 + p^2*sqrt(m^2 + p^2)*cospi(t)^3))/tmp^5 / ew^4

        else=#

            tmp1 = exp(w)*(p+sqrt(m^2+p^2)*cospi(t))/2 + exp(-w)*(p-sqrt(p^2+m^2)*cospi(t))/2
            tmp2 = exp(w)*(sqrt(m^2+p^2)+p*cospi(t))/2 + exp(-w)*(sqrt(p^2+m^2)-p*cospi(t))/2

            DF = p^2*tmp1
            DF /= (-m^2+(tmp2)^2)^(3/2)

        #end

    end

    if isinf(DF)
        println("DF = $DF")
        println("w = $w")
        println("p = $p")
        println("t = $t")
        error("DF inf")
    end

    pv[5] = abs(DF)

end

"""
    DopplerFactor!(pv::Vector{Float64},m::Float64,w::Float64)

Computes the doppler factor aka cosine angle jacobian for a particle of mass `m` between the CM frame and lab frame, for a boost with rapididty `w` in the z direction (w must be +ve). The result is stored in the fifth entry of the vector `pv`. Which is expected to have components [p,cos(theta),phi,theta,doppler].
"""
function DopplerFactor!(pv::Vector{Float64},m::Float64,w::Float64)

    p::Float64 = pv[1]
    t::Float64 = pv[4]

    if m == 0e0

        ew = exp(w)

        DF = 1/(ew*cospi(t/2)^2+sinpi(t/2)^2/ew)^2

    elseif (z=p/m) <= 1e-4 # small p approximation 
        # valid to at least second order in p/m

        # COULD be more accurate
        DF = cospi(t)/sinh(w)^2 * z^2
        DF -= (4-3*sinpi(t)^2)*cosh(w)/(2*sinh(w)^3) * z^3

    elseif z >= 1e4 # large p approximation
        # valid to at least second order in m/p

        # COULD be more accurate
        ew = exp(w)
        # tmp is cosh(w)+cospi(t)*sinh(w)
        tmp = ew*cospi(t/2)^2+sinpi(t/2)^2/ew
        # tmp2 is cosh(2w)-cospi(t)*sinh(2w)
        tmp2 = exp(2*w)*sinpi(t/2)^2+cospi(t/2)^2*exp(-2*w)
        DF = 1/tmp^2
        DF += (-1+2*tmp2-tmp^2)/2/tmp^4 / z^2


    else # no approximation

        #=if w >= 8e0

            ew = exp(w)
            tmp = (sqrt(p^2+m^2)+p*cospi(t))

            DF = 4*p^2*(p+sqrt(p^2+m^2)*cospi(t))/tmp^3 / ew^2
            DF += (8*p^2*(2*m^2*p - p^3 + (m^2 - p^2)*sqrt(m^2 + p^2)*cospi(t) + (-m^2*p + p^3)*cospi(t)^2 + p^2*sqrt(m^2 + p^2)*cospi(t)^3))/tmp^5 / ew^4

        else=#

            tmp1 = exp(w)*(p+sqrt(m^2+p^2)*cospi(t))/2 + exp(-w)*(p-sqrt(p^2+m^2)*cospi(t))/2
            tmp2 = exp(w)*(sqrt(m^2+p^2)+p*cospi(t))/2 + exp(-w)*(sqrt(p^2+m^2)-p*cospi(t))/2

            DF = p^2*tmp1
            DF /= (-m^2+(tmp2)^2)^(3/2)

        #end

    end

    if isinf(DF)
        println("DF = $DF")
        println("w = $w")
        println("p = $p")
        println("t = $t")
        error("DF inf")
    end

    pv[5] = abs(DF)

end

"""
    RotateToCentre(pv::Vector{Float64},t::Float64,h::Float64)

    Rotates the momentum vector `pv` from the lab spherical coordinates to the spherical coordinates aligned with the centre of momentum velocity direction, The rotation angle is given by `t` theta from the lab z axis and `h` phi about the lab z axis.
"""
function RotateToCentre!(pv::Vector{Float64},tβ::Float64,hβ::Float64)

    tv = pv[4]
    hv = pv[3]

    if abs(tv-tβ) < 1e-7 && abs(hv-hβ) < 1e-7

        stv = sinpi(tv)
        stβ = sinpi(tβ)

        t = sqrt((tv-tβ)^2+(hv-hβ)^2*stv*stβ) # /pi not needed as already factored into approximation
        h = atan((hv-hβ)*stv,tv-tβ)/pi
        h = mod(h,2)
        #println("here, $(tv/tβ), $(hv/hβ)")

    else

        (stβ,ctβ) = sincospi(tβ)
        (stv,ctv) = sincospi(tv)
        (shvhβ, chvhβ) = sincospi(hv-hβ)
        
        # theta 
        t = acos(ctv*ctβ+chvhβ*stv*stβ)/pi # chvhβ is inaccurate when pv[3] is approx hβ i.e. closely aligned
        # phi
        x = -ctv*stβ+chvhβ*stv*ctβ
        y = shvhβ*stv
        h = mod(atan(y,x)/pi,2)
        #println("there,$(tv/tβ), $(hv/hβ)")

    end

    pv[4] = t
    pv[3] = h
    pv[2] = cospi(pv[4])

    return nothing

end

"""
    RotateToLab(pv::Vector{Float64},t::Float64,h::Float64)

    Rotates the momentum vector `pv` from the centre of momentum vector aligned spherical coordinates to the spherical coordinates aligned with the lab frame. This is the inverse rotation from `RotateToCentre` The rotation angle is given by `t` theta from the lab z axis and `h` phi about the lab z axis.
"""
function RotateToLab!(pv::Vector{Float64},t::Float64,h::Float64)

    (st,ct) = sincospi(t)
    (stv,ctv) = sincospi(pv[4])
    (shv,chv) = sincospi(pv[3])
    (sh,ch) = sincospi(h)
    
    # theta 
    pv[4] = acos(ctv*ct - stv*st*chv)/pi
    pv[2] = cospi(pv[4])
    # phi
    x = -sh*shv*stv+ch*(chv*ct*stv+ctv*st)
    y = ch*shv*stv+sh*(chv*ct*stv+ctv*st)
    pv[3] = mod2pi(atan(y,x))/pi

    return nothing

end

function RotateToLab!(p3v::Vector{Float64},p4v::Vector{Float64},t::Float64,h::Float64)

    (st,ct) = sincospi(t)
    (sh,ch) = sincospi(h)
    (st3,ct3) = sincospi(p3v[4])
    (sh3,ch3) = sincospi(p3v[3])
    (st4,ct4) = sincospi(p4v[4])
    (sh4,ch4) = sincospi(p4v[3])

    # theta3 
    p3v[4] = acos(ct3*ct - st3*st*ch3)/pi
    p3v[2] = cospi(p3v[4])
    # phi3
    x = -sh*sh3*st3+ch*(ch3*ct*st3+ct3*st)
    y = ch*sh3*st3+sh*(ch3*ct*st3+ct3*st)
    p3v[3] = mod2pi(atan(y,x))/pi

    # theta4 
    p4v[4] = acos(ct4*ct - st4*st*ch4)/pi
    p4v[2] = cospi(p4v[4])
    # phi4
    x = -sh*sh4*st4+ch*(ch4*ct*st4+ct4*st)
    y = ch*sh4*st4+sh*(ch4*ct*st4+ct4*st)
    p4v[3] = mod2pi(atan(y,x))/pi

    return nothing

end

