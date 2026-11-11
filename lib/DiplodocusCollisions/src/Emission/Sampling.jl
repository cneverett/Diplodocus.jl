function ImportanceSamplingSync!(p1v::Vector{Float64},p3v::Vector{Float64},localGainTally3::Array{UInt32,3},localGainTotal3::Array{Float64,3},SmallParameters,WeightFactors::Tuple{Float64,Float64,Float64})

    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,m1,m2,m3,z1,z2,z3,BMag) = SmallParameters

    (w,t,h) = WeightFactors

    prob = RPointSphereWeighted!(p3v,w) # sample angles aligned to p1v
    RotateToLab!(p3v,t,h)   # rotate to z aligned
    RPointLogMomentum!(p3v,p3_low,p3_up,p3_num)

    # calculate S value
    Sval = SyncKernel(p3v,p1v,m1,z1,BMag)

    # find S array location 
    p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)

    localGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
    localGainTotal3[p3loc,u3loc,h3loc] += Sval/prob

    return nothing

end