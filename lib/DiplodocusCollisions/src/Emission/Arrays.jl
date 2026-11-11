"""
    MonteCarloArraysEmission(Parameters)

Generates arrays for Monte Carlo sampling for emissive (1->23) interactions.
"""
function MonteCarloArraysEmission(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64})

    (name1,name2,name3,type,mu1,mu2,mu3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num) = Parameters

    GainTotal2::Array{Float64,6} = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num); 
    GainTotal3::Array{Float64,6} = zeros(Float64,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num); 
    LossTotal1::Array{Float64,3} = zeros(Float64,p1_num,u1_num,h1_num);

    GainTallyN2::Array{UInt32,6} = zeros(UInt32,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
    GainTallyN3::Array{UInt32,6} = zeros(UInt32,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
    LossTallyN1::Array{UInt32,3} = zeros(UInt32,p1_num,u1_num,h1_num);

    GainTallyK2::Array{UInt32,6} = zeros(UInt32,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
    GainTallyK3::Array{UInt32,6} = zeros(UInt32,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
    LossTallyK1::Array{UInt32,3} = zeros(UInt32,p1_num,u1_num,h1_num);

    GainMatrix2::Array{Float64,6} = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
    GainMatrix3::Array{Float64,6} = zeros(Float64,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
    LossMatrix1::Array{Float64,3} = zeros(Float64,p1_num,u1_num,h1_num);


    return (GainTotal2,GainTotal3,LossTotal1,GainTallyN2,GainTallyK2,GainTallyN3,GainTallyK3,LossTallyN1,LossTallyK1,GainMatrix2,GainMatrix3,LossMatrix1)

end


"""
    OldMonteCarloArraysEmission(Parameters)

Load/Generates arrays for Monte Carlo sampling for emissive (1->23) interactions for if there is/is not a previously saved sample.
"""
function OldMonteCarloArraysEmission(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64, Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64,Float64},filePath::String)

    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext) = Parameters
    fileExist = isfile(filePath)

    if fileExist

        println("Previous Monte Carlo arrays found, loading from file.")
        f = jldopen(filePath,"r+");
        OldGainTallyK2 = f["GainTallyK2"];
        OldGainTallyN2 = f["GainTallyN2"];
        OldGainMatrix2 = f["GainMatrix2"];
        OldGainTallyK3 = f["GainTallyK3"];
        OldGainTallyN3 = f["GainTallyN3"];
        OldGainMatrix3 = f["GainMatrix3"];
        OldLossTallyK1 = f["LossTallyK1"];
        OldLossTallyN1 = f["LossTallyN1"];
        OldLossMatrix1 = f["LossMatrix1"];
        close(f)

    else

        OldGainTallyK2::Array{UInt32,6} = zeros(UInt32,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
        OldGainTallyK3::Array{UInt32,6} = zeros(UInt32,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
        OldLossTallyK1::Array{UInt32,3} = zeros(UInt32,p1_num,u1_num,h1_num);

        OldGainTallyN2::Array{UInt32,6} = zeros(UInt32,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
        OldGainTallyN3::Array{UInt32,6} = zeros(UInt32,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
        OldLossTallyN1::Array{UInt32,3} = zeros(UInt32,p1_num,u1_num,h1_num);

        OldGainMatrix2::Array{Float64,6} = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
        OldGainMatrix3::Array{Float64,6} = zeros(Float64,p3_num,u3_num,h3_num,p1_num,u1_num,h1_num);
        OldLossMatrix1::Array{Float64,3} = zeros(Float64,p1_num,u1_num,h1_num);

    end

    return (OldGainTallyK2,OldGainTallyK3,OldLossTallyK1,OldGainTallyN2,OldGainTallyN3,OldLossTallyN1,OldGainMatrix2,OldGainMatrix3,OldLossMatrix1)

end