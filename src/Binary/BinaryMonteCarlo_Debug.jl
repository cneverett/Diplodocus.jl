
#= 
This module provides functions for MonteCarlo Integration of S and T Matrices
=#

"""
    BinaryMonteCarlo_Debug(Arrays,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,scale,prog,thread_id)

Debug version of `BinaryMonteCarlo` which is only run when `numThreads=1`, useful for testing. Function performs the Monte-Carlo sampling of incoming and outgoing particle states for binary interactions.
"""
function BinaryMonteCarlo_Debug!(OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix,CorrectedGainMatrix3,CorrectedGainMatrix4,CorrectedLossMatrix,sigma::Function,dsigmadt::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},numLoss::Int64,numGain::Int64,indices::Vector{CartesianIndex{2}},scale::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},prog::Progress,thread_id::Int64)

    # Set Parameters
    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,h1_grid_st,h1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,h2_grid_st,h2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,h3_grid_st,h3_num,p4_low,p4_up,p4_grid_st,p4_num,u4_grid_st,u4_num,h4_grid_st,h4_num) = Parameters

    # allocate arrays for each thread
    p1v::Vector{Float64} = zeros(Float64,4)
    p2v::Vector{Float64} = zeros(Float64,4)
    p3v::Vector{Float64} = zeros(Float64,4)
    p3pv::Vector{Float64} = zeros(Float64,4)
    p4v::Vector{Float64} = zeros(Float64,4)
    p4pv::Vector{Float64} = zeros(Float64,4)
    GainVal::Float64 = 0e0
    GainValp::Float64 = 0e0
    LossVal::Float64 = 0e0
    p_physical::Bool = true
    pp_physical::Bool = true
    NumStates::Int64 = 2
    prob3::Float64 = 0e0
    prob4::Float64 = 0e0
    w3::Float64 = 0e0
    w4::Float64 = 0e0
    t::Float64 = 0e0
    h::Float64 = 0e0

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0
    h1loc::Int64 = 0
    h2loc::Int64 = 0
    p3loc::Int64 = 0
    u3loc::Int64 = 0
    h3loc::Int64 = 0
    p3ploc::Int64 = 0
    u3ploc::Int64 = 0
    h3ploc::Int64 = 0
    p4loc::Int64 = 0
    u4loc::Int64 = 0
    h4loc::Int64 = 0
    p4ploc::Int64 = 0
    u4ploc::Int64 = 0
    h4ploc::Int64 = 0
    loc12::CartesianIndex{4} = CartesianIndex(0,0,0,0)

    p1_grid::GridType = Grid_String_to_Type(p1_grid_st)
    p2_grid::GridType = Grid_String_to_Type(p2_grid_st)
    p3_grid::GridType = Grid_String_to_Type(p3_grid_st)
    p4_grid::GridType = Grid_String_to_Type(p4_grid_st)
    u1_grid::GridType = Grid_String_to_Type(u1_grid_st)
    u2_grid::GridType = Grid_String_to_Type(u2_grid_st)
    u3_grid::GridType = Grid_String_to_Type(u3_grid_st)
    u4_grid::GridType = Grid_String_to_Type(u4_grid_st)
    h1_grid::GridType = Grid_String_to_Type(h1_grid_st)
    h2_grid::GridType = Grid_String_to_Type(h2_grid_st)
    h3_grid::GridType = Grid_String_to_Type(h3_grid_st)
    h4_grid::GridType = Grid_String_to_Type(h4_grid_st)

    symmetric_grid::Bool = (u1_grid == "u" || u1_grid == "b") && (u2_grid == "u" || u2_grid == "b") && (u3_grid == "u" || u3_grid == "b") && (u4_grid == "u" || u4_grid == "b")

    Indistinguishable_12::Bool = name1 == name2
    Indistinguishable_34::Bool = name3 == name4

    u1_r::Vector{Float64} = bounds(u_low,u_up,u1_num,u1_grid_st)
    h1_r::Vector{Float64} = bounds(h_low,h_up,h1_num,h1_grid_st)
    u2_r::Vector{Float64} = bounds(u_low,u_up,u2_num,u2_grid_st)
    h2_r::Vector{Float64} = bounds(h_low,h_up,h2_num,h2_grid_st)

    u1_up::Float64 = 0.0
    u1_low::Float64 = 0.0
    h1_up::Float64 = 0.0
    h1_low::Float64 = 0.0
    u2_up::Float64 = 0.0
    u2_low::Float64 = 0.0
    h2_up::Float64 = 0.0
    h2_low::Float64 = 0.0

    # Angle / Momentum Ranges
    u3val = bounds(u_low,u_up,u3_num,u3_grid_st)
    u4val = bounds(u_low,u_up,u4_num,u4_grid_st)
    h3val = bounds(h_low,h_up,h3_num,h3_grid_st).*pi
    h4val = bounds(h_low,h_up,h4_num,h4_grid_st).*pi

    # Struct for gain correction
    GainCorrectionTmp = GainCorrectionStruct(Parameters)
    
    # local arrays are size of each chunk of Zarr that store data
    # for gain arrays this is p3_num+2,u3_num,h3_num,u1_num,h1_num,u2_num,h2_num
    # as p1 and p2 are what define the chunks
    # for loss arrays this is u1_num,h1_num,u2_num,h2_num as p1 and p2 are what define the chunks
    # these local arrays are in-memory and per thread
    # tallies for Gain terms are [k_underflow,k1,k2,k3,...,kn,k_overflow,N]
    ChunkGainTotal3::Array{Float64,7} = zeros(Float64,(p3_num+2),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkGainTally3::Array{UInt32,7} = zeros(UInt32,(p3_num+3),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkGainMatrix3::Array{Float64,7} = zeros(Float64,(p3_num+2),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkLossTotal::Array{Float64,4} = zeros(Float64,u1_num,h1_num,u2_num,h2_num)
    ChunkLossTally::Array{UInt32,4} = zeros(UInt32,u1_num,h1_num,u2_num,h2_num)
    if m3 != m4
        ChunkGainTotal4::Array{Float64,7} = zeros(Float64,(p4_num+2),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
        ChunkGainTally4::Array{UInt32,7} = zeros(UInt32,(p4_num+3),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
        ChunkGainMatrix4::Array{Float64,7} = zeros(Float64,(p4_num+2),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
    end
    ChunkLossMatrix::Array{Float64,4} = zeros(Float64,u1_num,h1_num,u2_num,h2_num)

    # N values are last element of the tally array
    ChunkGainTally3_N = @view(ChunkGainTally3[end,:,:,:,:,:,:])
    # K value are all but last element of the tally array
    ChunkGainTally3_K = @view(ChunkGainTally3[1:end-1,:,:,:,:,:,:])
    if m3 != m4 
        ChunkGainTally4_N = @view(ChunkGainTally4[end,:,:,:,:,:,:])
        ChunkGainTally4_K = @view(ChunkGainTally4[1:end-1,:,:,:,:,:,:])
    end 

    # local arrays for faster memory access that will be written to chunk arrays 
    LocalGainTotal3::Array{Float64,3} = zeros(Float64,size(ChunkGainTotal3)[1:3])
    LocalGainTally3::Array{UInt32,3} = zeros(UInt32,size(ChunkGainTally3)[1:3])
    LocalGainTally3N = @view(LocalGainTally3[end,:,:])
    if m3 != m4
        LocalGainTotal4::Array{Float64,3} = zeros(Float64,size(ChunkGainTotal4)[1:3])
        LocalGainTally4::Array{UInt32,3} = zeros(UInt32,size(ChunkGainTally4)[1:3])
        LocalGainTally4N = @view(LocalGainTally4[end,:,:])
    end

    # old chunk arrays 
    OldChunkGainMatrix3Full::Array{Float64,9} = zeros(Float64,(p3_num+2),u3_num,h3_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkGainMatrix4Full::Array{Float64,9} = zeros(Float64,(p4_num+2),u4_num,h4_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkLossMatrixFull::Array{Float64,6} = zeros(Float64,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkGainWeights3Full::Array{Float64,9} = zeros(Float64,(p3_num+2),u3_num,h3_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkGainWeights4Full::Array{Float64,9} = zeros(Float64,(p4_num+2),u4_num,h4_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkLossTallyFull::Array{UInt32,6} = zeros(UInt32,1,u1_num,h1_num,1,u2_num,h2_num)

    OldChunkGainMatrix3 = @view(OldChunkGainMatrix3Full[:,:,:,1,:,:,1,:,:])
    OldChunkGainMatrix4 = @view(OldChunkGainMatrix4Full[:,:,:,1,:,:,1,:,:])
    OldChunkLossMatrix = @view(OldChunkLossMatrixFull[1,:,:,1,:,:])
    OldChunkGainWeights3 = @view(OldChunkGainWeights3Full[:,:,:,1,:,:,1,:,:])
    OldChunkGainWeights4 = @view(OldChunkGainWeights4Full[:,:,:,1,:,:,1,:,:])
    OldChunkLossTally = @view(OldChunkLossTallyFull[1,:,:,1,:,:])

    for index in eachindex(indices)

        # define chunk indices of Zarr
        p1loc = indices[index][1]
        p2loc = indices[index][2]

        gain3loc = CartesianIndex(1,1,1,p1loc,1,1,p2loc,1,1)
        gain4loc = CartesianIndex(1,1,1,p1loc,1,1,p2loc,1,1)
        lossloc = CartesianIndex(p1loc,1,1,p2loc,1,1)

        # Load old chunk arrays from Zarr
        readblock_singlechunk!(OldChunkGainMatrix3Full,OldGainMatrix3,gain3loc)
        readblock_singlechunk!(OldChunkGainMatrix4Full,OldGainMatrix4,gain4loc)
        readblock_singlechunk!(OldChunkLossMatrixFull,OldLossMatrix,lossloc)
        readblock_singlechunk!(OldChunkGainWeights3Full,OldGainWeights3,gain3loc)
        readblock_singlechunk!(OldChunkGainWeights4Full,OldGainWeights4,gain4loc)
        readblock_singlechunk!(OldChunkLossTallyFull,OldLossTally,lossloc)

        # reset in-memory local chunk arrays to zero 
        fill!(ChunkGainTotal3,Float64(0))
        fill!(ChunkGainTally3,UInt32(0))
        fill!(ChunkLossTotal,Float64(0))
        fill!(ChunkLossTally,UInt32(0))
        if m3 != m4
            fill!(ChunkGainTotal4,Float64(0))
            fill!(ChunkGainTally4,UInt32(0))
        end

        for (idx_scale,scale_val) in enumerate(scale)

        for u1loc in 1:u1_num, h1loc in 1:h1_num, u2loc in 1:u2_num, h2loc in 1:h2_num

            u1_up = u1_r[u1loc+1]
            u1_low = u1_r[u1loc]
            h1_up = h1_r[h1loc+1]
            h1_low = h1_r[h1loc]
            u2_up = u2_r[u2loc+1]
            u2_low = u2_r[u2loc]
            h2_up = h2_r[h2loc+1]
            h2_low = h2_r[h2loc]
            loc12 = CartesianIndex(u1loc,h1loc,u2loc,h2loc)

            ChunkGainTally3View = @view(ChunkGainTally3[:,:,:,loc12])
            ChunkGainTotal3View = @view(ChunkGainTotal3[:,:,:,loc12])
            if m3 != m4
                ChunkGainTally4View = @view(ChunkGainTally4[:,:,:,loc12])
                ChunkGainTotal4View = @view(ChunkGainTotal4[:,:,:,loc12])
            end

        for _ in 1:numLoss # sample incoming sates
        
            # generate p1 and p2 vectors initially as to not have to re-calculate
            RPointSphereCosThetaPhiBounds!(p1v,u1_low,u1_up,h1_low,h1_up)
            RPointSphereCosThetaPhiBounds!(p2v,u2_low,u2_up,h2_low,h2_up)
            RPointLogMomentum!(p1v,p1_up,p1_low,p1_num,p1loc)
            RPointLogMomentum!(p2v,p2_up,p2_low,p2_num,p2loc)

            # reset local array tallies 
            fill!(LocalGainTally3,UInt32(0))
            if m3 != m4
                fill!(LocalGainTally4,UInt32(0))
            end

            # LossVal
            (LossVal,sBig,sSmol) = LossValue(p1v,p2v,sigma,m1,m2,m3,m4)

            if LossVal != 0e0 # i.e. it is a valid interaction state

                (w3,w4,t,h) = WeightedFactors(p1v,p2v,m1,m2,m3,m4,sBig,sSmol,scale_val)

                # reset local arrays for gain totals 
                fill!(LocalGainTotal3,Float64(0))
                if m3 != m4
                    fill!(LocalGainTotal4,Float64(0))
                end
                    
                for _ in 1:(numGain*p3_num*u3_num*h3_num)

                    prob3 = RPointSphereWeighted!(p3v,w3)  
                    RotateToLab!(p3v,t,h)
                    @. p3pv = p3v

                    # Calculate p3 value
                    (p_physical,pp_physical,NumStates) = MomentumValue!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4,p3_low,p3_up)

                    # Gain Array Tallies
                    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states. NOTE: This has been removed due to difference in sampling probability not being accounted for  
                    #if NumStates != 0
                        u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
                        LocalGainTally3[end,u3loc,h3loc] += UInt32(1)
                    #end

                    # Calculate Gain Array totals
                    if NumStates == 1
                        if p_physical
                            p3loc = locationUnderOver(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                            GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                        end
                    end

                    if NumStates == 2
                        if p_physical
                            p3loc = locationUnderOver(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                            GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                        end
                        if pp_physical
                            u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                            p3ploc = locationUnderOver(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                            GainValp = GainValue3(p3pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3ploc,u3ploc,h3ploc] += GainValp/prob3
                            LocalGainTally3[p3ploc,u3ploc,h3ploc] += UInt32(1)
                        end
                    end

                end # p3 outgoing loop

                if m3 != m4

                    for _ in 1:(numGain*p4_num*u4_num*h4_num)

                        prob4 = RPointSphereWeighted!(p4v,w4)
                        RotateToLab!(p4v,t,h)
                        @. p4pv = p4v

                        # Calculate p4 value
                        (p_physical,pp_physical,NumStates) = MomentumValue!(p4v,p4pv,p2v,p1v,m2,m1,m4,m3,p4_low,p4_up)

                        # S Array Tallies
                        # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states. NOTE: This has been removed due to difference in sampling probability not being accounted for 
                        u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                        h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
                        #u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                        #h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
                        LocalGainTally4[end,u4loc,h4loc] += UInt32(1)
                        #LocalGainTally4[end,u4locMirror,h4locMirror] += UInt32(1)

                        # Calculate S Array totals
                        if NumStates == 1
                            if p_physical
                                p4loc = locationUnderOver(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                                GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                                LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                            end
                        end

                        if NumStates == 2
                            if p_physical
                                p4loc = locationUnderOver(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                                GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                                LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                            end
                            if pp_physical
                                u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                                h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                                if u4ploc != u4loc || h4ploc != h4loc
                                println("u and h loc don't match up")
                                end
                                println("here")
                                p4ploc = locationUnderOver(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                                GainValp = GainValue4(p4pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4ploc,u4ploc,h4ploc] += GainValp/prob4
                                LocalGainTally4[p4ploc,u4ploc,h4ploc] += UInt32(1)
                            end
                        end

                    end # p4 outgoing loop
            
                end

            else # no valid interaction state
                # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
                LocalGainTally3N .+= UInt32(1)
                if m3 != m4
                    LocalGainTally4N .+= UInt32(1)
                end
            end

            # assign local arrays to chunk arrays
            ChunkLossTotal[loc12] += LossVal
            ChunkLossTally[loc12] += UInt32(1)
            ChunkGainTally3View .+= LocalGainTally3
            if m3 != m4 
                ChunkGainTally4View .+= LocalGainTally4
            end
            if LossVal != 0e0
                ChunkGainTotal3View .+= LocalGainTotal3
                if m3 != m4
                    ChunkGainTotal4View .+= LocalGainTotal4
                end
            end

        end # numLoss loop

        end # u1,h1,u2,h2 loop

        # === Update Gain and Loss Matrices === #

            # calculate the gain and loss matrices
            @inbounds for i in axes(ChunkGainTotal3,1)
                tot = @view(ChunkGainTotal3[i,:,:,:,:,:,:])
                out = @view(ChunkGainMatrix3[i,:,:,:,:,:,:])
                for lin in eachindex(ChunkGainTally3_N)
                    c = ChunkGainTally3_N[lin]
                    out[lin] = c == 0 ? 0.0 : tot[lin] / c
                end
            end
            if m3 != m4
                @inbounds for i in axes(ChunkGainTotal4,1)
                    tot = @view(ChunkGainTotal4[i,:,:,:,:,:,:])
                    out = @view(ChunkGainMatrix4[i,:,:,:,:,:,:])
                    for lin in eachindex(ChunkGainTally4_N)
                        c = ChunkGainTally4_N[lin]
                        out[lin] = c == 0 ? 0.0 : tot[lin] / c
                    end
                end
            end
            @inbounds for lin in eachindex(ChunkLossTally)
                c = ChunkLossTally[lin]
                ChunkLossMatrix[lin] = c == 0 ? 0.0 : ChunkLossTotal[lin] / c
            end

            # Momentum space volume elements
            if m3 == m4
                MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3,u3val,h3val,Indistinguishable_12)
            else
                MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3,ChunkGainMatrix4,u3val,h3val,u4val,h4val,Indistinguishable_12)
            end

            # perform weighted average of old and new gain matrices and loss matrices
            # old arrays are modified in this process
            if m3 == m4
                WeightedAverageGainBinaryChunk!(ChunkGainMatrix3,OldChunkGainMatrix3,ChunkGainTally3_K,ChunkGainTally3_N,OldChunkGainWeights3)
                if Indistinguishable_34 == false # particles are distinguishable 
                    @. OldChunkGainMatrix4 = OldChunkGainMatrix3
                    @. OldChunkGainWeights4 = OldChunkGainWeights3
                end
            else
                WeightedAverageGainBinaryChunk!(ChunkGainMatrix3,OldChunkGainMatrix3,ChunkGainTally3_K,ChunkGainTally3_N,OldChunkGainWeights3,ChunkGainMatrix4,OldChunkGainMatrix4,ChunkGainTally4_K,ChunkGainTally4_N,OldChunkGainWeights4)
            end
            WeightedAverageLossBinaryChunk!(ChunkLossMatrix,OldChunkLossMatrix,ChunkLossTally,OldChunkLossTally)

        end # scale loop 

        # ===== Saving Unsymmetrised/Uncorrected Arrays ===== #

            Zarr.write_singlechunk_fastpath!(OldGainMatrix3,OldChunkGainMatrix3Full,gain3loc)
            Zarr.write_singlechunk_fastpath!(OldGainMatrix4,OldChunkGainMatrix4Full,gain4loc)
            Zarr.write_singlechunk_fastpath!(OldLossMatrix,OldChunkLossMatrixFull,lossloc)
            Zarr.write_singlechunk_fastpath!(OldGainWeights3,OldChunkGainWeights3Full,gain3loc)
            Zarr.write_singlechunk_fastpath!(OldGainWeights4,OldChunkGainWeights4Full,gain4loc)
            Zarr.write_singlechunk_fastpath!(OldLossTally,OldChunkLossTallyFull,lossloc)
        
        # ========= Apply Symmetries ========== # 

            # Apply Symmetries to the Gain and Loss Matrices, this does not affect the weighting of the average, which is already done in the previous step
            # TODO: improve this averaging to account for the weights and tallies to do a weighted average symmetrising
            GainLossPolarSymmetryMatrixBinaryChunk!(OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix,OldChunkGainWeights3,OldChunkGainWeights4,OldChunkLossTally,m1,m2,m3,m4,symmetric_grid)
            GainLossAzimuthalSymmetryMatrixBinaryChunk!(OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix,OldChunkGainWeights3,OldChunkGainWeights4,OldChunkLossTally,m1,m2,m3,m4)

        # ===== Generate Corrected Arrays ===== #

            GainCorrectionChunk!(Parameters,OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix,p1loc,p2loc,GainCorrectionTmp)

        # ========== Save Symmetrised/Corrected Chunks to Zarr ============== #

            Zarr.write_singlechunk_fastpath!(CorrectedGainMatrix3,OldChunkGainMatrix3Full,gain3loc)
            Zarr.write_singlechunk_fastpath!(CorrectedGainMatrix4,OldChunkGainMatrix4Full,gain4loc)
            Zarr.write_singlechunk_fastpath!(CorrectedLossMatrix,OldChunkLossMatrixFull,lossloc)

            println(stdout,"Completed MC loop for p1loc=$p1loc, p2loc=$p2loc")
            flush(stdout)
            
            # Update progress 
            next!(prog)

    end # indices loop

    return nothing

end # function 
