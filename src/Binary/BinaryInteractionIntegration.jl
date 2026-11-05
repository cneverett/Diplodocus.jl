#= Script for running the ST integration and returning data arrays =#

#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidance of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

"""
    BinaryInteractionIntegration(Setup)

Function to build/load the collision matrices, then run the Monte Carlo integration of the Gain and Loss arrays in a multi-threaded environment. The function will run the Monte Carlo integration in parallel across the number of threads specified in the global variable `numThreads``. The function will then evaluate the Gain and Loss matrices using a weighted averaging method and save the results to a file specified by `fileLocation` and `fileName` in the input `Setup`.
"""
function BinaryInteractionIntegration(Setup::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},Tuple{Int64,Int64,Int64,Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

    (Parameters,bins,scale,numLoss,numGain,numThreads,fileLocation,fileName) = Setup
    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters
    (p1loc_low,p1loc_up,p2loc_low,p2loc_up) = bins 

    # ====================================== #

    # ========= u Symmetric Grids? =============== #
        
        symmetric_grid::Bool = (u1_grid == "u" || u1_grid == "b") && (u2_grid == "u" || u2_grid == "b") && (u3_grid == "u" || u3_grid == "b") && (u4_grid == "u" || u4_grid == "b")

    # ====================================== #

    # ===== Are states Distinguishable ===== #

        Indistinguishable_12::Bool = name1 == name2
        Indistinguishable_34::Bool = name3 == name4

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ========= #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma = getfield(DiplodocusCollisions,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt = getfield(DiplodocusCollisions,name_dsigmadt)

    # ============================================================================ #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p1_num] 

    # ===================================== #

    # =========== Load Old and Generate New Arrays ========= #

        println(stdout,"Loading Old and Allocating New Sampling Arrays")
        flush(stdout)

        if isdir(fileLocation) == false
            mkpath(fileLocation)
        end
                
        filePath = joinpath(fileLocation,fileName)

        (OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix,CorrectedChunkGainMatrix3,CorrectedChunkGainMatrix4,CorrectedChunkLossMatrix) = OldMonteCarloArraysBinary(Parameters,filePath)

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        println(stdout,"Running Monte Carlo Integration")
        flush(stdout)

        indices::Vector{CartesianIndex{2}} = CartesianIndices((p1loc_low:p1loc_up,p2loc_low:p2loc_up))[1:end]
        shuffle!(indices) # better balances workload between threads
        length_indices::Int64 = length(indices)

        indices_chunks = [indices[r] for r in index_chunks(indices;n=numThreads)] # split indices into chunks for each thread

        @time if numThreads == 1
            numProgress = length_indices
            prog = Progress(numProgress)
            # Run in serial if only one thread, easier to use for debugging
            BinaryMonteCarlo_Debug!(OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix,CorrectedChunkGainMatrix3,CorrectedChunkGainMatrix4,CorrectedChunkLossMatrix,sigma,dsigmadt,Parameters,numLoss,numGain,indices_chunks[1],scale,prog,1)
            finish!(prog)
        else
            numProgress = length_indices
            prog = Progress(numProgress)
            workers = [BinaryMonteCarlo!(OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix,CorrectedChunkGainMatrix3,CorrectedChunkGainMatrix4,CorrectedChunkLossMatrix,sigma,dsigmadt,Parameters,numLoss,numGain,indices_chunks[thread],scale,prog,thread) for thread in 1:numThreads]
            wait.(workers) # Allow all workers to finish
            finish!(prog)
        end

    # ===================================== #

    # ============= Error Estimates ======= # 

        println(stdout,"Calculating Error Estimates")
        flush(stdout)

        ErrorOutput =  DoesConserve((Parameters,OldGainMatrix3,OldGainMatrix4,OldLossMatrix))

    # ===================================== #

    # ========== Save Error estimate and Parameters============== #

        jldopen(joinpath(filePath,"data.jld2"),"w") do f # creates file and overwrites previous file if one existed
            write(f,"Parameters",Parameters)
            write(f,"ErrorEstimates",ErrorOutput)
        end

    # ===================================== #

    GC.gc()

    return nothing

end # function