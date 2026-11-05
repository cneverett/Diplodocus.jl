"""
    EmissionInteractionIntegration(userSyncInputMultiThread)

Function to run the Monte Carlo integration of the S array in a serial environment. 
"""
function EmissionInteractionIntegration(Setup::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64},Tuple{Int64,Int64,Int64,Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},Int64,Int64,Int64,String,String})

    # ======= Load User Parameters ======= #

    (Parameters,bins,scale,numLoss,numGain,numThreads,fileLocation,fileName) = Setup;
    #(name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag) = Parameters;
    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext) = Parameters
    (p1loc_low,p1loc_up,p3loc_low,p3loc_up) = bins 

    # ==================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

    # ====================================== #

    # ======= Define Kernel Functions Based on Emission Type  ========= #

        name_kernel = Symbol(type*"Kernel")
        EmissionKernel = getfield(DiplodocusCollisions,name_kernel)

    # ============================================================================ #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p2_num]    

    # ===================================== #

    # ======== Load Old Arrays ========= #

        println("Loading Old and Allocating New Sampling Arrays")

        if isdir(fileLocation) == false
            mkpath(fileLocation)
        end
        
        filePath = joinpath(fileLocation,fileName)

        (OldGainTallyK2,OldGainTallyK3,OldLossTallyK1,OldGainTallyN2,OldGainTallyN3,OldLossTallyN1,OldGainMatrix2,OldGainMatrix3,OldLossMatrix1) = OldMonteCarloArraysEmission(Parameters,filePath)

        (GainTotal2,GainTotal3,LossTotal1,GainTallyN2,GainTallyK2,GainTallyN3,GainTallyK3,LossTallyN1,LossTallyK1,GainMatrix2,GainMatrix3,LossMatrix1) = MonteCarloArraysEmission(Parameters)

    # ================================= #

    # ===== Run MonteCarlo Integration ==== #

        println("Running Monte Carlo Integration")
        println("Ext Value = $Ext")

        for (ii,scale_val) in enumerate(scale)

            println("")
            println("scale = $scale_val, iteration = $ii out of $(length(scale))")
            println("")

            indices = CartesianIndices((p1loc_low:p1loc_up,p3loc_low:p3loc_up))
            length_indices::Int64 = length(indices)

            if length_indices/numThreads > 1.0
                length_div_threads::Int64 = ceil(Int64,length_indices/numThreads)
                index_range = Vector(0:length_div_threads:numThreads*length_div_threads)
                index_range[end] = length_indices
            else 
                index_range = Vector(0:length_indices)
            end

            # reset new arrays
            fill!(LossTotal1,Float64(0))
            fill!(GainTotal2,Float64(0))
            fill!(GainTotal3,Float64(0))
            fill!(LossTallyN1,UInt32(0))
            fill!(GainTallyN2,UInt32(0))
            fill!(GainTallyN3,UInt32(0))
            fill!(LossTallyK1,UInt32(0))
            fill!(GainTallyK2,UInt32(0))
            fill!(GainTallyK3,UInt32(0))

            if numThreads == 1
                numProgress = numLoss*index_range[end]*u1_num*h1_num
                prog = Progress(numProgress)
                EmissionMonteCarlo_Debug!(GainTotal2,GainTallyN2,GainTallyK2,GainTotal3,GainTallyN3,GainTallyK3,LossTotal1,LossTallyN1,LossTallyK1,ArrayOfLocks,EmissionKernel,Parameters,numLoss,numGain,indices[1:end],scale_val,prog,1)
                finish!(prog)
            else 
                numProgress = numLoss*index_range[1+1]*u1_num*h1_num
                prog = Progress(numProgress)
                workers  = [EmissionMonteCarlo!(GainTotal2,GainTallyN2,GainTallyK2,GainTotal3,GainTallyN3,GainTallyK3,LossTotal1,LossTallyN1,LossTallyK1,ArrayOfLocks,EmissionKernel,Parameters,numLoss,numGain,indices[index_range[thread]+1:index_range[thread+1]],scale_val,prog,thread) for thread in 1:(length(index_range)-1)]
                wait.(workers) # Allow all workers to finish#
                finish!(prog)
            end

    # ===================================== #


    # ===== Update Gain and Loss Matrices === #

            # Calculate Gain and Loss matrix
            GainMatrix2 = GainTotal2 ./ GainTallyN2
            GainMatrix3 = GainTotal3 ./ GainTallyN3
            LossMatrix1 = LossTotal1 ./ LossTallyN1

            replace!(GainMatrix2,NaN=>0e0)
            replace!(GainMatrix3,NaN=>0e0)
            replace!(LossMatrix1,NaN=>0e0)

            println("")
            println("Applying Momentum Space Factors")

            # Angle / Momentum Ranges
            p3val = bounds(p3_low,p3_up,p3_num,p3_grid)
            u3val = bounds(u_low,u_up,u3_num,u3_grid)
            h3val = bounds(h_low,h_up,h3_num,h3_grid).*pi

            # Apply Momentum space volume elements
            MomentumSpaceFactorsEmission!(LossMatrix1,GainMatrix2,GainMatrix3,p3val,u3val,h3val)

            println("Weighting average of New and Old Sampling Arrays")

            # old arrays are modified in this process 
            WeightedAverageGainEmission!(GainMatrix2,OldGainMatrix2,GainTallyN2,OldGainTallyN2,GainMatrix3,OldGainMatrix3,GainTallyN3,OldGainTallyN3)
            WeightedAverageLossEmission!(LossMatrix1,OldLossMatrix1,LossTallyN1,OldLossTallyN1)

        end # scale loop

        println("Applying Symmetries to Gain and Loss Matrices")

        # Apply symmetries 
        GainLossPolarSymmetryEmission!(OldGainMatrix2,OldGainMatrix3,OldLossMatrix1,OldGainTallyN2,OldGainTallyN3,OldLossTallyN1)
        GainLossAzimuthalSymmetryEmission!(OldGainMatrix2,OldGainMatrix3,OldLossMatrix1,OldGainTallyN2,OldGainTallyN3,OldLossTallyN1)

    # ===================================== #

    # ========== Save Arrays ============== #

        println("")
        println("Saving Arrays")

        jldopen(filePath,"w";compress=true) do f # creates file and overwrites previous file if one existed
            write(f,"GainTallyK2",OldGainTallyK2)
            write(f,"GainTallyN2",OldGainTallyN2)
            write(f,"GainMatrix2",OldGainMatrix2)

            write(f,"GainTallyK3",OldGainTallyK3)
            write(f,"GainTallyN3",OldGainTallyN3)
            write(f,"GainMatrix3",OldGainMatrix3)
            
            write(f,"LossTallyK1",OldLossTallyK1)
            write(f,"LossTallyN1",OldLossTallyN1)
            write(f,"LossMatrix1",OldLossMatrix1)

            write(f,"Parameters",Parameters)
        end

    # ===================================== #

end # function 