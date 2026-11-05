# ============== Chunked Versions ============= # 

"""
    GainLossPolarSymmetryMatrixBinaryChunk!(ChunkGainMatrix3,ChunkGainMatrix4,ChunkLossMatrix,ChunkGainWeights3,ChunkGainWeights4,ChunkLossTally,m1,m2,m3,m4,symmetric_grid::Bool)

Applies various physical polar angle symmetries to the Gain and Loss Matrices and Weights/Tally for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossPolarSymmetryMatrixBinaryChunk!(GainMatrix3::AbstractArray{Float32,7},GainMatrix4::AbstractArray{Float32,7},LossMatrix::AbstractArray{Float32,4},GainWeights3::AbstractArray{Float32,7},GainWeights4::AbstractArray{Float32,7},LossTally::AbstractArray{UInt32,4},m1::Float64,m2::Float64,m3::Float64,m4::Float64,symmetric_grid::Bool)

    # The Gain and Loss matrices are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    #   (This FIRST condition only applies if the grids are symmetric in u, e.g. uniform or binary grids)
    # SECOND: if the incident masses are equal (m1==m2) then Gain and Loss are symmetric to swapping the incident particles

    nu1::Int64 = size(GainMatrix3)[4]
    nu2::Int64 = size(GainMatrix3)[6]
    nu3::Int64 = size(GainMatrix3)[2]
    nu4::Int64 = size(GainMatrix4)[2]

    # Particle 3 Gain terms
    @inbounds for h2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,6), h1 in axes(GainMatrix3,5), u1 in axes(GainMatrix3,4), h3 in axes(GainMatrix3,3)

        if m1 == m2 && u1 > u2
            continue  # Skip entire inner loop for redundant swap case
        end
        
        u1_mir = nu1 - u1 + 1
        u2_mir = nu2 - u2 + 1

        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2)
            # Only process the "lower" half to avoid redundant computation
            if (symmetric_grid && u3 > nu3 - u3 + 1)
                continue  # Skip u3_mir equivalents; they were already processed
            end

            u3_mir = nu3 - u3 + 1

            if m1 == m2 && symmetric_grid # Both first and second Symmetry true

                # 4 element averaging two are u angles mirrors and two are particle 12 swaps 
                w_sum = (GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainWeights3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] + GainWeights3[p3,u3,h3,u2,h2,u1,h1] + GainWeights3[p3,u3_mir,h3,u2_mir,h2,u1_mir,h1])
                val_sum = (GainMatrix3[p3,u3,h3,u1,h1,u2,h2] * GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainMatrix3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] * GainWeights3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] + GainMatrix3[p3,u3,h3,u2,h2,u1,h1] * GainWeights3[p3,u3,h3,u2,h2,u1,h1] + GainMatrix3[p3,u3_mir,h3,u2_mir,h2,u1_mir,h1] * GainWeights3[p3,u3_mir,h3,u2_mir,h2,u1_mir,h1])

                avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                GainMatrix3[p3,u3,h3,u1,h1,u2,h2] = avg
                GainMatrix3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] = avg
                GainMatrix3[p3,u3,h3,u2,h2,u1,h1] = avg
                GainMatrix3[p3,u3_mir,h3,u2_mir,h2,u1_mir,h1] = avg

                if m3 == m4 # no need to recalculate average if same mass
                    GainMatrix4[p3,u3,h3,u1,h1,u2,h2] = avg
                    GainMatrix4[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] = avg
                    GainMatrix4[p3,u3,h3,u2,h2,u1,h1] = avg
                    GainMatrix4[p3,u3_mir,h3,u2_mir,h2,u1_mir,h1] = avg
                end

            elseif m1 == m2 # only second symmetry true

                # 2 element averaging two are particle 12 swaps
                w_sum = (GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainWeights3[p3,u3,h3,u2,h2,u1,h1])
                val_sum = (GainMatrix3[p3,u3,h3,u1,h1,u2,h2] * GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainMatrix3[p3,u3,h3,u2,h2,u1,h1] * GainWeights3[p3,u3,h3,u2,h2,u1,h1])

                avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                GainMatrix3[p3,u3,h3,u1,h1,u2,h2] = avg
                GainMatrix3[p3,u3,h3,u2,h2,u1,h1] = avg

                if m3 == m4 # no need to recalculate average if same mass
                    GainMatrix4[p3,u3,h3,u1,h1,u2,h2] = avg
                    GainMatrix4[p3,u3,h3,u2,h2,u1,h1] = avg
                end

            elseif symmetric_grid # only first symmetry true

                # 2 element averaging two are u angles mirrors
                w_sum = (GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainWeights3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2])
                val_sum = (GainMatrix3[p3,u3,h3,u1,h1,u2,h2] * GainWeights3[p3,u3,h3,u1,h1,u2,h2] + GainMatrix3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] * GainWeights3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2])

                avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                GainMatrix3[p3,u3,h3,u1,h1,u2,h2] = avg
                GainMatrix3[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] = avg

                if m3 == m4 # no need to recalculate average if same mass
                    GainMatrix4[p3,u3,h3,u1,h1,u2,h2] = avg
                    GainMatrix4[p3,u3_mir,h3,u1_mir,h1,u2_mir,h2] = avg
                end

            end

        end

    end

    if m3 != m4

        # Particle 4 Gain terms
        @inbounds for h2 in axes(GainMatrix4,7), u2 in axes(GainMatrix4,6), h1 in axes(GainMatrix4,5), u1 in axes(GainMatrix4,4), h4 in axes(GainMatrix4,3)

            if m1 == m2 && u1 > u2
                continue  # Skip entire inner loop for redundant swap case
            end
    
            u1_mir = nu1 - u1 + 1
            u2_mir = nu2 - u2 + 1

            for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2)
                # Only process the "lower" half to avoid redundant computation
                if (symmetric_grid && u4 > nu4 - u4 + 1)
                    continue  # Skip u4_mir equivalents; they were already processed
                end

                u4_mir = nu4 - u4 + 1

                if m1 == m2 && symmetric_grid # Both first and second Symmetry true

                    # 4 element averaging two are u angles mirrors and two are particle 12 swaps 
                    w_sum = (GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainWeights4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] + GainWeights4[p4,u4,h4,u2,h2,u1,h1] + GainWeights4[p4,u4_mir,h4,u2_mir,h2,u1_mir,h1])
                    val_sum = (GainMatrix4[p4,u4,h4,u1,h1,u2,h2] * GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainMatrix4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] * GainWeights4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] + GainMatrix4[p4,u4,h4,u2,h2,u1,h1] * GainWeights4[p4,u4,h4,u2,h2,u1,h1] + GainMatrix4[p4,u4_mir,h4,u2_mir,h2,u1_mir,h1] * GainWeights4[p4,u4_mir,h4,u2_mir,h2,u1_mir,h1])

                    avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                    GainMatrix4[p4,u4,h4,u1,h1,u2,h2] = avg
                    GainMatrix4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] = avg
                    GainMatrix4[p4,u4,h4,u2,h2,u1,h1] = avg
                    GainMatrix4[p4,u4_mir,h4,u2_mir,h2,u1_mir,h1] = avg

                elseif m1 == m2  # only second symmetry true

                    # 2 element averaging two are particle 12 swaps
                    w_sum = (GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainWeights4[p4,u4,h4,u2,h2,u1,h1])
                    val_sum = (GainMatrix4[p4,u4,h4,u1,h1,u2,h2] * GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainMatrix4[p4,u4,h4,u2,h2,u1,h1] * GainWeights4[p4,u4,h4,u2,h2,u1,h1])

                    avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                    GainMatrix4[p4,u4,h4,u1,h1,u2,h2] = avg
                    GainMatrix4[p4,u4,h4,u2,h2,u1,h1] = avg

                elseif symmetric_grid # only first symmetry true

                    # 2 element averaging two are u angles mirrors
                    w_sum = (GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainWeights4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2])
                    val_sum = (GainMatrix4[p4,u4,h4,u1,h1,u2,h2] * GainWeights4[p4,u4,h4,u1,h1,u2,h2] + GainMatrix4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] * GainWeights4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2])

                    avg = w_sum == 0f0 ? 0f0 : Float32(val_sum / w_sum) # avoid div by zero

                    GainMatrix4[p4,u4,h4,u1,h1,u2,h2] = avg
                    GainMatrix4[p4,u4_mir,h4,u1_mir,h1,u2_mir,h2] = avg

                end

            end

        end

    end

    # Particle Loss Terms
    @inbounds for h2 in axes(LossMatrix,4), u2 in axes(LossMatrix,3), h1 in axes(LossMatrix,2), u1 in axes(LossMatrix,1)

        if (symmetric_grid && u1 > nu1 - u1 + 1) || (m1 == m2 && (u1 > u2 || (u1 == u2 && h1 > h2)))
            continue
        end
        
        u1_mir = nu1 - u1 + 1
        u2_mir = nu2 - u2 + 1

        if m1 == m2 && symmetric_grid # Both first and second Symmetry true

            # 4 element averaging two are u angles mirrors and two are particle 12 swaps
            w_sum = (LossTally[u1,h1,u2,h2] + LossTally[u2,h2,u1,h1] + LossTally[u1_mir,h1,u2_mir,h2] + LossTally[u2_mir,h2,u1_mir,h1])
            v_sum = (LossMatrix[u1,h1,u2,h2] * LossTally[u1,h1,u2,h2]) + (LossMatrix[u2,h2,u1,h1] * LossTally[u2,h2,u1,h1]) + (LossMatrix[u1_mir,h1,u2_mir,h2] * LossTally[u1_mir,h1,u2_mir,h2]) + (LossMatrix[u2_mir,h2,u1_mir,h1] * LossTally[u2_mir,h2,u1_mir,h1])

            avg = w_sum == 0f0 ? 0f0 : Float32(v_sum / w_sum) # avoid div by zero

            LossMatrix[u1,h1,u2,h2] = avg
            LossMatrix[u2,h2,u1,h1] = avg
            LossMatrix[u1_mir,h1,u2_mir,h2] = avg
            LossMatrix[u2_mir,h2,u1_mir,h1] = avg

        elseif m1 == m2 # only second symmetry true

            # 2 element averaging two are particle 12 swaps
            w_sum = (LossTally[u1,h1,u2,h2] + LossTally[u2,h2,u1,h1])
            v_sum = (LossMatrix[u1,h1,u2,h2] * LossTally[u1,h1,u2,h2]) + (LossMatrix[u2,h2,u1,h1] * LossTally[u2,h2,u1,h1])

            avg = w_sum == 0f0 ? 0f0 : Float32(v_sum / w_sum) # avoid div by zero

            LossMatrix[u1,h1,u2,h2] = avg
            LossMatrix[u2,h2,u1,h1] = avg

        elseif symmetric_grid # only first symmetry true

            # 2 element averaging two are u angles mirrors
            w_sum = (LossTally[u1,h1,u2,h2] + LossTally[u1_mir,h1,u2_mir,h2])
            v_sum = (LossMatrix[u1,h1,u2,h2] * LossTally[u1,h1,u2,h2]) + (LossMatrix[u1_mir,h1,u2_mir,h2] * LossTally[u1_mir,h1,u2_mir,h2])

            avg = w_sum == 0f0 ? 0f0 : Float32(v_sum / w_sum) # avoid div by zero

            LossMatrix[u1,h1,u2,h2] = avg
            LossMatrix[u1_mir,h1,u2_mir,h2] = avg

        end

    end 

    return nothing

end

"""
    GainLossAzimuthalSymmetryMatrixBinaryChunk!(GainMatrix3,GainMatrix4,LossMatrix,GainWeights3,GainWeights4,LossTally,m1,m2,m3,m4)

Applies various physical azimuthal angle symmetries to the Gain and Loss Matrices and Weights/Tally for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossAzimuthalSymmetryMatrixBinaryChunk!(GainMatrix3::AbstractArray{Float32,7},GainMatrix4::AbstractArray{Float32,7},LossMatrix::AbstractArray{Float32,4},GainWeights3::AbstractArray{Float32,7},GainWeights4::AbstractArray{Float32,7},LossTally::AbstractArray{UInt32,4},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # The Gain and Loss matrices are symmetric in with respect to rotations of the azimuthal angle. 
    # If the azimuthal grid is uniform then we can apply this symmetry to all variations of the azimuthal bins.

    num_h1 = size(GainMatrix3,5)
    num_h2 = size(GainMatrix3,7)
    num_h3 = size(GainMatrix3,3)
    num_h4 = size(GainMatrix4,3)
    num_sections12 = lcm(num_h1,num_h2) # number of sections to divide the azimuthal bins into for averaging
    num_sections123 = lcm(num_h1,num_h2,num_h3) # number of sections to divide the azimuthal bins into for averaging
    num_sections124 = lcm(num_h1,num_h2,num_h4) # number of sections to divide the azimuthal bins into for averaging

    num_sec123divnum_h1 = num_sections123 / num_h1
    num_sec123divnum_h2 = num_sections123 / num_h2
    num_sec123divnum_h3 = num_sections123 / num_h3

    num_sec124divnum_h1 = num_sections124 / num_h1
    num_sec124divnum_h2 = num_sections124 / num_h2
    num_sec124divnum_h4 = num_sections124 / num_h4

    num_sec12divnum_h1 = num_sections12 / num_h1
    num_sec12divnum_h2 = num_sections12 / num_h2

    # Particle 3 Gain terms
    @inbounds for u2 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,4), u3 in axes(GainMatrix3,2), p3 in axes(GainMatrix3,1)

        for off2 in 0:num_sections123-1, off3 in 0:num_sections123-1 # loop over the maximum number of azimuthal bins for the particles
            
            tmp_matrix = zero(Float32)
            tmp_weight = zero(Float32)

            for h in 1:num_sections123

                h1 = mod(floor(Int64, h / (num_sec123divnum_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sec123divnum_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sec123divnum_h3)),num_h3) + 1

                tmp_matrix += GainMatrix3[p3,u3,h3,u1,h1,u2,h2] * GainWeights3[p3,u3,h3,u1,h1,u2,h2]
                tmp_weight += GainWeights3[p3,u3,h3,u1,h1,u2,h2]

            end

            val = tmp_weight == 0f0 ? 0f0 : Float32(tmp_matrix / tmp_weight) # avoid div by zero

            for h in 1:num_sections123

                h1 = mod(floor(Int64, h / (num_sec123divnum_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sec123divnum_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sec123divnum_h3)),num_h3) + 1

                GainMatrix3[p3,u3,h3,u1,h1,u2,h2] = val # average over number of sections/rotations times weight of each section
                if m3 == m4 # no need to recalculate average if same mass
                    GainMatrix4[p3,u3,h3,u1,h1,u2,h2] = val # average over number of sections/rotations times weight of each section
                end

            end

        end

    end

    if m3 != m4
        # Particle 4 Gain terms
        @inbounds for u2 in axes(GainMatrix4,6), u1 in axes(GainMatrix4,4), u4 in axes(GainMatrix4,2), p4 in axes(GainMatrix4,1)

            for off2 in 0:num_sections124-1, off4 in 0:num_sections124-1 # loop over the maximum number of azimuthal bins for the particles

                tmp_matrix = zero(Float32)
                tmp_weight = zero(Float32)

                for h in 1:num_sections124

                    h1 = mod(floor(Int64, h / (num_sec124divnum_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sec124divnum_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sec124divnum_h4)),num_h4) + 1

                    tmp_matrix += GainMatrix4[p4,u4,h4,u1,h1,u2,h2] * GainWeights4[p4,u4,h4,u1,h1,u2,h2]
                    tmp_weight += GainWeights4[p4,u4,h4,u1,h1,u2,h2]

                end

                val = tmp_weight == 0f0 ? 0f0 : Float32(tmp_matrix / tmp_weight) # avoid div by zero

                for h in 1:num_sections124 # loop over the maximum number of azimuthal bins for the particles

                    h1 = mod(floor(Int64, h / (num_sec124divnum_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sec124divnum_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sec124divnum_h4)),num_h4) + 1

                    GainMatrix4[p4,u4,h4,u1,h1,u2,h2] = val # average over number of sections/rotations times weight of each section

                end

            end

        end

    end

    # Particle Loss Terms
    @inbounds for u2 in axes(LossMatrix,3) , u1 in axes(LossMatrix,2)

        for off2 in 0:num_sections12-1 # loop over the maximum number of azimuthal bins for the particles

            tmp_matrix = zero(Float32)
            tmp_weight = zero(UInt32)

            for h in 1:num_sections12

                h1 = mod(floor(Int64, h / (num_sec12divnum_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sec12divnum_h2)),num_h2) + 1

                tmp_matrix += LossMatrix[u1,h1,u2,h2] * LossTally[u1,h1,u2,h2] # LossMatrix2 is just a permutation of LossMatrix1 so we only need to sum over one of them
                tmp_weight += LossTally[u1,h1,u2,h2]

            end

            val = tmp_weight == 0 ? 0f0 : Float32(tmp_matrix / tmp_weight) # avoid div by zero

            for h in 1:num_sections12 # loop over the maximum number of azimuthal bins for the particles

                h1 = mod(floor(Int64, h / (num_sec12divnum_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sec12divnum_h2)),num_h2) + 1

                LossMatrix[u1,h1,u2,h2] = val # average over number of sections/rotations times weight of each section

            end

        end

    end 

    return nothing

end


###### LEGACY VERSIONS BELOW ######

#=

"""
    GainLossPolarSymmetryBinary(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,Indistinguishable_34,m1,m2,m3,m4)

Applies various physical polar angle symmetries to the Gain and Loss terms for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossPolarSymmetryBinary!(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)
    # if the outgoing particles are indistinguishable or have identical masses (both conditions give m3==m4) then the Gain terms are identical assuming they have the same discretisation
    #if m3 == m4
    #    @. GainTotal3 = GainTotal3 + GainTotal4
    #    @. GainTotal4 = GainTotal3
    #    @. GainTally3 = GainTally3 + GainTally4
    #    @. GainTally4 = GainTally3
    #end

    # The Gain and Loss matrices are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    #   (This FIRST condition only applies if the grids are symmetric in u, e.g. uniform or binary grids)
    # SECOND: if the incident masses are equal (m1==m2) then Gain and Loss are symmetric to swapping the incident particles
    
    tmp_total3 = zeros(Float64,size(GainTotal3)[1:2])
    tmp_tally3 = zeros(Int64,size(GainTally3)[1:2])

    # Particle 3 Gain terms
    for h2 in axes(GainTotal3,9), u2 in axes(GainTotal3,8), p2 in axes(GainTotal3,7), h1 in axes(GainTotal3,6), u1 in axes(GainTotal3,5), p1 in axes(GainTotal3,4), h3 in axes(GainTotal3,3)

        nu1 = size(GainTotal3)[5]
        nu2 = size(GainTotal3)[8]

        ViewGainTotal3 = @view(GainTotal3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        ViewGainTotal3Mirror123 = @view(GainTotal3[:,end:-1:1,h3,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        ViewGainTally3 = @view(GainTally3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        ViewGainTally3Mirror123 = @view(GainTally3[:,end:-1:1,h3,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        if m1 == m2 # Both first and second Symmetry true

            ViewGainTotal3Swap12 = @view(GainTotal3[:,1:end,h3,p2,u2,h2,p1,u1,h1])
            ViewGainTally3Swap12 = @view(GainTally3[:,1:end,h3,p2,u2,h2,p1,u1,h1])

            if symmetric_grid

                ViewGainTotal3Swap12Mirror123 = @view(GainTotal3[:,end:-1:1,h3,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])
                ViewGainTally3Swap12Mirror123 = @view(GainTally3[:,end:-1:1,h3,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

                @. tmp_total3 = ViewGainTotal3 + ViewGainTotal3Mirror123 + ViewGainTotal3Swap12 + ViewGainTotal3Swap12Mirror123
                @. tmp_tally3 = ViewGainTally3 + ViewGainTally3Mirror123 + ViewGainTally3Swap12 + ViewGainTally3Swap12Mirror123

                @. ViewGainTotal3 = tmp_total3
                @. ViewGainTotal3Mirror123 = tmp_total3
                @. ViewGainTotal3Swap12 = tmp_total3
                @. ViewGainTotal3Swap12Mirror123 = tmp_total3

                @. ViewGainTally3 = tmp_tally3
                @. ViewGainTally3Mirror123 = tmp_tally3
                @. ViewGainTally3Swap12 = tmp_tally3
                @. ViewGainTally3Swap12Mirror123 = tmp_tally3

            else

                @. tmp_total3 = ViewGainTotal3 + ViewGainTotal3Swap12 
                @. tmp_tally3 = ViewGainTally3 + ViewGainTally3Swap12

                @. ViewGainTotal3 = tmp_total3
                @. ViewGainTotal3Swap12 = tmp_total3

                @. ViewGainTally3 = tmp_tally3
                @. ViewGainTally3Swap12 = tmp_tally3

            end

        elseif symmetric_grid # only first symmetry true

            tmp_total3 = ViewGainTotal3 + ViewGainTotal3Mirror123
            tmp_tally3 = ViewGainTally3 + ViewGainTally3Mirror123

            @. ViewGainTotal3 = tmp_total3
            @. ViewGainTotal3Mirror123 = tmp_total3

            @. ViewGainTally3 = tmp_tally3
            @. ViewGainTally3Mirror123 = tmp_tally3

        end

    end

    if m3 != m4

        tmp_total4 = zeros(Float64,size(GainTotal4)[1:2])
        tmp_tally4 = zeros(Int64,size(GainTally4)[1:2])

        # Particle 4 Gain terms
        for h2 in axes(GainTotal4,9), u2 in axes(GainTotal4,8), p2 in axes(GainTotal4,7), h1 in axes(GainTotal4,6), u1 in axes(GainTotal4,5), p1 in axes(GainTotal4,4), h4 in axes(GainTotal4,3)

            nu1 = size(GainTotal4)[5]
            nu2 = size(GainTotal4)[8]

            ViewGainTotal4 = @view(GainTotal4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
            ViewGainTotal4Mirror123 = @view(GainTotal4[:,end:-1:1,h4,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

            ViewGainTally4 = @view(GainTally4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
            ViewGainTally4Mirror123 = @view(GainTally4[:,end:-1:1,h4,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

            if m1 == m2 # Both first and second Symmetry true

                ViewGainTotal4Swap12 = @view(GainTotal4[:,1:end,h4,p2,u2,h2,p1,u1,h1])
                ViewGainTally4Swap12 = @view(GainTally4[:,1:end,h4,p2,u2,h2,p1,u1,h1])

                if symmetric_grid

                    ViewGainTotal4Swap12Mirror123 = @view(GainTotal4[:,end:-1:1,h4,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])
                    ViewGainTally4Swap12Mirror123 = @view(GainTally4[:,end:-1:1,h4,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

                    @. tmp_total4 = ViewGainTotal4 + ViewGainTotal4Mirror123 + ViewGainTotal4Swap12 + ViewGainTotal4Swap12Mirror123
                    @. tmp_tally4 = ViewGainTally4 + ViewGainTally4Mirror123 + ViewGainTally4Swap12 + ViewGainTally4Swap12Mirror123

                    @. ViewGainTotal4 = tmp_total4
                    @. ViewGainTotal4Mirror123 = tmp_total4
                    @. ViewGainTotal4Swap12 = tmp_total4
                    @. ViewGainTotal4Swap12Mirror123 = tmp_total4

                    @. ViewGainTally4 = tmp_tally4
                    @. ViewGainTally4Mirror123 = tmp_tally4
                    @. ViewGainTally4Swap12 = tmp_tally4
                    @. ViewGainTally4Swap12Mirror123 = tmp_tally4

                else

                    @. tmp_total4 = ViewGainTotal4 + ViewGainTotal4Swap12
                    @. tmp_tally4 = ViewGainTally4 + ViewGainTally4Swap12

                    @. ViewGainTotal4 = tmp_total4
                    @. ViewGainTotal4Swap12 = tmp_total4

                    @. ViewGainTally4 = tmp_tally4
                    @. ViewGainTally4Swap12 = tmp_tally4

                end

            elseif symmetric_grid # only first symmetry true

                @. tmp_total4 = ViewGainTotal4 + ViewGainTotal4Mirror123
                @. tmp_tally4 = ViewGainTally4 + ViewGainTally4Mirror123

                @. ViewGainTotal4 = tmp_total4
                @. ViewGainTotal4Mirror123 = tmp_total4

                @. ViewGainTally4 = tmp_tally4
                @. ViewGainTally4Mirror123 = tmp_tally4

            end

        end

    end

    # Particle Loss Terms
    for h2 in axes(LossTotal,6), u2 in axes(LossTotal,5), p2 in axes(LossTotal,4), h1 in axes(LossTotal,3), u1 in axes(LossTotal,2), p1 in axes(LossTotal,1)

        nu1 = size(LossTotal)[2]
        nu2 = size(LossTotal)[5]        

        if m1==m2 # Both first and second Symmetry true

            if symmetric_grid

                tmp_total = LossTotal[p1,u1,h1,p2,u2,h2] + LossTotal[p2,u2,h2,p1,u1,h1] + LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] + LossTotal[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1]
                tmp_tally = LossTally[p1,u1,h1,p2,u2,h2] + LossTally[p2,u2,h2,p1,u1,h1] + LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] + LossTally[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1]

                LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
                LossTotal[p2,u2,h2,p1,u1,h1] = tmp_total
                LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total
                LossTotal[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1] = tmp_total

                LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally
                LossTally[p2,u2,h2,p1,u1,h1] = tmp_tally
                LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_tally
                LossTally[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1] = tmp_tally

            else

                tmp_total = LossTotal[p1,u1,h1,p2,u2,h2] + LossTotal[p2,u2,h2,p1,u1,h1]
                tmp_tally = LossTally[p1,u1,h1,p2,u2,h2] + LossTally[p2,u2,h2,p1,u1,h1] 

                LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
                LossTotal[p2,u2,h2,p1,u1,h1] = tmp_total

                LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally
                LossTally[p2,u2,h2,p1,u1,h1] = tmp_tally

            end

        elseif symmetric_grid # only first symmetry true

            tmp_total = LossTotal[p1,u1,h1,p2,u2,h2] + LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2]
            tmp_tally = LossTally[p1,u1,h1,p2,u2,h2] + LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2]

            LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
            LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total

            LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally
            LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_tally

        end
    end 

end

"""
    GainLossPolarSymmetryMatrixBinary!(GainMatrix3,GainMatrix4,LossMatrix,GainWeights3,GainWeights4,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)

Applies various physical polar angle symmetries to the Gain and Loss Matrices and Weights/Tally for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossPolarSymmetryMatrixBinary!(GainMatrix3,GainMatrix4,LossMatrix,GainWeights3,GainWeights4,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)

    # The Gain and Loss matrices are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    #   (This FIRST condition only applies if the grids are symmetric in u, e.g. uniform or binary grids)
    # SECOND: if the incident masses are equal (m1==m2) then Gain and Loss are symmetric to swapping the incident particles
    
    tmp_total3 = zeros(Float64,size(GainMatrix3)[1:2])

    # Particle 3 Gain terms
    for h2 in axes(GainMatrix3,9), u2 in axes(GainMatrix3,8), p2 in axes(GainMatrix3,7), h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4), h3 in axes(GainMatrix3,3)

        nu1 = size(GainMatrix3)[5]
        nu2 = size(GainMatrix3)[8]

        ViewGainMatrix3 = @view(GainMatrix3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        ViewGainMatrix3Mirror123 = @view(GainMatrix3[:,end:-1:1,h3,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        if m1 == m2 # Both first and second Symmetry true

            ViewGainMatrix3Swap12 = @view(GainMatrix3[:,1:end,h3,p2,u2,h2,p1,u1,h1])

            if symmetric_grid

                ViewGainMatrix3Swap12Mirror123 = @view(GainMatrix3[:,end:-1:1,h3,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

                @. tmp_total3 = ViewGainMatrix3 + ViewGainMatrix3Mirror123 + ViewGainMatrix3Swap12 + ViewGainMatrix3Swap12Mirror123
                @. tmp_total3 /= 4.0 # average over the 4 symmetric terms

                @. ViewGainMatrix3 = tmp_total3
                @. ViewGainMatrix3Mirror123 = tmp_total3
                @. ViewGainMatrix3Swap12 = tmp_total3
                @. ViewGainMatrix3Swap12Mirror123 = tmp_total3

            else

                @. tmp_total3 = ViewGainMatrix3 + ViewGainMatrix3Swap12 
                @. tmp_total3 /= 2.0 # average over the 2 symmetric terms

                @. ViewGainMatrix3 = tmp_total3
                @. ViewGainMatrix3Swap12 = tmp_total3

            end

        elseif symmetric_grid # only first symmetry true

            @. tmp_total3 = ViewGainMatrix3 + ViewGainMatrix3Mirror123
            @. tmp_total3 /= 2.0 # average over the 2 symmetric terms

            @. ViewGainMatrix3 = tmp_total3
            @. ViewGainMatrix3Mirror123 = tmp_total3

        end

    end

    if m3 != m4

        tmp_total4 = zeros(Float64,size(GainMatrix4)[1:2])

        # Particle 4 Gain terms
        for h2 in axes(GainMatrix4,9), u2 in axes(GainMatrix4,8), p2 in axes(GainMatrix4,7), h1 in axes(GainMatrix4,6), u1 in axes(GainMatrix4,5), p1 in axes(GainMatrix4,4), h4 in axes(GainMatrix4,3)

            nu1 = size(GainMatrix4)[5]
            nu2 = size(GainMatrix4)[8]

            ViewGainMatrix4 = @view(GainMatrix4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
            ViewGainMatrix4Mirror123 = @view(GainMatrix4[:,end:-1:1,h4,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

            if m1 == m2 # Both first and second Symmetry true

                ViewGainMatrix4Swap12 = @view(GainMatrix4[:,1:end,h4,p2,u2,h2,p1,u1,h1])

                if symmetric_grid

                    ViewGainMatrix4Swap12Mirror123 = @view(GainMatrix4[:,end:-1:1,h4,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

                    @. tmp_total4 = ViewGainMatrix4 + ViewGainMatrix4Mirror123 + ViewGainMatrix4Swap12 + ViewGainMatrix4Swap12Mirror123
                    @. tmp_total4 /= 4.0 # average over the 4 symmetric terms

                    @. ViewGainMatrix4 = tmp_total4
                    @. ViewGainMatrix4Mirror123 = tmp_total4
                    @. ViewGainMatrix4Swap12 = tmp_total4
                    @. ViewGainMatrix4Swap12Mirror123 = tmp_total4

                else

                    @. tmp_total4 = ViewGainMatrix4 + ViewGainMatrix4Swap12
                    @. tmp_total4 /= 2.0 # average over the 2 symmetric terms

                    @. ViewGainMatrix4 = tmp_total4
                    @. ViewGainMatrix4Swap12 = tmp_total4

                end

            elseif symmetric_grid # only first symmetry true

                @. tmp_total4 = ViewGainMatrix4 + ViewGainMatrix4Mirror123
                @. tmp_total4 /= 2.0 # average over the 2 symmetric terms

                @. ViewGainMatrix4 = tmp_total4
                @. ViewGainMatrix4Mirror123 = tmp_total4

            end

        end

    end

    # Particle Loss Terms
    for h2 in axes(LossMatrix,6), u2 in axes(LossMatrix,5), p2 in axes(LossMatrix,4), h1 in axes(LossMatrix,3), u1 in axes(LossMatrix,2), p1 in axes(LossMatrix,1)

        nu1 = size(LossMatrix)[2]
        nu2 = size(LossMatrix)[5]        

        if m1==m2 # Both first and second Symmetry true

            if symmetric_grid

                tmp_total = LossMatrix[p1,u1,h1,p2,u2,h2] + LossMatrix[p2,u2,h2,p1,u1,h1] + LossMatrix[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] + LossMatrix[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1]
                tmp_total /= 4.0 # average over the 4 symmetric terms

                LossMatrix[p1,u1,h1,p2,u2,h2] = tmp_total
                LossMatrix[p2,u2,h2,p1,u1,h1] = tmp_total
                LossMatrix[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total
                LossMatrix[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1] = tmp_total

            else

                tmp_total = LossMatrix[p1,u1,h1,p2,u2,h2] + LossMatrix[p2,u2,h2,p1,u1,h1]
                tmp_total /= 2.0 # average over the 2 symmetric terms 

                LossMatrix[p1,u1,h1,p2,u2,h2] = tmp_total
                LossMatrix[p2,u2,h2,p1,u1,h1] = tmp_total

            end

        elseif symmetric_grid # only first symmetry true

            tmp_total = LossMatrix[p1,u1,h1,p2,u2,h2] + LossMatrix[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2]
            tmp_total /= 2.0 # average over the 2 symmetric terms

            LossMatrix[p1,u1,h1,p2,u2,h2] = tmp_total
            LossMatrix[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total

        end
    end 

    return nothing

end

"""
    GainLossAzimuthalSymmetryBinary(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,Indistinguishable_34,m1,m2,m3,m4)

Applies various physical azimuthal angle symmetries to the Gain and Loss terms for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossAzimuthalSymmetryBinary!(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)

    # The Gain and Loss matrices are symmetric in with respect to rotations of the azimuthal angle. 
    # If the azimuthal grid is uniform then we can apply this symmetry to all variations of the azimuthal bins.

    num_h1 = size(GainTotal3,6)
    num_h2 = size(GainTotal3,9)
    num_h3 = size(GainTotal3,3)
    num_h4 = size(GainTotal4,3)
    num_sections12 = lcm(num_h1,num_h2) # number of sections to divide the azimuthal bins into for averaging
    num_sections123 = lcm(num_h1,num_h2,num_h3) # number of sections to divide the azimuthal bins into for averaging
    num_sections124 = lcm(num_h1,num_h2,num_h4) # number of sections to divide the azimuthal bins into for averaging

    # Particle 3 Gain terms
    Threads.@threads for idx in CartesianIndices((axes(GainTotal3,8),axes(GainTotal3,7),axes(GainTotal3,5),axes(GainTotal3,4),axes(GainTotal3,2),axes(GainTotal3,1)))

        u2, p2, u1, p1, u3, p3 = Tuple(idx)
    #Threads.@threads for u2 in axes(GainTotal3,8), p2 in axes(GainTotal3,7), u1 in axes(GainTotal3,5), p1 in axes(GainTotal3,4), u3 in axes(GainTotal3,2), p3 in axes(GainTotal3,1)

        for off2 in 0:num_sections123-1, off3 in 0:num_sections123-1 # loop over the maximum number of azimuthal bins for the particles

            tmp_total = zero(Float64)
            tmp_tally = zero(Int64)

            for h in 1:num_sections123, 

                h1 = mod(floor(Int64, h / (num_sections123 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections123 / num_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections123 / num_h3)),num_h3) + 1

                if Threads.threadid() == 1
                println("Particle 3 - p3: $p3, u3: $u3, p1: $p1, u1: $u1, p2: $p2, u2: $u2, h: $h, off2: $off2, off3: $off3, h1: $h1, h2: $h2, h3: $h3")
                end

                tmp_total += GainTotal3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
                tmp_tally += GainTally3[p3,u3,h3,p1,u1,h1,p2,u2,h2]

            end

            for h in 1:num_sections123

                h1 = mod(floor(Int64, h / (num_sections123 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections123 / num_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections123 / num_h3)),num_h3) + 1

                GainTotal3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = tmp_total
                GainTally3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = tmp_tally

            end

        end

    end

    if m3 != m4
        # Particle 4 Gain terms
        Threads.@threads for idx in CartesianIndices((axes(GainTotal4,8),axes(GainTotal4,7),axes(GainTotal4,5),axes(GainTotal4,4),axes(GainTotal4,2),axes(GainTotal4,1)))

            u2, p2, u1, p1, u4, p4 = Tuple(idx)
        #Threads.@threads for u2 in axes(GainTotal4,8), p2 in axes(GainTotal4,7), u1 in axes(GainTotal4,5), p1 in axes(GainTotal4,4), u4 in axes(GainTotal4,2), p4 in axes(GainTotal4,1)

            for off2 in 0:num_sections124-1, off4 in 0:num_sections124-1 # loop over the maximum number of azimuthal bins for the particles

                tmp_total = zero(Float64)
                tmp_tally = zero(Int64)

                for h in 1:num_sections124

                    h1 = mod(floor(Int64, h / (num_sections124 / num_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sections124 / num_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sections124 / num_h4)),num_h4) + 1

                    tmp_total += GainTotal4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
                    tmp_tally += GainTally4[p4,u4,h4,p1,u1,h1,p2,u2,h2]

                end

                for h in 1:num_sections124 # loop over the maximum number of azimuthal bins for the particles

                    h1 = mod(floor(Int64, h / (num_sections124 / num_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sections124 / num_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sections124 / num_h4)),num_h4) + 1

                    GainTotal4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = tmp_total
                    GainTally4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = tmp_tally

                end

            end

        end

    end

    # Particle Loss Terms
    Threads.@threads for idx in CartesianIndices((axes(LossTotal,6),axes(LossTotal,5),axes(LossTotal,4),axes(LossTotal,3),axes(LossTotal,2),axes(LossTotal,1)))

        u2, h2, p2, h1, u1, p1 = Tuple(idx)
    #Threads.@threads for u2 in axes(LossTotal,5), p2 in axes(LossTotal,4), u1 in axes(LossTotal,2), p1 in axes(LossTotal,1)

        for off2 in 0:num_sections12-1 # loop over the maximum number of azimuthal bins for the particles

            tmp_total = zero(Float64)
            tmp_tally = zero(Int64)

            for h in 1:num_sections12

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                tmp_total += LossTotal[p1,u1,h1,p2,u2,h2]
                tmp_tally += LossTally[p1,u1,h1,p2,u2,h2]
                
            end

            for h in 1:num_sections12 # loop over the maximum number of azimuthal bins for the particles

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
                LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally

            end

        end

    end 

    return nothing

end

"""
    GainLossAzimuthalSymmetryMatrixBinary(GainMatrix3,GainMatrix4,LossMatrix,GainWeights3,GainWeights4,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)

Applies various physical azimuthal angle symmetries to the Gain and Loss Matrices and Weights/Tally for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossAzimuthalSymmetryMatrixBinary!(GainMatrix3,GainMatrix4,LossMatrix,GainWeights3,GainWeights4,LossTally,m1,m2,m3,m4,symmetric_grid::Bool)

    # The Gain and Loss matrices are symmetric in with respect to rotations of the azimuthal angle. 
    # If the azimuthal grid is uniform then we can apply this symmetry to all variations of the azimuthal bins.

    num_h1 = size(GainMatrix3,6)
    num_h2 = size(GainMatrix3,9)
    num_h3 = size(GainMatrix3,3)
    num_h4 = size(GainMatrix4,3)
    num_sections12 = lcm(num_h1,num_h2) # number of sections to divide the azimuthal bins into for averaging
    num_sections123 = lcm(num_h1,num_h2,num_h3) # number of sections to divide the azimuthal bins into for averaging
    num_sections124 = lcm(num_h1,num_h2,num_h4) # number of sections to divide the azimuthal bins into for averaging

    # Particle 3 Gain terms
    Threads.@threads for idx in CartesianIndices((axes(GainMatrix3,8),axes(GainMatrix3,7),axes(GainMatrix3,5),axes(GainMatrix3,4),axes(GainMatrix3,2),axes(GainMatrix3,1)))

        u2, p2, u1, p1, u3, p3 = Tuple(idx)

        for off2 in 0:num_sections123-1, off3 in 0:num_sections123-1 # loop over the maximum number of azimuthal bins for the particles
            
            tmp_matrix = zero(Float64)

            for h in 1:num_sections123

                h1 = mod(floor(Int64, h / (num_sections123 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections123 / num_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections123 / num_h3)),num_h3) + 1

                if Threads.threadid() == 1
                println("Particle 3 - p3: $p3, u3: $u3, p1: $p1, u1: $u1, p2: $p2, u2: $u2, h: $h, off2: $off2, off3: $off3, h1: $h1, h2: $h2, h3: $h3")
                end

                tmp_matrix += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]

            end

            for h in 1:num_sections123

                h1 = mod(floor(Int64, h / (num_sections123 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections123 / num_h2)),num_h2) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections123 / num_h3)),num_h3) + 1

                GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = tmp_matrix / num_sections123 # average over number of sections/rotations

            end

        end

    end

    if m3 != m4
        # Particle 4 Gain terms
        Threads.@threads for idx in CartesianIndices((axes(GainMatrix4,8),axes(GainMatrix4,7),axes(GainMatrix4,5),axes(GainMatrix4,4),axes(GainMatrix4,2),axes(GainMatrix4,1)))

            u2, p2, u1, p1, u4, p4 = Tuple(idx)
        #Threads.@threads for u2 in axes(GainMatrix4,8), p2 in axes(GainMatrix4,7), u1 in axes(GainMatrix4,5), p1 in axes(GainMatrix4,4), u4 in axes(GainMatrix4,2), p4 in axes(GainMatrix4,1)

            for off2 in 0:num_sections124-1, off4 in 0:num_sections124-1 # loop over the maximum number of azimuthal bins for the particles

                tmp_matrix = zero(Float64)

                for h in 1:num_sections124

                    h1 = mod(floor(Int64, h / (num_sections124 / num_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sections124 / num_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sections124 / num_h4)),num_h4) + 1

                    tmp_matrix += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]

                end

                for h in 1:num_sections124 # loop over the maximum number of azimuthal bins for the particles

                    h1 = mod(floor(Int64, h / (num_sections124 / num_h1)),num_h1) + 1
                    h2 = mod(floor(Int64, (h+off2) / (num_sections124 / num_h2)),num_h2) + 1
                    h4 = mod(floor(Int64, (h+off4) / (num_sections124 / num_h4)),num_h4) + 1

                    GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = tmp_matrix / num_sections124 # average over number of sections/rotations

                end

            end

        end

    end

    # Particle Loss Terms
    Threads.@threads for idx in CartesianIndices((axes(LossMatrix1,5),axes(LossMatrix1,4),axes(LossMatrix1,2),axes(LossMatrix1,1)))

        u2, p2, u1, p1 = Tuple(idx)

        for off2 in 0:num_sections12-1 # loop over the maximum number of azimuthal bins for the particles

            tmp_matrix = zero(Float64)

            for h in 1:num_sections12

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                tmp_matrix += LossMatrix[p1,u1,h1,p2,u2,h2] # LossMatrix2 is just a permutation of LossMatrix1 so we only need to sum over one of them
                
            end

            for h in 1:num_sections12 # loop over the maximum number of azimuthal bins for the particles

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                LossMatrix[p1,u1,h1,p2,u2,h2] = tmp_matrix / num_sections12 # average over number of sections/rotations

            end

        end

    end 

    return nothing

end

=#
