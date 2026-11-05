"""
    WeightedAverageGainBinary!(GainMatrix3,OldGainMatrix3,GainTally3_K,OldGainWeights3,GainMatrix4,OldGainMatrix4,GainTally4_K,OldGainWeights4_K)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
```math
    I = (I1w1 + I2w2)/(w1 + w2)
```
where `I1` and `I1` are the old and new gain matrix element estimates and `w1` and `w2` are the corresponding weights. Here the weights are taken to be `w=k`
"""
function WeightedAverageGainBinary!(GainMatrix3::Array{Float64,9},OldGainMatrix3::Array{Float64,9},GainTally3_K::AbstractArray{UInt32,9},GainTally3_N::AbstractArray{UInt32,8},OldGainWeights3::Array{Float64,9},GainMatrix4::Array{Float64,9},OldGainMatrix4::Array{Float64,9},GainTally4_K::AbstractArray{UInt32,9},GainTally4_N::AbstractArray{UInt32,8},OldGainWeights4::Array{Float64,9})

    # new weights k^2/N
    NewGainWeights = similar(OldGainWeights3)
    for i in axes(GainTally3_K,1)
        @view(NewGainWeights[i,:,:,:,:,:,:,:,:]) .= @view(GainTally3_K[i,:,:,:,:,:,:,:,:]) #./ GainTally3_N
    end
    replace!(NewGainWeights,NaN=>0e0)
    #weighted average
    @. OldGainMatrix3 = (GainMatrix3*NewGainWeights+OldGainMatrix3*OldGainWeights3)/(NewGainWeights+OldGainWeights3)
    replace!(OldGainMatrix3,NaN=>0e0)
    # adjust weights for next integration step
    @. OldGainWeights3 += NewGainWeights

    NewGainWeights = 0.0
    GC.gc()
    
    # repeat above for 4 to save memory 
    NewGainWeights = similar(OldGainWeights4)
    for i in axes(GainTally4_K,1)
        @view(NewGainWeights[i,:,:,:,:,:,:,:,:]) .= @view(GainTally4_K[i,:,:,:,:,:,:,:,:]) #./ GainTally4_N
    end
    replace!(NewGainWeights,NaN=>0e0)
    @. OldGainMatrix4 = (GainMatrix4*NewGainWeights+OldGainMatrix4*OldGainWeights4)/(NewGainWeights+OldGainWeights4)
    replace!(OldGainMatrix4,NaN=>0e0)
    @. OldGainWeights4 += NewGainWeights

end

function WeightedAverageGainBinary!(GainMatrix3::Array{Float64,9},OldGainMatrix3::Array{Float64,9},GainTally3_K::AbstractArray{UInt32,9},GainTally3_N::AbstractArray{UInt32,8},OldGainWeights3::Array{Float64,9})
    # Version for if mu3 == mu4

    # new weights k^2/N
    NewGainWeights = similar(OldGainWeights3)
    for i in axes(GainTally3_K,1)
        @view(NewGainWeights[i,:,:,:,:,:,:,:,:]) .= @view(GainTally3_K[i,:,:,:,:,:,:,:,:]) #./ GainTally3_N
    end
    replace!(NewGainWeights,NaN=>0e0)
    #weighted average
    @. OldGainMatrix3 = (GainMatrix3*NewGainWeights+OldGainMatrix3*OldGainWeights3)/(NewGainWeights+OldGainWeights3)
    replace!(OldGainMatrix3,NaN=>0e0)
    # adjust weights for next integration step
    @. OldGainWeights3 += NewGainWeights

end

"""
    WeightedAverageLossBinary!(LossMatrix,OldLossMatrix,LossTally,OldLossTally)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageLossBinary!(LossMatrix::Array{Float64,6},OldLossMatrix::Array{Float64,6},LossTally::Array{UInt32,6},OldLossTally::Array{UInt32,6})

    # weighted average 
    @. OldLossMatrix = (LossMatrix*LossTally+OldLossMatrix*OldLossTally)/(LossTally+OldLossTally)

    replace!(OldLossMatrix,NaN=>0e0)

    # adding tallies
    @. OldLossTally += LossTally

end



"""
    WeightedAverageGainBinaryChunk!(ChunkGainMatrix3,OldChunkGainMatrix3,ChunkGainTally3_K,OldChunkGainWeights3,ChunkGainMatrix4,OldChunkGainMatrix4,ChunkGainTally4_K,OldChunkGainWeights4_K)

Computes the integral estimate by weighted average of the old and new chunk gain matrices. Mutating the old gain and tally terms.
```math
    I = (I1w1 + I2w2)/(w1 + w2)
```
where `I1` and `I1` are the old and new gain matrix element estimates and `w1` and `w2` are the corresponding weights. Here the weights are taken to be `w=k`
"""
function WeightedAverageGainBinaryChunk!(ChunkGainMatrix3::AbstractArray{Float32,7},OldChunkGainMatrix3::AbstractArray{Float32,7},ChunkGainTally3_K::AbstractArray{UInt32,7},GainTally3_N::AbstractArray{UInt32,6},OldChunkGainWeights3::AbstractArray{Float32,7},ChunkGainMatrix4::AbstractArray{Float32,7},OldChunkGainMatrix4::AbstractArray{Float32,7},ChunkGainTally4_K::AbstractArray{UInt32,7},GainTally4_N::AbstractArray{UInt32,6},OldChunkGainWeights4::AbstractArray{Float32,7})

    # new weights k
    @inbounds @simd for I in eachindex(OldChunkGainMatrix3)
        neww = Float64(ChunkGainTally3_K[I]) 
        oldw = OldChunkGainWeights3[I]
        denom = neww + oldw
        oldm = OldChunkGainMatrix3[I]
        newm = ChunkGainMatrix3[I]
        OldChunkGainMatrix3[I] = denom == 0f0 ? 0f0 : Float32(muladd(newm, neww, oldm * oldw) / denom)
        OldChunkGainWeights3[I] = Float32(denom)
    end

    @inbounds @simd for I in eachindex(OldChunkGainMatrix4)
        neww = Float64(ChunkGainTally4_K[I]) 
        oldw = OldChunkGainWeights4[I]
        denom = neww + oldw
        oldm = OldChunkGainMatrix4[I]
        newm = ChunkGainMatrix4[I]
        OldChunkGainMatrix4[I] = denom == 0f0 ? 0f0 : Float32(muladd(newm, neww, oldm * oldw) / denom)
        OldChunkGainWeights4[I] = Float32(denom)
    end

end

function WeightedAverageGainBinaryChunk!(ChunkGainMatrix3::AbstractArray{Float32,7},OldChunkGainMatrix3::AbstractArray{Float32,7},ChunkGainTally3_K::AbstractArray{UInt32,7},ChunkGainTally3_N::AbstractArray{UInt32,6},OldChunkGainWeights3::AbstractArray{Float32,7})

    # Version for if mu3 == mu4

    # new weights k
    @inbounds @simd for I in eachindex(OldChunkGainMatrix3)
        neww = Float64(ChunkGainTally3_K[I]) 
        oldw = OldChunkGainWeights3[I]
        denom = neww + oldw
        oldm = OldChunkGainMatrix3[I]
        newm = ChunkGainMatrix3[I]
        OldChunkGainMatrix3[I] = denom == 0f0 ? 0f0 : Float32(muladd(newm, neww, oldm * oldw) / denom)
        OldChunkGainWeights3[I] = Float32(denom)
    end

end

"""
    WeightedAverageLossBinaryChunk!(ChunkLossMatrix,OldChunkLossMatrix,ChunkLossTally,OldChunkLossTally)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageLossBinaryChunk!(ChunkLossMatrix::AbstractArray{Float32,4},OldChunkLossMatrix::AbstractArray{Float32,4},ChunkLossTally::AbstractArray{UInt32,4},OldChunkLossTally::AbstractArray{UInt32,4})

    # weighted average 
    @. OldChunkLossMatrix = (ChunkLossMatrix*ChunkLossTally+OldChunkLossMatrix*OldChunkLossTally)/(ChunkLossTally+OldChunkLossTally)

    # TODO: replace with a more memory efficient method
    replace!(OldChunkLossMatrix,NaN=>0f0)

    # adding tallies
    @. OldChunkLossTally += ChunkLossTally

end