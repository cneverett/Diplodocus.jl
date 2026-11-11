"""
    WeightedAverageGainEmission!(GainMatrix2,OldGainMatrix2,GainTallyK2,OldGainTallyK2,GainMatrix3,OldGainMatrix3,GainTallyK3,OldGainTallyK3)

Computes the integral estimate by weighted average of the old and new gain matrices, weighted by the number of non-zero samples obtained. Mutating the old gain and tally terms.
"""
function WeightedAverageGainEmission!(GainMatrix2::AbstractArray{Float64,6},OldGainMatrix2::AbstractArray{Float64,6},GainTallyK2::AbstractArray{UInt32,6},OldGainTallyK2::AbstractArray{UInt32,6},GainMatrix3::AbstractArray{Float64,6},OldGainMatrix3::AbstractArray{Float64,6},GainTallyK3::AbstractArray{UInt32,6},OldGainTallyK3::AbstractArray{UInt32,6})

    # weighted average
    @. OldGainMatrix2 = (GainMatrix2*GainTallyK2+OldGainMatrix2*OldGainTallyK2)/(GainTallyK2+OldGainTallyK2) 
    @. OldGainMatrix3 = (GainMatrix3*GainTallyK3+OldGainMatrix3*OldGainTallyK3)/(GainTallyK3+OldGainTallyK3)

    replace!(OldGainMatrix2,NaN=>0e0)
    replace!(OldGainMatrix3,NaN=>0e0)

    # adding tallies 
    @. OldGainTallyK2 += GainTallyK2
    @. OldGainTallyK3 += GainTallyK3

end

"""
    WeightedAverageLossEmission!(LossMatrix1,OldLossMatrix1,LossTallyK1,OldLossTallyK1)

Computes the integral estimate by weighted average of the old and new gain matrices, weighted by the number of non-zero samples obtained. Mutating the old gain and tally terms.
"""
function WeightedAverageLossEmission!(LossMatrix1::AbstractArray{Float64,3},OldLossMatrix1::AbstractArray{Float64,3},LossTallyK1::AbstractArray{UInt32,3},OldLossTallyK1::AbstractArray{UInt32,3})

    # weighted average 
    @. OldLossMatrix1 = (LossMatrix1*LossTallyK1+OldLossMatrix1*OldLossTallyK1)/(LossTallyK1+OldLossTallyK1)

    replace!(OldLossMatrix1,NaN=>0e0)

    # adding tallies
    @. OldLossTallyK1 += LossTallyK1

end