
"""
    MomentumSpaceFactorsBinary!(GainMatrix3,GainMatrix4,u3_bounds,u4_bounds,h3_bounds,h4_bounds,Indistinguishable_12)

Applies momentum space volume element to the Gain and Loss matrices such that they have the correct dimensions for use in particle transport, i.e. they have units of length^3/time (normalised by σT*c)
"""
function MomentumSpaceFactorsBinary!(GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Momentum space volume elements
    for h2 in axes(GainMatrix3,9), u2 in axes(GainMatrix3,8), p2 in axes(GainMatrix3,7), h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4) # common axes
        for h3 in axes(GainMatrix3,3), u3 in axes(GainMatrix3,2), p3 in 1:size(GainMatrix3,1)
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3]) # du3dh3
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
        for h4 in axes(GainMatrix4,3), u4 in axes(GainMatrix4,2), p4 in 1:size(GainMatrix4,1)
            GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4]) # du4dh4
            GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
    end

    return nothing

end

function MomentumSpaceFactorsBinary!(GainMatrix3::Array{Float64,9},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Momentum space volume elements
    for h2 in axes(GainMatrix3,9), u2 in axes(GainMatrix3,8), p2 in axes(GainMatrix3,7), h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4) # common axes
        for h3 in axes(GainMatrix3,3), u3 in axes(GainMatrix3,2), p3 in 1:size(GainMatrix3,1)
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3]) # du3dh3
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
    end

    return nothing

end


# ==================== Chunked versions ==================== #


"""
    MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3,ChunkGainMatrix4,u3_bounds,u4_bounds,h3_bounds,h4_bounds,Indistinguishable_12)

Applies momentum space volume element to the Chunked Gain and Loss matrices such that they have the correct dimensions for use in particle transport, i.e. they have units of length^3/time (normalised by σT*c)
"""
function MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3::Array{Float32,7},ChunkGainMatrix4::Array{Float32,7},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Momentum space volume elements
    for h2 in axes(ChunkGainMatrix3,7), u2 in axes(ChunkGainMatrix3,6), h1 in axes(ChunkGainMatrix3,5), u1 in axes(ChunkGainMatrix3,4) # common axes
        for h3 in axes(ChunkGainMatrix3,3), u3 in axes(ChunkGainMatrix3,2), p3 in 1:size(ChunkGainMatrix3,1)
            ChunkGainMatrix3[p3,u3,h3,u1,h1,u2,h2] *= Float32((u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])) # du3dh3
            ChunkGainMatrix3[p3,u3,h3,u1,h1,u2,h2] /= Float32(1+Indistinguishable_12)
        end
        for h4 in axes(ChunkGainMatrix4,3), u4 in axes(ChunkGainMatrix4,2), p4 in 1:size(ChunkGainMatrix4,1)
            ChunkGainMatrix4[p4,u4,h4,u1,h1,u2,h2] *= Float32((u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])) # du4dh4
            ChunkGainMatrix4[p4,u4,h4,u1,h1,u2,h2] /= Float32(1+Indistinguishable_12)
        end
    end

    return nothing

end

function MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3::Array{Float32,7},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Momentum space volume elements
    for h2 in axes(ChunkGainMatrix3,7), u2 in axes(ChunkGainMatrix3,6), h1 in axes(ChunkGainMatrix3,5), u1 in axes(ChunkGainMatrix3,4) # common axes
        for h3 in axes(ChunkGainMatrix3,3), u3 in axes(ChunkGainMatrix3,2), p3 in 1:size(ChunkGainMatrix3,1)
            ChunkGainMatrix3[p3,u3,h3,u1,h1,u2,h2] *= Float32((u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])) # du3dh3
            ChunkGainMatrix3[p3,u3,h3,u1,h1,u2,h2] /= Float32(1+Indistinguishable_12)
        end
    end

    return nothing

end
