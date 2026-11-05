
function MomentumSpaceFactorsEmission!(LossMatrix1,GainMatrix2,GainMatrix3::Array{Float64,6},p3val::Vector{Float64},u3val::Vector{Float64},h3val::Vector{Float64})

    for h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4), h3 in axes(GainMatrix3,3), u3 in axes(GainMatrix3,2), p3 in axes(GainMatrix3,1)
        GainMatrix3[p3,u3,h3,p1,u1,h1] *= (u3val[u3+1]-u3val[u3])*(p3val[p3+1]-p3val[p3])*(h3val[h3+1]-h3val[h3]) #dp3du3dh3
    end

    # TODO: ADD GAINMATRIX2 and LOSSMATRIX1

    return nothing

end

function GainLossPolarSymmetryEmission!(GainMatrix2::Array{Float64,6},GainMatrix3::Array{Float64,6},LossMatrix1::Array{Float64,3},GainTally2::Array{UInt32,6},GainTally3::Array{UInt32,6},LossTally1::Array{UInt32,3})

    GainMatrix2Mirror = @view(GainMatrix2[:,end:-1:1,:,:,end:-1:1,:])
    GainMatrix3Mirror = @view(GainMatrix3[:,end:-1:1,:,:,end:-1:1,:])
    LossMatrix1Mirror = @view(LossMatrix1[:,end:-1:1,:,:,end:-1:1,:])

    GainTally2Mirror = @view(GainTally2[:,end:-1:1,:,:,end:-1:1,:])
    GainTally3Mirror = @view(GainTally3[:,end:-1:1,:,:,end:-1:1,:])
    LossTally1Mirror = @view(LossTally1[:,end:-1:1,:,:,end:-1:1,:])

    @. GainMatrix2 = (GainMatrix2Mirror + GainMatrix2) / 2
    @. GainMatrix3 = (GainMatrix3Mirror + GainMatrix3) / 2
    @. LossMatrix1 = (LossMatrix1Mirror + LossMatrix1) / 2

    @. GainTally2 = round(UInt32, (GainTally2Mirror + GainTally2) / 2)
    @. GainTally3 = round(UInt32, (GainTally3Mirror + GainTally3) / 2)
    @. LossTally1 = round(UInt32, (LossTally1Mirror + LossTally1) / 2)

    return nothing

end # function

function GainLossAzimuthalSymmetryEmission!(GainMatrix2::Array{Float64,6},GainMatrix3::Array{Float64,6},LossMatrix1::Array{Float64,3},GainTally2::Array{UInt32,6},GainTally3::Array{UInt32,6},LossTally1::Array{UInt32,3})

    num_h1 = size(GainMatrix3,6)
    num_h2 = size(GainMatrix2,3)
    num_h3 = size(GainMatrix3,3)
    num_sections13 = lcm(num_h1,num_h3) # number of sections to divide the azimuthal bins into for averaging
    num_sections12 = lcm(num_h1,num_h2) # number of sections to divide the azimuthal bins into for averaging
    num_sections1 = num_h1

    # Particle 3 Gain terms
    Threads.@threads for idx in CartesianIndices((axes(GainMatrix3,5),axes(GainMatrix3,4),axes(GainMatrix3,2),axes(GainMatrix3,1),axes(GainMatrix3,2)))

        u1, p1, u3, p3 = Tuple(idx)

        for off3 in 0:num_sections13-1 # loop over the maximum number of azimuthal bins for the particles
            
            tmp_total = zero(Float64)
            tmp_tally = zero(UInt32)

            for h in 1:num_sections13

                h1 = mod(floor(Int64, h / (num_sections13 / num_h1)),num_h1) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections13 / num_h3)),num_h3) + 1

                tmp_total += GainMatrix3[p3,u3,h3,p1,u1,h1]
                tmp_tally += GainTally3[p3,u3,h3,p1,u1,h1]
                
            end

            for h in 1:num_sections13

                h1 = mod(floor(Int64, h / (num_sections13 / num_h1)),num_h1) + 1
                h3 = mod(floor(Int64, (h+off3) / (num_sections13 / num_h3)),num_h3) + 1

                GainMatrix3[p3,u3,h3,p1,u1,h1] = tmp_total / num_sections13 # average over number of sections/rotations
                GainTally3[p3,u3,h3,p1,u1,h1] = round(UInt32, tmp_tally / num_sections13)
            end

        end

    end

    # Particle 2 Gain terms
    Threads.@threads for idx in CartesianIndices((axes(GainMatrix2,5),axes(GainMatrix2,4),axes(GainMatrix2,2),axes(GainMatrix2,1),axes(GainMatrix2,2)))

        u1, p1, u2, p2 = Tuple(idx)

        for off2 in 0:num_sections12-1 # loop over the maximum number of azimuthal bins for the particles
            
            tmp_total = zero(Float64)
            tmp_tally = zero(UInt32)

            for h in 1:num_sections12

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                tmp_total += GainMatrix2[p2,u2,h2,p1,u1,h1]
                tmp_tally += GainTally2[p2,u2,h2,p1,u1,h1]

            end

            for h in 1:num_sections12

                h1 = mod(floor(Int64, h / (num_sections12 / num_h1)),num_h1) + 1
                h2 = mod(floor(Int64, (h+off2) / (num_sections12 / num_h2)),num_h2) + 1

                GainMatrix2[p2,u2,h2,p1,u1,h1] = tmp_total / num_sections12 # average over number of sections/rotations
                GainTally2[p2,u2,h2,p1,u1,h1] = round(UInt32, tmp_tally / num_sections12)

            end

        end

    end

    # Particle 1 Loss Terms
    Threads.@threads for idx in CartesianIndices((axes(LossMatrix1,2),axes(LossMatrix1,1)))

        u1, p1 = Tuple(idx)

        for off1 in 0:num_sections1-1 # loop over the maximum number of azimuthal bins for the particles

            tmp_total = zero(Float64)
            tmp_tally = zero(UInt32)

            for h in 1:num_sections1

                h1 = mod(floor(Int64, h / (num_sections1 / num_h1)),num_h1) + 1

                tmp_total += LossMatrix1[p1,u1,h1] # LossMatrix2 is just a permutation of LossMatrix1 so we only need to sum over one of them
                tmp_tally += LossTally1[p1,u1,h1] # LossMatrix2 is just a permutation of LossMatrix1 so we only need to sum over one of them
                
            end

            for h in 1:num_sections1 # loop over the maximum number of azimuthal bins for the particles

                h1 = mod(floor(Int64, h / (num_sections1 / num_h1)),num_h1) + 1
   
                LossMatrix1[p1,u1,h1] = tmp_total / num_sections1 # average over number of sections/rotations
                LossTally1[p1,u1,h1] = round(UInt32, tmp_tally / num_sections1)

            end

        end

    end 

    return nothing

end




## ============= TO Be Removed ================== ##

#=
"""
    PhaseSpaceFactorsSync1!(SMatrix,p1val,t1val,p2val,t2val)

Applies phase space volume element factors for 'SMatrix' terms in order to correctly apply 'SyncSymmetry' corrections. 
"""
function PhaseSpaceFactorsSync1!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64},p2val::Vector{Float64},t2val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] *= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
        SMatrix[ll,kk,jj,ii] *= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) #dp2dmu2
    end

end # function

"""
    PhaseSpaceFactorsSync2!(SMatrix,p1val,t1val)

To follow 'PhaseSpaceFactorsSync1' and 'SyncSymmetry'. Correct phase space factors on 'SMatrix' for use in kinetic codes. 
"""
function PhaseSpaceFactorsSync2!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) #dp1dmu1 
    end

end # function 


"""
    SyncSymmetry!(SMatrix)

To follow 'PhaseSpaceFactorsSync1'. Synchrotron emission has a symmetry with respect to cos(theta) -> -cos(theta) for both initial particle and photon momenta.
"""
function SymmetryEmission!(TMatrix1::Array{Float64,3},SMatrix2::Array{Float64,6},SMatrix3::Array{Float64,6})

    avgS = zeros(Float64,size(SMatrix2))
    @. avgS = (SMatrix2[:,end:-1:1,:,:,end:-1:1,:] + SMatrix2)/2
    SMatrix2 .= avgS

    avgS = zeros(Float64,size(SMatrix3))
    @. avgS = (SMatrix3[:,end:-1:1,:,:,end:-1:1,:] + SMatrix3)/2
    SMatrix3 .= avgS

end # function
=#