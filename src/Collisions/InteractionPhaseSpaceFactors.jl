function PhaseSpaceFactors_Emi_Undo!(SMatrix::Array{Float64,4},pxr2::Vector{Float64},pyr2::Vector{Float64},pxr1::Vector{Float64},pyr1::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (pyr1[ii+1]-pyr1[ii])*(pxr1[jj+1]-pxr1[jj]) #dp1dmu1
        SMatrix[ll,kk,jj,ii] *= (pyr2[kk+1]-pyr2[kk])*(pxr2[ll+1]-pxr2[ll]) #dp2dmu2
    end

end # function

function PhaseSpaceFactors_Binary_Undo!(p3_bounds::Vector{Float64},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},p4_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},p1_bounds::Vector{Float64},u1_bounds::Vector{Float64},h1_bounds::Vector{Float64},p2_bounds::Vector{Float64},u2_bounds::Vector{Float64},h2_bounds::Vector{Float64};SMatrix3=0f0,SMatrix4=0f0,TMatrix1=0f0,TMatrix2=0f0)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.

    if SMatrix3 == 0f0 && SMatrix4 == 0f0 && TMatrix1 == 0f0 && TMatrix2 == 0f0
        error("No Matrices to apply phase space factors to.")
    end

    # === SMatrix3 === #
    # Momentum space volume elements
    if SMatrix3 != 0f0
        for  h2 in axes(SMatrix3,9), u2 in axes(SMatrix3,8), p2 in axes(SMatrix3,7), h1 in axes(SMatrix3,6), u1 in axes(SMatrix3,5), p1 in axes(SMatrix3,4) # common axes
            for h3 in axes(SMatrix3,3), u3 in axes(SMatrix3,2), p3 in axes(SMatrix3,1)
                if p3 == 1 # must account for underflow values increasing bin size 
                    SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[p3+1])# dp3du3dh3
                elseif p3 == size(SMatrix3,1) # overflow bin size assumed to be 1*maximum p3_bounds 
                    SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[end]) # dp3du3dh3
                else
                    SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3])*(p3_bounds[p3+1]-p3_bounds[p3]) # dp3du3dh3 
                end
                SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (u1_bounds[u1+1]-u1_bounds[u1])*(h1_bounds[h1+1]-h1_bounds[h1])*(p1_bounds[p1+1]-p1_bounds[p1]) # dp1du1dh1
                SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (u2_bounds[u2+1]-u2_bounds[u2])*(h2_bounds[h2+1]-h2_bounds[h2])*(p2_bounds[p2+1]-p2_bounds[p2]) # dp2du2dh2
            end
            if TMatrix1 != 0f0
                TMatrix1[p1,u1,h1,p2,u2,h2] /= (u2_bounds[u2+1]-u2_bounds[u2])*(h2_bounds[h2+1]-h2_bounds[h2])*(p2_bounds[p2+1]-p2_bounds[p2]) # dp2du2dh2   
            end
            if TMatrix2 != 0f0
                TMatrix2[p2,u2,h2,p1,u1,h1] /= (u1_bounds[u1+1]-u1_bounds[u1])*(h1_bounds[h1+1]-h1_bounds[h1])*(p1_bounds[p1+1]-p1_bounds[p1]) # dp1du1dh1
            end
        end
    end

    # === SMatrix4 === #
    # Momentum space volume elements
    if SMatrix4 != 0f0
        for h2 in axes(SMatrix4,9), u2 in axes(SMatrix4,8), p2 in axes(SMatrix4,7), h1 in axes(SMatrix4,6), u1 in axes(SMatrix4,5), p1 in axes(SMatrix4,4) # common axes
            for h4 in axes(SMatrix4,3), u4 in axes(SMatrix4,2), p4 in axes(SMatrix4,1)
                if p4 == 1 # must account for underflow values increasing bin size 
                    SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[p4+1])# dp4du4dh4 
                elseif p4 == size(SMatrix4,1) # overflow bin size assumed to be 1*maximum p4_bounds 
                    SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[end]) # dp4du4dh4
                else
                    SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4])*(p4_bounds[p4+1]-p4_bounds[p4]) # dp4du4dh4
                end 
                SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (u1_bounds[u1+1]-u1_bounds[u1])*(h1_bounds[h1+1]-h1_bounds[h1])*(p1_bounds[p1+1]-p1_bounds[p1]) # dp1du1dh1
                SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (u2_bounds[u2+1]-u2_bounds[u2])*(h2_bounds[h2+1]-h2_bounds[h2])*(p2_bounds[p2+1]-p2_bounds[p2]) # dp2du2dh2
            end
        end
    end

    return nothing

end
