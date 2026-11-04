"""
    BuildBigMatrices(PhaseSpace,DataDirectory;loading_check)

Function that builds the big matrices associated with binary and emissive interactions. If there are such interactions, first space is allocated for the arrays, then data is loaded into these arrays from the desired `DataDirectory` location and finally the big matrices are returned as an immutable `BinaryMatricesStruct`.
"""
function BuildBinaryMatrices(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;symmetric::Bool=false,loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false,gain_loss::Union{Nothing,String}=nothing)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    n = PhaseSpace.Grids.n_momentum
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            M_Bin = spzeros(Precision,0,0)
        else
            M_Bin = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            M_Bin_I::Vector{Int32} = Int32[]
            M_Bin_J::Vector{Int32} = Int32[]
            M_Bin_V::Vector{Precision} = Precision[]
            LoadMatrices_Binary(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;symmetric=symmetric,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V,gain_loss=gain_loss)
            println("Building sparse M_Bin")
            println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
            println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
            M_Bin = SparseArrays.sparse!(M_Bin_I,M_Bin_J,M_Bin_V,m,n)::SparseMatrixCSC{Precision,Int32}
   
            GC.gc()
        else
            M_Bin = zeros(Precision,m,n)::Matrix{Precision}
            LoadMatrices_Binary(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;symmetric=symmetric,M_Bin=M_Bin,gain_loss=gain_loss)
        end
    end

    size = Base.summarysize(M_Bin)

    if size > 1e9
        println("M_Bin is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("M_Bin is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("M_Bin is approx. $(size/1e3) KB in memory")
    else
        println("M_Bin is approx. $size bytes in memory")
    end

    return BinaryMatricesStruct{Precision}(M_Bin,Binary_list,Domain)

end


"""
    Fill_M_Bin!(...,name_locs,PhaseSpace;GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

Fills the big matrix `M_Bin` directly if dense or the vectors of rows, columns and values `M_Bin_I`, `M_Bin_J``, `M_Bin_V` if sparse, with the interaction rates for a specific binary interactions given by `name_locs` and the collision arrays `GainMatrix3`, `GainMatrix4`, `LossMatrix1`, `LossMatrix2`.
"""
function Fill_M_Bin!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::ZArray,GainMatrix4::ZArray,LossMatrix::ZArray,GainScale::Float64,LossScale::Float64,Indistinguishable_12::Bool;mode::AbstractMode=Ani(),symmetric::Bool=false,M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int32}}=nothing,M_Bin_J::Union{Nothing,Vector{Int32}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    n_momentum = PhaseSpace.Grids.n_momentum
    m = n_momentum*n_momentum

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]

    if !isempty(M_Bin_I) && !isempty(M_Bin_J) && !isempty(M_Bin_V)
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
    end

    if GainScale != 0.0
        GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,GainScale;symmetric,M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

        # before moving onto the next interaction lets combine the duplicate entries in the sparse matrix representation of M_Bin
        println("before combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))
        if isnothing(M_Bin)
            dedup_triplets!(M_Bin_I,M_Bin_J,M_Bin_V)
        end
        println("after combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))
        GC.gc()

        GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,GainScale;symmetric,M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

        # before moving onto the next interaction lets combine the duplicate entries in the sparse matrix representation of M_Bin
        println("before combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))
        if isnothing(M_Bin)
            dedup_triplets!(M_Bin_I,M_Bin_J,M_Bin_V)
        end
        println("after combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))
        GC.gc()
    end

    if LossScale != 0.0
        LossMatrix_to_M_Bin!(PhaseSpace,LossMatrix,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,LossScale,Indistinguishable_12;symmetric,M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

        # before moving onto the next interaction lets combine the duplicate entries in the sparse matrix representation of M_Bin
        println("before combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))
        if isnothing(M_Bin)
            dedup_triplets!(M_Bin_I,M_Bin_J,M_Bin_V)
        end
        println("after combining duplicates")
        println("max I: ",maximum(M_Bin_I)," max J: ",maximum(M_Bin_J)," max V: ",maximum(M_Bin_V))
        println("min I: ",minimum(M_Bin_I)," min J: ",minimum(M_Bin_J)," min V: ",minimum(M_Bin_V))
        println("length I: ",length(M_Bin_I)," length J: ",length(M_Bin_J)," length V: ",length(M_Bin_V))

    end

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix = nothing

    GC.gc()

    return M_Bin_I, M_Bin_J, M_Bin_V

end

function GainMatrix_to_M_Bin!(PhaseSpace::PhaseSpaceStruct,GainMatrix::ZArray,offset3::Int64,offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},dpy3::Vector{Float64},dpz3::Vector{Float64},GainScale::Float64;symmetric::Bool=false,M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int32}}=nothing,M_Bin_J::Union{Nothing,Vector{Int32}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = PhaseSpace.Grids.n_momentum
    println("N = $N")
    #println("$offset2")

    #=E = zeros(Float64,N)
    for species in eachindex(PhaseSpace.name_list)
        px_num = PhaseSpace.Momentum.px_num_list[species]
        py_num = PhaseSpace.Momentum.py_num_list[species]
        pz_num = PhaseSpace.Momentum.pz_num_list[species]
        dE = PhaseSpace.Grids.dE_list[species]
        for px in 1:px_num
            for py in 1:py_num
                for pz in 1:pz_num
                    idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                    E[idx] = dE[px]
                end
            end
        end
    end=#

    is_sparse = isnothing(M_Bin)

    ChunkGainMatrix::Array{Float64,9} = zeros(Float64,px3_num+2,py3_num,pz3_num,1,py1_num,pz1_num,1,py2_num,pz2_num)

    GainMatrixLocal = @view(ChunkGainMatrix[:,:,:,1,:,:,1,:,:])

    wIso = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3) * sum(dpy1) * sum(dpy2) * sum(dpy3))
    wAxi = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3))
    wAni = 1.0

    for px2 in 1:px2_num, px1 in 1:px1_num

        r = CartesianIndices((1:(px3_num+2),1:py3_num,1:pz3_num,px1:px1,1:py1_num,1:pz1_num,px2:px2,1:py2_num,1:pz2_num))
        Zarr.readblock!(ChunkGainMatrix,GainMatrix,r) # load chunk of Zarray

        #GainMatrixLocal .= GainMatrix[:,:,:,px1,:,:,px2,:,:] # load chunk of Zarray
        
    for px3 in 1:px3_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso

            val = 0.0 
            w = wIso

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in 1:py1_num, pz1 in 1:pz1_num, py2 in 1:py2_num, pz2 in 1:pz2_num, py3 in 1:py3_num, pz3 in 1:pz3_num
                val += GainMatrixLocal[px3+1,py3,pz3,py1,pz1,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3] * dpy1[py1] * dpy2[py2] * dpy3[py3]
            end

        end

        for py1 in 1:py1_num, py2 in 1:py2_num, py3 in 1:py3_num

            if mode isa Axi

                val = 0.0 
                w = wAxi

                # average over incoming and outgoing phi angles (pz)
                for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                    val += GainMatrixLocal[px3+1,py3,pz3,py1,pz1,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3]
                end

            end

            for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                
                if mode isa Ani

                    #=if is_sparse
                        #GainMax = maximum(@view(GainMatrix[:,py3,pz3,px1,py1,pz1,px2,py2,pz2]))
                        GainMax = maximum(@view(GainMatrixLocal[:,:,:,py1,pz1,py2,pz2]))
                    end=#

                    val = GainMatrixLocal[px3+1,py3,pz3,py1,pz1,py2,pz2]
                    w = wAni

                    #if is_sparse && val*w < eps(GainMax) # skips values smaller than this value to reduce memory usage
                    #    continue
                    #end

                end

                if val == 0.0 || !isfinite(val)
                    if !isfinite(val)
                        println("val = $val, px1 = $px1, py1 = $py1, pz1 = $pz1, px2 = $px2, py2 = $py2, pz2 = $pz2, px3 = $px3, py3 = $py3, pz3 = $pz3")
                    end
                    continue
                end

                a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
                b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                if is_sparse
                    if symmetric 
                        # symmetric in jk, non M-Matrix structure but good for Jacobian 
                        if b == c 
                            push!(M_Bin_I,Int32((b-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(c))
                            push!(M_Bin_V,convert(F,val*w * GainScale #=* E[a] / E[b] / E[c]=#))
                        else
                            push!(M_Bin_I,Int32((b-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(c))
                            push!(M_Bin_V,convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#))
                            push!(M_Bin_I,Int32((c-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(b))
                            push!(M_Bin_V,convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#))
                        end
                    else
                        # M_Bin terms allocated symmetrically except for if offset3==offset1 or offset3==offset2 but not both then assigned to the ij diagonal to match M-Matrix structure 
                        if offset3 == offset1 && offset3 != offset2
                            push!(M_Bin_I,Int32((b-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(c))
                            push!(M_Bin_V,convert(F,val*w * GainScale #=* E[a] / E[b] / E[c]=#))
                        elseif offset3 == offset2 && offset3 != offset1 
                            push!(M_Bin_I,Int32((c-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(b))
                            push!(M_Bin_V,convert(F,val*w * GainScale #=* E[a] / E[b] / E[c]=#))
                        else # asign symmetrically in jk
                            push!(M_Bin_I,Int32((b-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(c))
                            push!(M_Bin_V,convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#))
                            push!(M_Bin_I,Int32((c-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(b))
                            push!(M_Bin_V,convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#))
                        end
                    end
                else
                    if symmetric 
                        # symmetric in jk, non M-Matrix structure but good for Jacobian 
                        M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#)
                        M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#)
                    else
                        # M_Bin terms allocated symmetrically except for if offset3==offset1 or offset3==offset2 but not both then assigned to the ij diagonal to match M-Matrix structure 
                        if offset3 == offset1 && offset3 != offset2
                            M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w * GainScale #=* E[a] / E[b] / E[c]=#)
                        elseif offset3 == offset2 && offset3 != offset1 
                            M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w * GainScale #=* E[a] / E[b] / E[c]=#)
                        else # asign symmetrically in jk
                            M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#)
                            M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2 * GainScale #=* E[a] / E[b] / E[c]=#)
                        end
                    end
                end

            end # pz loop

        end # py loop

    end # px3 loop

    end # px1 and px2 loop

end

function LossMatrix_to_M_Bin!(PhaseSpace::PhaseSpaceStruct,LossMatrix::ZArray,offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},LossScale::Float64,Indistinguishable_12::Bool;symmetric::Bool=false,M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int32}}=nothing,M_Bin_J::Union{Nothing,Vector{Int32}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = PhaseSpace.Grids.n_momentum
    println("N = $N")

    is_sparse = isnothing(M_Bin)

    #=E = zeros(Float64,N)
    for species in eachindex(PhaseSpace.name_list)
        px_num = PhaseSpace.Momentum.px_num_list[species]
        py_num = PhaseSpace.Momentum.py_num_list[species]
        pz_num = PhaseSpace.Momentum.pz_num_list[species]
        dE = PhaseSpace.Grids.dE_list[species]
        for px in 1:px_num
            for py in 1:py_num
                for pz in 1:pz_num
                    idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                    E[idx] = dE[px]
                end
            end
        end
    end=#

    ChunkLossMatrix::Array{Float64,6} = zeros(Float64,1,py1_num,pz1_num,1,py2_num,pz2_num)
    LossMatrixLocal = @view(ChunkLossMatrix[1,:,:,1,:,:])

    wIso = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))
    wAxi = 1.0 / (sum(dpz1) * sum(dpz2))
    wAni = 1.0

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        r = CartesianIndices((px1:px1, 1:py1_num, 1:pz1_num, px2:px2, 1:py2_num, 1:pz2_num))
        Zarr.readblock!(ChunkLossMatrix,LossMatrix,r) # load chunk of Zarray

        #LossMatrixLocal .= LossMatrix[px1,:,:,px2,:,:] # load chunk of Zarray

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = wIso

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrixLocal[py1,pz1,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = wAxi

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrixLocal[py1,pz1,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrixLocal[py1,pz1,py2,pz2]
                    w = wAni

                end

                if val == 0.0 || !isfinite(val)
                    if !isfinite(val)
                        println("val = $val, px1 = $px1, py1 = $py1, pz1 = $pz1, px2 = $px2, py2 = $py2, pz2 = $pz2, px3 = $px3, py3 = $py3, pz3 = $pz3")
                    end
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                b = a
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2
                d = c

                # M_Bin terms allocated symmetrically
                if is_sparse
                    if symmetric 
                        # symmetric in jk, non M-Matrix structure but good for Jacobian 
                        push!(M_Bin_I,Int32((b-1)*N+(a-1)+1))
                        push!(M_Bin_J,Int32(c))
                        push!(M_Bin_V,-convert(F,val*w/2 * LossScale))
                        push!(M_Bin_I,Int32((c-1)*N+(a-1)+1))
                        push!(M_Bin_J,Int32(b))
                        push!(M_Bin_V,-convert(F,val*w/2 * LossScale))
                        if !Indistinguishable_12 # need to add Loss terms corresponding to particle 2 being different from particle 1 (swap a with c and b with d)
                            push!(M_Bin_I,Int32((d-1)*N+(c-1)+1))
                            push!(M_Bin_J,Int32(a))
                            push!(M_Bin_V,-convert(F,val*w/2 * LossScale))
                            push!(M_Bin_I,Int32((a-1)*N+(c-1)+1))
                            push!(M_Bin_J,Int32(d))
                            push!(M_Bin_V,-convert(F,val*w/2 * LossScale))
                        end
                    else
                        # diagonal in ij entries so that total matrix has an M-Matrix structure
                        push!(M_Bin_I,Int32((a-1)*N+(a-1)+1))
                        push!(M_Bin_J,Int32(c))
                        push!(M_Bin_V,-convert(F,val*w/2 * LossScale #=* E[a] / E[b] / E[c]=#))
                        push!(M_Bin_I,Int32((c-1)*N+(c-1)+1))
                        push!(M_Bin_J,Int32(b))
                        push!(M_Bin_V,-convert(F,val*w/2 * LossScale #=* E[a] / E[b] / E[c]=#))
                        if !Indistinguishable_12 # need to add Loss terms corresponding to particle 2 being different from particle 1 (swap a with c and b with d)
                            push!(M_Bin_I,Int32((c-1)*N+(c-1)+1))
                            push!(M_Bin_J,Int32(a))
                            push!(M_Bin_V,-convert(F,val*w/2 * LossScale #=* E[a] / E[b] / E[c]=#))
                            push!(M_Bin_I,Int32((a-1)*N+(a-1)+1))
                            push!(M_Bin_J,Int32(d))
                            push!(M_Bin_V,-convert(F,val*w/2 * LossScale #=* E[a] / E[b] / E[c]=#))
                        end
                    end
                else
                    if symmetric 
                        # symmetric in jk, non M-Matrix structure but good for Jacobian 
                        M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2 * LossScale)
                        M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2 * LossScale)
                        if !Indistinguishable_12 # need to add Loss terms corresponding to particle 2 being different from particle 1 (swap a with c and b with d)
                            M_Bin[(d-1)*N+(c-1)+1,a] -= convert(F,val*w/2 * LossScale)
                            M_Bin[(a-1)*N+(c-1)+1,d] -= convert(F,val*w/2 * LossScale)
                        end
                    else
                        # diagonal in ij entries so that total matrix has an M-Matrix structure
                        M_Bin[(a-1)*N+(a-1)+1,c] -= convert(F,val*w/2 * LossScale)
                        M_Bin[(c-1)*N+(c-1)+1,b] -= convert(F,val*w/2 * LossScale)
                        if !Indistinguishable_12 # need to add Loss terms corresponding to particle 2 being different from particle 1 (swap a with c and b with d)
                            M_Bin[(c-1)*N+(c-1)+1,a] -= convert(F,val*w/2 * LossScale)
                            M_Bin[(a-1)*N+(a-1)+1,d] -= convert(F,val*w/2 * LossScale)
                        end
                    end
                end

            end # pz loop

        end # py loop

    end # px loop

end

@inline function keylt(rows, cols, i::Int, rj, cj)
    ri = rows[i]
    ci = cols[i]
    ri < rj || (ri == rj && ci < cj)
end

@inline function swap3!(rows, cols, vals, i::Int, j::Int)
    i == j && return
    rows[i], rows[j] = rows[j], rows[i]
    cols[i], cols[j] = cols[j], cols[i]
    vals[i], vals[j] = vals[j], vals[i]
end

function quicksort3!(rows, cols, vals, lo::Int, hi::Int)
    lo >= hi && return

    i = lo
    j = hi

    # pivot value, not pivot index
    mid = (lo + hi) >>> 1
    pr = rows[mid]
    pc = cols[mid]

    while i <= j
        while keylt(rows, cols, i, pr, pc)
            i += 1
        end
        while (pr < rows[j]) || (pr == rows[j] && pc < cols[j])
            j -= 1
        end

        if i <= j
            swap3!(rows, cols, vals, i, j)
            i += 1
            j -= 1
        end
    end

    if lo < j
        quicksort3!(rows, cols, vals, lo, j)
    end
    if i < hi
        quicksort3!(rows, cols, vals, i, hi)
    end

    return nothing
end

function combine_sorted_duplicates!(rows, cols, vals; op=+)
    n = length(rows)
    n <= 1 && return n

    k = 1
    for i in 2:n
        if rows[i] == rows[k] && cols[i] == cols[k]
            vals[k] = op(vals[k], vals[i])
        else
            k += 1
            if k != i
                rows[k] = rows[i]
                cols[k] = cols[i]
                vals[k] = vals[i]
            end
        end
    end

    resize!(rows, k)
    resize!(cols, k)
    resize!(vals, k)

    return nothing
end

function dedup_triplets!(rows::Vector{I}, cols::Vector{I}, vals::Vector{T}; op=+) where {I,T}
    @assert length(rows) == length(cols) == length(vals)

    n = length(rows)
    n <= 1 && return n

    quicksort3!(rows, cols, vals, 1, n)
    combine_sorted_duplicates!(rows, cols, vals; op=op)

    return nothing
end

#=
function BuildBinaryMatricesPatankar(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            Gijk = spzeros(Precision,0,0)
            Liij = spzeros(Precision,0,0)
        else
            Gijk = zeros(Precision,0,0)
            Liij = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            Gijk_I::Vector{Int64} = Int64[]
            Gijk_J::Vector{Int64} = Int64[]
            Gijk_V::Vector{Precision} = Precision[]
            
            Lij_I::Vector{Int64} = Int64[]
            Lij_J::Vector{Int64} = Int64[]
            Lij_V::Vector{Precision} = Precision[]
            LoadMatrices_BinaryPatankar(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk_I=Gijk_I,Gijk_J=Gijk_J,Gijk_V=Gijk_V,Lij_I=Lij_I,Lij_J=Lij_J,Lij_V=Lij_V)
            println("Building sparse Gijk and Liij")
            Gijk = sparse(Gijk_I,Gijk_J,Gijk_V,m,n)::SparseMatrixCSC{Precision,Int64}
            Lij = sparse(Lij_I,Lij_J,Lij_V,n,n)::SparseMatrixCSC{Precision,Int64}
   
            GC.gc()
        else
            Gijk = zeros(Precision,m,n)::Matrix{Precision}
            Lij = zeros(Precision,n,n)::Matrix{Precision}
            LoadMatrices_BinaryPatankar(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk=Gijk,Lij=Lij)
        end
    end

    size = Base.summarysize(Gijk)

    if size > 1e9
        println("Gijk is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("Gijk is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("Gijk is approx. $(size/1e3) KB in memory")
    else
        println("Gijk is approx. $size bytes in memory")
    end

    return BinaryMatricesStructPatankar{Precision}(Gijk,Lij,Binary_list,Domain)

end

function BuildBinaryMatricesPatankarSymmetric(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            Gijk = spzeros(Precision,0,0)
            Lijk = spzeros(Precision,0,0)
            Liij = spzeros(Precision,0,0)
            Liik = spzeros(Precision,0,0)
        else
            Gijk = zeros(Precision,0,0)
            Lijk = zeros(Precision,0,0)
            Liij = zeros(Precision,0,0)
            Liik = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            Gijk_I::Vector{Int64} = Int64[]
            Gijk_J::Vector{Int64} = Int64[]
            Gijk_V::Vector{Precision} = Precision[]
            
            Lijk_I::Vector{Int64} = Int64[]
            Lijk_J::Vector{Int64} = Int64[]
            Lijk_V::Vector{Precision} = Precision[]

            Liij_I::Vector{Int64} = Int64[]
            Liij_J::Vector{Int64} = Int64[]
            Liij_V::Vector{Precision} = Precision[]

            Liik_I::Vector{Int64} = Int64[]
            Liik_J::Vector{Int64} = Int64[]
            Liik_V::Vector{Precision} = Precision[]

            LoadMatrices_BinaryPatankarSymmetric(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk_I=Gijk_I,Gijk_J=Gijk_J,Gijk_V=Gijk_V,Lijk_I=Lijk_I,Lijk_J=Lijk_J,Lijk_V=Lijk_V,Liij_I=Liij_I,Liij_J=Liij_J,Liij_V=Liij_V,Liik_I=Liik_I,Liik_J=Liik_J,Liik_V=Liik_V)
            println("Building sparse Gijk and Ljk")
            Gijk = sparse(Gijk_I,Gijk_J,Gijk_V,m,n)::SparseMatrixCSC{Precision,Int64}
            Lijk = sparse(Lijk_I,Lijk_J,Lijk_V,m,n)::SparseMatrixCSC{Precision,Int64}

            Liij = sparse(Liij_I,Liij_J,Liij_V,m,n)::SparseMatrixCSC{Precision,Int64}
            Liik = sparse(Liik_I,Liik_J,Liik_V,m,n)::SparseMatrixCSC{Precision,Int64}
   
            GC.gc()
        else
            Gijk = zeros(Precision,m,n)::Matrix{Precision}
            Lijk = zeros(Precision,m,n)::Matrix{Precision}
            Liij = zeros(Precision,m,n)::Matrix{Precision}
            Liik = zeros(Precision,m,n)::Matrix{Precision}
            LoadMatrices_BinaryPatankarSymmetric(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk=Gijk,Lijk=Lijk,Liij=Liij,Liik=Liik)
        end
    end

    size = Base.summarysize(Gijk)

    if size > 1e9
        println("Gijk is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("Gijk is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("Gijk is approx. $(size/1e3) KB in memory")
    else
        println("Gijk is approx. $size bytes in memory")
    end

    return BinaryMatricesStructPatankarSymmetric{Precision}(Gijk,Lijk,Liij,Liik,Binary_list,Domain)

end

function BuildBinaryMatricesGraphLaplacian(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            M_Bin = spzeros(Precision,0,0)
        else
            M_Bin = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            M_Bin_I::Vector{Int64} = Int64[]
            M_Bin_J::Vector{Int64} = Int64[]
            M_Bin_V::Vector{Precision} = Precision[]
            LoadMatrices_BinaryGraphLaplacian(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
            println("Building sparse M_Bin")
            M_Bin = sparse(M_Bin_I,M_Bin_J,M_Bin_V,m,n)::SparseMatrixCSC{Precision,Int64}
   
            GC.gc()
        else
            M_Bin = zeros(Precision,m,n)::Matrix{Precision}
            LoadMatrices_BinaryGraphLaplacian(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;M_Bin=M_Bin)
        end
    end

    size = Base.summarysize(M_Bin)

    if size > 1e9
        println("M_Bin is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("M_Bin is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("M_Bin is approx. $(size/1e3) KB in memory")
    else
        println("M_Bin is approx. $size bytes in memory")
    end

    return BinaryMatricesStruct{Precision}(M_Bin,Binary_list,Domain)

end



function Fill_M_BinPatankar!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6},n_momentum::Int64;mode::AbstractMode=Ani(),Gijk::Union{Nothing,Matrix{F}}=nothing,Gijk_I::Union{Nothing,Vector{Int64}}=nothing,Gijk_J::Union{Nothing,Vector{Int64}}=nothing,Gijk_V::Union{Nothing,Vector{F}}=nothing,Lij::Union{Nothing,Matrix{F}}=nothing,Lij_I::Union{Nothing,Vector{Int64}}=nothing,Lij_J::Union{Nothing,Vector{Int64}}=nothing,Lij_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]

    GainMatrix_to_M_BinPatankar!(GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)
    GainMatrix_to_M_BinPatankar!(GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)
    LossMatrix_to_M_BinPatankar!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=Lij,M_Bin_I=Lij_I,M_Bin_J=Lij_J,M_Bin_V=Lij_V)
    LossMatrix_to_M_BinPatankar!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=Lij,M_Bin_I=Lij_I,M_Bin_J=Lij_J,M_Bin_V=Lij_V)

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    return nothing

end

function Fill_M_BinPatankarSymmetric!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6},n_momentum::Int64;mode::AbstractMode=Ani(),Gijk::Union{Nothing,Matrix{F}}=nothing,Gijk_I::Union{Nothing,Vector{Int64}}=nothing,Gijk_J::Union{Nothing,Vector{Int64}}=nothing,Gijk_V::Union{Nothing,Vector{F}}=nothing,Lijk::Union{Nothing,Matrix{F}}=nothing,Lijk_I::Union{Nothing,Vector{Int64}}=nothing,Lijk_J::Union{Nothing,Vector{Int64}}=nothing,Lijk_V::Union{Nothing,Vector{F}}=nothing,Liij::Union{Nothing,Matrix{F}}=nothing,Liij_I::Union{Nothing,Vector{Int64}}=nothing,Liij_J::Union{Nothing,Vector{Int64}}=nothing,Liij_V::Union{Nothing,Vector{F}}=nothing,Liik::Union{Nothing,Matrix{F}}=nothing,Liik_I::Union{Nothing,Vector{Int64}}=nothing,Liik_J::Union{Nothing,Vector{Int64}}=nothing,Liik_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]

    GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)
    GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)

    LossMatrix_to_M_BinPatankarjk!(PhaseSpace,LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=Lijk,M_Bin_I=Lijk_I,M_Bin_J=Lijk_J,M_Bin_V=Lijk_V)
    LossMatrix_to_M_BinPatankarjk!(PhaseSpace,LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=Lijk,M_Bin_I=Lijk_I,M_Bin_J=Lijk_J,M_Bin_V=Lijk_V)
    LossMatrix_to_M_BinPatankarij!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=Liij,M_Bin_I=Liij_I,M_Bin_J=Liij_J,M_Bin_V=Liij_V)
    LossMatrix_to_M_BinPatankarij!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=Liij,M_Bin_I=Liij_I,M_Bin_J=Liij_J,M_Bin_V=Liij_V)

    LossMatrix_to_M_BinPatankarik!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=Liik,M_Bin_I=Liik_I,M_Bin_J=Liik_J,M_Bin_V=Liik_V)
    LossMatrix_to_M_BinPatankarik!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=Liik,M_Bin_I=Liik_I,M_Bin_J=Liik_J,M_Bin_V=Liik_V)

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    return nothing

end

function Fill_M_BinGraphLaplacian!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6},n_momentum::Int64;mode::AbstractMode=Ani(),M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]


        GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
        GainMatrix_to_M_Bin!(PhaseSpace,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
        #LossMatrix_to_M_Bin!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
        #LossMatrix_to_M_Bin!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

        # Apply redistribution of gain spectra to form the graph Laplacian form of the M_Bin matrix
        M_BinGraphLaplacian!(PhaseSpace,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)


    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    return nothing

end

function M_BinGraphLaplacian!(PhaseSpace::PhaseSpaceStruct,n_momentum::Int64;M_Bin::Union{Nothing,Matrix{T}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{T}}=nothing) where T<:Union{Float32,Float64}

    N = n_momentum
    tol = T(1e-16)
    Precision::DataType = getfield(Main,Symbol("Precision"))

    E = zeros(Float64,N)
    for species in eachindex(PhaseSpace.name_list)
        px_num = PhaseSpace.Momentum.px_num_list[species]
        py_num = PhaseSpace.Momentum.py_num_list[species]
        pz_num = PhaseSpace.Momentum.pz_num_list[species]
        dE = PhaseSpace.Grids.dE_list[species]
        for px in 1:px_num
            for py in 1:py_num
                for pz in 1:pz_num
                    idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                    E[idx] = dE[px]
                end
            end
        end
    end

    is_sparse = isnothing(M_Bin)
    # build Bijk, a version of M_Bin where the gain spectra are distributed equally among jk and loss terms are on the ij diagonal
    if is_sparse
        B = reshape(sparse(M_Bin_I,M_Bin_J,M_Bin_V,N*N,N),N,N,N)
        # empty the triplet vectors ready to be rebuilt
        empty!(M_Bin_I)
        empty!(M_Bin_J)
        empty!(M_Bin_V)
    else
        B = reshape(M_Bin,N,N,N)
        # M_Bin to build M_Bin to overwrite it at the end
        M_Bin_I::Vector{Int64} = Int64[]
        M_Bin_J::Vector{Int64} = Int64[]
        M_Bin_V::Vector{Precision} = Precision[]
    end

    # Clear output triplets


    # Small channel tolerance based on energy scale
    Etol = T(100) * sqrt(eps(T))

    # Precompute energy sort once
    order = zeros(Int64, N)

    # delta_mass_target
    delta_mass_target = 0.0

    for j in 1:N
        for k in j:N

            # Symmetric gain slice for unordered pair {j,k}
            if j == k
                g_real = collect(@view B[:, j, j])
            else
                g_real = collect(@view B[:, j, k]) .+ collect(@view B[:, k, j])
            end

            if j == 1 || j == 288 || j == 888 || k == 1 || k == 288 || k == 888
                # TODO: this prevents states which energy are not in the convex hull from be affected
                continue
            end

            if !all(isfinite, g_real)
                # no interaction for this channel, skip
                continue
            end

            M_real = sum(g_real)
            if M_real == T(0.0)
                # no interaction for this channel, skip
                continue
            end

            M = M_real / M_real
            g = g_real ./ M_real

            # energy check 
            energy_check = sum(g_real .* E) - 0.5 * sum(g_real) * (E[j] + E[k])
            println("Energy check for channel ($j,$k): $energy_check, Egain = $(sum(g_real .* E)), Eloss = $(0.5 * sum(g_real) * (E[j] + E[k])))")

            Ebar = dot(E, g)
            if !isfinite(Ebar)
                @warn "Non-finite Ebar"
                continue
            end

            if j == k
                # Single-donor case
                if abs(Ebar - E[j]) > Etol
                    @warn "Skipping self-channel ($j,$k): energy mismatch" Ebar E[j]
                    continue
                end

                # Flattened tensor row for (i,j)
                @inbounds for i in 1:N
                    a = g_real[i]
                    if a != 0
                        push!(M_Bin_I, flatidx(i, j, N))
                        push!(M_Bin_J, j)
                        push!(M_Bin_V, a)
                    end
                end

                # Diagonal loss: A[j,j,j] -= M_real
                push!(M_Bin_I, flatidx(j, j, N))
                push!(M_Bin_J, j)
                push!(M_Bin_V, -M_real)

                continue
            end

            # Two-donor split
            midpoint = T(0.5) * M * (E[j] + E[k])

            # Mass split target (equal donors)
            mass_target   = T(0.5) * M
            energy_target = T(0.5) * M * E[j]

            if abs(Ebar - midpoint) > T(1e-4)
                @warn "Channel ($j,$k) midpoint mismatch; using symmetric half split, Ebar=$Ebar, midpoint=$midpoint, E[j]=$(E[j]), E[k]=$(E[k])"
                xj = T(0.5) .* g
                xk = g .- xj
                ok = true
            else
                xj, xk, ok, delta_mass_target = ordered_greedy_two_donor_split(g, E, mass_target, energy_target,j,k)
                if ok
                    break
                end

                if !ok
                    @warn "Greedy split failed for channel ($j,$k); using best result"
                end
            end

            # Store donor-j column for partner k
            @inbounds for i in 1:N
                a = xj[i] * M_real
                if a != 0
                    push!(M_Bin_I, flatidx(i, j, N))
                    push!(M_Bin_J, k)
                    push!(M_Bin_V, a)
                end
            end
            push!(M_Bin_I, flatidx(j, j, N))
            push!(M_Bin_J, k)
            push!(M_Bin_V, - (mass_target+delta_mass_target) * M_real)

            # Store donor-k column for partner j
            @inbounds for i in 1:N
                a = xk[i] * M_real
                if a != 0
                    push!(M_Bin_I, flatidx(i, k, N))
                    push!(M_Bin_J, j)
                    push!(M_Bin_V, a)
                end
            end
            push!(M_Bin_I, flatidx(k, k, N))
            push!(M_Bin_J, j)
            push!(M_Bin_V, -mass_target * M_real)
        end
    end

    if !is_sparse
        # Build dense M_Bin from triplets
        fill!(M_Bin, zero(T))
        for idx in 1:length(M_Bin_I)
            i = M_Bin_I[idx]
            j = M_Bin_J[idx]
            v = M_Bin_V[idx]
            M_Bin[i, j] = v
        end
        M_Bin_I = nothing
        M_Bin_J = nothing
        M_Bin_V = nothing
    end

    return nothing
end


@inline flatidx(i::Int, j::Int, N::Int) = i + (j - 1) * N

"""
Greedy two-donor split of a nonnegative gain spectrum g.

We seek xj and xk = g - xj such that
    sum(xj)      = mass_target
    dot(E, xj)   = energy_target

The routine starts from a proportional (or symmetric) split and then
greedily swaps mass between low-energy and high-energy bins to match energy.

Returns (xj, xk, ok).
"""
function greedy_two_donor_split(g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T,order::AbstractVector{Int},j,k;tol::T = T(sqrt(eps(T))),maxiter::Int = 10_000) where {T<:Real}
    N = length(g)
    @assert length(E) == N
    @assert length(order) == N

    M = sum(g)
    if !isfinite(M) || M <= tol
        return zeros(T, N), zeros(T, N), false
    end

    λ = clamp(mass_target / M, zero(T), one(T))
    xj = λ .* g
    err = energy_target - dot(E, xj)

    lo = 1
    hi = N

    δ_needed = 0.0
    cap = 0.0

    for _ in 1:maxiter
        if abs(err) <= tol
            break
        end
        #println("$j, $k, $hi, $lo, energy error: $err, xjlo: $(xj[order[lo]]), xjhi: $(xj[order[hi]]), ghi: $(g[order[hi]]), glo: $(g[order[lo]]), Elo: $(E[order[lo]]), Ehi: $(E[order[hi]]), δneeded: $δ_needed, cap: $cap")


        # Advance lo to first bin with mass left in xj
        while lo <= N && xj[order[lo]] <= 0.0
            lo += 1
        end

        # Advance hi to first bin with room left in xj
        while hi >= 1 && xj[order[hi]] >= g[order[hi]]
            hi -= 1
        end

        if lo >= hi
            break
        end

        iL = order[lo]
        iH = order[hi]

        # Need a genuine energy gap
        if E[iH] <= E[iL]
            hi -= 1
            continue
        end

        dE = E[iH] - E[iL]

        if err > 0
            # Need MORE energy in xj:
            # move mass from low-energy bin -> high-energy bin
            cap = min(xj[iL], g[iH] - xj[iH])

            if cap <= 0
                # Exhausted pair: advance whichever side is exhausted
                if xj[iL] <= 0
                    lo += 1
                end
                if xj[iH] >= g[iH]
                    hi -= 1
                end
                continue
            end

            δ_needed = err / dE
            δ = min(δ_needed, cap)

            if δ <= tol
                # No effective move possible
                if xj[iL] <= tol
                    lo += 1
                end
                if xj[iH] >= g[iH] - tol
                    hi -= 1
                end
                continue
            end

            xj[iL] -= δ
            xj[iH] += δ

        else
            # Need LESS energy in xj:
            # move mass from high-energy bin -> low-energy bin
            cap = min(g[iL] - xj[iL], xj[iH])

            if cap <= tol
                if xj[iL] >= g[iL] - tol
                    lo += 1
                end
                if xj[iH] <= tol
                    hi -= 1
                end
                continue
            end

            δ_needed = (-err) / dE
            δ = min(δ_needed, cap)

            if δ <= tol
                if xj[iL] >= g[iL] - tol
                    lo += 1
                end
                if xj[iH] <= tol
                    hi -= 1
                end
                continue
            end

            xj[iL] += δ
            xj[iH] -= δ
        end

        # Recompute residual exactly (prevents sign drift)
        err = energy_target - dot(E, xj)

        # If the bins we just used are now exhausted, advance pointers
        if xj[iL] <= 0
            lo += 1
        end
        if xj[iH] >= g[iH]
            hi -= 1
        end
    end

    xj = clamp.(xj, zero(T), g)
    xk = g .- xj

    ok = abs(sum(xj) - mass_target) <= 10tol && abs(dot(E, xj) - energy_target) <= 10tol && all(xk .>= zero(T))

    if !ok
        println("$(abs(sum(xj) - mass_target) <= 10tol), $(abs(dot(E, xj) - energy_target) <= 10tol), $(all(xk .>= zero(T)))")
        println("Greedy split did not converge to target within tolerance, mass_target: $mass_target, sum(xj): $(sum(xj)), sum(xk): $(sum(xk)), energy_target: $energy_target, dot(E, xj): $(dot(E, xj)), energy error $(dot(E, xj) - energy_target), hi: $hi, lo: $lo, delta needed $δ_needed, cap: $cap")
    end

    return xj, xk, ok
end

function greedy_two_donor_split_bestpair(g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T,j,k;tol::T = 100T(sqrt(eps(T))),maxiter::Int = 10_000) where {T<:Real}

    N = length(g)
    M = sum(g)
    if !isfinite(M) || M <= tol
        return zeros(T, N), zeros(T, N), false
    end

    # Start from proportional split
    λ = clamp(mass_target / M, zero(T), one(T))
    xj = λ .* g

    for _ in 1:maxiter
        err = energy_target - dot(E, xj)

        if abs(err) <= tol
            xk = g .- xj
            ok = abs(sum(xj) - mass_target) <= 10tol &&
                 abs(dot(E, xj) - energy_target) <= 10tol &&
                 all(xk .>= -tol)
            return xj, xk, ok
        end

        iL, iH, δ, dE, exact = best_swap_pair(xj, g, E, err; tol=tol)

        if iL == 0 || iH == 0 || δ <= tol || dE <= tol
            break
        end

        if err > 0
            # move low -> high
            xj[iL] -= δ
            xj[iH] += δ
        else
            # move high -> low
            xj[iH] -= δ
            xj[iL] += δ
        end

        err = energy_target - dot(E, xj)

        if exact && abs(err) <= tol
            break
        end
    end

    xj = clamp.(xj, zero(T), g)
    xk = g .- xj

    ok = abs(sum(xj) - mass_target) <= 10tol &&
         abs(dot(E, xj) - energy_target) <= 10tol &&
         all(xk .>= -tol)

         
    if !ok
        println("$(abs(sum(xj) - mass_target) <= 10tol), $(abs(dot(E, xj) - energy_target) <= 10tol), $(all(xk .>= zero(T)))")
        println("Greedy split did not converge to target within tolerance, mass_target: $mass_target, sum(xj): $(sum(xj)), sum(xk): $(sum(xk)), energy_target: $energy_target, dot(E, xj): $(dot(E, xj)), energy error $(dot(E, xj) - energy_target)")
    end

    return xj, xk, ok
end

function ordered_greedy_two_donor_split(g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T,j::Int,k::Int;tol::T = T(sqrt(eps(T))),frac::T=T(1.0),maxiter::Int = 10_000) where {T<:Real}

    N = length(g)
    @assert length(E) == N

    M = sum(g)
    if !isfinite(M) || M <= tol
        return zeros(T, N), zeros(T, N), false
    end

    # Start from proportional split
    λ = mass_target / M

    g_scaled = similar(g)
    xj = similar(g)
    xk = similar(g)

    iter = 0 
    iter2 = 0
    delta_mass_target = 0.0

    ok = false

    while !ok && iter2 < 3
        iter2 += 1

        mass_target_scaled = mass_target + delta_mass_target

        g_scaled .= g
        g_scaled[j] += delta_mass_target

        xj .= λ .* g
        xj[j] += delta_mass_target

        # Helper: recompute residual
        err = energy_target - dot(E, xj)

        for _ in 1:maxiter
            iter += 1
            if abs(err) <= tol
                break
            end

            # Build candidate orderings from the current state
            if err < 0
                # Too much energy in xj:
                # source = high current energy contribution
                # destination = low energy bins with room
                src_order = sortperm(E .* xj, rev = true)
                dst_order = sortperm(E)
            else
                # Too little energy in xj:
                # source = low current energy contribution
                # destination = high energy bins with room
                src_order = sortperm(E .* xj)
                dst_order = sortperm(E, rev = true)
            end

            moved = false

            if err < 0
                # Need to decrease energy in xj:
                # move mass from a high-energy source bin to a low-energy destination bin
                for h in src_order

                    if xj[h] <= 0
                        continue
                    end

                    for l in dst_order
                        if g[l] - xj[l] <= 0
                            continue
                        end
                        if E[h] <= E[l]
                            continue
                        end

                        dE = E[h] - E[l]
                        cap = min(xj[h], g_scaled[l] - xj[l])
                        if cap <= 0
                            continue
                        end

                        delta_needed = (-err) / dE

                        delta = frac * min(delta_needed, cap)

                        if delta <= 0
                            continue
                        end

                        # Apply the swap
                        xj[h] -= delta
                        xj[l] += delta

                        # Recompute residual exactly
                        err = energy_target - dot(E, xj)
                        moved = true
                        break
                    end
                    if moved
                        break
                    end
                end

            else
                # Need to increase energy in xj:
                # move mass from a low-energy source bin to a high-energy destination bin
                for l in src_order
                
                    if xj[l] <= 0
                        continue
                    end

                    for h in dst_order
                        if g[h] - xj[h] <= 0
                            continue
                        end
                        if E[h] <= E[l]
                            continue
                        end

                        dE = E[h] - E[l]
                        cap = min(xj[l], g_scaled[h] - xj[h])
                        if cap <= 0
                            continue
                        end

                        delta_needed = err / dE

                        delta = frac * min(delta_needed, cap)

                        if delta <= 0
                            continue
                        end

                        # Apply the swap
                        xj[l] -= delta
                        xj[h] += delta

                        # Recompute residual exactly
                        err = energy_target - dot(E, xj)
                        moved = true
                        break
                    end
                    if moved
                        break
                    end
                end
            end

            if !moved
                break
            end
        end

        # Clean up
        xj = clamp.(xj, zero(T), g_scaled)
        xk = g_scaled .- xj

        ok = abs(sum(xj) - mass_target_scaled) <= tol && abs(dot(E, xj) - energy_target)/energy_target <= tol && all(xk .>= -tol)
        notok = false

        if !ok
            println("ok: $ok: $(abs(sum(xj) - mass_target_scaled) <= tol), $(abs(dot(E, xj) - energy_target)/energy_target <= tol), $(all(xk .>= zero(T)))")
            println("Greedy split did not converge to target within tolerance, mass_target: $mass_target_scaled, sum(xj): $(sum(xj)), sum(xk): $(sum(xk)), energy_target_j: $energy_target, dot(E, xj): $(dot(E, xj)), energy_target_k: $(T(0.5) * E[k]), dot(E, xk): $(dot(E, xk)), energy error $(dot(E, xj) - energy_target), total energy gain: $(dot(E,g)), iter: $iter, iter2: $iter2, delta_mass_target: $delta_mass_target")
            #println("Final xj:")
            #display(reshape(xj,37,24))
            #println("Final Exj:")
            #display(reshape(E .* xj,37,24))
            #println("Final xk:")
            #display(reshape(xk,37,24))
            #println("Final Exk:")
            #display(reshape(E .* xk,37,24))

            delta_mass_target = err/E[j]

            notok = true

        else
            if notok
                println("Greedy split did converge to target within tolerance, mass_target: $mass_target_scaled, sum(xj): $(sum(xj)), sum(xk): $(sum(xk)), energy_target_j: $energy_target, dot(E, xj): $(dot(E, xj)), energy_target_k: $(T(0.5) * E[k]), dot(E, xk): $(dot(E, xk)), energy error $(dot(E, xj) - energy_target), total energy gain: $(dot(E,g)), iter: $iter, iter2: $iter2, delta_mass_target: $delta_mass_target")
            end
        end

    end

    return xj, xk, ok, delta_mass_target
end

function best_swap_pair(xj::AbstractVector{T},g::AbstractVector{T},E::AbstractVector{T},err::T;tol::T = T(1e-10)) where {T<:Real}

    N = length(xj)
    best_iL = 0
    best_iH = 0
    best_delta = zero(T)
    best_dE = zero(T)

    target = abs(err)

    # Track:
    # 1) exact-fix candidates: d = cap*dE >= target
    #    choose the one with smallest d (least disturbance)
    # 2) if no exact fix exists, choose the one with largest d
    exact_found = false
    best_exact_d = typemax(T)

    best_partial_d = zero(T)

    if err > 0
        # Need to increase energy in xj:
        # move mass from low-E donor bin -> high-E receiver bin.
        for iL in 1:N
            if xj[iL] <= 0
                continue
            end
            for iH in 1:N
                if g[iH] - xj[iH] <= 0
                    continue
                end
                if E[iH] <= E[iL]
                    continue
                end

                cap = min(xj[iL], g[iH] - xj[iH])
                dE  = E[iH] - E[iL]
                d   = cap * dE

                if d >= target
                    # exact fix possible
                    if !exact_found || d < best_exact_d
                        exact_found = true
                        best_exact_d = d
                        best_iL = iL
                        best_iH = iH
                        best_delta = target / dE
                        best_dE = dE
                    end
                else
                    # partial fix only; maximize reduction
                    if !exact_found && d > best_partial_d
                        best_partial_d = d
                        best_iL = iL
                        best_iH = iH
                        best_delta = cap
                        best_dE = dE
                    end
                end
            end
        end
    else
        # Need to decrease energy in xj:
        # move mass from high-E donor bin -> low-E receiver bin.
        for iH in 1:N
            if xj[iH] <= 0
                continue
            end
            for iL in 1:N
                if g[iL] - xj[iL] <= 0
                    continue
                end
                if E[iH] <= E[iL]
                    continue
                end

                cap = min(xj[iH], g[iL] - xj[iL])
                dE  = E[iH] - E[iL]
                d   = cap * dE

                if d >= target
                    if !exact_found || d < best_exact_d
                        exact_found = true
                        best_exact_d = d
                        best_iL = iL
                        best_iH = iH
                        best_delta = target / dE
                        best_dE = dE
                    end
                else
                    if !exact_found && d > best_partial_d
                        best_partial_d = d
                        best_iL = iL
                        best_iH = iH
                        best_delta = cap
                        best_dE = dE
                    end
                end
            end
        end
    end

    return best_iL, best_iH, best_delta, best_dE, exact_found
end


"""
Project a nonnegative vector g onto the set

    0 <= x <= g
    sum(x) = mass_target
    dot(E, x) = energy_target

starting from a proportional split. This is a small 2-moment box-constrained
projection. It returns an exact feasible split when one exists.
"""
function project_two_moment_box(y::AbstractVector{T},g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T;tol::T = T(1e-16),maxiter::Int = 50) where {T<:Real}

    α = zero(T)
    β = zero(T)

    for _ in 1:maxiter
        x = clamp.(y .- α .- β .* E, zero(T), g)

        r1 = sum(x) - mass_target
        r2 = dot(E, x) - energy_target

        if abs(r1) < sqrt(tol) && abs(r2) < sqrt(tol)
            return x
        end

        active = (x .> tol) .& (x .< g .- tol)
        if !any(active)
            # No feasible split found for this channel
            println("r1: $(abs(r1)), r2: $(abs(r2)). No active set for projection. Returning zero split.")
            return zeros(T, length(g))
        end

        Ea = E[active]
        n1 = length(Ea)
        sE = sum(Ea)
        sE2 = sum(Ea .* Ea)

        # Jacobian of residual wrt (α, β)
        J = T[-n1   -sE; -sE   -sE2]

        δ = J \ T[-r1, -r2]

        α += δ[1]
        β += δ[2]
    end

    error("2-moment projection did not converge.")
end

"""
Find x in [0, g] such that sum(x)=mass_target and dot(E,x)=energy_target.

Smooth bounded parameterization:
    x_i = g_i * sigmoid(logit(y_i/g_i) + α + β * Ehat_i)

where Ehat_i = (E_i - Eref)/Escale.
"""
function smooth_two_moment_box(y::AbstractVector{T},g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T;Eref::T = sum(E) / length(E),Escale::T = max(maximum(abs.(E .- Eref)), one(T)),tol::T = T(1e-12),maxiter::Int = 50,eps::T = T(1e-14)) where {T<:Real}

    N = length(g)
    @assert length(y) == N && length(E) == N

    Ehat = (E .- Eref) ./ Escale

    # Initial logit state from y/g
    z0 = similar(g, T)
    @inbounds for i in 1:N
        if g[i] > zero(T)
            p = clamp(y[i] / g[i], eps, one(T) - eps)
            z0[i] = logit(p)
        else
            z0[i] = -T(100)   # effectively zero after sigmoid
        end
    end

    α = zero(T)
    β = zero(T)
    x = similar(g, T)

    for iter in 1:maxiter
        # Current x
        @inbounds for i in 1:N
            if g[i] == zero(T)
                x[i] = zero(T)
            else
                z = z0[i] + α + β * Ehat[i]
                x[i] = g[i] * sigmoid(z)
            end
        end

        r1 = sum(x) - mass_target
        r2 = dot(E, x) - energy_target

        if abs(r1) < tol && abs(r2) < tol
            return x
        end

        # Derivatives: dx/dz = x * (1 - x/g)
        w = similar(x, T)
        @inbounds for i in 1:N
            if g[i] == zero(T)
                w[i] = zero(T)
            else
                w[i] = x[i] * (one(T) - x[i] / g[i])
            end
        end

        # 2x2 Jacobian
        J11 = sum(w)
        J12 = sum(Ehat .* w)
        J21 = sum(E .* w)
        J22 = sum(E .* Ehat .* w)

        J = T[J11 J12; J21 J22]
        δ = J \ T[-r1, -r2]

        α += δ[1]
        β += δ[2]

        # Backtracking on the 2 residuals
        step = one(T)
        #=while step > T(1e-12)
            αt = α + step * δ[1]
            βt = β + step * δ[2]

            r1t = zero(T)
            r2t = zero(T)
            @inbounds for i in 1:N
                if g[i] == zero(T)
                    xi = zero(T)
                else
                    xi = g[i] * sigmoid(z0[i] + αt + βt * Ehat[i])
                end
                r1t += xi
                r2t += E[i] * xi
            end
            r1t -= mass_target
            r2t -= energy_target

            if norm(T[r1t, r2t]) <= (one(T) - T(1e-4) * step) * norm(T[r1, r2])
                α = αt
                β = βt
                break
            end

            step *= T(0.5)
        end=#

        if step <= T(1e-12)
            # No feasible split found for this channel
            println("r1: $(abs(r1)), r2: $(abs(r2)). No active set for projection. Returning zero split.")
            return zeros(T, length(g))
        end
    end

    error("smooth_two_moment_box: did not converge.")
end

function project_two_moment_exp_clamp(y::AbstractVector{T},g::AbstractVector{T},E::AbstractVector{T},mass_target::T,energy_target::T;tol::T = T(1e-12),maxiter::Int = 50) where {T<:Real}

    N = length(g)
    @assert length(y) == N && length(E) == N

    # Here E is already expected to be scaled if you pass Ehat.
    Ehat = E

    α = zero(T)
    β = zero(T)
    x = similar(g, T)
    r1 = zero(T)
    r2 = zero(T)

    J = zeros(T, 2, 2)

    for iter in 1:maxiter
        @inbounds for i in 1:N
            if g[i] <= zero(T) || y[i] <= zero(T)
                x[i] = zero(T)
            else
                z = α + β * Ehat[i]
                z = clamp(z, T(-40), T(40))   # avoid overflow
                xi = y[i] * exp(z)
                x[i] = clamp(xi, zero(T), g[i])
            end
        end

        r1 = sum(x) - mass_target
        r2 = dot(Ehat, x) - energy_target

        if abs(r1) < tol && abs(r2) < tol
            return x
        end

        active = (x .> T(0.0)) .& (x .< g)
        if !any(active)
            println("r1: $r1, r2: $r2. No active set for projection. Returning zero split. alpha: $α, beta: $β, cond(J): $(cond(J))")
            return nothing
        end

        w  = x[active]
        Eh = Ehat[active]

        J11 = sum(w)
        J12 = sum(Eh .* w)
        J21 = J12
        J22 = sum(Eh .* Eh .* w)

        J .= T[J11 J12;
              J21 J22]

        δ = J \ T[-r1, -r2]
        α += δ[1]
        β += δ[2]
    end

    println("r1: $r1, r2: $r2. Max iterations reached. Returning zero split. alpha: $α, beta: $β, cond(J): $(cond(J))")

    return nothing
end

# Numerically stable sigmoid
@inline function sigmoid(z::T) where {T<:Real}
    if z ≥ zero(T)
        return one(T) / (one(T) + exp(-z))
    else
        ez = exp(z)
        return ez / (one(T) + ez)
    end
end

@inline logit(p::T) where {T<:Real} = log(p) - log1p(-p)



function GainMatrix_to_M_BinPatankar!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},dpy3::Vector{Float64},dpz3::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = n_momentum
    #println("N = $N")
    #println("$offset2")

    is_sparse = isnothing(M_Bin)

    for px2 in 1:px2_num, px1 in 1:px1_num, px3 in 1:px3_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3) * sum(dpy1) * sum(dpy2) * sum(dpy3))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in 1:py1_num, pz1 in 1:pz1_num, py2 in 1:py2_num, pz2 in 1:pz2_num, py3 in 1:py3_num, pz3 in 1:pz3_num
                val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3] * dpy1[py1] * dpy2[py2] * dpy3[py3]
            end

        end

        for py1 in 1:py1_num, py2 in 1:py2_num, py3 in 1:py3_num

            if mode isa Axi

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3))

                # average over incoming and outgoing phi angles (pz)
                for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                    val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3]
                end

            end

            for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                
                if mode isa Ani

                    if is_sparse
                        #GainMax = maximum(@view(GainMatrix[:,py3,pz3,px1,py1,pz1,px2,py2,pz2]))
                        GainMax = maximum(@view(GainMatrix[:,:,:,px1,py1,pz1,px2,py2,pz2]))
                    end

                    val = GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                    #if is_sparse && val*w < eps(GainMax) # skips values smaller than this value to reduce memory usage
                    #    continue
                    #end

                end

                if val == 0.0
                    continue
                end

                a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
                b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                    if is_sparse
                    # M_Bin terms allocated symmetrically
                    if a != b
                        push!(M_Bin_I,(b-1)*N+(a-1)+1)
                        push!(M_Bin_J,c)
                        push!(M_Bin_V,convert(F,val*w/2))
                    end
                    if a != c
                        push!(M_Bin_I,(c-1)*N+(a-1)+1)
                        push!(M_Bin_J,b)
                        push!(M_Bin_V,convert(F,val*w/2))
                    end
                    # M_Bin terms allocated asymmetrically (save memory and avoids implicit coupling)
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,convert(F,val*w))
                else
                    # M_Bin terms allocated symmetrically
                    if a != b
                        M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                    end
                    if a != c
                        M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2)
                    end
                    # M_Bin terms allocated asymmetrically
                    #M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w)
                end

            end # pz loop

        end # py loop

    end # px loop

end



function LossMatrix_to_M_BinPatankar!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = n_momentum

    is_sparse = isnothing(M_Bin)

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrix[px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                end

                if val == 0.0
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                # M_Bin terms allocated symmetrically
                if is_sparse
                    # symmetric in b and c
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                    # symmetric in b and c
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/2 * (a==b)))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/2 * (a==c)))
                    # asymmetric in b and c (save memory and works better with MPE solver)
                    push!(M_Bin_I,a)
                    push!(M_Bin_J,c)
                    push!(M_Bin_V,-convert(F,val*w))
                    #push!(M_Bin_I,(c-1)*N+(c-1)+1)
                    #push!(M_Bin_J,a)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                else
                    # symmetric in b and c
                    #M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                    #M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
                    # asymmetric in b and c
                    M_Bin[a,c] -= convert(F,val*w)
                end

            end # pz loop

        end # py loop

    end # px loop

end

function LossMatrix_to_M_BinPatankarjk!(PhaseSpace::PhaseSpaceStruct,LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = n_momentum

        E = zeros(Float64,N)
        for species in eachindex(PhaseSpace.name_list)
            px_num = PhaseSpace.Momentum.px_num_list[species]
            py_num = PhaseSpace.Momentum.py_num_list[species]
            pz_num = PhaseSpace.Momentum.pz_num_list[species]
            dE = PhaseSpace.Grids.dE_list[species]
            for px in 1:px_num
                for py in 1:py_num
                    for pz in 1:pz_num
                        idx = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,species)
                        E[idx] = dE[px]
                    end
                end
            end
        end

    is_sparse = isnothing(M_Bin)

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrix[px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                end

                if val == 0.0
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                # M_Bin terms allocated symmetrically
                if is_sparse
                    # symmetric in b and c
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/4))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/4))
                    #push!(M_Bin_I,(b-1)*N+(c-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/4))
                    #push!(M_Bin_I,(c-1)*N+(c-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/4))

                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#))

                    push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    push!(M_Bin_J,c)
                    push!(M_Bin_V,-convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#))
                    push!(M_Bin_I,(c-1)*N+(c-1)+1)
                    push!(M_Bin_J,b)
                    push!(M_Bin_V,-convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#))
                else
                    # symmetric in b and c
                    M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#)
                    M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2 #=* E[a] / E[b] / E[c]=#)
                end

            end # pz loop

        end # py loop

    end # px loop

end

function LossMatrix_to_M_BinPatankarij!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = n_momentum

    is_sparse = isnothing(M_Bin)

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrix[px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                end

                if val == 0.0
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                # M_Bin terms as Likj where k == i
                if is_sparse
                    push!(M_Bin_I,(c-1)*N+(c-1)+1)
                    push!(M_Bin_J,b)
                    push!(M_Bin_V,-convert(F,val*w))
                else
                    if c == a
                    M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w)
                    end
                end

            end # pz loop

        end # py loop

    end # px loop

end

function LossMatrix_to_M_BinPatankarik!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = n_momentum

    is_sparse = isnothing(M_Bin)

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrix[px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                end

                if val == 0.0
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                b = a
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                # M_Bin terms as Lijk where j == i
                if is_sparse
                    # symmetric in b and c
                    push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    push!(M_Bin_J,c)
                    push!(M_Bin_V,-convert(F,val*w))
                else
                    # symmetric in b and c
                    if b == a
                    M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w)
                    end
                end

            end # pz loop

        end # py loop

    end # px loop

end

=#

#===== OLD CODE =======# 
#=
"""
    Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct)

Allocates a big matrix `M_Bin` which stores the interaction rates for binary interactions between all particles in the simulation. If `n` is the size of the momentum domain, then `M_Bin` is an `n^2 x n` matrix. The size of `M_Bin` in memory is printed to the console upon allocation.
"""
function Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType,Bin_sparse::Bool)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if loading_check

        if size > 1e9
            println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Bin will be approx. $size bytes in memory if dense")
        end

        #println("Do you want to proceed? (y/n)")
        #proceed = readline()
        #if proceed != "y"
        #    error("Error: User aborted binary interaction matrix loading")
        #end

    else

        println("Building binary interaction matrix")
        if size > 1e9
            println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Bin will be approx. $size bytes in memory if dense")
        end

    end

    if Bin_sparse
        M_Bin = spzeros(Precision,m,n)
    else 
        M_Bin = zeros(Precision,m,n)
    end

    return M_Bin

end

function GainMatrix_to_M_Bin_Axi!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpz1,dpz2,dpz3;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = size(M_Bin,2)

    is_sparse = isnothing(M_Bin)

    for py2 in axes(GainMatrix,8), px2 in axes(GainMatrix,7), py1 in axes(GainMatrix,5), px1 in axes(GainMatrix,4), py3 in axes(GainMatrix,2), px3 in 1:px3_num

        val = 0.0 # pz averaged GainMatrix term
        for pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6), pz3 in axes(GainMatrix,3)
            val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3]
        end

        for pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6), pz3 in axes(GainMatrix,3)

            a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            # weight for the size of the incoming and outgoing phase space
            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpz3)
        
            # M_Bin terms allocated symmetrically
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 
            end

        end

    end

end

function GainMatrix_to_M_Bin_Iso!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = size(M_Bin,2)

    for px2 in axes(GainMatrix,7), px1 in axes(GainMatrix,4), px3 in 1:px3_num

        val = 0.0 # py pz averaged GainMatrix term
        for py2 in axes(GainMatrix,8), pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6),py1 in axes(GainMatrix,5), pz3 in axes(GainMatrix,3), py3 in axes(GainMatrix,2)
            val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpy1[py1] * dpz2[pz2] * dpy2[py2] * dpz3[pz3] * dpy3[py3]
        end

        for py2 in axes(GainMatrix,8), pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6),py1 in axes(GainMatrix,5), pz3 in axes(GainMatrix,3), py3 in axes(GainMatrix,2)

            a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            # weight for the size of the incoming and outgoing phase space
            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpz3) / sum(dpy1) / sum(dpy2) / sum(dpy3)
        
            # M_Bin terms allocated symmetrically
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 
            end

        end

    end

end

function LossMatrix_to_M_Bin_Axi!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpz1,dpz2;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F <:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    is_sparse = isnothing(M_Bin)

    for py2 in axes(LossMatrix,5), px2 in axes(LossMatrix,4), py1 in axes(LossMatrix,2), px1 in axes(GainMatrix,1)

        val = 0.0
        for pz2 in axes(GainMatrix,6), pz1 in axes(GainMatrix,3)
            val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
        end

        for pz2 in axes(GainMatrix,6), pz1 in axes(GainMatrix,3)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            w = 1 / sum(dpz1) / sum(dpz2)
        
            # M_Bin terms allocated symmetrically
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,-convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,-convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
            end

        end

    end

end

function LossMatrix_to_M_Bin_Iso!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F <:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    is_sparse = isnothing(M_Bin)

    for px2 in axes(LossMatrix,4), px1 in axes(LossMatrix,1)

        val = 0.0
        for pz2 in axes(LossMatrix,6), py2 in axes(LossMatrix,5), pz1 in axes(LossMatrix,3),py1 in axes(LossMatrix,2)
            val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
        end

        for pz2 in axes(LossMatrix,6), py2 in axes(LossMatrix,5), pz1 in axes(LossMatrix,3),py1 in axes(LossMatrix,2)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpy1) / sum(dpy2)
        
            # M_Bin terms allocated symmetrically
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,-convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,-convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
            end

        end

    end

end

=#
