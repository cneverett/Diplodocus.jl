"""
    BinaryFileLoad_All(fileLocation,fileName;corrected=false)

Loads all the data stored in `fileName` stored at `fileLocation`. If `corrected` is set to `true` then the corrected gain and loss matrices will be loaded, otherwise the uncorrected matrices will be loaded.

# Example
```julia-repl
    (Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,GainTally3,GainTally4,LossTally) = BinaryFileLoad_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `GainTally3` : A 8D matrix of tallies of the number of emission spectrum values sampled for 12->34 interaction.
- `GainTally4` : A 8D matrix of tallies of the number of emission spectrum values sampled for 12->43 interaction.
- `LossTally` : A 6D matrix of tallies of the number of absorption spectrum values sampled.
- `GainMatrix3` : A 9D matrix of the emission spectrum for 12->34 interaction.
- `GainMatrix4` : A 9D matrix of the emission spectrum for 12->43 interaction.
- `LossMatrix1` : A 6D matrix of the absorption spectrum for 12->34 interaction.
- `LossMatrix2` : A 6D matrix of the absorption spectrum for 21->34 interaction i.e. by permutation of GainMatrix1 and correct application of phase space factors if species 1 != species 2.
"""
function BinaryFileLoad_All(fileLocation::String,fileName::String;corrected::Bool=false)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isdir(filePath)

    if fileExist
        f = zopen(filePath,"r");
        Parameters = f["Parameters"]
        if corrected
            GainMatrix3 = f["CorrectedGainMatrix3"];
            GainMatrix4 = f["CorrectedGainMatrix4"];
            LossMatrix1 = f["CorrectedLossMatrix1"];
            LossMatrix2 = f["CorrectedLossMatrix2"];
        else
            GainMatrix3 = f["GainMatrix3"];
            GainMatrix4 = f["GainMatrix4"];
            LossMatrix1 = f["LossMatrix1"];
            LossMatrix2 = f["LossMatrix2"]; 
        end
        GainWeights3 = f["GainWeights3"];
        GainWeights4 = f["GainWeights4"];
        LossTally = f["LossTally"];

        j = jldopen(filePath*"/data.jld2","r") 
        Parameters = j["Parameters"]
        close(j)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,GainWeights3,GainWeights4,LossTally)

end

"""
    BinaryFileLoad_Matrix(fileLocation,fileName;corrected=false)

Loads just the Gain and Loss Matrices stored in `fileName` stored at `fileLocation`. If `corrected` is set to `true` then the corrected gain and loss matrices will be loaded, otherwise the uncorrected matrices will be loaded.

# Example
```julia-repl
    (Parameters,GainMatrix3,GainMatrix4,LossMatrix) = BinaryFileLoad_Matrix(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `GainMatrix3` : A 9D matrix of the emission spectrum for 12->34 interaction.
- `GainMatrix4` : A 9D matrix of the emission spectrum for 12->43 interaction.
- `LossMatrix` : A 6D matrix of the absorption spectrum for 12->34 interaction.
"""
function BinaryFileLoad_Matrix(fileLocation::String,fileName::String;corrected::Bool=false)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isdir(filePath)

    if fileExist
        f = zopen(filePath,"r");
        if corrected
            GainMatrix3 = f["CorrectedGainMatrix3"];
            GainMatrix4 = f["CorrectedGainMatrix4"];
            LossMatrix = f["CorrectedLossMatrix"];
        else
            GainMatrix3 = f["GainMatrix3"];
            GainMatrix4 = f["GainMatrix4"];
            LossMatrix = f["LossMatrix"];
        end 
        j = jldopen(filePath*"/data.jld2","r") 
        Parameters = j["Parameters"]
        close(j) 
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3,GainMatrix4,LossMatrix)

end

"""
    fload_Matrix_ISO(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation` first converting them to an isotropic form by summing over angles. (The dimensions of the matrices stay the same i.e. 6D->6D with three dimensions having a size of 1)

# Example
```julia-repl
    Matrices = fload_All_ISO(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction.

If initial or final particles are identical then only one of the SMatrices or TMatrices will be returned for that state.

"""
function fload_Matrix_ISO(fileLocation::String,fileName::String)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        SMatrix3 = f["SMatrix3"];
        SMatrix4 = f["SMatrix4"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    if Parameters[1] == Parameters[2] && Parameters[3] == Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (Parameters,SMatrix3ISO,TMatrix1ISO)
    end

    if Parameters[1] == Parameters[2] && Parameters[3] != Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (Parameters,SMatrix3ISO,SMatrix4ISO,TMatrix1ISO)
    end

    if Parameters[1] != Parameters[2] && Parameters[3] == Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (Parameters,SMatrix3ISO,TMatrix1ISO,TMatrix2ISO)
    end

    if Parameters[1] != Parameters[2] && Parameters[3] != Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (SMatrix3ISO,SMatrix4ISO,TMatrix1ISO,TMatrix2ISO)
    end

end

"""
    DoesConserve(Output)

Function prints a set of measures to gauge the accuracy of the Monte-Carlo integration of the Gain and Loss terms. This includes the total number of particles gained and lost as well as the energy gained and lost and the ratio of these values. Argument `Output` is generated by either the `BinaryFileLoad_Matrix` or `BinaryFileLoad_All` functions. 
"""
function DoesConserve(Output::Tuple{Tuple,ZArray,ZArray,ZArray};Tuple_Output::Bool = false)

    # Output is tuple generated by BinaryFileLoad_Matrix

    println("Calculating Conservation Statistics for ZArrays")

    Parameters = Output[1] # always Output[1]

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    Indistinguishable_12::Bool = name1 == name2

    GainMatrix3 = Output[2]
    GainMatrix4 = Output[3]
    LossMatrix = Output[4]

    # underflow and overflow bins are taken to have size 0 -> p_r[1] and p_r[end] -> 2*p_r[end] respectively and are only needed for p3 and p4 the outgoing particles

    p1_r::Vector{Float64} = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d::Vector{Float64} = deltaVector(p1_r);
    #p1_d_full = [p1_r[1] ; p1_d; deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ::Vector{Float64} = deltaEVector(p1_r,mu1);
    #E1_Δ_full = [E1_Δ; deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    #E1_d_full = E1_Δ_full ./ p1_d_full;
    E1_d = E1_Δ ./ p1_d

    p2_r::Vector{Float64} = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d::Vector{Float64} = deltaVector(p2_r);
    #p2_d_full = [p2_d; deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ::Vector{Float64} = deltaEVector(p2_r,mu2);
    #E2_Δ_full = [E2_Δ; deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    #E2_d_full = E2_Δ_full ./ p2_d_full;
    E2_d = E2_Δ ./ p2_d

    p3_r::Vector{Float64} = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d::Vector{Float64} = deltaVector(p3_r);
    #p3_d_full = [p3_d; deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ::Vector{Float64} = deltaEVector(p3_r,mu3);
    #E3_Δ_full = [E3_Δ; deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    #E3_d_full = E3_Δ_full ./ p3_d_full
    E3_d = E3_Δ ./ p3_d

    p4_r::Vector{Float64} = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d::Vector{Float64} = deltaVector(p4_r);
    #p4_d_full = [p4_d; deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ::Vector{Float64} = deltaEVector(p4_r,mu4);
    #E4_Δ_full = [E4_Δ; deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    #E4_d_full = E4_Δ_full ./ p4_d_full
    E4_d = E4_Δ ./ p4_d

    #println("E1_d: $E1_d \n E2_d: $E2_d \n E3_d: $E3_d \n E4_d: $E4_d")

    SsumN3::Float64 = 0.0
    TsumN1::Float64 = 0.0
    SsumE3::Float64 = 0.0
    TsumE1::Float64 = 0.0

    SsumN4::Float64 = 0.0
    TsumN2::Float64 = 0.0
    SsumE4::Float64 = 0.0
    TsumE2::Float64 = 0.0

    NGainMatrix3::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    NLossMatrix1::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    EGainMatrix3::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    ELossMatrix1::Array{Float64,6} = zeros(Float64,size(LossMatrix))

    NGainMatrix4::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    NLossMatrix2::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    EGainMatrix4::Array{Float64,6} = zeros(Float64,size(LossMatrix))
    ELossMatrix2::Array{Float64,6} = zeros(Float64,size(LossMatrix))

    ChunkGainMatrix3Full::Array{Float64,9} = zeros(Float64,p3_num+2,u3_num,h3_num,1,u1_num,h1_num,1,u2_num,h2_num)
    ChunkGainMatrix4Full::Array{Float64,9} = zeros(Float64,p4_num+2,u4_num,h4_num,1,u1_num,h1_num,1,u2_num,h2_num)
    ChunkLossMatrixFull::Array{Float64,6} = zeros(Float64,1,u1_num,h1_num,1,u2_num,h2_num)

    ChunkGainMatrix3 = @view(ChunkGainMatrix3Full[:,:,:,1,:,:,1,:,:])
    ChunkGainMatrix4 = @view(ChunkGainMatrix4Full[:,:,:,1,:,:,1,:,:])
    ChunkLossMatrix = @view(ChunkLossMatrixFull[1,:,:,1,:,:])

    NErrMatrix = zeros(Float64,size(LossMatrix))
    EErrMatrix = zeros(Float64,size(LossMatrix))

    #=@inbounds=# for p1 in axes(GainMatrix3, 4), p2 in axes(GainMatrix3,7)

        r3 = CartesianIndices((1:(p3_num+2), 1:u3_num, 1:h3_num,p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))
        r4 = CartesianIndices((1:(p4_num+2), 1:u4_num, 1:h4_num,p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))
        r = CartesianIndices((p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))

        Zarr.ZarrCore.readblock!(ChunkGainMatrix3Full,GainMatrix3,r3)
        Zarr.ZarrCore.readblock!(ChunkGainMatrix4Full,GainMatrix4,r4)
        Zarr.ZarrCore.readblock!(ChunkLossMatrixFull,LossMatrix,r)

        #=@inbounds=# for u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
            tmp = ChunkGainMatrix3[p3,u3,h3,u1,h1,u2,h2]
            SsumN3 += tmp
            SsumE3 += tmp*E3_d[p3]
            NGainMatrix3[p1,u1,h1,p2,u2,h2] += tmp
            EGainMatrix3[p1,u1,h1,p2,u2,h2] += tmp*E3_d[p3]
            end
        end

        #=@inbounds=# for u1 in axes(GainMatrix4,5), h1 in axes(GainMatrix4,6), u2 in axes(GainMatrix4,8), h2 in axes(GainMatrix4,9)
            for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
            tmp = ChunkGainMatrix4[p4,u4,h4,u1,h1,u2,h2]
            SsumN4 += tmp
            SsumE4 += tmp*E4_d[p4]
            NGainMatrix4[p1,u1,h1,p2,u2,h2] += tmp
            EGainMatrix4[p1,u1,h1,p2,u2,h2] += tmp*E4_d[p4]
            end
        end

        #=@inbounds=# for u1 in axes(LossMatrix, 2), h1 in axes(LossMatrix,3), u2 in axes(LossMatrix,5), h2 in axes(LossMatrix,6)
            tmp1 = ChunkLossMatrix[u1,h1,u2,h2]
            if !Indistinguishable_12
                tmp2 = ChunkLossMatrix[u1,h1,u2,h2]
            else 
                tmp2 = 0.0
            end
            TsumN1 += tmp1
            TsumN2 += tmp2
            NLossMatrix1[p1,u1,h1,p2,u2,h2] += tmp1
            NLossMatrix2[p1,u1,h1,p2,u2,h2] += tmp2
            if !Indistinguishable_12 # so LossMatrix2 is non-zero
                TsumE1 += tmp1*E1_d[p1]
                TsumE2 += tmp2*E2_d[p2]
                ELossMatrix1[p1,u1,h1,p2,u2,h2] += tmp1*E1_d[p1]
                ELossMatrix2[p1,u1,h1,p2,u2,h2] += tmp2*E2_d[p2]
            else
                TsumE1 += tmp1*(E1_d[p1]+E2_d[p2])/2
                TsumE2 += 0.0
                ELossMatrix1[p1,u1,h1,p2,u2,h2] += tmp1*(E1_d[p1]+E2_d[p2])/2
                ELossMatrix2[p1,u1,h1,p2,u2,h2] += 0.0#LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            end
        end

    end


    NErrMatrix = (NGainMatrix3 .+ NGainMatrix4 .- NLossMatrix1 .- NLossMatrix2) ./ (NLossMatrix1 .+ NLossMatrix2)
    NErrList = filter(isfinite, NErrMatrix)
    println("sum gain 3 Matrices:",sum(NGainMatrix3))
    println("sum gain 4 Matrices:",sum(NGainMatrix4))
    println("sum loss 1 Matrices:",sum(NLossMatrix1))
    println("sum loss 2 Matrices:",sum(NLossMatrix2))
    println("sum gain Matrices:",sum(NGainMatrix3 .+ NGainMatrix4))
    println("sum gain-loss Matrices:",sum(NGainMatrix3 .+ NGainMatrix4 .- NLossMatrix1 .- NLossMatrix2))
    println("sum loss Matrices:",sum(NLossMatrix1 .+ NLossMatrix2))
    println("sum NErrMatrix:",sum(NErrMatrix))
    println("sum NErrList:",sum(NErrList))
    println("length NErrList:",length(NErrList))
    meanNErr = sum(NErrList) / length(NErrList)
    stdN = sqrt(sum((NErrList .- meanNErr).^2)/length(NErrList))

    EErrMatrix = (EGainMatrix3 .+ EGainMatrix4 .- ELossMatrix1 .- ELossMatrix2) ./ (ELossMatrix1 .+ ELossMatrix2)
    EErrList = filter(isfinite, EErrMatrix)
    meanEErr = sum(EErrList) / length(EErrList)
    stdE = sqrt(sum((EErrList .- meanEErr).^2)/length(EErrList))

    println("sumGainN3 = "*string(SsumN3))
    println("sumGainN4 = "*string(SsumN4))
    println("sumLossN1 = "*string(TsumN1))    
    println("sumLossN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumGainN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumLossN = "*string(TsumN))

    println("#")

    println("sumGainE3 = "*string(SsumE3))
    println("sumGainE4 = "*string(SsumE4))
    println("sumLossE1 = "*string(TsumE1))  
    println("sumLossE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumGainE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumLossE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string(SsumE-TsumE))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

    println("#")
    println("mean error in N = $meanNErr")
    println("std of error in  N = $stdN")
    println("#")
    println("mean error in E = $meanEErr")
    println("std of error in E = $stdE")

    if Tuple_Output == true
        return NGainMatrix3, NLossMatrix1, NErrMatrix,NGainMatrix4, NLossMatrix2
    else
        return nothing
    end
    

end


####### LEGACY CODE BELOW, NOT USED IN CURRENT VERSION OF THE PACKAGE, KEPT FOR REFERENCE #####

#=

function readblock_singlechunk!(aout::AbstractArray{<:Any,N}, z::ZArray{<:Any, N}, r::CartesianIndex{N}) where {N}

    # Fast path: single-chunk full-read decodes directly into `aout`, skipping the readtask channel and scratch buffer.
    chunk_compressed = Zarr.store_readchunk(z.storage, z.path, r, z.metadata.chunk_key_encoding)
    Zarr.uncompress_raw!(aout, z, chunk_compressed)

end


function DoesConserve(Output::Tuple{Tuple,Array{Float64,9},Array{Float64,9},Array{Float64,6},Array{Float64,6}};Tuple_Output::Bool = false)

    # Output is tuple generated by BinaryFileLoad_Matrix

    Parameters = Output[1] # always Output[1]

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    GainMatrix3 = Output[2]
    GainMatrix4 = Output[3]
    LossMatrix1 = Output[4]
    LossMatrix2 = Output[5]

    # underflow and overflow bins are taken to have size 0 -> p_r[1] and p_r[end] -> 2*p_r[end] respectively and are only needed for p3 and p4 the outgoing particles

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    #p1_d_full = [p1_r[1] ; p1_d; deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ = deltaEVector(p1_r,mu1);
    #E1_Δ_full = [E1_Δ; deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    #E1_d_full = E1_Δ_full ./ p1_d_full;
    E1_d = E1_Δ ./ p1_d

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    #p2_d_full = [p2_d; deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ = deltaEVector(p2_r,mu2);
    #E2_Δ_full = [E2_Δ; deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    #E2_d_full = E2_Δ_full ./ p2_d_full;
    E2_d = E2_Δ ./ p2_d

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    #p3_d_full = [p3_d; deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ = deltaEVector(p3_r,mu3);
    #E3_Δ_full = [E3_Δ; deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    #E3_d_full = E3_Δ_full ./ p3_d_full
    E3_d = E3_Δ ./ p3_d

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    #p4_d_full = [p4_d; deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ = deltaEVector(p4_r,mu4);
    #E4_Δ_full = [E4_Δ; deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    #E4_d_full = E4_Δ_full ./ p4_d_full
    E4_d = E4_Δ ./ p4_d

    #println("E1_d: $E1_d \n E2_d: $E2_d \n E3_d: $E3_d \n E4_d: $E4_d")

    SsumN3 = 0
    TsumN1 = 0
    SsumE3 = 0
    TsumE1 = 0

    SsumN4 = 0
    TsumN2 = 0
    SsumE4 = 0
    TsumE2 = 0

    NGainMatrix3 = zeros(Float64,size(LossMatrix1))
    NLossMatrix1 = zeros(Float64,size(LossMatrix1))
    EGainMatrix3 = zeros(Float64,size(LossMatrix1))
    ELossMatrix1 = zeros(Float64,size(LossMatrix1))

    NGainMatrix4 = zeros(Float64,size(LossMatrix1))
    NLossMatrix2 = zeros(Float64,size(LossMatrix1))
    EGainMatrix4 = zeros(Float64,size(LossMatrix1))
    ELossMatrix2 = zeros(Float64,size(LossMatrix1))

    @inbounds for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
        SsumN3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        SsumE3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
        NGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        EGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
        end
    end

    @inbounds for p1 in axes(GainMatrix4, 4), u1 in axes(GainMatrix4,5), h1 in axes(GainMatrix4,6), p2 in axes(GainMatrix4,7), u2 in axes(GainMatrix4,8), h2 in axes(GainMatrix4,9)
        for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
        SsumN4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        SsumE4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
        NGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        EGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
        end
    end

    @inbounds for p1 in axes(LossMatrix1,1), u1 in axes(LossMatrix1, 2), h1 in axes(LossMatrix1,3), p2 in axes(LossMatrix1,4), u2 in axes(LossMatrix1,5), h2 in axes(LossMatrix1,6)
        TsumN1 += LossMatrix1[p1,u1,h1,p2,u2,h2]
        TsumN2 += LossMatrix2[p2,u2,h2,p1,u1,h1]
        NLossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]
        NLossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]
        if name1 != name2 # Indistinguishable_12 == false so LossMatrix2 is non-zero
            TsumE1 += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
            TsumE2 += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
            ELossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
            ELossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
        else
            TsumE1 += LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            TsumE2 += 0.0
            ELossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            ELossMatrix2[p1,u1,h1,p2,u2,h2] += 0.0#LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
        end
    end

    NErrMatrix = (NGainMatrix3 .+ NGainMatrix4 .- NLossMatrix1 .- NLossMatrix2) ./ (NLossMatrix1 .+ NLossMatrix2)
    NErrList = filter(!isnan, NErrMatrix)
    #meanNErr = sum(abs.(NErrMatrix)) / length(NErrMatrix)
    meanNErr = sum(NErrList) / length(NErrList)
    #stdN = sqrt(sum((NErrMatrix .- meanNErr).^2)/length(NLossMatrix1))
    stdN = sqrt(sum((NErrList .- meanNErr).^2)/length(NErrList))

    EErrMatrix = (EGainMatrix3 .+ EGainMatrix4 .- ELossMatrix1 .- ELossMatrix2) ./ (ELossMatrix1 .+ ELossMatrix2)
    EErrList = filter(!isnan, EErrMatrix)
    #meanEErr = sum(abs.(EErrMatrix)) / length(NErrMatrix)
    meanEErr = sum(EErrList) / length(EErrList)
    #stdE = sqrt(sum((EErrMatrix .- meanEErr).^2)/length(NLossMatrix1))
    stdE = sqrt(sum((EErrList .- meanEErr).^2)/length(EErrList))

    println("sumGainN3 = "*string(SsumN3))
    println("sumGainN4 = "*string(SsumN4))
    println("sumLossN1 = "*string(TsumN1))    
    println("sumLossN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumGainN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumLossN = "*string(TsumN))

    println("#")

    println("sumGainE3 = "*string(SsumE3))
    println("sumGainE4 = "*string(SsumE4))
    println("sumLossE1 = "*string(TsumE1))  
    println("sumLossE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumGainE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumLossE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string(SsumE-TsumE))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

    println("#")
    println("mean error in N = $meanNErr")
    println("std of error in  N = $stdN")
    #println("mean error in E = $meanEErr1")
    #println("std of error in E = $stdE1")
    #println("#")
    #println("mean error in N2 = $meanNErr2")
    #println("std of error in  N2 = $stdN2")
    #println("mean error in E = $meanEErr2")
    #println("std of error in E = $stdE2")
    println("#")
    println("mean error in E = $meanEErr")
    println("std of error in E = $stdE")

    if Tuple_Output == true
        return NGainMatrix3, NLossMatrix1, NErrMatrix,NGainMatrix4, NLossMatrix2
    else
        return nothing
    end
    

end

function DoesConserve(Output::Tuple{Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, Array{Float64, 9}, Array{Float64, 9}, Array{Float64, 6}, Array{Float64, 6}, Array{UInt32, 9}, Array{UInt32, 9}, Array{UInt32, 6}};Tuple_Output::Bool = false)

    # Output is tuple generated by BinaryFileLoad_All

    Parameters = Output[1] # always Output[1]

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    GainMatrix3 = Output[2]
    GainMatrix4 = Output[3]
    LossMatrix1 = Output[4]
    LossMatrix2 = Output[5]
    GainSamples3 = Output[6]
    GainSamples4 = Output[7]
    LossSamples1 = Output[8]
    perm = [4,5,6,1,2,3]
    LossSamples2 = permutedims(Output[8],perm)

    # underflow and overflow bins are taken to have size 0 -> p_r[1] and p_r[end] -> 2*p_r[end] respectively and are only needed for p3 and p4 the outgoing particles

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    #p1_d_full = [p1_d; deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ = deltaEVector(p1_r,mu1);
    #E1_Δ_full = [E1_Δ; deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    #E1_d_full = E1_Δ_full ./ p1_d_full;
    E1_d = E1_Δ ./ p1_d

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    #p2_d_full = [p2_d; deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ = deltaEVector(p2_r,mu2);
    #E2_Δ_full = [E2_Δ; deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    #E2_d_full = E2_Δ_full ./ p2_d_full;
    E2_d = E2_Δ ./ p2_d


    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    #p3_d_full = [p3_d; deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ = deltaEVector(p3_r,mu3);
    #E3_Δ_full = [E3_Δ; deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    #E3_d_full = E3_Δ_full ./ p3_d_full
    E3_d = E3_Δ ./ p3_d

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    #p4_d_full = [p4_d; deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ = deltaEVector(p4_r,mu4);
    #E4_Δ_full = [E4_Δ; deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    #E4_d_full = E4_Δ_full ./ p4_d_full
    E4_d = E4_Δ ./ p4_d

    #println("E1_d: $E1_d, E2_d: $E2_d, E3_d: $E3_d, E4_d: $E4_d")

    SsumN3 = 0
    TsumN1 = 0
    SsumE3 = 0
    TsumE1 = 0

    SsumN4 = 0
    TsumN2 = 0
    SsumE4 = 0
    TsumE2 = 0

    NGainMatrix3 = zeros(Float64,size(LossMatrix1))
    NLossMatrix1 = zeros(Float64,size(LossMatrix1))
    EGainMatrix3 = zeros(Float64,size(LossMatrix1))
    ELossMatrix1 = zeros(Float64,size(LossMatrix1))

    NGainMatrix4 = zeros(Float64,size(LossMatrix1))
    NLossMatrix2 = zeros(Float64,size(LossMatrix1))
    EGainMatrix4 = zeros(Float64,size(LossMatrix1))
    ELossMatrix2 = zeros(Float64,size(LossMatrix1))

    @inbounds for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
        SsumN3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        SsumE3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
        NGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        EGainMatrix3[p1,u1,h1,p2,u2,h2] += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
        end
    end

    @inbounds for p1 in axes(GainMatrix4, 4), u1 in axes(GainMatrix4,5), h1 in axes(GainMatrix4,6), p2 in axes(GainMatrix4,7), u2 in axes(GainMatrix4,8), h2 in axes(GainMatrix4,9)
        for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
        SsumN4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        SsumE4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
        NGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        EGainMatrix4[p1,u1,h1,p2,u2,h2] += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
        end
    end

    @inbounds for p1 in axes(LossMatrix1,1), u1 in axes(LossMatrix1, 2), h1 in axes(LossMatrix1,3), p2 in axes(LossMatrix1,4), u2 in axes(LossMatrix1,5), h2 in axes(LossMatrix1,6)
        TsumN1 += LossMatrix1[p1,u1,h1,p2,u2,h2]
        TsumE1 += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
        TsumN2 += LossMatrix2[p2,u2,h2,p1,u1,h1]
        TsumE2 += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
        NLossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]
        NLossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]
        ELossMatrix1[p1,u1,h1,p2,u2,h2] += LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
        ELossMatrix2[p1,u1,h1,p2,u2,h2] += LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
    end

    NErrMatrix = (NGainMatrix3 .+ NGainMatrix4 .- NLossMatrix1 .- NLossMatrix2) ./ (NLossMatrix1 .+ NLossMatrix2)
    #EErrMatrix1 = (EGainMatrix3 .- ELossMatrix1) ./ ELossMatrix1
    meanNErr = sum(abs.(NErrMatrix)) / length(NLossMatrix1)
    #meanEErr1 = sum(abs.(EErrMatrix1)) / length(NLossMatrix1)
    stdN = sqrt(sum((NErrMatrix .- meanNErr).^2)/length(NLossMatrix1))
    #stdE1 = sqrt(sum((EErrMatrix1 .- meanEErr1).^2)/length(NLossMatrix1))

    #NErrMatrix2 = (NGainMatrix4 .- NLossMatrix2) ./ NLossMatrix2
    #EErrMatrix2 = (EGainMatrix4 .- ELossMatrix2) ./ ELossMatrix2
    #meanNErr2 = sum(abs.(NErrMatrix2)) / length(NLossMatrix2)
    #meanEErr2 = sum(abs.(EErrMatrix2)) / length(NLossMatrix2)
    #stdN2 = sqrt(sum((NErrMatrix2 .- meanNErr2).^2)/length(NLossMatrix2))
    #stdE2 = sqrt(sum((EErrMatrix2 .- meanEErr2).^2)/length(NLossMatrix2))

    EErrMatrix = (EGainMatrix3 .+ EGainMatrix4 .- ELossMatrix1 .- ELossMatrix2) ./ (ELossMatrix1 .+ ELossMatrix2)
    meanEErr = sum(abs.(EErrMatrix)) / length(NLossMatrix1)
    stdE = sqrt(sum((EErrMatrix .- meanEErr).^2)/length(NLossMatrix1))

    println("sumGainN3 = "*string(SsumN3))
    println("sumGainN4 = "*string(SsumN4))
    println("sumLossN1 = "*string(TsumN1))    
    println("sumLossN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumGainN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumLossN = "*string(TsumN))

    println("#")

    println("sumGainE3 = "*string(SsumE3))
    println("sumGainE4 = "*string(SsumE4))
    println("sumLossE1 = "*string(TsumE1))  
    println("sumLossE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumGainE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumLossE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string(SsumE-TsumE))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

    println("#")
    println("mean error in N = $meanNErr")
    println("std of error in  N = $stdN")
    #println("mean error in E = $meanEErr1")
    #println("std of error in E = $stdE1")
    #println("#")
    #println("mean error in N2 = $meanNErr2")
    #println("std of error in  N2 = $stdN2")
    #println("mean error in E = $meanEErr2")
    #println("std of error in E = $stdE2")
    println("#")
    println("mean error in E = $meanEErr")
    println("std of error in E = $stdE")

    if Tuple_Output == true
        return NGainMatrix3, NLossMatrix1, NErrMatrix,NGainMatrix4, NLossMatrix2, GainSamples3, GainSamples4, LossSamples1, LossSamples2
    else
        return nothing
    end 
    
end

=#