function readblock_singlechunk!(aout::AbstractArray{<:Any,N}, z::ZArray{<:Any, N}, r::CartesianIndex{N}) where {N}

    # Fast path: single-chunk full-read decodes directly into `aout`, skipping the readtask channel and scratch buffer.
    chunk_compressed = Zarr.store_readchunk(z.storage, z.path, r, z.metadata.chunk_key_encoding)
    Zarr.uncompress_raw!(aout, z, chunk_compressed)

end


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
    meanNErr_sum = 0.0
    meanNErr2_sum = 0.0
    meanNErr_count = 0
    meanEErr_sum = 0.0
    meanEErr2_sum = 0.0
    meanEErr_count = 0

    #=@inbounds=# for p1 in axes(GainMatrix3, 4), p2 in axes(GainMatrix3,7)

        if p1 == 1 && p2 == 24 
        #    continue 
        end

        #gain3loc = CartesianIndex(1,1,1,p1,1,1,p2,1,1)
        #gain4loc = CartesianIndex(1,1,1,p1,1,1,p2,1,1)
        #lossloc = CartesianIndex(p1,1,1,p2,1,1)

        r3 = CartesianIndices((1:(p3_num+2), 1:u3_num, 1:h3_num,p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))
        r4 = CartesianIndices((1:(p4_num+2), 1:u4_num, 1:h4_num,p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))
        r = CartesianIndices((p1:p1, 1:u1_num, 1:h1_num,p2:p2, 1:u2_num, 1:h2_num))

        Zarr.ZarrCore.readblock!(ChunkGainMatrix3Full,GainMatrix3,r3)
        Zarr.ZarrCore.readblock!(ChunkGainMatrix4Full,GainMatrix4,r4)
        Zarr.ZarrCore.readblock!(ChunkLossMatrixFull,LossMatrix,r)

        # does not work, do not know why
        #readblock_singlechunk!(ChunkGainMatrix3Full,GainMatrix3,gain3loc)
        #readblock_singlechunk!(ChunkGainMatrix4Full,GainMatrix4,gain4loc)
        #readblock_singlechunk!(ChunkLossMatrixFull,LossMatrix,lossloc)

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

## TODO: remove old correction functions

#=
"""
    GainCorrection(Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

MC sampling introduces noise that can lead to poor number and energy conservation. `GainCorrection` provides a corrective step to ensure number and energy conservation to numerical precision. If there is no GainMatrix element then the value of the LossMatrix is applied to the same bin as the input state (if they are identical particles); if not identical particles if there is no GainMatrix element then the value of the LossMatrix is set to zero to ensure particle conservation (with good MC sampling this should rarely occur).
"""
function GainCorrection(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, GainMatrix3::Array{Float64, 9}, GainMatrix4::Array{Float64, 9}, LossMatrix1::Array{Float64, 6}, LossMatrix2::Array{Float64, 6})

    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    CorrectedGainMatrix3 = similar(GainMatrix3)
    CorrectedGainMatrix4 = similar(GainMatrix4)
    CorrectedLossMatrix1 = similar(LossMatrix1)
    CorrectedLossMatrix2 = similar(LossMatrix2)
    fill!(CorrectedGainMatrix3,Float64(0))
    fill!(CorrectedGainMatrix4,Float64(0))
    CorrectedLossMatrix1 .= LossMatrix1
    CorrectedLossMatrix2 .= LossMatrix2

    # underflow and overflow bins are taken to have size 0 -> p1_r[1] and p1_r[end] -> 2*p1_r[end] respectively

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    u1_r = bounds(u_low,u_up,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    u2_r = bounds(u_low,u_up,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    u3_r = bounds(u_low,u_up,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    u4_r = bounds(u_low,u_up,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);

    for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        
        GainSumN3 = zero(Float64)
        GainSumN4 = zero(Float64)
        LossSumN1 = LossMatrix1[p1,u1,h1,p2,u2,h2]
        LossSumN2 = LossMatrix2[p2,u2,h2,p1,u1,h1]
        
        for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
            GainSumN3 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
        end
        for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
            GainSumN4 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
        end

        if name1 == name3 # particle 1 and particle 3 are identical
            if GainSumN3 != 0e0
                Correction = LossSumN1/GainSumN3
                @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])
            else
                CorrectedGainMatrix3[p1,u1,h1,p1,u1,h1,p2,u2,h2] += LossSumN1
            end
        end
        if name2 == name4 # particle 2 and 4 are identical
            if GainSumN4 != 0e0
                Correction = LossSumN2/GainSumN4
                @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])
            else
                CorrectedGainMatrix4[p2,u2,h2,p1,u1,h1,p2,u2,h2] += LossSumN2
            end
        end    
        if (name1 != name3) && (name2 != name4) && m1 == m2 && m3 == m4 # incoming and outgoing states are symmetric but not the same particles
            if (GainSumN3 + GainSumN4) != 0e0
                Correction = (LossSumN1+LossSumN2)/(GainSumN3+GainSumN4)
                #println(GainSumN3)
                #println(GainSumN4)
                #println(LossSumN1)
                #println(LossSumN2)
                #println(Correction)
                @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])
                @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])
            elseif (LossSumN1+LossSumN2) != 0e0
                # no gain term but there is a loss term, as outgoing states are not identical to incoming there is no way to correct for this so have to set loss terms to 0.0
                CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0
            end
        end

    end

    return CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2

end


function GainCorrection2(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, GainMatrix3::Array{Float64, 9}, GainMatrix4::Array{Float64, 9}, LossMatrix1::Array{Float64, 6}, LossMatrix2::Array{Float64, 6})

    # Function that applies the correct phase space factors to GainMatrix and LossMatrix derived from Stotal and Ttotal arrays

    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    CorrectedGainMatrix3 = similar(GainMatrix3)
    CorrectedGainMatrix4 = similar(GainMatrix4)
    CorrectedLossMatrix1 = similar(LossMatrix1)
    CorrectedLossMatrix2 = similar(LossMatrix2)
    fill!(CorrectedGainMatrix3,Float64(0))
    fill!(CorrectedGainMatrix4,Float64(0))
    CorrectedLossMatrix1 .= LossMatrix1
    CorrectedLossMatrix2 .= LossMatrix2

    # underflow and overflow bins are taken to have size 0 -> p1_r[1] and p1_r[end] -> 2*p1_r[end] respectively

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    u1_r = bounds(u_low,u_up,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    u2_r = bounds(u_low,u_up,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    u3_r = bounds(u_low,u_up,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    u4_r = bounds(u_low,u_up,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);

    
    E1_Δ = deltaEVector(p1_r,m1);
    E1_d = E1_Δ ./ p1_d

    E2_Δ = deltaEVector(p2_r,m2);
    E2_d = E2_Δ ./ p2_d

    E3_Δ = deltaEVector(p3_r,m3);
    E3_d = E3_Δ ./ p3_d

    E4_Δ = deltaEVector(p4_r,m4);
    E4_d = E4_Δ ./ p4_d

    num_wrong = 0
    num_right = 0

    for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        
        LossSumN1 = LossMatrix1[p1,u1,h1,p2,u2,h2]
        LossSumN2 = LossMatrix2[p2,u2,h2,p1,u1,h1]
        if name1 != name2 # Indistinguishable_12==false so LossMatrix2 is non-zero
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
            LossSumE2 = LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
        else
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            LossSumE2 = 0.0#LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
        end

        
        wrong = true
        nonzero_gain = true
        max_high_bins = 1

        while wrong && nonzero_gain

            GainSumN31 = zero(Float64)
            GainSumN32 = zero(Float64)
            GainSumN41 = zero(Float64)
            GainSumN42 = zero(Float64)
            GainSumE31 = zero(Float64)
            GainSumE32 = zero(Float64)
            GainSumE41 = zero(Float64)
            GainSumE42 = zero(Float64)

            numbeta = 0
            
            for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                high_bins = 0
                for p3 in reverse(axes(GainMatrix3,1))
                    if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] != 0e0 && high_bins < max_high_bins # number of non-zero element starting from high p3 values
                        GainSumN32 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
                        GainSumE32 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
                        high_bins += 1
                        numbeta +=1
                    else
                        GainSumN31 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
                        GainSumE31 += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
                    end
                end
            end
            for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                high_bins = 0
                for p4 in reverse(axes(GainMatrix4,1))
                    if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] != 0e0 && high_bins < max_high_bins # number of non-zero element starting from high p3 values
                        GainSumN42 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
                        GainSumE42 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
                        high_bins += 1
                        numbeta +=1
                    else
                        GainSumN41 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
                        GainSumE41 += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
                    end
                end
            end

            if (name1==name2==name3==name4) || (name1!=name2!=name3!=name4) # all particles identical or all particles different
                GainN1 = GainSumN31 + GainSumN41
                GainE1 = GainSumE31 + GainSumE41
                GainN2 = GainSumN32 + GainSumN42
                GainE2 = GainSumE32 + GainSumE42
                LossN = LossSumN1 + LossSumN2
                LossE = LossSumE1 + LossSumE2

                a = GainN1
                if a == 0e0 # loop has not produced any gain terms
                    nonzero_gain = false
                    #println("No gain terms for p1=$p1,p2=$p2")
                    num_wrong += 1
                    continue
                end
                ac = GainE1
                c = GainE1/GainN1

                b = GainN2
                bd = GainE2
                d = GainE2/GainN2

                e = a+b-LossN
                f = ac+bd-LossE

                alpha = (e*d-f)/(ac-a*d)+1
                beta = (e*c-f)/(bd-b*c)+1

                if alpha < 0e0 || beta < 0e0
                    #println("LossN = $LossN, LossE = $LossE, LossN1 = $LossSumN1, LossN2 = $LossSumN2, LossE1 = $LossSumE1, LossE2 = $LossSumE2")
                    #println("alpha = $alpha, beta = $beta, a+b-L = f = $(a+b-LossN), ac+bd-LE = e = $(ac+bd-LossE), ac = $ac, bd = $bd, c = $c, d = $d, a = $a, b = $b")
                    #println("Old method N Correct = $(LossN/(GainN1+GainN2))")
                    #println("Old method E Correct = $((LossE)/(GainE1+GainE2))")
                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins")
                    wrong = false
                    num_right += 1
                end

                for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                    high_bins = 0
                    for p3 in reverse(axes(GainMatrix3,1))
                        if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] != 0e0 && high_bins < max_high_bins # number of non-zero element starting from high p3 values
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] =  GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] * beta
                            high_bins += 1
                        else
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] * alpha
                        end
                    end
                end

                for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                    high_bins = 0
                    for p4 in reverse(axes(GainMatrix4,1))
                        if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] != 0e0 && high_bins < max_high_bins # number of non-zero element starting from high p3 values
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] * beta
                            high_bins += 1
                        else
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] * alpha
                        end
                    end
                end

                #println("N: $(sum(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) + sum(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) vs $(LossMatrix1[p1,u1,h1,p2,u2,h2] + LossMatrix2[p2,u2,h2,p1,u1,h1])")

                #println("E: $((GainSumE31+GainSumE41)*alpha + (GainSumE32+GainSumE42+beta*((GainSumE32+GainSumE42)/(GainSumN32+GainSumN42)))) vs $LossE")
                #println("N: $((GainSumN31+GainSumN41)*alpha + (GainSumN32+GainSumN42+beta)) vs $LossN")

            end

        end # while

        #=if name1 == name3 # particle 1 and particle 3 are identical
            if GainSumN3 != 0e0
                Correction = LossSumN1/GainSumN3
                @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])
            else
                CorrectedGainMatrix3[p1,u1,h1,p1,u1,h1,p2,u2,h2] += LossSumN1
            end
        end
        if name2 == name4 # particle 2 and 4 are identical
            if GainSumN4 != 0e0
                Correction = LossSumN2/GainSumN4
                @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])
            else
                CorrectedGainMatrix4[p2,u2,h2,p1,u1,h1,p2,u2,h2] += LossSumN2
            end
        end    
        if (name1 != name3) && (name2 != name4) && m1 == m2 && m3 == m4 # incoming and outgoing states are symmetric but not the same particles
            if (GainSumN3 + GainSumN4) != 0e0
                Correction = (LossSumN1+LossSumN2)/(GainSumN3+GainSumN4)``
                @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])
                @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]) .= Correction * @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])
            elseif (LossSumN1+LossSumN2) != 0e0
                # no gain term but there is a loss term, as outgoing states are not identical to incoming there is no way to correct for this so have to set loss terms to 0.0
                CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0
            end
        end=#

    end
    println("wrong = $num_wrong")
    println("right = $num_right")

    #println(sum(E3_d[1:2:end]))
    #println(sum(E3_d[2:2:end]))

    return CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2

end

function GainCorrection3(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, GainMatrix3::Array{Float64, 9}, GainMatrix4::Array{Float64, 9}, LossMatrix1::Array{Float64, 6}, LossMatrix2::Array{Float64, 6})

    
    #= Different possible combinations of identical and different particles that affect how to apply conservation corrections:

        1. n1=n3 != n2=n4 (1 and 3 identical, 2 and 4 identical)
            3 parameters needed:
                alpha: consv. num for 1 and 3
                beta: consv. num for 2 and 4
                gamma: consv. energy for all
            In this setup gain of n3 or n4 could be zero from MC, so have to be careful when applying correction    
        2. all other combinations:
            2 parameters needed:
                alpha: consv. num for all
                beta: consv. energy for all

    =#

    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    tol = sqrt(eps(Float64)) # tolerance for how low the gain terms can be compared to the loss terms, if less than this tolerance then the gain terms are set to zero and not included in the correction calculation.

    if (name1 == name3) && (name2 == name4) && (name1 != name2)
        CorrType = 1
    else
        CorrType = 2
    end

    CorrectedGainMatrix3 = similar(GainMatrix3)
    CorrectedGainMatrix4 = similar(GainMatrix4)
    CorrectedLossMatrix1 = similar(LossMatrix1)
    CorrectedLossMatrix2 = similar(LossMatrix2)
    fill!(CorrectedGainMatrix3,Float64(0))
    fill!(CorrectedGainMatrix4,Float64(0))
    CorrectedLossMatrix1 .= LossMatrix1
    CorrectedLossMatrix2 .= LossMatrix2

    # underflow and overflow bins are taken to have size 0 -> p1_r[1] and p1_r[end] -> 2*p1_r[end] respectively

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    u1_r = bounds(-1.0,1.0,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    u2_r = bounds(-1.0,1.0,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    u3_r = bounds(-1.0,1.0,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    u4_r = bounds(-1.0,1.0,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);

    E1_Δ = deltaEVector(p1_r,m1);
    E1_d = E1_Δ ./ p1_d

    E2_Δ = deltaEVector(p2_r,m2);
    E2_d = E2_Δ ./ p2_d

    E3_Δ = deltaEVector(p3_r,m3);
    E3_d = E3_Δ ./ p3_d

    E4_Δ = deltaEVector(p4_r,m4);
    E4_d = E4_Δ ./ p4_d

    num_wrong = 0
    num_right = 0

    for p1 in axes(GainMatrix3, 4), u1 in axes(GainMatrix3,5), h1 in axes(GainMatrix3,6), p2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,8), h2 in axes(GainMatrix3,9)
        
        LossSumN1 = LossMatrix1[p1,u1,h1,p2,u2,h2]
        LossSumN2 = LossMatrix2[p2,u2,h2,p1,u1,h1]
        if name1 != name2 # Indistinguishable_12==false so LossMatrix2 is non-zero
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
            LossSumE2 = LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
        else
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            LossSumE2 = 0.0#LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
        end

        if LossSumN1 + LossSumN2 == 0e0
            # no loss term so no interaction at these incoming states

            num_right += 1

            continue
        end
        
        wrong = true
        nonzero_gain = true
        max_high_bins = 0
        Gain3False = false
        Gain4False = false

        alpha1 = 0.0
        alpha2 = 0.0
        beta = 0.0

        a1 = 0.0
        b1 = 0.0
        a2 = 0.0
        b2 = 0.0

        p1Big = false
        p2Big = false

        p3_offset = length(axes(GainMatrix3,1))
        p4_offset = length(axes(GainMatrix4,1))

        if E1_d[p1] > 1e0*E2_d[p2] && m1 > m2
            p1Big = true
        elseif E2_d[p2] > 1e0*E1_d[p1] && m2 > m1
            p2Big = true
        end

        while wrong && nonzero_gain

            GainSumN31 = zero(Float64)
            GainSumN32 = zero(Float64)
            GainSumN41 = zero(Float64)
            GainSumN42 = zero(Float64)
            GainSumE31 = zero(Float64)
            GainSumE32 = zero(Float64)
            GainSumE41 = zero(Float64)
            GainSumE42 = zero(Float64)
            
            high_bins = 0
            searching = true
            p3_offset = length(axes(GainMatrix3,1))
            for p3 in reverse(axes(GainMatrix3,1))
                tmpN = 0.0
                tmpE = 0.0
                for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                    if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1 # only consider gain terms that are above a certain tolerance compared to loss term to avoid numerical issues when applying correction
                        tmpN += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]
                        tmpE += GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_d[p3]
                    else 
                        continue
                    end
                end
                if searching && tmpN == 0e0 #&& p3_offset > 0
                    p3_offset -= 1
                    # keep looking at a lower p3 bin
                elseif high_bins <= max_high_bins && !p2Big
                    GainSumN32 += tmpN
                    GainSumE32 += tmpE
                    high_bins += 1
                    searching = false
                else
                    GainSumN31 += tmpN
                    GainSumE31 += tmpE
                    searching = false
                end
            end

            high_bins = 0
            searching = true
            p4_offset = length(axes(GainMatrix4,1))
            for p4 in reverse(axes(GainMatrix4,1))
                tmpN = 0.0
                tmpE = 0.0
                for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                    if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2 # only consider gain terms that are above a certain tolerance compared to loss term to avoid numerical issues when applying correction
                        tmpN += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]
                        tmpE += GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_d[p4]
                    else
                        continue
                    end
                end
                if searching && tmpN == 0e0 #&& p4_offset > 0
                    p4_offset -= 1
                    # keep looking at a lower p4 bin
                elseif high_bins <= max_high_bins && !p1Big
                    GainSumN42 += tmpN
                    GainSumE42 += tmpE
                    high_bins += 1
                    searching = false
                else
                    GainSumN41 += tmpN
                    GainSumE41 += tmpE
                    searching = false
                end
            end

            if CorrType == 1 

                a1 = GainSumN31
                b1 = GainSumN32

                a2 = GainSumN41
                b2 = GainSumN42

                ac1 = GainSumE31
                c1 = ac1 != 0.0 ? GainSumE31/GainSumN31 : 0.0
                bd1 = GainSumE32
                d1 = bd1 != 0.0 ? GainSumE32/GainSumN32 : 0.0

                ac2 = GainSumE41
                c2 = ac2 != 0.0 ? GainSumE41/GainSumN41 : 0.0
                bd2 = GainSumE42
                d2 = bd2 != 0.0 ? GainSumE42/GainSumN42 : 0.0
                
                if (a1+b1 == 0e0 || (max_high_bins==0 && a1==0.0 && b1!=0.0))# There is not gain term produced by MC or only a single bin
                    #c1 = p2Big ? E3_d[p1] : E3_d[p1-1]
                    #d1 = p2Big ? E3_d[p1+1] : E3_d[p1]
                    # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
                    c1 = p2Big ? E3_d[p1+1] : E3_d[p1]
                    d1 = p2Big ? E3_d[p1+2] : E3_d[p1+1]
                    L = LossSumN1
                    LE = LossSumE1+LossSumE2
                    e2 = a2 + b2 - LossSumN2
                    #c2 = GainSumE41/GainSumN41

                    a1 = ((ac2-c2*e2)+(d1*L-LE))/(d1-c1)
                    b1 = L-a1
                    ac1 = a1*c1
                    bd1 = b1*d1

                    Gain3False = true
                end

                if (a2+b2 == 0e0 || (max_high_bins==0 && a2==0.0 && b2!=0.0))# There is not gain term produced by MC or only a single bin
                    #c2 = p1Big ? E4_d[p2] : E4_d[p2-1]
                    #d2 = p1Big ? E4_d[p2+1] : E4_d[p2]
                    # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
                    c2 = p1Big ? E4_d[p2+1] : E4_d[p2]
                    d2 = p1Big ? E4_d[p2+2] : E4_d[p2+1]
                    L = LossSumN2
                    LE = LossSumE2 + LossSumE1
                    e1 = a1 + b1 - LossSumN1
                    #c1 = GainSumE31/GainSumN31

                    a2 = ((ac1-c1*e1)+(d2*L-LE))/(d2-c2)
                    b2 = L-a2
                    ac2 = a2*c2
                    bd2 = b2*d2
                    Gain4False = true
                end

                e1 = a1 + b1 - LossSumN1
                e2 = a2 + b2 - LossSumN2
                f  = ac1 + bd1 + ac2 + bd2 - LossSumE1 - LossSumE2

                alpha1 = (b2*(d2-c2)*e1+b1*(d1*e1+c2*e2-f))/(a1*(b1*(c1-d1)+b2*(c2-d2))) + 1
                alpha2 = (b1*(d1-c1)*e2+b2*(d2*e2+c1*e1-f))/(a2*(b1*(c1-d1)+b2*(c2-d2))) + 1
                beta = (f-c1*e1-c2*e2)/(b1*(c1-d1)+b2*(c2-d2)) + 1

                if a1 <= 0e0 || a2 <= 0e0 || max_high_bins > max(size(GainMatrix3,1),size(GainMatrix4,1))# loop has not produced any corrected gain terms

                    nonzero_gain = false
                    #println("No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    #println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42, LossSumN1: $LossSumN1, LossSumN2: $LossSumN2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
                    alpha1 = 0.0
                    alpha2 = 0.0
                    beta = 0.0
                    CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                    CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0

                    num_wrong += 1

                    continue

                end

                if isnan(alpha1) || isnan(alpha2) || isnan(beta)
                    println("NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42")
                    #println("$(a+b-LossN), $(ac+bd-LossE)")
                    #println("LossE/LossN: $(LossE/LossN), bd/b: $(bd/b), ac/a: $(ac/a)")
                end

                if  alpha1 < 0e0 || alpha2 < 0e0 || beta < 0e0

                    #println("p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    #println("new gain1 = $(a1*alpha1+b1*beta), new gain2 = $(a2*alpha2+b2*beta), new gain = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)), Loss = $(LossSumN1+LossSumN2), err = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)- (LossSumN1+LossSumN2)), L1 = $LossSumN1, L2 = $LossSumN2")
                    wrong = false
                    num_right += 1
                end

            elseif CorrType == 2

                GainN1 = GainSumN31 + GainSumN41
                GainE1 = GainSumE31 + GainSumE41
                GainN2 = GainSumN32 + GainSumN42
                GainE2 = GainSumE32 + GainSumE42
                LossN = LossSumN1 + LossSumN2
                LossE = LossSumE1 + LossSumE2

                a = GainN1
                if a == 0e0 # loop has not produced any gain terms
                    nonzero_gain = false
                    println("No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")

                    alpha1 = 0.0
                    alpha2 = 0.0
                    beta = 0.0
                    CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                    CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0

                    num_wrong += 1
                    
                    continue
                end

                ac = GainE1
                c = GainE1/GainN1

                b = GainN2
                bd = GainE2
                d = GainE2/GainN2

                e = a+b-LossN
                f = ac+bd-LossE

                if isapprox(c,d)
                    # In this case deltaE is the approx the same for all bins involved, so energy is approx conserved if number is conserved so only one equation and one parameter needed which we take scales b
                    #alpha = -e/(a+b) + 1
                    #beta = alpha
                    beta = -e/b +1
                    alpha = 1.0
                else
                    alpha = (e*d-f)/(ac-a*d)+1
                    beta = (e*c-f)/(bd-b*c)+1
                end

                if isnan(alpha) || isnan(beta)
                    println("NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    println("a: $a, b: $b, ac: $ac, bd: $bd, e: $e, f: $f, c: $c, d: $d")
                    println("$(a+b-LossN), $(ac+bd-LossE)")
                    println("LossE/LossN: $(LossE/LossN), bd/b: $(bd/b), ac/a: $(ac/a)")
                end

                if alpha < 0e0 || beta < 0e0

                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    wrong = false
                    num_right += 1
                    alpha1 = alpha
                    alpha2 = alpha
                end

            end # if CorrType

        end # while

        if Gain3False 
            # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
            CorrectedGainMatrix3[p1+1+1,u1,h1,p1,u1,h1,p2,u2,h2] = p2Big * b1 * beta
            CorrectedGainMatrix3[p1+1,u1,h1,p1,u1,h1,p2,u2,h2] = p2Big ? a1 * alpha1 : b1 * beta
            CorrectedGainMatrix3[p1-1+1,u1,h1,p1,u1,h1,p2,u2,h2] = !p2Big * a1 * alpha1  
        else
            for p3 in axes(GainMatrix3,1)
                for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                    if p3 >= p3_offset-(max_high_bins) && !p2Big
                        if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] * beta
                        else
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        end
                    else
                        if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] * alpha1
                        else
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        end
                    end
                end
            end
        end

        if Gain4False 
            # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
            CorrectedGainMatrix4[p2+1+1,u2,h2,p1,u1,h1,p2,u2,h2] = p1Big * b2 * beta
            CorrectedGainMatrix4[p2+1,u2,h2,p1,u1,h1,p2,u2,h2] = p1Big ? a2 * alpha2 : b2 * beta
            CorrectedGainMatrix4[p2-1+1,u2,h2,p1,u1,h1,p2,u2,h2] = !p1Big * a2 * alpha2
        else
            for p4 in axes(GainMatrix4,1)
                for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                    if p4 >= p4_offset-(max_high_bins) && !p1Big
                        if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] * beta
                        else
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        end
                    else
                        if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] * alpha2
                        else
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        end
                    end
                end
            end
        end

        # check correction has worked for this bin

        num_check = sum( @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) + sum(@view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) - CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] - CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1]

        eng_check = sum(E3_d .* @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) + sum(E4_d .* @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) - E1_d[p1] * CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] - E2_d[p2] * CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1]

        if !isapprox(num_check,0e0,atol=sqrt(eps(Float64))) || !isapprox(eng_check,0e0,atol=sqrt(eps(Float64)))
            println("Correction failed for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2, Number check: $num_check, Energy check: $eng_check")
            #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
            #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
        end

        #println("p1=$p1,p2=$p2,u1=$u1,u2=$u2,h1=$h1, h2=$h2")

        #println("$(sum(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) , $(sum(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) , $(CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2]) , $(CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1])")

    end
    println("Number of bins that couldn't be corrected = $num_wrong")
    println("Number of bins corrected = $num_right")

    return CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2

end


=#


