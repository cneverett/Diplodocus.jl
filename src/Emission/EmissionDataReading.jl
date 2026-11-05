"""
    EmissionFileLoad_All(fileLocation,fileName)

Loads all the data stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    (Run_Parameters,GainMatrix3,GainTally3) = EmissionFileLoad_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `GainMatrix3` : A 6D matrix of the emission spectrum for 1->23 interaction.
- `GainTally3` : A 6D matrix of the tallies of the number of emission spectrum values sampled for 1->23 interaction.

"""
function EmissionFileLoad_All(fileLocation::String,fileName::String)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r");
        Parameters = f["Parameters"]

        GainTally3 = f["GainTally3"];
        GainMatrix3 = f["GainMatrix3"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3,GainTally3);

end

"""
    EmissionFileLoad_Matrix(fileLocation,fileName)

Loads just the Gain and Loss Matrices stored in `fileName` stored at `fileLocation`. 

# Example
```julia-repl
    (Parameters,GainMatrix3) = EmissionFileLoad_Matrix(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `GainMatrix3` : A 6D matrix of the emission spectrum for 1->23 interaction.

"""
function EmissionFileLoad_Matrix(fileLocation::String,fileName::String)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r");
        Parameters = f["Parameters"]
        GainMatrix3 = f["GainMatrix3"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3)

end

"""
    fload_Matrix_SyncISO(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation` first converting them to an isotropic form by summing over angles. (The dimensions of the matrices stay the same i.e. 6D->6D with three dimensions having a size of 1)

# Example
```julia-repl
    Matrices = fload_Matrix_SyncISO(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `SMatrix` : A 4D matrix of the emission spectrum for Synchrotron.

"""
function fload_Matrix_SyncISO(fileLocation::String,fileName::String)
        
    filePath = joinpath(fileLocation,fileName)
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r");
        Parameters = f["Parameters"]
        SMatrix = f["SMatrix"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    SMatrixISO = sum(SMatrix,dims=(2,4))
    return (Parameters,SMatrixISO)

end