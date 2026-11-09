"""
    SolutionFileLoad(fileLocation::String, fileName::String)

Loads the solution and PhaseSpace struct from a file at the specified `fileLocation` with the given `fileName`.
"""
function SolutionFileLoad(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r");

        PhaseSpace = f["PhaseSpace"];
        sol = f["sol"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (PhaseSpace, sol)

end

"""
    SolutionFileSave(sol,output,fileLocation,fileName)

Saves the solution `sol` and the `PhaseSpace` struct to a file at the specified `fileLocation` with the given `fileName`.
"""
function SolutionFileSave(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,fileLocation::String,fileName::String)
        
    filePath = joinpath(fileLocation,fileName)

    f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
    
    write(f,"PhaseSpace",PhaseSpace)
    write(f,"sol",sol)

    close(f)

end