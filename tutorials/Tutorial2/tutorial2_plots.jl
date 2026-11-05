using Diplodocus

fileName = "tutorial2_output.jld2"
fileLocation = joinpath(pwd(),"Data")
(PhaseSpace,sol) = SolutionFileLoad(fileLocation,fileName);
