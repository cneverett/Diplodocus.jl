using Diplodocus

name1 = "Sph"
name2 = "Sph"
name3 = "Sph"
name4 = "Sph"

p_low_Sph = -5.0
p_up_Sph = 4.0
p_grid_Sph = "l"
p_num_Sph = 72

u_grid_Sph = "u"
u_num_Sph = 9

h_grid_Sph = "u"
h_num_Sph = 1

numLoss = 16
numGain = 16
numThreads = Threads.nthreads()

scale = 0.0:0.1:0.2

fileLocation = joinpath(pwd(),"Data")

(Setup,fileName) = UserBinaryParameters()

BinaryInteractionIntegration(Setup)

Output = BinaryFileLoad_Matrix(fileLocation,fileName);
InteractiveBinaryGainLossPlot(Output)

DoesConserve(Output)

Output = BinaryFileLoad_Matrix(fileLocation,fileName,corrected=true);
DoesConserve(Output)
