using Diplodocus

name1 = "Ele"
name2 = "Pho"
name3 = "Ele"
name4 = "Pho"

    p_low_Ele = -3.0
    p_up_Ele = 9.0
    p_grid_Ele = "l"
    p_num_Ele = 288
    u_grid_Ele = "u"
    u_num_Ele = 1
    h_grid_Ele = "u"
    h_num_Ele = 1

    p_low_Pho = -16.0
    p_up_Pho = 9.0
    p_grid_Pho = "l"
    p_num_Pho = 600
    u_grid_Pho = "u"
    u_num_Pho = 1
    h_grid_Pho = "u"
    h_num_Pho = 1

numLoss = 16
numGain = 16
numThreads = Threads.nthreads()

scale = 0.0:0.1:1.0

fileLocation = joinpath(pwd(),"Data")

(Setup,fileName) = UserBinaryParameters()

BinaryInteractionIntegration(Setup)

Output = BinaryFileLoad_Matrix(fileLocation,fileName);
InteractiveBinaryGainLossPlot(Output)

DoesConserve(Output)

Output = BinaryFileLoad_Matrix(fileLocation,fileName,corrected=true);
DoesConserve(Output)
