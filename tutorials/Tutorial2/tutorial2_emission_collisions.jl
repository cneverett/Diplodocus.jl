using Diplodocus

name1 = "Ele"
name2 = "Ele"
name3 = "Pho"
type = "Sync"

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

    Ext::Vector{Float64} = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1];

    numLoss = 256
    numGain = 256
    numThreads = Threads.nthreads()

scale = 0.0:0.1:0.0

fileLocation = joinpath(pwd(),"Data")

(Setup,fileName) = UserEmissiveParameters()

EmissiveInteractionIntegration(Setup)

Output = EmissionFileLoad_Matrix(fileLocation,fileName);
InteractiveEmissiveGainLossPlot(Output)
