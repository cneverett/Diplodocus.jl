    using Diplodocus

# ========= Generate Collision Matrices ========== #

    # Particle definitions
    name1::String = "Ele";
    name2::String = "Ele";
    name3::String = "Pho";
    type::String = "Sync";

    # Momentum space grids
    p_low_Ele = -3e0
    p_up_Ele = 7e0
    p_grid_Ele = "l"
    p_num_Ele = 80
    u_grid_Ele = "u"
    u_num_Ele = 15
    h_grid_Ele = "u"
    h_num_Ele = 1

    p_low_Pho = -20e0
    p_up_Pho = 0e0
    p_grid_Pho = "l"
    p_num_Pho = 80
    u_grid_Pho = "u"
    u_num_Pho = 15
    h_grid_Pho = "u"
    h_num_Pho = 1

    # Define the External parameters. For synchrotron Ext[1] = B field in Tesla 
    Ext::Vector{Float64} = [1e-4,];

    # number of points to sample
    numLoss = 128*p_num_Pho*u_num_Pho*h_num_Pho
    numGain = 1024*p_num_Pho*u_num_Pho*h_num_Pho
    numThreads = 18
    # scale factor range
    scale = 0.0:0.1:0.0

    # file location and setup
    fileLocation = pwd()*"\\examples\\Synchrotron\\Data";
    (Setup,fileName) = UserEmissionParameters()

    fileName = "SyncEleElePho#-3.0-7.0l80#u15#u1#-3.0-7.0l80#u15#u1#-20.0-0.0l80#u15#u1#0.0001 copy 3.jld2"

    #(Parameters,scale, numLoss, numGain, numThreads, fileLocation, fileName) = Setup

    #fileName = "weightedtest.jld2"

    #Setup = (Parameters,scale, numLoss, numGain, numThreads, fileLocation, fileName)

    # run the integration
    EmissionInteractionIntegration(Setup)

# ===== Check Accuracy of  Collision Matrices ===== #

    Output = EmissionFileLoad_Matrix(fileLocation,fileName); 
    InteractiveEmissionGainLossPlot(Output)