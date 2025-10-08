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
    u_num_Ele = 9
    h_grid_Ele = "u"
    h_num_Ele = 1

    p_low_Pho = -15e0
    p_up_Pho = 7e0
    p_grid_Pho = "l"
    p_num_Pho = 88
    u_grid_Pho = "u"
    u_num_Pho = 9
    h_grid_Pho = "u"
    h_num_Pho = 1

    # Define the External parameters. For synchrotron Ext[1] = B field in Tesla 
    Ext::Vector{Float64} = [1e-4,];

    # number of points to sample
    numLoss = 256
    numGain = 256
    numThreads = Threads.nthreads()
    
    # scale factor range
    scale = 0.0:0.1:0.0

    # file location and setup
    fileLocation = pwd()*"\\examples\\Data";
    (Setup,fileName) = UserEmissionParameters()

    # run the integration
    EmissionInteractionIntegration(Setup)

# ===== Check Accuracy of  Collision Matrices ===== #

    Output = EmissionFileLoad_Matrix(fileLocation,fileName); 
    InteractiveEmissionGainLossPlot(Output)