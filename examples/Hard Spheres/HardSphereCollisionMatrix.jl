
using Diplodocus

# ========= Generate Collision Matrices ========== #  

    # Particle definitions
    name1 = "Sph"
    name2 = "Sph"
    name3 = "Sph"
    name4 = "Sph" 

    # Momentum space grids
    p_low_Sph = -5.0
    p_up_Sph = 4.0
    p_grid_Sph = "l"
    p_num_Sph = 72

    u_grid_Sph = "u"
    u_num_Sph = 8

    h_grid_Sph = "u"
    h_num_Sph = 1

    numLoss = p_num_Sph^2*u_num_Sph^2*h_num_Sph^2
    numGain = p_num_Sph*u_num_Sph*h_num_Sph
    numThreads = 2

    scale = 0.0:0.1:0.0

    fileLocation = pwd()*"\\examples\\Hard Spheres\\Data"

    (Setup,fileName) = UserBinaryParameters()
    
    BinaryInteractionIntegration(Setup)

# ===== Check Accuracy of  Collision Matrices ===== #

    Output = BinaryFileLoad_Matrix(fileLocation,fileName); 
    InteractiveBinaryGainLossPlot(Output)

    DoesConserve(Output)

    Output = BinaryFileLoad_Matrix(fileLocation,fileName,corrected=true); 
    DoesConserve(Output)


