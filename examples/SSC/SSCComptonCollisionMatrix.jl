#=
This script is provided as an example of how to run the BinaryInteractionSpectra.jl module.
=#

    using Diplodocus

# First select the binary interaction (12->34) you want to evaluate the collision integral for by setting the names of particles 1 through 4. The pairs of names (12) and (34) should be in alphabetical order.

    name1::String = "Ele";
    name2::String = "Pho";
    name3::String = "Ele";
    name4::String = "Pho";

# Define the Momentum space discretisation for each particle species named. This includes the upper and lower momentum bounds, grid type "l", "u", or "b", and the number of bins. Must be of the format `p_low_name`, p_low_name`, `p_grid_name`, `p_num_name`, `u_grid_name` and `u_num_name` where `name` is the abbreviated three letter name of the particle species. 


    p_low_Ele = -3e0
    p_up_Ele = 7e0
    p_grid_Ele = "l"
    p_num_Ele = 80
    u_num_Ele = 9
    u_grid_Ele = "u"
    h_num_Ele = 1
    h_grid_Ele = "u"

    p_low_Pho = -15e0
    p_up_Pho = 7e0
    p_grid_Pho = "l"
    p_num_Pho = 88
    u_num_Pho = 9
    u_grid_Pho = "u"
    h_num_Pho = 1
    h_grid_Pho = "u"

    numLoss = 16
    numGain = 16
    numThreads = 15

    scale = 0.0:0.1:1.0

    fileLocation = pwd()*"\\examples\\Data"
    (Setup,fileName) = UserBinaryParameters()

    BinaryInteractionIntegration(Setup)

    Output = BinaryFileLoad_Matrix(fileLocation,fileName,corrected=true); 
    InteractiveBinaryGainLossPlot(Output)

    DoesConserve(Output,Tuple_Output=false)

    #=

    NGain3, NLoss1, NErr, NGain4, NLoss2= DoesConserve(Output,Tuple_Output=true);
    NGain3, NLoss1, NErr1, NGain4, NLoss2, NErr2, GainS3,GainS4,LossS1,LossS2 = DiplodocusCollisions.DoesConserve2(OutputAll, Tuple_Output=true);

    test1 = dropdims(sum(NErr,dims=(2,3,5,6))./length(@view(NErr[1,:,:,1,:,:])),dims=(2,3,5,6))
    test2 = dropdims(sum(NGain3,dims=(2,3,5,6))./length(@view(NErr[1,:,:,1,:,:])),dims=(2,3,5,6))
    test3 = dropdims(sum(NLoss1,dims=(2,3,5,6))./length(@view(NErr[1,:,:,1,:,:])),dims=(2,3,5,6))

    =#


