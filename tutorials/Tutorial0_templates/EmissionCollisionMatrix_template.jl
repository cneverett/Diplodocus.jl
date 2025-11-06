using Diplodocus

# Particle definitions
name1::String = "name1";
name2::String = "name2";
name3::String = "name3";
type::String = "type";

# Momentum space grids
p_low_name1 = -5.0
p_up_name1 = 4.0
p_grid_name1 = "l"
p_num_name1 = 72
u_grid_name1 = "u"
u_num_name1 = 8
h_grid_name1 = "u"
h_num_name1 = 1

p_low_name2 = -5.0
p_up_name2 = 4.0
p_grid_name2 = "l"
p_num_name2 = 72
u_grid_name2 = "u"
u_num_name2 = 8
h_grid_name2 = "u"
h_num_name2 = 1

p_low_name3 = -5.0
p_up_name3 = 4.0
p_grid_name3 = "l"
p_num_name3 = 72
u_grid_name3 = "u"
u_num_name3 = 8
h_grid_name3 = "u"
h_num_name3 = 1


# Define the External parameters. For synchrotron Ext[1] = B field in Tesla 
Ext::Vector{Float64} = [1e-4,];

# number of points to sample
numLoss = 256
numGain = 256
numThreads = Threads.nthreads()

# scale factor range
scale = 0.0:0.1:0.0

# file location and setup
fileLocation = pwd()*"\\Data";
(Setup,fileName) = UserEmissionParameters()

# run the integration
EmissionInteractionIntegration(Setup)
