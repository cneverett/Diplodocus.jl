
using Diplodocus

# Particle definitions
name1 = "name1"
name2 = "name2"
name3 = "name3"
name4 = "name4" 

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

p_low_name4 = -5.0
p_up_name4 = 4.0
p_grid_name4 = "l"
p_num_name4 = 72
u_grid_name4 = "u"
u_num_name4 = 8
h_grid_name4 = "u"
h_num_name4 = 1

# number of points to sample
numLoss = 16
numGain = 16
numThreads = Threads.nthreads()
# scale factor range
scale = 0.0:0.1:1.0

# file location and setup
fileLocation = pwd()*"\\Data"
(Setup,fileName) = UserBinaryParameters()

# run the integration
BinaryInteractionIntegration(Setup)



