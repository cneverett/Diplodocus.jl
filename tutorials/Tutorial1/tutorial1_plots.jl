using Diplodocus

fileName = "tutorial1_output.jld2"
fileLocation = joinpath(pwd(),"Data")
(PhaseSpace,sol) = SolutionFileLoad(fileLocation,fileName);

PAndUDisPlot = MomentumAndPolarAngleDistributionPlot(Static(),PhaseSpace,sol,["Sph"],(1,11,1001))

save("Tutorial1_PAndUDisPlot.svg",PAndUDisPlot)

PDisPlot = MomentumDistributionPlot0D(Static(),PhaseSpace,sol,["Sph"];step=10,thermal=true,plot_limits=(-1.4,2.1,-3.4,0.9))
save("Tutorial1_PDisPlot.svg",PDisPlot)

TandIPlot = IsThermalAndIsotropicPlot0D(Static(),PhaseSpace,sol,"Sph")
save("Tutorial1_TAndIPlot.svg",TandIPlot)

NumPlot = NumberDensityPlot0D(Static(),PhaseSpace,sol,species="Sph";frac=true)
save("Tutorial1_NumPlot.svg",NumPlot)
EngPlot = EnergyDensityPlot0D(Static(),PhaseSpace,sol,species="Sph";frac=true,perparticle=true)
save("Tutorial1_EngPlot.svg",EngPlot)

AniPlot = MomentumComboAnimation0D(Animated(),PhaseSpace,sol,["Sph"];plot_limits_momentum=(-0.2,1.9,-2.1,0.8),filename="Tutorial1_Animation.mp4",thermal=true,framerate=24)
