module Diplodocus

    using DiplodocusCollisions
    using DiplodocusTransport
    using DiplodocusPlots
    using RecursiveArrayTools

    # Collision Exports
    export UserBinaryParameters
    export BinaryInteractionIntegration
    export BinaryFileLoad_Matrix, BinaryFileLoad_All, DoesConserve

    # Transport Exports
    export LoadMatrices, BigMatrices, FluxMatrices
    export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct, OutputStruct
    export BinaryStruct, EmiStruct, ForceType
    export Cylindrical, Spherical, Cartesian, Axi
    export SyncRadReact
    export BuildBigMatrices, BuildFluxMatrices
    export Initial_Constant, Initial_MaxwellJuttner, Initial_PowerLaw
    export Solve, EulerStruct
    export SolutionFileLoad

    # Plot Exports 
    export DiplodocusDark, DiplodocusLight
    export MomentumDistributionPlot, MomentumAndPolarAngleDistributionPlot
    export InteractiveBinaryGainLossPlot
    export FracNumberDensityPlot, NumberDensityPlot
    export FracEnergyDensityPlot, EnergyDensityPlot
    export IsThermalPlot, IsIsotropicPlot, IsThermalAndIsotropicPlot

    # RecursiveArrayTools
    export ArrayPartition

end
