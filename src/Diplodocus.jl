module Diplodocus

    using DiplodocusCollisions
    using DiplodocusTransport
    using DiplodocusPlots

    # Collision Exports
    export BinaryInteractionIntegration, BinaryFileLoad_All, DoesConserve, BinaryFileLoad_Matrix
    export EmissionInteractionIntegration, EmissionFileLoad_Matrix, EmissionFileLoad_All
    export UserBinaryParameters, UserEmissionParameters

    # Transport Exports
    export LoadMatrices, BigMatrices, FluxMatrices
    export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct, OutputStruct, CharacteristicStruct, GridsStruct, ElectroMagneticFieldStruct
    export BinaryStruct, EmiStruct, ForceType
    export BackendType
    export CoordinateType, Cylindrical, Spherical, Cartesian
    export ModeType, Ani, Axi, Iso
    export BoundaryType, Periodic, Open, Closed, Reflective
    export CoordinateForce, SyncRadReact, GradBInvZDecay
    export BuildBinaryMatrices, BuildEmissionMatrices, BuildFluxMatrices
    export Initialise_Initial_Condition, Location_Species_To_StateVector, Initial_Constant!, Initial_MaxwellJuttner!, Initial_PowerLaw!, Initial_BoostedPowerLaw!, Initial_BlackBody!
    export Initialise_Injection_Condition, Injection_Constant!, Injection_MaxwellJuttner!, Injection_PowerLaw!, Injection_BoostedPowerLaw!, Injection_BlackBody!
    export ElectroMagneticField_Constant, ElectroMagneticField_InvZDecay
    export CollisionDomain
    export Solve, ForwardEulerStruct, ForwardSymplecticEulerStruct
    export SolutionFileLoad
    export MaxwellJuttner_Distribution, BlackBody_Distribution, Wien_Distribution
    export CUDABackend, CPUBackend

    # Plot Exports 
    export DiplodocusDark, DiplodocusLight
    export Static, Animated, Interactive
    export InteractiveBinaryGainLossPlot, InteractiveEmissionGainLossPlot
    export MomentumDistributionPlot0D, MomentumAndPolarAngleDistributionPlot0D
    export NumberDensityPlot0D
    export EnergyDensityPlot0D
    export IsThermalAndIsotropicPlot0D

end
