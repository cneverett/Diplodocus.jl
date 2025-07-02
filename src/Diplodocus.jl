module Diplodocus

    # Transport Exports
    export LoadMatrices, BigMatrices, FluxMatrices
    export PhaseSpaceStruct, MomentumStruct, SpaceStruct, TimeStruct, OutputStruct
    export BinaryStruct, EmiStruct, ForceType
    export Cylindrical, Spherical, Cartesian, Axi
    export SyncRadReact
    export BuildBigMatrices, BuildFluxMatrices
    export Initial_Constant, Initial_MaxwellJuttner, Initial_PowerLaw
    export Solve, EulerStruct

    # Plot Exports 
    export DiplodocusDark, DiplodocusLight
    export MomentumDistributionPlot

    # RecursiveArrayTools
    export ArrayPartition

    using DiplodocusCollisions
    using DiplodocusTransport
    using DiplodocusPlots
    using RecursiveArrayTools

end
