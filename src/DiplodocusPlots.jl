module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export Static, Animated, Interactive
    export InteractiveBinaryGainLossPlot, InteractiveEmissionGainLossPlot

    # OD plots
    export MomentumDistributionPlot0D, MomentumAndPolarAngleDistributionPlot0D, MomentumComboAnimation0D
    # export AngleDistributionPlot, MomentumAndAzimuthalAngleDistributionPlot, AzimuthalAngleDistributionPlot
    #export TimeScalePlot
    export NumberDensityPlot0D
    export EnergyDensityPlot0D
    export IsThermalAndIsotropicPlot0D

    # 1D plots
    export MomentumDistributionPlot1D, NumberDensityPlot1D, NumberPlot1D

    # 2D plots
    export NumberDensityPlot2D, NumberPlot2D, BulkGammaPlot2D, EnergyDensityPlot2D

    # Utilities 
    #export IsThermalPlot, IsIsotropicPlot, 
    export CodeToSIUnitsTime, CodeToCodeUnitsTime, SIToCodeUnitsTime, SyncToCodeUnitsTime, CodeToSyncUnitsTime
    # export ObserverFluxPlot

    using CairoMakie
    using GLMakie
    using StaticArrays
    using LinearAlgebra
    using SparseArrays

    import DiplodocusCollisions as DC
    using DiplodocusTransport

    include("Themes.jl")
    include("PlotTypes.jl")
    include("TimeUnits.jl")
    include("CollisionPlots/Interactive/AnalyticalSpectra.jl/InverseCompton.jl")
    include("CollisionPlots/Interactive/GainLossAnalysis.jl")

    # Transport Plots
    #Interactive Plots
    #include("TransportPlots/Interactive/ObserverFlux.jl")
    # Static Plots
    include("TransportPlots/Static/ObserverFlux.jl")
    include("TransportPlots/Static/TimeScale.jl")

    # 0D (space) plots
    include("TransportPlots/Plots_0D/Distribution0D.jl")
    include("TransportPlots/Plots_0D/NumberDensity0D.jl")
    include("TransportPlots/Plots_0D/EnergyDensity0D.jl")
    include("TransportPlots/Plots_0D/Number0D.jl")
    include("TransportPlots/Plots_0D/ThermalAndIsotropic0D.jl")

    # 1D (space) plots
    include("TransportPlots/Static/Distribution1D.jl")
    include("TransportPlots/Static/NumberDensity1D.jl")
    include("TransportPlots/Static/Number1D.jl")
    include("TransportPlots/Static/Hydro1D.jl")

    # 2D (space) plots
    include("TransportPlots/Plots_2D/2DPlotUtilities.jl")
    include("TransportPlots/Static/NumberDensity2D.jl")
    include("TransportPlots/Static/Hydro2D.jl")
    include("TransportPlots/Static/Number2D.jl")
    include("TransportPlots/Static/EnergyDensity2D.jl")
    include("TransportPlots/Static/Magnetisation2D.jl")
    include("TransportPlots/Static/ElectroMagneticEnergyDensity2D.jl")
    include("TransportPlots/Static/ElectroMagneticField2D.jl")

end
