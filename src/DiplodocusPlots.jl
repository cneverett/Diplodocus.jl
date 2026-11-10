module DiplodocusPlots

    export DiplodocusDark, DiplodocusLight
    export Static, Animated, Interactive
    export InteractiveBinaryGainLossPlot, InteractiveEmissionGainLossPlot

    # OD plots
    export MomentumDistributionPlot0D, MomentumAndPolarAngleDistributionPlot0D, MomentumComboAnimation0D, MomentumAndAzimuthalAngleDistributionPlot0D, AzimuthalAngleDistributionPlot0D
    #export AngleDistributionPlot,
    export TimeScalePlot
    export NumberDensityPlot0D
    export EnergyDensityPlot0D
    export IsThermalAndIsotropicPlot0D

    # 1D plots
    export VelocityZPlot1D, PressurePlot1D, TemperaturePlot1D
    export NumberPlot1D, NumberDensityPlot1D
    #export MomentumDistributionPlot1D, NumberDensityPlot1D, NumberPlot1D

    # 2D plots
    export NumberDensityPlot2D, NumberPlot2D, BulkGammaPlot2D, EnergyDensityPlot2D

    # Utilities 
    #export IsThermalPlot, IsIsotropicPlot, 
    export CodeToSIUnitsTime, CodeToCodeUnitsTime, SIToCodeUnitsTime, SyncToCodeUnitsTime, CodeToSyncUnitsTime
    # export ObserverFluxPlot

    # Not currently exported plots from papers etc.
    # export EquatorialNumberDensityPlot1D
    # export EquatorialTemperaturePlot1D

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
    
    # Static Plots
    include("TransportPlots/Experimental/ObserverFlux.jl")
    include("TransportPlots/StaticAnimated/TimeScale.jl")

    # 0D (space) plots
    include("TransportPlots/StaticAnimated/Distribution0D.jl")
    include("TransportPlots/StaticAnimated/NumberDensity0D.jl")
    include("TransportPlots/StaticAnimated/EnergyDensity0D.jl")
    include("TransportPlots/StaticAnimated/Number0D.jl")
    include("TransportPlots/StaticAnimated/ThermalAndIsotropic0D.jl")

    # 1D (space) plots
    include("TransportPlots/StaticAnimated/Distribution1D.jl")
    include("TransportPlots/StaticAnimated/NumberDensity1D.jl")
    include("TransportPlots/StaticAnimated/Number1D.jl")
    include("TransportPlots/StaticAnimated/Hydro1D.jl")

    # 2D (space) plots
    include("TransportPlots/StaticAnimated/2DPlotUtilities.jl")
    include("TransportPlots/StaticAnimated/NumberDensity2D.jl")
    include("TransportPlots/StaticAnimated/Hydro2D.jl")
    include("TransportPlots/StaticAnimated/Number2D.jl")
    include("TransportPlots/StaticAnimated/EnergyDensity2D.jl")
    include("TransportPlots/StaticAnimated/Magnetisation2D.jl")
    include("TransportPlots/StaticAnimated/ElectroMagneticEnergyDensity2D.jl")
    include("TransportPlots/StaticAnimated/ElectroMagneticField2D.jl")

    # paper specific plots are not exported 
    include("TransportPlots/PaperSpecificPlots.jl")

end
