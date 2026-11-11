"""
    MagnetisationPlot2D(type,PhaseSpace,sol,species,coordinates,surface,time_idx;...)

Returns a density plot of the magnetisation of a `species` at a given `time_idx` on a 2D surface of spatial `coordinates` (These coordinates may or may not be the same coordinates used in a simulation). The surface is selected using the tuple `surface`, where the first entry indicates the surface e.g. `"xy"`, `"yz"` or `"xz"` and the second entry is an index defining the coordinate index of the coordinate remaining constant.

Arguments:
- `type::PlotType` determines the type of plot to generate, it can be either `Static` or `Animated`.
- `PhaseSpace::PhaseSpaceStruct` is the structure containing the phase space information of the simulation.
- `sol::OutputStruct` is the solution object containing the distribution function of all particles and time stepping information of the simulation. 
- `species::String` is the abbreviated three letter name of the particle species to plot.
- `coordinates::AbstractCoordinates` is the coordinate system for the plot to use (may or may not be the same as the simulation coordinates e.g. a simulation run in paraboloidal coordinates may be plotted in cylindrical coordinates).
- `surface::Tuple{String,Int64}` the first entry is a string that defines which 2D coordinates are being plotted while the third defines the index of the third coordinate that remains constant on this surface.

Common keyword arguments: 
- `theme=DiplodocusDark()`: the colour theme to use for the plot, can be either `DiplodocusDark()` for dark mode or `DiplodocusLight()` for light mode.
- `plot_limits=(nothing,nothing)`: the limits of the x and y axes, default sets the limits automatically based on the data.
- `wide=false`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), if `false` (single column 4:3 aspect ratio).
- `x_units::String="\\text{Code Units}"`: LaTeX string to be plotted as the x axis units
- `y_units::String="\\text{Code Units}"`: LaTeX string to be plotted as the x axis units
- `log_sigma=false`: If `true` the Log10 of the magnetisation will be plotted. 
"""
function MagnetisationPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,fig=nothing,theme=DiplodocusDark(),title=nothing,wide=false,plot_limits=(nothing,nothing),log_sigma=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    surface_string = surface[1]
    surface_index = surface[2]

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Spacetime = PhaseSpace.Spacetime
    Grids = PhaseSpace.Grids
    metric = Spacetime.metric
    sim_coordinates = Spacetime.coordinates
    tetrad = Spacetime.tetrad
    mx_sim = Grids.mx
    my_sim = Grids.my
    mz_sim = Grids.mz
    xr_sim = Grids.xr
    yr_sim = Grids.yr
    zr_sim = Grids.zr
    mass_list = PhaseSpace.Grids.mass_list
    B_field = PhaseSpace.Grids.B_field
    CHAR_mass_density = PhaseSpace.Characteristic.CHAR_mass_density
    CHAR_B = PhaseSpace.Characteristic.CHAR_magnetic_field
    μ0 = DiplodocusTransport.CONST_μ0
    c = DiplodocusTransport.CONST_c

    if wide
        fig = Figure(size=figsize) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    sigma = zeros(Float64,nx_dir,ny_dir)
    x_points = zeros(Float64,nx_dir,ny_dir)
    y_points = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer

    t = sol.t[time_idx]

    pos_sim[1] = t

    for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

        (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

        # TODO: add sub sampling of plot points
        pos_sim[2] = mx_sim[x]
        pos_sim[3] = my_sim[y]
        pos_sim[4] = mz_sim[z]

        B = Float64.(B_field[x,y,z])

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))

        num_den = ScalarNumberDensity(Nᵇ,Uₐ)

        sigma_local = B^2 / (2μ0 * num_den * mass_list[species_index] * c^2) 
        sigma_local *= CHAR_B^2 / (CHAR_mass_density) # convert to SI units

        sigma[x_plot_idx,y_plot_idx] += sigma_local
        if surface_string == "xy"
            x_points[x_plot_idx,y_plot_idx] = pos_plot[2]
            y_points[x_plot_idx,y_plot_idx] = pos_plot[3]
        elseif surface_string == "yz"
            x_points[x_plot_idx,y_plot_idx] = pos_plot[3]
            y_points[x_plot_idx,y_plot_idx] = pos_plot[4]
        elseif surface_string == "xz"
            x_points[x_plot_idx,y_plot_idx] = pos_plot[2]
            y_points[x_plot_idx,y_plot_idx] = pos_plot[4]
        end

    end

    replace!(sigma,0.0=>NaN)

    if log_sigma
        sigma = log10.(sigma)
        replace!(sigma, -Inf => NaN)
    end

    println("sigma = $(sigma)")

    #hm = heatmap!(ax,xr_dir,yr_dir,sigma,colormap=theme.colormap_var) 
    hm = voronoiplot!(ax,x_points,y_points,sigma,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       

    z_lab = log_sigma ? L"\log_{10}\left(\sigma\right)" : L"\sigma"

    Colorbar(fig[1,2],hm,label=z_lab)

    return fig

    end # with_theme

end