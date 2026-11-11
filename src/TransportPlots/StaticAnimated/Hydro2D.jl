"""
    BulkGammaPlot2D(type,PhaseSpace,sol,species,coordinates,surface,time_idx;...)

Returns a density plot of the bulk gamma of a `species` at a given `time_idx` on a 2D surface of spatial `coordinates` (These coordinates may or may not be the same coordinates used in a simulation). The surface is selected using the tuple `surface`, where the first entry indicates the surface e.g. `"xy"`, `"yz"` or `"xz"` and the second entry is an index defining the coordinate index of the coordinate remaining constant.

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
- `gamma_sign::Bool=true`: If `true`, the output values will be signed based on the sign of the local momentum space `z` velocity.
- `gamma_z::Bool=false`: If `true`, only the component of the bulk gamma with respect to the global coordinate `z` direction will be plotted. 
- `obs_frame::Bool=false`: If `true` the anisotropy of the local momentum distribution is ignored and the bulk velocity of the local observer frame itself if plotted.
"""
function BulkGammaPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64},time_idx::Int64;fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing),gamma_sign::Bool=true,gamma_z::Bool=false,obs_frame::Bool=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    surface_string = surface[1]
    surface_index = surface[2]

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

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

    if wide
        fig = Figure(size=figsize) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    Gamma = zeros(Float64,nx_dir,ny_dir)
    x_points = zeros(Float64,nx_dir,ny_dir)
    y_points = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)
    g::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer

    t = sol.t[time_idx]

    pos_sim[1] = t

    for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

        (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

        # TODO: add sub sampling of plot points
        pos_sim[2] = mx_sim[x]
        pos_sim[3] = my_sim[y]
        pos_sim[4] = mz_sim[z]

        MetricComponents!(pos_sim,g,metric,sim_coordinates)

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

        InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
        TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
        TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
        InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

        Nᵃ_sim = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
        if obs_frame     
            Nᵃ_sim = Float64[1.0,0.0,0.0,0.0]  # for testing
        end
        n = ScalarNumberDensity(Nᵃ_sim,Uₐ)

        if n != 0.0
            Vᵇ = HydroFourVelocity(Nᵃ_sim)
            CoordinateBasisTransform!(e_sim' * Vᵇ,Vᵇ,pos_sim,sim_coordinates,coordinates)
            Vᵃ = inve_plot' * Vᵇ # transform to stationary tetrad frame
            if gamma_z # static observed z component of bulk gamma
                if gamma_sign
                    βz = Vᵃ[4]/Vᵃ[1]
                    Γ = 1/sqrt(1 - βz^2) * sign(Vᵃ[4]) # bulk gamma with sign of velocity to distinguish inflow vs outflow
                else
                    βz = Vᵃ[4]/Vᵃ[1]
                    Γ = 1/sqrt(1 - βz^2) # bulk gamma without sign, so all flow contributes positively
                end
            else
                if gamma_sign
                    Γ = (-transpose(Uₐ) * Vᵃ) * sign(Vᵃ[4]) # bulk gamma with sign of velocity to distinguish inflow vs outflow
                else
                    Γ = -transpose(Uₐ) * Vᵃ # bulk gamma without sign, so all flow contributes positively
                end
            end
        else
            Γ = NaN
        end

        println("Γ = $Γ at (x,y,z) = ($x,$y,$z)")

        Gamma[x_plot_idx,y_plot_idx] += Γ
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

    pos_sim_min = [t,xr_sim[1],yr_sim[1],zr_sim[1]]
    pos_sim_max = [t,xr_sim[end],yr_sim[end],zr_sim[end]]
    pos_plot_min = zeros(Float64,4)
    pos_plot_max = zeros(Float64,4)
    CoordinateTransform!(pos_sim_min,pos_plot_min,sim_coordinates,coordinates)
    CoordinateTransform!(pos_sim_max,pos_plot_max,sim_coordinates,coordinates)

    #=if surface_string == "xy"
        bb = BBox(pos_plot_min[2],pos_plot_max[2],pos_plot_min[3],pos_plot_max[3])
    elseif surface_string == "yz"
        bb = BBox(pos_plot_min[3],pos_plot_max[3],pos_plot_min[4],pos_plot_max[4])
    elseif surface_string == "xz"
        bb = BBox(pos_plot_min[2],pos_plot_max[2],pos_plot_min[4],pos_plot_max[4])
    end=#

    #hm = heatmap!(ax,xr_dir,yr_dir,Gamma,colormap=theme.colormap_var)  
    hm = voronoiplot!(ax,x_points,y_points,Gamma,colormap=theme.colormap_var,show_generators=false,strokewidth=0.0)       

    t_unit_string = TimeUnits()

    Colorbar(fig[1,2],hm,label= gamma_z ? L"$\Gamma_z\,$" : L"$\Gamma\,$")

    end # with_theme

    return fig

end

"""
    HydroZVelocityPlot2D(type,sol,PhaseSpace,surface,species,time_idx;...)

Returns a density plot of the local momentum space z component of the bulk velocity of a `species` at a given `time_idx` on a 2D surface of spatial `coordinates` (These coordinates may or may not be the same coordinates used in a simulation). The surface is selected using the tuple `surface`, where the first entry indicates the surface e.g. `"xy"`, `"yz"` or `"xz"` and the second entry is an index defining the coordinate index of the coordinate remaining constant.

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
"""
function HydroZVelocityPlot2D(type::Static,sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,surface::String,species::String,time_idx::Int64;fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing))

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Space = PhaseSpace.Space
    Time = PhaseSpace.Time
    Grids = PhaseSpace.Grids
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    dx = Grids.dx
    dy = Grids.dy
    dz = Grids.dz
    mx = Grids.mx
    my = Grids.my
    mz = Grids.mz
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    if surface == "xy"
        xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        ylab = L"$\left(y\,[\text{Code Units}]\right)$"
        nx_dir = x_num
        ny_dir = y_num
        xgrid_type = Space.x_grid  
        ygrid_type = Space.y_grid
        if xgrid_type == "u"
            mxdir_plot = mx
            xr_dir = xr
            xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(mx)
            xr_dir = log10.(xr)
            xlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = my
            yr_dir = yr
            ylab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(my)
            yr_dir = log10.(yr)
            ylab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
    elseif surface == "yz"
        xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        nx_dir = y_num
        ny_dir = z_num
        xgrid_type = Space.y_grid
        ygrid_type = Space.z_grid
        if xgrid_type == "u"
            mxdir_plot = my
            xr_dir = yr
            xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(my)
            xr_dir = log10.(yr)
            xlab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = mz
            yr_dir = zr
            ylab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(mz)
            yr_dir = log10.(zr)
            ylab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    elseif surface == "xz"
        nx_dir = x_num
        ny_dir = z_num
        xgrid_type = Space.x_grid
        ygrid_type = Space.z_grid
        if xgrid_type == "u"
            mxdir_plot = mx
            xr_dir = xr
            xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif xgrid_type == "l"
            mxdir_plot = log10.(mx)
            xr_dir = log10.(xr)
            xlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
        if ygrid_type == "u"
            mydir_plot = mz
            yr_dir = zr
            ylab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif ygrid_type == "l"
            mydir_plot = log10.(mz)
            yr_dir = log10.(zr)
            ylab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    else
        error("Invalid surface specified. Must be one of 'xy', 'yz', or 'xz'.")
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    Vz = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    pr = Grids.pxr_list[species_index]
    ur = Grids.pyr_list[species_index]
    hr = Grids.pzr_list[species_index]
    mass = Grids.mass_list[species_index]

    f1D = zeros(Float32,p_num*u_num*h_num)

    t = sol.t[time_idx]

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        Nᵃ = DiplodocusTransport.FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)

        HydroFourVelocity = DiplodocusTransport.HydroFourVelocity(Nᵃ)

        Vz_local = HydroFourVelocity[4]/HydroFourVelocity[1]

        println("Vz_local = $Vz_local at (x,y,z) = ($x,$y,$z)")

        #println("n = $n at (x,y,z) = ($x,$y,$z)")

        if surface == "xy"
            Vz[x,y] += Vz_local * dz[z] / (sum(dz))
        elseif surface == "yz"
            Vz[y,z] += Vz_local * dx[x] / (sum(dx))
        elseif surface == "xz"
            Vz[x,z] += Vz_local * dy[y] / (sum(dy))
        end

    end

    replace!(Vz,0.0=>NaN)

    hm = heatmap!(ax,xr_dir,yr_dir,Vz,colormap=theme.colormap_var)       

    t_unit_string = TimeUnits()

    Colorbar(fig[1,2],hm,label=L"$v_z\,$")

    end # with_theme

    return fig

end

