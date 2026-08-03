"""
    NumberPlot2D(type,PhaseSpace,sol,species,surface;...)

Returns a density plot of the number of a `species` at a given `time_idx` on a specific spatial coordinate `surface`.
"""
function NumberPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing))

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

    if wide
        fig = Figure(size=(2*576,2*216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    num = zeros(Float64,nx_dir,ny_dir)
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

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

        TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
        InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
        TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
        InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
        CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
        Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
        num_den = ScalarNumberDensity(Nᵃ,Uₐ)
        Vol = VolFunction(PhaseSpace,1,x,y,z)

        if (y==1 && z==3) || (y==3 && z==1)
        println("n = $num_den at (x,y,z) = ($x,$y,$z), vol = $Vol")
        end

        num[x_plot_idx,y_plot_idx] += num_den * Vol
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

    replace!(num,0.0=>NaN)

    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var)   
    hm = voronoiplot!(ax,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1,nan_color=:transparent)      

    t_unit_string = TimeUnits()

    Colorbar(fig[1,2],hm,label=L"$\text{num}$")

    return fig

    end # with_theme

end

function NumberPlot2D(type::Animated,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing),framerate::Int64=12,filename="Number2D.mp4",figure=nothing)

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

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates)
    
    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(sol.t[$time_idx])
        if wide
            fig = Figure(size=(2*576,2*216))
        else
            fig = Figure()
        end
        #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab)
    end

    ax.limits = plot_limits

    num = zeros(Float64,nx_dir,ny_dir)
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

    t = sol.t[1]
    pos_sim[1] = t

    # setup grid of x, y points (only once since the grid is static)
    for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

        (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

        # TODO: add sub sampling of plot points
        pos_sim[2] = mx_sim[x]
        pos_sim[3] = my_sim[y]
        pos_sim[4] = mz_sim[z]

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

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

    # compute number density at each point
    num = @lift begin
        num = zeros(Float64,nx_dir,ny_dir)
        for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

            (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

            # TODO: add sub sampling of plot points
            pos_sim[2] = mx_sim[x]
            pos_sim[3] = my_sim[y]
            pos_sim[4] = mz_sim[z]

            CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

            TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
            InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
            InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

            Nᵇ = Float64.(FourFlow(sol.f[$time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
            CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
            Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
            num_den = ScalarNumberDensity(Nᵃ,Uₐ)
            Vol = VolFunction(PhaseSpace,1,x,y,z)

            num[x_plot_idx,y_plot_idx] += num_den * Vol

        end
        replace!(num,0.0=>NaN)
        num
    end

    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var)   
    hm = voronoiplot!(ax,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1,nan_color=:transparent)      

    t_unit_string = TimeUnits()

    Colorbar(fig[1,2],hm,label=L"$\text{num}$")

    if !isnothing(filename)
        # recording the animation
        time_idxs = 1:length(sol.t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            print("generating frame: $frame \r")
            time_idx[] = frame
        end
    end

    end # with_theme

end