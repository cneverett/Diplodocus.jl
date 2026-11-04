"""
    EnergyDensityPlot2D(type,PhaseSpace,sol,species,surface,time_idx;...)

Returns a density plot of the energy density of a `species` at a given `time_idx` on a specific spatial coordinate `surface`.
"""
function EnergyDensityPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,fig=nothing,theme=DiplodocusDark(),title=nothing,EngDenUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),wide=false,plot_limits=(nothing,nothing),log_eng_dens=false,perparticle::Bool=false,figsize::Tuple{Int64,Int64}=(576,216),comoving::Bool=false,x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

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

    engden_unit_string = EngDenUnits[2]
    engden_unit_scale = EngDenUnits[1]

    if wide
        fig = Figure(size=figsize) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    eng = zeros(Float64,nx_dir,ny_dir)
    x_points = zeros(Float64,nx_dir,ny_dir)
    y_points = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer
    Ub::Vector{Float64} = zeros(Float64,4) # static observer in simulation frame
    η::Matrix{Float64} = [-1.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0 ; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    t = sol.t[time_idx]

    pos_sim[1] = t

    for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

        (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

        # TODO: add sub sampling of plot points
        pos_sim[2] = mx_sim[x]
        pos_sim[3] = my_sim[y]
        pos_sim[4] = mz_sim[z]

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
        Tᵃᵇ = Float64.(StressEnergyTensor(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))

        if comoving

            num_den = ScalarNumberDensity(Nᵇ,Uₐ)
            eng_den = ScalarEnergyDensity(Tᵃᵇ,Uₐ,num_den;perparticle=perparticle)

        else
            InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
            InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())
            TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
            TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)

            # for number density it is easy to transform the four-flow vector to the stationary observer frame
            CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,pos_sim,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
            Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
            num_den = ScalarNumberDensity(Nᵃ,Uₐ)
            # for energy density at Tab is a rank-2 tensor it maybe easier to transform U_a to the simulation frame and then contract with Tab in the simulation frame, as energy density is an invariant scalar this will give the same result as transforming Tab to the stationary tetrad frame and contracting with U_a in the stationary tetrad frame
            CoordinateBasisTransform!(e_plot' * η * Uₐ,Ub,pos_plot,coordinates,sim_coordinates) # transform U_a to simulation frame
            Ub = η * inve_sim' * Ub

            #println("Ub = $(Ub) at (x,y,z) = ($x,$y,$z),")
            #display(e_sim)

            eng_den = ScalarEnergyDensity(Tᵃᵇ,Ub,num_den;perparticle=perparticle)
            eng_den *= engden_unit_scale # convert to specified units

        end

        eng[x_plot_idx,y_plot_idx] += eng_den
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

    replace!(eng,0.0=>NaN)

    if log_eng_dens
        eng = log10.(eng)
        replace!(eng, -Inf => NaN)
    end

    println("eng = $(eng)")

    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var) 
    hm = voronoiplot!(ax,x_points,y_points,eng,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       

    z_lab = log_eng_dens ? L"\log_{10}\left(e\, \left[ %$(engden_unit_string)\right] \right)" : L"e\, \left[ %$(engden_unit_string) \right]"

    Colorbar(fig[1,2],hm,label=z_lab)

    return fig

    end # with_theme

end

"""
    EquipartitionPlot2D(type,PhaseSpace,sol,species,surface,time_idx;...)

Returns a density plot of the ratio of particle energy density of a `species` ti that of the magnetic field as measured by a local observer at a given `time_idx` on a specific spatial coordinate `surface`.
"""
function EquipartitionPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,fig=nothing,theme=DiplodocusDark(),title=nothing,wide=false,plot_limits=(nothing,nothing),log_equi=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

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
    CHAR_mass = PhaseSpace.Characteristic.CHAR_mass
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

    equi = zeros(Float64,nx_dir,ny_dir)
    x_points = zeros(Float64,nx_dir,ny_dir)
    y_points = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer
    Ub::Vector{Float64} = zeros(Float64,4) # static observer in simulation frame
    η::Matrix{Float64} = [-1.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0 ; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    t = sol.t[time_idx]

    pos_sim[1] = t

    for x_plot_idx in 1:nx_dir, y_plot_idx in 1:ny_dir

        (x,y,z) = PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,sim_coordinates,coordinates)

        # TODO: add sub sampling of plot points
        pos_sim[2] = mx_sim[x]
        pos_sim[3] = my_sim[y]
        pos_sim[4] = mz_sim[z]

        CoordinateTransform!(pos_sim,pos_plot,sim_coordinates,coordinates)

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
        Tᵃᵇ = Float64.(StressEnergyTensor(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))

        num_den = ScalarNumberDensity(Nᵇ,Uₐ)
        eng_den = ScalarEnergyDensity(Tᵃᵇ,Uₐ,num_den)

        B = Float64(B_field[x,y,z])

        eng_den_B = (B * CHAR_B)^2 / (2μ0) 
        eng_den_B /= PhaseSpace.Characteristic.CHAR_energy_density # convert to code units

        equi_local = eng_den / eng_den_B

        equi[x_plot_idx,y_plot_idx] += equi_local
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

    replace!(equi,0.0=>NaN)

    if log_equi
        equi = log10.(equi)
        replace!(equi, -Inf => NaN)
    end

    println("equi = $(equi)")

    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var) 
    hm = voronoiplot!(ax,x_points,y_points,equi,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       


    z_lab = log_equi ? L"\log_{10}\left(e/e_B\right)" : L"e/e_B"
    Colorbar(fig[1,2],hm,label=z_lab)

    return fig

    end # with_theme

end

#=function NumberDensityPlot2D(type::Animated,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing),framerate::Int64=12,filename="NumberDensity2D.mp4",figure=nothing,log_num_dens=false,figsize::Tuple{Int64,Int64}=(576,216))

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
            fig = Figure(size=figsize)
        else
            fig = Figure()
        end
        #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure # use the provided figure and time index
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end

    ax.limits = plot_limits

    num = zeros(Float64,nx_dir,ny_dir)
    x_points = zeros(Float64,nx_dir,ny_dir)
    y_points = zeros(Float64,nx_dir,ny_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
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

            InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())

            Nᵇ = Float64.(FourFlow(sol.f[$time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
            CoordinateTransform!(inve_sim * Nᵇ,Nᵇ,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
            Nᵃ = e_plot * Nᵇ # transform to stationary tetrad frame
            num_den = ScalarNumberDensity(Nᵃ,Uₐ)

            num[x_plot_idx,y_plot_idx] += num_den

        end

        if log_num_dens
            num = log10.(num)
            replace!(num, -Inf => NaN)
        end 
        #replace!(num,0.0=>NaN)

        num
    end


    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var)   
    hm = voronoiplot!(ax,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1,nan_color=:transparent)      

    t_unit_string = TimeUnits()

    z_lab = log_num_dens ? L"\log_{10}(n)\,\, \text{code units}" : L"n\,\, \text{code units}"
    Colorbar(fig[1,2],hm,label=z_lab)

    if !isnothing(filename)
        # recording the animation
        time_idxs = 1:length(sol.t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            print("generating frame: $frame \r")
            time_idx[] = frame
        end
    end

    end # with_theme

end=#