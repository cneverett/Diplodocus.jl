"""
    NumberDensityPlot2D(type,PhaseSpace,sol,species,surface,time_idx;...)

Returns a density plot of the number density of a `species` at a given `time_idx` on a specific spatial coordinate `surface`.
"""
function NumberDensityPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,fig=nothing,theme=DiplodocusDark(),title=nothing,NumDenUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),wide=false,plot_limits=(nothing,nothing),log_num_dens=false,figsize::Tuple{Int64,Int64}=(576,216),comoving::Bool=false,x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

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

    numden_unit_string = NumDenUnits[2]
    numden_unit_scale = NumDenUnits[1]

    if wide
        fig = Figure(size=figsize) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
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

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))

        if comoving 

            num_den = ScalarNumberDensity(Nᵇ,Uₐ)

        else

            InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
            InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

            CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,pos_sim,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
            Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
            num_den = ScalarNumberDensity(Nᵃ,Uₐ) * numden_unit_scale # convert to specified units

        end

        num[x_plot_idx,y_plot_idx] += num_den
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

    if log_num_dens
        num = log10.(num)
        replace!(num, -Inf => NaN)
    end

    println("num = $(num)")

    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var) 
    hm = voronoiplot!(ax,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       

    z_lab = log_num_dens ? L"\log_{10}\left(n\, \left[ %$(numden_unit_string)\right] \right)" : L"n\, \left[ %$(numden_unit_string) \right]"

    Colorbar(fig[1,2],hm,label=z_lab,tellheight=false,tellwidth=false)

    #resize_to_layout!(fig)

    return fig

    end # with_theme

end

function NumberDensityPlot2D(type::Animated,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};theme=DiplodocusDark(),title=nothing,NumDenUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),wide=false,plot_limits=(nothing,nothing),framerate::Int64=12,filename="NumberDensity2D.mp4",figure=nothing,log_num_dens=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

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

    numden_unit_string = NumDenUnits[2]
    numden_unit_scale = NumDenUnits[1]

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

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

            InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
            TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
            InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

            Nᵇ = Float64.(FourFlow(sol.f[$time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))
            CoordinateTransform!(e_sim' * Nᵇ,Nᵇ,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
            Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
            num_den = ScalarNumberDensity(Nᵃ,Uₐ) * numden_unit_scale # convert to specified units

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

    z_lab = log_num_dens ? L"\log_{10}\left(n\, \left[ %$(numden_unit_string)\right] \right)" : L"n\, \left[ %$(numden_unit_string) \right]"
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

end

"""
    GJMultiplicityPlot2D(type,PhaseSpace,sol,species,surface,time_idx;...)

Returns a density plot of the ratio (multiplicity) between the number density of a `species` at a given `time_idx` on a specific spatial coordinate `surface` and the Goldreich-Julian number density (calculated assuming a paraboloidal force-free field structure).
"""
function GJMultiplicityPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,fig=nothing,theme=DiplodocusDark(),title=nothing,wide=false,plot_limits=(nothing,nothing),log_multi=false,figsize::Tuple{Int64,Int64}=(576,216),max_multi::Union{Nothing,Float64}=nothing,x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    @assert PhaseSpace.Spacetime.tetrad isa DiplodocusTransport.ParabolicForceFreeFieldTetrad "GJMultiplicityPlot2D only works with a ParabolicForceFreeFieldTetrad"

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
        fig = Figure(size=figsize) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    multi = zeros(Float64,nx_dir,ny_dir)
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

        Nᵇ = Float64.(FourFlow(sol.f[time_idx],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z))

        InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
        TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
        TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
        InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

        CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,pos_sim,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
        Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
        num = ScalarNumberDensity(Nᵃ,Uₐ)

        Ω = PhaseSpace.Spacetime.tetrad.Ω
        B0 = PhaseSpace.Spacetime.tetrad.B0
        # need to calculate Bϕ 
        u = pos_sim[3]
        v = pos_sim[4]
        if u^2-v^2 > 0 # z > 0
            ρ0 = v^2
            ΩF = Ω(ρ0)
            num_GJ = 2*DiplodocusTransport.CONST_ε0*B0*ΩF / ((u^2+v^2)*(1+(v^2*ΩF)^2)^(3/2)) * (1/PhaseSpace.Characteristic.CHAR_time) / DiplodocusTransport.CONST_q # Goldreich-Julian number density in SI units
        else
            ρ0 = u^2
            ΩF = Ω(ρ0)
            num_GJ = 2*DiplodocusTransport.CONST_ε0*B0*ΩF / ((u^2+v^2)*(1+(u^2*ΩF)^2)^(3/2)) * (1/PhaseSpace.Characteristic.CHAR_time) / DiplodocusTransport.CONST_q # Goldreich-Julian number density in SI units
        end
        num_GJ /= PhaseSpace.Characteristic.CHAR_number_density # convert to code units
        num_GJ = abs(num_GJ) # ensure positive number density

        println("num = $(num), num_GJ = $(num_GJ), ΩF = $(ΩF), $u, $v, $(u*v), $(PhaseSpace.Grids.B_field[x,y,z])")

        multi_local = num / num_GJ

        multi[x_plot_idx,y_plot_idx] += multi_local
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

    replace!(multi,0.0=>NaN)

    if log_multi
        multi = log10.(multi)
        replace!(multi, -Inf => NaN)
        max = maximum(x for x in multi if !isnan(x))
        min = minimum(x for x in multi if !isnan(x))
        if isnothing(max_multi)
            col_range = (min,max)
        else
            col_range = (min,log10(max_multi))
        end
    else 
        max = maximum(x for x in multi if !isnan(x))
        min = minimum(x for x in multi if !isnan(x))
        if isnothing(max_multi)
            col_range = (min,max)
        else
            col_range = (min,max_multi)
        end
    end

    println("multi = $(multi)")
    
    #hm = heatmap!(ax,xr_dir,yr_dir,num,colormap=theme.colormap_var)
    if isnothing(max_multi)
        hm = voronoiplot!(ax,x_points,y_points,multi,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1) 
    else 
        hm = voronoiplot!(ax,x_points,y_points,multi,colormap=theme.colormap_var,colorrange=col_range,show_generators=false,strokewidth=0.1,highclip=theme.colormap_var[][1.0])  
    end     

    z_lab = log_multi ? L"\log_{10}\left(n/n_\text{GJ}\right)" : L"n/n_\text{GJ}"

    Colorbar(fig[1,2],hm,label=z_lab)

    return fig

    end # with_theme

end




######### Just for Thesis ##########

# Jet loading plot 

function JetLoadingPlot2D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64};time_idx::Int64=1,theme=DiplodocusDark(),NumDenUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),EngDenUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),perparticle::Bool=false,plot_limits=(nothing,nothing),log_num_dens=false,log_eng_dens=false,figsize::Tuple{Int64,Int64}=(576,216),comoving::Bool=false,x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

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

    numden_unit_string = NumDenUnits[2]
    numden_unit_scale = NumDenUnits[1]

    engden_unit_string = EngDenUnits[2]
    engden_unit_scale = EngDenUnits[1]

    fig = Figure(size=figsize) 

    xlab, ylab, nx_dir, ny_dir = SimToPlotGrid(PhaseSpace,surface,coordinates,sim_coordinates;x_units=x_units,y_units=y_units) 

    axnumpho = Axis(fig[1,2],xlabel=xlab,ylabel=ylab,aspect=DataAspect(),xreversed=true)
    axnumpho.limits = plot_limits
    hidexdecorations!(axnumpho)

    axnumele = Axis(fig[1,3],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    axnumele.limits = plot_limits
    hidedecorations!(axnumele)

    axengpho = Axis(fig[2,2],xlabel=xlab,ylabel=ylab,aspect=DataAspect(),xreversed=true,yreversed=true)
    axengpho.limits = plot_limits
    hidedecorations!(axengpho)

    axengele = Axis(fig[2,3],xlabel=xlab,ylabel=ylab,aspect=DataAspect(),yreversed=true)
    axengele.limits = plot_limits
    hideydecorations!(axengele)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer
    Ub::Vector{Float64} = zeros(Float64,4) # static observer in simulation frame
    η::Matrix{Float64} = [-1.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0 ; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    t = sol.t[time_idx]

    pos_sim[1] = t

    # Photon plots

        num = zeros(Float64,nx_dir,ny_dir)
        eng = zeros(Float64,nx_dir,ny_dir)
        x_points = zeros(Float64,nx_dir,ny_dir)
        y_points = zeros(Float64,nx_dir,ny_dir)

        species_index = findfirst(x->x=="Pho",name_list)

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
                TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
                TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
                InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

                CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,pos_sim,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
                Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
                num_den = ScalarNumberDensity(Nᵃ,Uₐ) * numden_unit_scale # convert to specified units

                CoordinateBasisTransform!(e_plot' * η * Uₐ,Ub,pos_plot,coordinates,sim_coordinates) # transform U_a to simulation frame
                Ub = η * inve_sim' * Ub
                eng_den = ScalarEnergyDensity(Tᵃᵇ,Ub,num_den;perparticle=perparticle)
                eng_den *= engden_unit_scale # convert to specified units

            end

            num[x_plot_idx,y_plot_idx] += num_den
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

        replace!(num,0.0=>NaN)
        if log_num_dens
            num = log10.(num)
            replace!(num, -Inf => NaN)
        end

        replace!(eng,0.0=>NaN)
        if log_eng_dens
            eng = log10.(eng)
            replace!(eng, -Inf => NaN)
        end

        # number plot
        hm = voronoiplot!(axnumpho,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       

        z_lab = log_num_dens ? L"\log_{10}\left(n\, \left[ %$(numden_unit_string)\right] \right)" : L"n\, \left[ %$(numden_unit_string) \right]"

        Colorbar(fig[1,1],hm,label=z_lab,flipaxis=false,ticks = 18:1:21)

        # energy plot
        hm = voronoiplot!(axengpho,x_points,y_points,eng,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1) 
        
        z_lab = log_eng_dens ? L"\log_{10}\left(e\, \left[ %$(engden_unit_string)\right] \right)" : L"e\, \left[ %$(engden_unit_string) \right]"

        Colorbar(fig[2,1],hm,label=z_lab,flipaxis=false,ticks = -5:1:-1)

    # Electron plots

        num = zeros(Float64,nx_dir,ny_dir)
        eng = zeros(Float64,nx_dir,ny_dir)
        x_points = zeros(Float64,nx_dir,ny_dir)
        y_points = zeros(Float64,nx_dir,ny_dir)

        species_index = findfirst(x->x=="Ele",name_list)

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
                TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
                TetradComponents!(pos_plot,e_plot,metric,coordinates,StationaryObserverTetrad())
                InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,StationaryObserverTetrad())

                CoordinateBasisTransform!(e_sim' * Nᵇ,Nᵇ,pos_sim,sim_coordinates,coordinates) # transform vector from one coordinate basis to another 
                Nᵃ = inve_plot' * Nᵇ # transform to stationary tetrad frame
                num_den = ScalarNumberDensity(Nᵃ,Uₐ) * numden_unit_scale # convert to specified units

                CoordinateBasisTransform!(e_plot' * η * Uₐ,Ub,pos_plot,coordinates,sim_coordinates) # transform U_a to simulation frame
                Ub = η * inve_sim' * Ub
                eng_den = ScalarEnergyDensity(Tᵃᵇ,Ub,num_den;perparticle=perparticle)
                eng_den *= engden_unit_scale # convert to specified units

            end

            num[x_plot_idx,y_plot_idx] += num_den
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

        replace!(num,0.0=>NaN)
        if log_num_dens
            num = log10.(num)
            replace!(num, -Inf => NaN)
        end

        replace!(eng,0.0=>NaN)
        if log_eng_dens
            eng = log10.(eng)
            replace!(eng, -Inf => NaN)
        end

        # number plot
        hm = voronoiplot!(axnumele,x_points,y_points,num,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1)       

        z_lab = log_num_dens ? L"\log_{10}\left(n\, \left[ %$(numden_unit_string)\right] \right)" : L"n\, \left[ %$(numden_unit_string) \right]"

        Colorbar(fig[1,4],hm,label=z_lab,ticks = 4:2:10)

        # energy plot
        hm = voronoiplot!(axengele,x_points,y_points,eng,colormap=theme.colormap_var,show_generators=false,strokewidth=0.1) 
        
        z_lab = log_eng_dens ? L"\log_{10}\left(e\, \left[ %$(engden_unit_string)\right] \right)" : L"e\, \left[ %$(engden_unit_string) \right]"

        Colorbar(fig[2,4],hm,label=z_lab,ticks = -8:2:-2)

    # design axis layout 

    colgap!(fig.layout,1, Relative(0.02))
    colgap!(fig.layout,2, 0.0)
    colgap!(fig.layout,3, Relative(0.02))
    rowgap!(fig.layout,1, 0.0)

    #resize_to_layout!(fig)

    return fig

    end # with_theme

end