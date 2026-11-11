############## Thesis Plots #################
#############################################

    ################# Plots that only work with Paraboloidal coordinates ############### 

    """
        EquatorialNumberDensityPlot1D(type,PhaseSpace,sol,species;...)

    Returns a plot of the number density of a `species` along the equatorial plane (z=0 in cylindrical coordinates) as a function of time.
    """
    function EquatorialNumberDensityPlot1D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,wide=false,plot_limits=(nothing,nothing),legend::Bool=true)

        @assert PhaseSpace.Spacetime.coordinates isa DiplodocusTransport.Paraboloidal "This plot is only designed currently for Paraboloidal coordinates"

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        name_list = PhaseSpace.name_list
        Momentum = PhaseSpace.Momentum
        Spacetime = PhaseSpace.Spacetime
        Grids = PhaseSpace.Grids
        x_num = Spacetime.x_num
        y_num = Spacetime.y_num
        z_num = Spacetime.z_num
        dx = Grids.dx
        dy = Grids.dy
        dz = Grids.dz
        mx = Grids.mx
        my = Grids.my
        mz = Grids.mz

        CHAR_n = PhaseSpace.Characteristic.CHAR_number_density

        if wide
            fig = Figure(size=(576,234)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        ylab = L"$\log_{10}\left(n\, \left[\text{m}^{-3}\right]\right)$"

        xlab = L"$\log_{10}\left(\rho\, \left[r_\text{s}\right]\right)$"

        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
        ax.limits = plot_limits

        linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

        legend_elements = []
        line_labels = []

        num = zeros(Float64,y_num)
        pos = zeros(Float64,y_num)

        for (species_idx, species_name) in enumerate(species) 

            species_index = findfirst(x->x==species_name,name_list)

            p_num = Momentum.px_num_list[species_index]  
            u_num = Momentum.py_num_list[species_index]
            h_num = Momentum.pz_num_list[species_index]
            pr = Grids.pxr_list[species_index]
            ur = Grids.pyr_list[species_index]
            hr = Grids.pzr_list[species_index]
            mass = Grids.mass_list[species_index]

            f1D = zeros(Float32,p_num*u_num*h_num)

            for t_idx in [1,length(sol.t)] # only plot initial and final time steps

                if t_idx == 1 
                    linestyle_idx = 1
                elseif t_idx == length(sol.t)
                    linestyle_idx = 2
                end

                fill!(num,NaN)
                fill!(pos,NaN)

            for space_idx in 1:y_num

                ρ = my[space_idx] * mz[space_idx] # radial coordinate

                Nᵃ = Float64.(DiplodocusTransport.FourFlow(sol.f[t_idx],PhaseSpace,species_index;x_idx=1,y_idx=space_idx,z_idx=space_idx))
                #Ua = HydroFourVelocity(Na)
                Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                n = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ) * CHAR_n

                println("n = $n at space_idx = $space_idx, ρ = $ρ")

                num[space_idx] = n == 0.0 ? NaN : log10(n)
                pos[space_idx] = log10(ρ)

            end
            println(num)
            println(pos)

            replace!(num,0.0=>NaN)

            scatterlines!(ax,pos,num,color = theme.textcolor[],markersize=5.0,linestyle=linestyles[linestyle_idx])

            push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[linestyle_idx],linewidth = 2.0))
            push!(line_labels,L"t=%$(sol.t[t_idx])")
            
            end # time loop
        
        end # species loop 

        # analytic solution for comparison 

        # Fiducual RIAF parameters of the M87 accretion flow
        ne0 = 1e10 # m^-3
        Te0 = 1e11 # K
        αn = 1.5
        αT = 0.84
        num_Init = zeros(Float64,y_num)
        fill!(num_Init,NaN)
        for i in 3:y_num
            ρ = PhaseSpace.Grids.yr[i]*PhaseSpace.Grids.zr[i]
            num_Init[i] = log10(ne0 * (ρ/2)^(-αn))
        end
        scatterlines!(ax,pos,num_Init,color = theme.textcolor[],markersize=5.0,linestyle=linestyles[3])
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[3],linewidth = 2.0))
        push!(line_labels,L"\text{Eq.}\,5.8")

        
        if legend
            axislegend(ax,legend_elements,line_labels,position = :rt)
        end

        end # with_theme

        return fig

    end

    """
        EquatorialTemperaturePlot1D(type,PhaseSpace,sol,species;...)

    Returns a plot of the temperature of the momentum distribution of a `species` along the equatorial plane (z=0 in cylindrical coordinates) as a function of time.
    """
    function EquatorialTemperaturePlot1D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,wide=false,plot_limits=(nothing,nothing),legend::Bool=true)

        @assert PhaseSpace.Spacetime.coordinates isa DiplodocusTransport.Paraboloidal "This plot is only designed currently for Paraboloidal coordinates"

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        name_list = PhaseSpace.name_list
        Momentum = PhaseSpace.Momentum
        Spacetime = PhaseSpace.Spacetime
        Grids = PhaseSpace.Grids
        x_num = Spacetime.x_num
        y_num = Spacetime.y_num
        z_num = Spacetime.z_num
        dx = Grids.dx
        dy = Grids.dy
        dz = Grids.dz
        mx = Grids.mx
        my = Grids.my
        mz = Grids.mz

        CHAR_n = PhaseSpace.Characteristic.CHAR_number_density

        if wide
            fig = Figure(size=(576,234)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        ylab = L"$\log_{10}\left(T\, \left[\text{K}\right]\right)$"

        xlab = L"$\log_{10}\left(\rho\, \left[r_\text{s}\right]\right)$"

        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
        ax.limits = plot_limits

        linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

        legend_elements = []
        line_labels = []

        Temp = zeros(Float64,y_num)
        pos = zeros(Float64,y_num)

        for (species_idx, species_name) in enumerate(species) 

            species_index = findfirst(x->x==species_name,name_list)

            p_num = Momentum.px_num_list[species_index]  
            u_num = Momentum.py_num_list[species_index]
            h_num = Momentum.pz_num_list[species_index]
            pr = Grids.pxr_list[species_index]
            ur = Grids.pyr_list[species_index]
            hr = Grids.pzr_list[species_index]
            mass = Grids.mass_list[species_index]

            f1D = zeros(Float32,p_num*u_num*h_num)

            for t_idx in [1,length(sol.t)] # only plot initial and final time steps

                if t_idx == 1 
                    linestyle_idx = 1
                elseif t_idx == length(sol.t)
                    linestyle_idx = 2
                end

                fill!(Temp,NaN)
                fill!(pos,NaN)

            for space_idx in 1:y_num

                ρ = my[space_idx] * mz[space_idx] # radial coordinate

                Nᵃ = Float64.(DiplodocusTransport.FourFlow(sol.f[t_idx],PhaseSpace,species_index;x_idx=1,y_idx=space_idx,z_idx=space_idx))
                #Ua = HydroFourVelocity(Na)
                Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ) * CHAR_n
                Δab = ProjectionTensor(Uₐ)
                Tᵃᵇ = Float64.(StressEnergyTensor(sol.f[t_idx],PhaseSpace,species_index;x_idx=1,y_idx=space_idx,z_idx=space_idx))
                Pressure = ScalarPressure(Tᵃᵇ,Δab)
                Temperature = ScalarTemperature(Pressure,num)

                println("T = $Temperature at space_idx = $space_idx, ρ = $ρ")

                Temp[space_idx] = Temperature == 0.0 ? NaN : log10(Temperature)
                pos[space_idx] = log10(ρ)

            end
            println(Temp)
            println(pos)

            replace!(Temp,0.0=>NaN)

            scatterlines!(ax,pos,Temp,color = theme.textcolor[],markersize=5.0,linestyle=linestyles[linestyle_idx])

            push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[linestyle_idx],linewidth = 2.0))
            push!(line_labels,L"t=%$(sol.t[t_idx])")
            
            end # time loop
        
        end # species loop 

        # analytic solution for comparison 

        # Fiducual RIAF parameters of the M87 accretion flow
        ne0 = 1e10 # m^-3
        Te0 = 1e11 # K
        αn = 1.5
        αT = 0.84
        T_Init = zeros(Float64,y_num)
        fill!(T_Init,NaN)
        for i in 3:y_num
            ρ = PhaseSpace.Grids.yr[i]*PhaseSpace.Grids.zr[i]
            T_Init[i] = log10(Te0 * (ρ/2)^(-αT))
        end
        scatterlines!(ax,pos,T_Init,color = theme.textcolor[],markersize=5.0,linestyle=linestyles[3])
        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[3],linewidth = 2.0))
        push!(line_labels,L"\text{Eqs.}\,5.8\text{ and }5.9")

        
        if legend
            axislegend(ax,legend_elements,line_labels,position = :rt)
        end

        end # with_theme

        return fig

    end

    ####### 3 solution Bulk Gamma comparisons #######

    function BulkGamma3SolPlot2D(type::Static,PhaseSpaces::Tuple{PhaseSpaceStruct,PhaseSpaceStruct,PhaseSpaceStruct},sols::Tuple{OutputStruct,OutputStruct,OutputStruct},species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64},time_idx::Int64;fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing),gamma_sign::Bool=true,gamma_z::Bool=false,obs_frame::Bool=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    surface_string = surface[1]
    surface_index = surface[2]

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    PhaseSpace = PhaseSpaces[1]

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

    ax1 = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax1.limits = plot_limits
    hidexdecorations!(ax1)

    ax2 = Axis(fig[2,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax2.limits = plot_limits
    hidexdecorations!(ax2)
    hideydecorations!(ax2)

    ax3 = Axis(fig[3,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax3.limits = plot_limits
    hideydecorations!(ax3)

    Gamma1 = zeros(Float64,nx_dir,ny_dir)
    Gamma2 = zeros(Float64,nx_dir,ny_dir)
    Gamma3 = zeros(Float64,nx_dir,ny_dir)
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

    for sol_idx in 1:3

        PhaseSpace = PhaseSpaces[sol_idx]
        sol = sols[sol_idx]

        Grids = PhaseSpace.Grids
        Spacetime = PhaseSpace.Spacetime
        metric = Spacetime.metric
        sim_coordinates = Spacetime.coordinates
        tetrad = Spacetime.tetrad

        if sol_idx == 1
            Gamma = @view(Gamma1[:,:])
        elseif sol_idx == 2
            Gamma = @view(Gamma2[:,:])
        elseif sol_idx == 3
            Gamma = @view(Gamma3[:,:])
        end

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


    end # sols

    colrange = (minimum([minimum(x for x in Gamma1 if isfinite(x)),minimum(x for x in Gamma2 if isfinite(x)),minimum(x for x in Gamma3 if isfinite(x))]),maximum([maximum(x for x in Gamma1 if isfinite(x)),maximum(x for x in Gamma2 if isfinite(x)),maximum(x for x in Gamma3 if isfinite(x))]))

    println("colrange = $colrange")
    println("Gamma1 = $Gamma1")
    println("Gamma2 = $Gamma2")

    println("Gamma3 = $Gamma3")

    #hm = heatmap!(ax,xr_dir,yr_dir,Gamma,colormap=theme.colormap_var)  
    hm1 = voronoiplot!(ax1,x_points,y_points,Gamma1,colormap=theme.colormap_var,show_generators=false,strokewidth=0.0,colorrange=colrange)       
    hm2 = voronoiplot!(ax2,x_points,y_points,Gamma2,colormap=theme.colormap_var,show_generators=false,strokewidth=0.0,colorrange=colrange)       
    hm3 = voronoiplot!(ax3,x_points,y_points,Gamma3,colormap=theme.colormap_var,show_generators=false,strokewidth=0.0,colorrange=colrange)       


    Colorbar(fig[1,2],hm1,label= gamma_z ? L"$\Gamma_z\,$" : L"$\Gamma\,$")

    rowgap!(fig.layout,1,0.0)
    rowgap!(fig.layout,2,0.0)

    end # with_theme

    return fig

    end

    function BulkGammaDecomp3SolPlot2D(type::Static,PhaseSpaces::Tuple{PhaseSpaceStruct,PhaseSpaceStruct,PhaseSpaceStruct},sols::Tuple{OutputStruct,OutputStruct,OutputStruct},species::String,coordinates::AbstractCoordinates,surface::Tuple{String,Int64},time_idx::Int64;fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,wide=false,plot_limits=(nothing,nothing),gamma_sign::Bool=true,gamma_z::Bool=false,obs_frame::Bool=false,figsize::Tuple{Int64,Int64}=(576,216),x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    surface_string = surface[1]
    surface_index = surface[2]

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    PhaseSpace = PhaseSpaces[1]

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

    xlab = L"$z\,\left[r_\text{s}\right]$"
    ylab1 = L"$\text{Lorentz factor}$"
    ylab2 = L"$\xi$"

    ax1_j = Axis(fig[1,1],xlabel=xlab,ylabel=ylab1,ygridvisible=false)
    ax2_j = Axis(fig[1,1],xlabel=xlab,ylabel=ylab2,yaxisposition=:right,ygridvisible=false)
    hidexdecorations!(ax1_j,grid=false)
    hidexdecorations!(ax2_j)
    ax1_j.limits = ((0,80),(1,nothing))
    ax2_j.limits = ((0,80),(0,1))

    ax1_s = Axis(fig[2,1],xlabel=xlab,ylabel=ylab1,ygridvisible=false)
    ax2_s = Axis(fig[2,1],xlabel=xlab,ylabel=ylab2,yaxisposition=:right,ygridvisible=false)
    ax1_s.limits = ((0,80),(1,nothing))
    ax2_s.limits = ((0,80),(0,1))

    Gamma = zeros(Float64,nx_dir)
    gamma = zeros(Float64,nx_dir)
    xi = zeros(Float64,nx_dir)
    z_points = zeros(Float64,nx_dir)

    species_index = findfirst(x->x==species,name_list)

    e_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_plot::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    e_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    inve_sim::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
    pos_sim::MVector{4,Float64} = zeros(Float64,4)
    pos_plot::MVector{4,Float64} = zeros(Float64,4)
    g::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer

    legend_elements_color = []
    line_labels_color = []
    legend_elements_style = []
    line_labels_style = []

    a_vals = [0.0,0.5,0.9]

    for sol_idx in 1:3 

        color = theme.colormap[][(sol_idx-1)/(3-1)]
        pushfirst!(legend_elements_color,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        pushfirst!(line_labels_color,L"$a=%$(round(a_vals[sol_idx],sigdigits=2))$")

    end

    pushfirst!(legend_elements_style,LineElement(color = :black, linestyle = :solid,linewidth = 2.0))
    pushfirst!(line_labels_style,L"$\Gamma$")
    pushfirst!(legend_elements_style,LineElement(color = :black, linestyle = :dash,linewidth = 2.0))
    pushfirst!(line_labels_style,L"$\gamma$")
    pushfirst!(legend_elements_style,LineElement(color = :black, linestyle = :dot,linewidth = 2.0))
    pushfirst!(line_labels_style,L"$\xi$")

    for y_plot_idx in [2,3]
    for sol_idx in 1:3

        PhaseSpace = PhaseSpaces[sol_idx]
        sol = sols[sol_idx]

        Grids = PhaseSpace.Grids
        Spacetime = PhaseSpace.Spacetime
        metric = Spacetime.metric
        sim_coordinates = Spacetime.coordinates
        tetrad = Spacetime.tetrad

        color = theme.colormap[][(sol_idx-1)/(3-1)]

        t = sol.t[time_idx]

        pos_sim[1] = t

    for x_plot_idx in 1:nx_dir

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
        Nᵃ_sim_obs = Float64[1.0,0.0,0.0,0.0]  # for testing
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

        
        if n != 0.0
            Vᵇ_obs = HydroFourVelocity(Nᵃ_sim_obs)
            CoordinateBasisTransform!(e_sim' * Vᵇ_obs,Vᵇ_obs,pos_sim,sim_coordinates,coordinates)
            Vᵃ = inve_plot' * Vᵇ_obs # transform to stationary tetrad frame
            if gamma_z # static observed z component of bulk gamma
                if gamma_sign
                    βz = Vᵃ[4]/Vᵃ[1]
                    γ = 1/sqrt(1 - βz^2) * sign(Vᵃ[4]) # bulk gamma with sign of velocity to distinguish inflow vs outflow
                else
                    βz = Vᵃ[4]/Vᵃ[1]
                    γ = 1/sqrt(1 - βz^2) # bulk gamma without sign, so all flow contributes positively
                end
            else
                if gamma_sign
                    γ = (-transpose(Uₐ) * Vᵃ) * sign(Vᵃ[4]) # bulk gamma with sign of velocity to distinguish inflow vs outflow
                else
                    γ = -transpose(Uₐ) * Vᵃ # bulk gamma without sign, so all flow contributes positively
                end
            end
        else
            γ = NaN
        end

        println("γ = $γ at (x,y,z) = ($x,$y,$z)")

        Gamma[x_plot_idx] = Γ
        gamma[x_plot_idx] = γ
        xi[x_plot_idx] = sqrt(1-(γ/Γ)^2)

        z_points[x_plot_idx] = pos_plot[4]

    end

        if y_plot_idx == 2
            scatterlines!(ax1_j,z_points,Gamma,color=color,linewidth=2.0,markersize=0.0)
            scatterlines!(ax1_j,z_points,gamma,color=color,linewidth=2.0,markersize=0.0,linestyle=:dash)
            scatterlines!(ax2_j,z_points,xi,color=color,linewidth=2.0,markersize=0.0,linestyle=:dot)
        elseif y_plot_idx == 3
            scatterlines!(ax1_s,z_points,Gamma,color=color,linewidth=2.0,markersize=0.0)
            scatterlines!(ax1_s,z_points,gamma,color=color,linewidth=2.0,markersize=0.0,linestyle=:dash)
            scatterlines!(ax2_s,z_points,xi,color=color,linewidth=2.0,markersize=0.0,linestyle=:dot)
        end

    end # sols
    end # x_idx 

    axislegend(ax1_j,legend_elements_color,line_labels_color,position=:lt,framevisible=true,rowgap=-8,nbanks=1)
    axislegend(ax1_j,legend_elements_style,line_labels_style,position=(0.3,1.0),framevisible=true,rowgap=-8,nbanks=1)

    rowgap!(fig.layout,1,0.0)

    end # with_theme

    return fig

    end

    ################ Jet Loading Plot ########################
    ##########################################################

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