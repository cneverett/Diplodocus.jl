"""
    VelocityZPlot1D(type,sol,PhaseSpace,direction,species;...)

Returns a plot of the velocity in the local z-direction of all species as a function of time along a specific spatial coordinate `direction`.
"""
function VelocityZPlot1D(type::Static,sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,direction::String,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,step::Int=1,wide=false,logt::Bool=false,plot_limits=(nothing,nothing))

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

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    ylab = L"$v_z\, [\text{Code Units}]$"

    if direction == "x"
        xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        n_dir = x_num
        grid_type = Spacetime.x_grid
        if grid_type == "u"
            mdir_plot = mx
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
        end
    elseif direction == "y"
        xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        n_dir = y_num
        grid_type = Spacetime.y_grid
        if grid_type == "u"
            mdir_plot = my
        elseif grid_type == "l"
            mdir_plot = log10.(my)
        end
    elseif direction == "z"
        xlab = L"$\left(z\,[\text{Code Units}]\right)$"
        n_dir = z_num
        grid_type = Spacetime.z_grid
        if grid_type == "u"
            mdir_plot = mz
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
        end
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    vz = zeros(Float64,n_dir)

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

        t_min = logt ? sol.t[2]/10 : sol.t[1]
        t_max = sol.t[end]

        for i in 1:t_save

            if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

                t = sol.t[i]
                #println("t=$(CodeToSIUnitsTime(t))")
                if logt
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                else
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                end

                # sum over other 2 spatial directions to get average number density along desired direction
                fill!(vz,0.0)

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    Nᵃ = FourFlow(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)
                    HydroFourVelocity = HydroFourVelocity(Nᵃ)

                    Vz_local = HydroFourVelocity[4]/HydroFourVelocity[1]

                    #println("n = $n at (x,y,z) = ($x,$y,$z)")

                    if direction == "x"
                        vz[x] += Vz_local * dy[y] * dz[z] / (sum(dy)*sum(dz))
                        #if species== "All"
                        #    num_total[x] += num[x]
                        #end
                    elseif direction == "y"
                        vz[y] += Vz_local * dx[x] * dz[z] / (sum(dx)*sum(dz))
                        #if species== "All"
                        #    num_total[y] += num[y]
                        #end
                    elseif direction == "z"
                        vz[z] += Vz_local * dx[x] * dy[y] / (sum(dx)*sum(dy))
                        #if species== "All"
                        #    num_total[z] += num[z]
                        #end
                    end

                end

                replace!(vz,0.0=>NaN)

                scatterlines!(ax,mdir_plot,vz,color = color,markersize=5.0)

            end # if times should be plotted
        
        end # time loop
       
    end # species loop 

    t_unit_string = TimeUnits()

    if logt
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(TimeUnits(sol.t[1]),sigdigits=5)),log10(round(TimeUnits(sol.t[end]),sigdigits=5))),label=L"$\log_{10}\left(t\,%$t_unit_string \right)$")
    else
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t\,$ $%$t_unit_string$")
    end
    
    #axislegend(ax,position = :rb)

    end # with_theme

    return fig

end


"""
    PressurePlot1D(type,sol,PhaseSpace,direction,species;...)

Returns a plot of the pressure of all `species` as a function of time along a specific spatial coordinate `direction`.
"""
function PressurePlot1D(type::Static,sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,direction::String,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,step::Int=1,wide=false,logt::Bool=false,plot_limits=(nothing,nothing))

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

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    ylab = L"$P\, [\text{Code Units}]$"

    if direction == "x"
        xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        n_dir = x_num
        grid_type = Space.x_grid
        if grid_type == "u"
            mdir_plot = mx
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
        end
    elseif direction == "y"
        xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        n_dir = y_num
        grid_type = Space.y_grid
        if grid_type == "u"
            mdir_plot = my
        elseif grid_type == "l"
            mdir_plot = log10.(my)
        end
    elseif direction == "z"
        xlab = L"$\left(z\,[\text{Code Units}]\right)$"
        n_dir = z_num
        grid_type = Space.z_grid
        if grid_type == "u"
            mdir_plot = mz
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
        end
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    P = zeros(Float64,n_dir)

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

        t_min = logt ? sol.t[2]/10 : sol.t[1]
        t_max = sol.t[end]

        for i in 1:t_save

            if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

                t = sol.t[i]
                #println("t=$(CodeToSIUnitsTime(t))")
                if logt
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                else
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                end

                # sum over other 2 spatial directions to get average number density along desired direction
                fill!(P,0.0)

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                    Δab = DiplodocusTransport.ProjectionTensor(Uₐ)

                    Tab = DiplodocusTransport.StressEnergyTensor(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)

                    P_local = DiplodocusTransport.ScalarPressure(Tab,Δab)

                    #println("n = $n at (x,y,z) = ($x,$y,$z)")

                    if direction == "x"
                        P[x] += P_local * dy[y] * dz[z] / (sum(dy)*sum(dz))
                        #if species== "All"
                        #    num_total[x] += num[x]
                        #end
                    elseif direction == "y"
                        P[y] += P_local * dx[x] * dz[z] / (sum(dx)*sum(dz))
                        #if species== "All"
                        #    num_total[y] += num[y]
                        #end
                    elseif direction == "z"
                        P[z] += P_local * dx[x] * dy[y] / (sum(dx)*sum(dy))
                        #if species== "All"
                        #    num_total[z] += num[z]
                        #end
                    end

                end

                replace!(P,0.0=>NaN)

                scatterlines!(ax,mdir_plot,P,color = color,markersize=5.0)

            end # if times should be plotted
        
        end # time loop
       
    end # species loop 

    t_unit_string = TimeUnits()

    if logt
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(TimeUnits(sol.t[1]),sigdigits=5)),log10(round(TimeUnits(sol.t[end]),sigdigits=5))),label=L"$\log_{10}\left(t\,%$t_unit_string \right)$")
    else
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t\,$ $%$t_unit_string$")
    end
    
    #axislegend(ax,position = :rb)

    end # with_theme

    return fig

end

"""
    TemperaturePlot1D(type,sol,PhaseSpace,direction,species;...)

Returns a plot of the temperature of all `species` as a function of time along a specific spatial coordinate `direction`.
"""
function TemperaturePlot1D(type::Static,sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,direction::String,species::Vector{String};fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,step::Int=1,wide=false,logt::Bool=false,plot_limits=(nothing,nothing))

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

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    ylab = L"$T\, [\text{K}]$"

    if direction == "x"
        xlab = L"$\left(x\,[\text{Code Units}]\right)$"
        n_dir = x_num
        grid_type = Space.x_grid
        if grid_type == "u"
            mdir_plot = mx
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
        end
    elseif direction == "y"
        xlab = L"$\left(y\,[\text{Code Units}]\right)$"
        n_dir = y_num
        grid_type = Space.y_grid
        if grid_type == "u"
            mdir_plot = my
        elseif grid_type == "l"
            mdir_plot = log10.(my)
        end
    elseif direction == "z"
        xlab = L"$\left(z\,[\text{Code Units}]\right)$"
        n_dir = z_num
        grid_type = Space.z_grid
        if grid_type == "u"
            mdir_plot = mz
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
        end
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    T = zeros(Float64,n_dir)

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

        t_min = logt ? sol.t[2]/10 : sol.t[1]
        t_max = sol.t[end]

        for i in 1:t_save

            if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

                t = sol.t[i]
                #println("t=$(CodeToSIUnitsTime(t))")
                if logt
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                else
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                end

                # sum over other 2 spatial directions to get average number density along desired direction
                fill!(T,0.0)

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                    Δab = DiplodocusTransport.ProjectionTensor(Uₐ)

                    FourFlow = DiplodocusTransport.FourFlow(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)
                    Tab = DiplodocusTransport.StressEnergyTensor(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)

                    n = DiplodocusTransport.ScalarNumberDensity(FourFlow,Uₐ)
                    P = DiplodocusTransport.ScalarPressure(Tab,Δab)
                    T_local = DiplodocusTransport.ScalarTemperature(P,n)

                    #println("n = $n at (x,y,z) = ($x,$y,$z)")

                    if direction == "x"
                        T[x] += T_local * dy[y] * dz[z] / (sum(dy)*sum(dz))
                        #if species== "All"
                        #    num_total[x] += num[x]
                        #end
                    elseif direction == "y"
                        T[y] += T_local * dx[x] * dz[z] / (sum(dx)*sum(dz))
                        #if species== "All"
                        #    num_total[y] += num[y]
                        #end
                    elseif direction == "z"
                        T[z] += T_local * dx[x] * dy[y] / (sum(dx)*sum(dy))
                        #if species== "All"
                        #    num_total[z] += num[z]
                        #end
                    end

                end

                replace!(T,0.0=>NaN)

                scatterlines!(ax,mdir_plot,T,color = color,markersize=5.0)

            end # if times should be plotted
        
        end # time loop
       
    end # species loop 

    t_unit_string = TimeUnits()

    if logt
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(TimeUnits(sol.t[1]),sigdigits=5)),log10(round(TimeUnits(sol.t[end]),sigdigits=5))),label=L"$\log_{10}\left(t\,%$t_unit_string \right)$")
    else
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(sol.t[1]),TimeUnits(sol.t[end])),label=L"$t\,$ $%$t_unit_string$")
    end
    
    #axislegend(ax,position = :rb)

    end # with_theme

    return fig

end