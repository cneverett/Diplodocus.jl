"""
    NumberPlot0D(type,PhaseSpace,sol,species;...)

Returns a plot of the number of a `species` as a function of time.

Arguments:
- `type::PlotType` determines the type of plot to generate, it can be either `Static` or `Animated`.
- `PhaseSpace::PhaseSpaceStruct` is the structure containing the phase space information of the simulation.
- `sol::OutputStruct` is the solution object containing the distribution function of all particles and time stepping information of the simulation. 
- `species::String` is the abbreviated three letter name of the particle species to plot. Can be set to `"All"` to plot all species.
- `domain::String` is the spatial domain to plot. This can either be "All" to plot the total number of particles over the entire spatial domain, or "Cell" to plot the number of particles in a single spatial cell defined by the `x_idx`, `y_idx` and `z_idx` keyword arguments.

Common keyword arguments:
- `theme=DiplodocusDark()`: the colour theme to use for the plot, can be either `DiplodocusDark()` for dark mode or `DiplodocusLight()` for light mode.
- `TimeUnits::Tuple{Float64,String}=(1.0, "\\text{Code Units}")`: a tuple that converts the time given in code units to the desired units for plotting. The first entry is the conversion factor and the second is a string that will be converted into a LaTeX string for the time label.
- `frac=false`: Whether to plot the fractional difference of number density between time-steps.
- `logt=false`: If `true` time will be converted to `log10` format for plotting.
- `title=nothing`: A string to place as the figure title, no title by default.
- `legend=true`: Whether to show the legend, `true` by default.
- `legend_pos=:rb`: The position of the legend, default is `:rb` (right bottom).
- `fig=nothing` is the figure to plot on, if `nothing` a new figure will be created.
- `x_idx::Int64=1`, the x index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `y_idx::Int64=1`, the y index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `z_idx::Int64=1`, the z index of the spatial grid cell that you want to plot the distribution for, default is 1.

"""
function NumberPlot0D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String,domain::String;fig=nothing,theme=DiplodocusDark(),title=nothing,legend=false,legend_pos::Symbol=:rb,frac::Bool=false,logt::Bool=false,TimeUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),x_idx=1,y_idx=1,z_idx=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    name_list_plot = similar(name_list)
    for (idx, name) in enumerate(name_list)
        if name == "Ele"
            name = "Electron"
        elseif name == "Pho"
            name = "Photon"
        elseif name == "Pos"
            name = "Positron"
        else
            name = name
        end
        name_list_plot[idx] = name
    end

    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Spacetime = PhaseSpace.Spacetime
    Characteristic = PhaseSpace.Characteristic
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    CHAR_number_density = Characteristic.CHAR_number_density
    CHAR_length = Characteristic.CHAR_length

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    num = zeros(Float64,length(sol.t))
    num_total = zeros(Float64,length(sol.t))
    val = 0.0

    t_unit_string = TimeUnits[2]
    t_unit_scale = TimeUnits[1]

    if logt
        xlab = L"\log_{10}($t\, [%$t_unit_string]$)"
    else
        xlab = L"$t\, [%$t_unit_string]$"    
    end
    if frac
        ylab = L"\text{Num. Frac. Change}"
    else
        ylab = L"$Num.$"
    end


    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    Uₐ::Vector{Float64} = [-1.0,0.0,0.0,0.0] # static observer

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            if domain == "All"
                val_prev = val
                val = 0.0
                for x in 1:x_num, y in 1:y_num, z in 1:z_num
                    Nᵃ = Float64.(FourFlow(sol.f[i],PhaseSpace,j;x_idx=x,y_idx=y,z_idx=z))
                    num_den = ScalarNumberDensity(Nᵃ,Uₐ)
                    Vol = VolFunction(PhaseSpace,1,x,y,z)
                    println("$x,$y,$z, Vol = $Vol, num_den = $num_den")
                    val += num_den * Vol * CHAR_number_density * CHAR_length^3
                end
            elseif domain == "Cell"
                Nᵃ = Float64.(FourFlow(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
                Vol = VolFunction(PhaseSpace,1,x_idx,y_idx,z_idx)
                val_prev = val
                val = ScalarNumberDensity(Nᵃ,Uₐ) * CHAR_number_density * Vol * CHAR_length^3
            else
                error("Invalid domain argument, must be either \"All\" or \"Cell\"")
            end

            if frac 
                if i == 1
                    num[i] = 0.0 # initial value
                else
                    num[i] = val/val_prev - 1.0
                end 
            else
                num[i] = val
            end

            if species== "All"
                num_total[i] += num[i]
            end
        
        end

        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,log10(t_unit_scale*sol.t[1]),log10(t_unit_scale*sol.t[end]))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,t_unit_scale*sol.t[1],t_unit_scale*sol.t[end])
        end

    end

    if species == "All"
        if frac
            for i in length(num_total):-1:2
                num_total[i] = num_total[i]/num_total[i-1] - 1.0
            end
            num_total[1] = 0.0
        end
        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(t_unit_scale*sol.t[1]),log10(t_unit_scale*sol.t[end]))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,num_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,t_unit_scale*sol.t[1],t_unit_scale*sol.t[end])
        end
    end
    
    if legend
        axislegend(ax,position = Legend_pos)
    end

    end # with_theme

    return fig

end

#="""
    FracNumberDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing)

Returns a plot of the fractional change in number density of all species between time setups as a function of time.
"""
function FracNumberDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel="Num. Den. Frac. Change",xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel="Num. Den. Frac. Change",xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    frac_num = zeros(Float64,length(sol.t))
    num = 0.0

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0]

            if i == 1
                frac_num[i] = 0.0 # initial value
            else
                frac_num[i] = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)/num - 1.0
            end

            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_num,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position = :rb)

    end # with_theme

    return fig
end
=#
