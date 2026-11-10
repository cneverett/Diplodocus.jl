"""
    EnergyDensityPlot0D(type,PhaseSpace,sol,species;...)

Returns a plot of the energy density of a `species`` as a function of time.

Arguments:
- `type::PlotType` determines the type of plot to generate, it can be either `Static` or `Animated`.
- `PhaseSpace::PhaseSpaceStruct` is the structure containing the phase space information of the simulation.
- `sol::OutputStruct` is the solution object containing the distribution function of all particles and time stepping information of the simulation. 
- `species::String` is the abbreviated three letter name of the particle species to plot. Can be set to `"All"` to plot all species.

Common keyword arguments:
- `theme=DiplodocusDark()`: the colour theme to use for the plot, can be either `DiplodocusDark()` for dark mode or `DiplodocusLight()` for light mode.
- `TimeUnits::Tuple{Float64,String}=(1.0, "\\text{Code Units}")`: a tuple that converts the time given in code units to the desired units for plotting. The first entry is the conversion factor and the second is a string that will be converted into a LaTeX string for the time label.
- `frac=false`: Whether to plot the fractional difference of number density between time-steps.
- `logt=false`: If `true` time will be converted to `log10` format for plotting.
- `loge=false`: If `true` energy density will be converted to `log10` format for plotting.
- `title=nothing`: A string to place as the figure title, no title by default.
- `legend=true`: Whether to show the legend, `true` by default.
- `legend_pos=:rb`: The position of the legend, default is `:rb` (right bottom).
- `fig=nothing` is the figure to plot on, if `nothing` a new figure will be created.
- `x_idx::Int64=1`, the x index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `y_idx::Int64=1`, the y index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `z_idx::Int64=1`, the z index of the spatial grid cell that you want to plot the distribution for, default is 1.
"""
function EnergyDensityPlot0D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String;fig=nothing,theme=DiplodocusDark(),title=nothing,legend::Bool=true,legend_pos::Symbol=:rb,plot_limits=(nothing,nothing),frac::Bool=false,only_all::Bool=false,perparticle::Bool=false,logt::Bool=false,loge::Bool=false,TimeUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),x_idx=1,y_idx=1,z_idx=1)

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

    Characteristic = PhaseSpace.Characteristic

    CHAR_number_density = Characteristic.CHAR_number_density

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    eng = zeros(Float64,length(sol.t))
    eng_total = zeros(Float64,length(sol.t))
    val = 0.0
    val_prev = 0.0

    t_unit_string = TimeUnits[2]
    t_unit_scale = TimeUnits[1]

    if logt
        xlab = L"\log_{10}($t\, [%$t_unit_string]$)"
    else
        xlab = L"$t\, [%$t_unit_string$]"
    end

    if loge && frac
        error("loge and frac cannot both be true, please choose one or the other.")
    end

    if perparticle
        if frac
            ylab = L"$\overline{p^0}\text{ Frac. Change}$"
        elseif loge
            ylab = L"\log_{10}\left(\overline{p^0} $[m_ec]$\right)"
        else
            ylab = L"$\overline{p^0}\, [m_ec]$"
        end
    else
        if frac
            ylab = L"$\text{Eng. Den. Frac. Change}$"
        elseif loge
            ylab = L"\log_{10}\left(e/c $[m_ec~\mathrm{m}^{-3}]$\right) "
        else
            ylab = L"$e/c$ $[m_ec\,\mathrm{m}^{-3}]$"
        end
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    end

    if !isnothing(title)
        ax.title = title
    end

    Uₐ::Vector{eltype(sol.t[1])} = [-1.0,0.0,0.0,0.0] # static observer

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            Nᵃ = FourFlow(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            num = ScalarNumberDensity(Nᵃ,Uₐ)

            Tᵃᵇ = StressEnergyTensor(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx) 

            val_prev = val
            val = ScalarEnergyDensity(Tᵃᵇ,Uₐ,num,perparticle=perparticle)
            val *= perparticle ? 1.0 : CHAR_number_density

            if frac
                if i == 1
                    eng[i] = 0.0 # initial value
                else
                    eng[i] = val/val_prev - 1.0
                end
            else
                if loge && (val == 0.0)
                    val = NaN
                end
                eng[i] = val
            end

            if species== "All"
                eng_total[i] += eng[i] == NaN ? 0.0 : eng[i]
            end
        
        end

        if loge 
            eng = log10.(eng)
        end

        if !only_all
        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,log10(t_unit_scale*sol.t[1]),log10(t_unit_scale*sol.t[end]))
            
        else
            scatterlines!(ax,sol.t .* t_unit_scale,eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,t_unit_scale*sol.t[1],t_unit_scale*sol.t[end])
        end
        end

    end

    if loge
        eng_total = log10.(eng_total)
    end

    if species == "All"
        if frac
            for i in length(eng_total):-1:2
                eng_total[i] = eng_total[i]/eng_total[i-1] - 1.0
            end
            eng_total[1] = 0.0
        end
        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(t_unit_scale*sol.t[1]),log10(t_unit_scale*sol.t[end]))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,t_unit_scale*sol.t[1],t_unit_scale*sol.t[end])
        end
    end

    if legend
        axislegend(ax,position = legend_pos)
    end

    end # with_theme

    return fig

end



#="""
    FracEnergyDensityPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark(),title=nothing,only_all=false)

Returns a plot of the fractional change in energy density of all species between time setups as a function of time.
"""
function FracEnergyDensityPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,only_all=false,TimeUnits::Function=CodeToCodeUnitsTime,logt::Bool=false,legend::Bool=true,perparticle=false,legend_pos::Symbol=:lc,plot_limits=(nothing,nothing))

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

    if t_grid == "l" || logt
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    elseif t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    end

    if perparticle
        ylab = "Mean Eng. Per Particle. Frac. Change"
    else
        ylab = "Eng. Den. Frac. Change"
    end

    if isnothing(fig)
        fig = Figure()
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false) # check units
    end

    if !isnothing(title)
        ax.title = title
    end

    ax.limits = plot_limits

    #j = findfirst(x->x==species,name_list)

    frac_eng = zeros(Float64,length(sol.t))
    frac_eng_total = zeros(Float64,length(sol.t))
    eng = 0.0
    eng_total = zeros(Float64,length(sol.t))

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f1D = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j])
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f1D,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass_list[j]) 

            if i == 1
                frac_eng[i] = 0.0 # initial value
            else
                frac_eng[i] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num,perparticle=perparticle)/eng - 1.0
            end

            eng = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num,perparticle=perparticle)

            if species== "All"
                eng_total[i] += eng
            end
        
        end

        if !only_all
            if t_grid == "u"
                t_plot = TimeUnits.(sol.t)
                if logt 
                    t_plot[1] = t_plot[2] / 10
                    t_plot = log10.(t_plot)
                end
                scatterlines!(ax,t_plot,frac_eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
                xlims!(ax,t_plot[1],t_plot[end])
            elseif t_grid == "l"
                scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_eng,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
                xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
            end
        end

    end

    if species == "All"
        for i in eachindex(sol.t)
            if i == 1
                frac_eng_total[i] = 0.0 # initial value
            else
                frac_eng_total[i] = eng_total[i]/eng_total[i-1] - 1.0
            end
        end
        if t_grid == "u"
            t_plot = copy(TimeUnits.(sol.t))
            if logt 
                t_plot[1] = t_plot[2] /10
                t_plot = log10.(t_plot)
            end
            scatterlines!(ax,t_plot,frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,t_plot[1],t_plot[end])
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=:dash,label="All")
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
    end

    #fig[1,2] = Legend(fig,ax,"Particles")
    if legend
        axislegend(ax,position = legend_pos)
    end

    end # with_theme

    return fig
end

function MulitSolEngDenPlot(sols::Vector{OutputStruct},species::Vector{String},PhaseSpaces::Vector{PhaseSpaceStruct};theme=DiplodocusDark(),TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true)

    with_theme(theme) do 

    t_unit_string = TimeUnits()

    t_grid = PhaseSpaces[1].Time.t_grid

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
    end

    ylab1 = L"$e/c$ $[m_ec\mathrm{m}^{-3}]$"
    ylab2 = L"\text{Eng. Den. Frac. Change}"

    fig = Figure(size=(616,216))
    ax1 = Axis(fig[1:2,1],xlabel=xlab,ylabel=ylab1,xgridvisible=false,ygridvisible=false) # check units
    ax2 = Axis(fig[1:2,2],xlabel=xlab,ylabel=ylab2,xgridvisible=false,ygridvisible=false) # check units

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    for sol_idx in eachindex(sols)

        sol = sols[sol_idx]
        PhaseSpace = PhaseSpaces[sol_idx]

        frac_eng = zeros(Float64,length(sol.t))
        frac_eng_total = zeros(Float64,length(sol.t))
        eng = zeros(Float64,length(sol.t))
        eng_total = zeros(Float64,length(sol.t))

        name_list = PhaseSpace.name_list
        name_list_plot = [L"\text{Electron}",L"\text{Photon}"]
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

        for name_idx in eachindex(name_list)

            for t_idx in eachindex(sol.t)

            f1D = copy(LocationSpeciesToStateVector(sol.f[t_idx],PhaseSpace,species_index=name_idx))

            Nᵃ = DiplodocusTransport.FourFlow(f1D,p_num_list[name_idx],u_num_list[name_idx],h_num_list[name_idx],pr_list[name_idx],ur_list[name_idx],hr_list[name_idx],mass_list[name_idx])
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)

            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f1D,p_num_list[name_idx],u_num_list[name_idx],h_num_list[name_idx],pr_list[name_idx],ur_list[name_idx],hr_list[name_idx],mass_list[name_idx]) 

            eng[t_idx] = DiplodocusTransport.ScalarEnergyDensity(Tᵃᵇ,Uₐ,num,perparticle=false)

            eng_total[t_idx] += eng[t_idx]

            end # t loop

            if t_grid == "u"
                scatterlines!(ax1,TimeUnits.(sol.t),eng,marker = :circle,color=theme.palette.color[][mod(2*name_idx-1,7)+1],markersize=0.0,label=name_list_plot[name_idx],linestyle = linestyles[sol_idx])
                xlims!(ax1,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
            elseif t_grid == "l"
                scatterlines!(ax1,log10.(TimeUnits.(sol.t)),eng,marker = :circle,color=theme.palette.color[][mod(2*name_idx-1,7)+1],markersize=0.0,label=name_list_plot[name_idx],linestyle = linestyles[sol_idx])
                xlims!(ax1,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
            end

        end # name loop

        if t_grid == "u"
            scatterlines!(ax1,TimeUnits.(sol.t),eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,label="All",linestyle=linestyles[sol_idx])
            xlims!(ax1,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax1,log10.(TimeUnits.(sol.t)),eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,label="All",linestyle=linestyles[sol_idx])
            xlims!(ax1,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end
        
        if sol_idx == 1

            Legend(fig[1,3],ax1,L"\text{Species}",padding=(6.0f0, 6.0f0, 6.0f0, 6.0f0))

        end

        for i in eachindex(sol.t)
            if i == 1
                frac_eng_total[i] = 0.0 # initial value
            else
                frac_eng_total[i] = eng_total[i]/eng_total[i-1] - 1.0
            end
        end

        grid_label = L"%$(2^(sol_idx-1))\times"

        if t_grid == "u"
            scatterlines!(ax2,TimeUnits.(sol.t),frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[sol_idx],label=grid_label)
            xlims!(ax2,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax2,log10.(TimeUnits.(sol.t)),frac_eng_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[sol_idx],label=grid_label)
            xlims!(ax2,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end # sols loop

    Legend(fig[2,3],ax2,L"\text{Grid Resolution}",nbanks=2,padding=(6.0f0, 6.0f0, 6.0f0, 6.0f0))

    rowsize!(fig.layout,1,Relative(0.6))

    return fig
        
    end # theme

end=#