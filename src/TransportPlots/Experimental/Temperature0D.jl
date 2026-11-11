"""
    function TemperaturePlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the temperature associated with each particles distribution function as a function of time.  
"""
function ThermalPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    t_grid = Time.t_grid

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    Temperature = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()

    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"$T$ $[\mathrm{K}]$"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"$T$ $[\mathrm{K}]$"
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

    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            name = name_list[j]
            mass = mass_list[j]

            f = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j))
            Nᵃ = DiplodocusTransport.FourFlow(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = DiplodocusTransport.ProjectionTensor(Uₐ)
            Tᵃᵇ = DiplodocusTransport.StressEnergyTensor(f,p_num_list[j],u_num_list[j],h_num_list[j],pr_list[j],ur_list[j],hr_list[j],mass)
            Pressure = DiplodocusTransport.ScalarPressure(Tᵃᵇ,Δab)
            Temperature[i] = DiplodocusTransport.ScalarTemperature(Pressure,num)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),Temperature,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),Temperature,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list_plot[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    ylims!(ax,0.0,nothing)

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,position=:lb)

    end # with_theme

    return fig

end