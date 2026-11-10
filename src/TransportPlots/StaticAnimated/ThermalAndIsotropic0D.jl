"""
    function IsThermalAndIsotropicPlot0D(type,PhaseSpace,sol,species;...)

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected thermal distribution and an isotropic distribution based on the current distribution of that species as a function of time.  

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
- `title=nothing`: A string to place as the figure title, no title by default.
- `fig=nothing` is the figure to plot on, if `nothing` a new figure will be created.
- `x_idx::Int64=1`, the x index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `y_idx::Int64=1`, the y index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `z_idx::Int64=1`, the z index of the spatial grid cell that you want to plot the distribution for, default is 1.
"""
function IsThermalAndIsotropicPlot0D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::String;fig=nothing,theme=DiplodocusDark(),title=nothing,logt::Bool=false,TimeUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),x_idx=1,y_idx=1,z_idx=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits[2]
    t_unit_scale = TimeUnits[1]

    if logt
        xlab = L"\log_{10}($t\, [%$t_unit_string]$)"
        ylab = L"\log_{10}(SSR)" 
    else
        xlab = L"$t\,[%$t_unit_string]$"
        ylab = L"\log_{10}(SSR)" 
    end

    if isnothing(fig)
        fig = Figure(size=(288,216))
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    else
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,xgridvisible=false,ygridvisible=false)
    end

    if !isnothing(title)
        ax.title = title
    end

    # Is Thermal?
    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
            Nᵃ = FourFlow(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = ProjectionTensor(Uₐ)
            Tᵃᵇ = StressEnergyTensor(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            Pressure = ScalarPressure(Tᵃᵇ,Δab)
            Temperature = ScalarTemperature(Pressure,num)

            if mass_list[j] != 0.0
                Ther = MaxwellJuttner_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            else
                Ther = BlackBody_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            end

            f1D = dropdims(sum(reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j])),dims=(2,3)),dims=(2,3))

            residuals = (f1D .- Ther).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j]*" Thermal?")
            xlims!(ax,log10(sol.t[1] * t_unit_scale),log10(sol.t[end] * t_unit_scale))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j]*" Thermal?")
            xlims!(ax,sol.t[1] * t_unit_scale,sol.t[end] * t_unit_scale)
        end

    end

    # Is Isotropic?
    for j in (species != "All" ? findfirst(x->x==species,name_list) : eachindex(name_list))

        for i in eachindex(sol.t)

            f = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
            f3D = reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j]))
            # sum over angles and divide by number of bins
            f_avg = dropdims(sum(f3D,dims=(2,3)),dims=(2,3)) / (u_num_list[j]*h_num_list[j])
 
            residuals = (f3D .- f_avg).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,linestyle=:dash,label=name_list[j]*" Isotropic?")
            xlims!(ax,log10(sol.t[1] * t_unit_scale),log10(sol.t[end] * t_unit_scale))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,linestyle=:dash,label=name_list[j]*" Isotropic?")
            xlims!(ax,sol.t[1] * t_unit_scale,sol.t[end] * t_unit_scale)
        end

    end

    #ylims!(ax,0.0,nothing)
    
    axislegend(ax,position=:lb)

    end # with_theme

    return fig

end


"""
    function IsThermalPlot0D(type,PhaseSpace,sol,species;...)

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected thermal distribution based on the current distribution of that species as a function of time.  

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
- `title=nothing`: A string to place as the figure title, no title by default.
- `fig=nothing` is the figure to plot on, if `nothing` a new figure will be created.
- `x_idx::Int64=1`, the x index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `y_idx::Int64=1`, the y index of the spatial grid cell that you want to plot the distribution for, default is 1.
- `z_idx::Int64=1`, the z index of the spatial grid cell that you want to plot the distribution for, default is 1.
"""
function IsThermalPlot0D(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct,species::String;fig=nothing,theme=DiplodocusDark(),title=nothing,logt::Bool=false,TimeUnits::Tuple{Float64,String}=(1.0,"\\text{Code Units}"),x_idx=1,y_idx=1,z_idx=1)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    pr_list = Grids.pxr_list
    ur_list = Grids.pyr_list
    hr_list = Grids.pzr_list

    mass_list = Grids.mass_list

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits[2]
    t_unit_scale = TimeUnits[1]

    if logt
        xlab = L"\log_{10}($t\, [%$t_unit_string]$)"
        ylab = L"\log_{10}(SSR)" 
    else
        xlab = L"$t\,[%$t_unit_string]$"
        ylab = L"SSR"
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

            f = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
            Nᵃ = FourFlow(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            Uₐ = [-1.0,0.0,0.0,0.0] # static observer
            num = ScalarNumberDensity(Nᵃ,Uₐ)
            Δab = ProjectionTensor(Uₐ)
            Tᵃᵇ = StressEnergyTensor(sol.f[i],PhaseSpace,j;x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            Pressure = ScalarPressure(Tᵃᵇ,Δab)
            Temperature = ScalarTemperature(Pressure,num)

            if mass_list[j] != 0.0
                Ther = MaxwellJuttner_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            else
                Ther = BlackBody_Distribution(PhaseSpace,name_list[j],Temperature;n=num)
            end

            f1D = dropdims(sum(reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j])),dims=(2,3)),dims=(2,3))

            residuals = (f1D .- Ther).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if logt
            scatterlines!(ax,log10.(sol.t .* t_unit_scale),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(sol.t[1] * t_unit_scale),log10(sol.t[end] * t_unit_scale))
        else
            scatterlines!(ax,sol.t .* t_unit_scale,SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,sol.t[1] * t_unit_scale,sol.t[end] * t_unit_scale)
        end

    end

    ylims!(ax,0.0,nothing)

    #fig[1,2] = Legend(fig,ax,"Particles")
    axislegend(ax,"Particles")

    end # with_theme

    return fig

end

"""
    function IsIsotopicPlot(sol,PhaseSpace;species="All",fig=nothing,theme=DiplodocusDark())

Returns a plot of the sum of squared residuals between the distribution function for each species and an expected isotropic distribution based on the current distribution of that species as a function of time.  
"""
function IsIsotropicPlot(sol::OutputStruct,PhaseSpace::PhaseSpaceStruct;species::String="All",fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime)

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

    SumSquaredResiduals = zeros(Float64,length(sol.t))

    t_unit_string = TimeUnits()
    if t_grid == "u"
        xlab = L"$t$ $%$t_unit_string$"
        ylab = L"SSR"
    elseif t_grid == "l"
        xlab = L"\log_{10}($t$ $%$t_unit_string$)"
        ylab = L"\log_{10}(SSR)"
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

            f = copy(LocationSpeciesToStateVector(sol.f[i],PhaseSpace,species_index=j))
            f3D = reshape(f,(p_num_list[j],u_num_list[j],h_num_list[j]))
            # sum over angles and divide by number of bins
            f_avg = dropdims(sum(f3D,dims=(2,3)),dims=(2,3)) / (u_num_list[j]*h_num_list[j])
 
            residuals = (f3D .- f_avg).^2
            residuals = residuals[isfinite.(residuals)]

            SumSquaredResiduals[i] = sum(residuals)
        
        end

        if t_grid == "u"
            scatterlines!(ax,TimeUnits.(sol.t),SumSquaredResiduals,marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,TimeUnits(sol.t[1]),TimeUnits(sol.t[end]))
        elseif t_grid == "l"
            scatterlines!(ax,log10.(TimeUnits.(sol.t)),log10.(SumSquaredResiduals),marker = :circle,color=theme.palette.color[][mod(2*j-1,7)+1],markersize=0.0,label=name_list[j])
            xlims!(ax,log10(TimeUnits(sol.t[1])),log10(TimeUnits(sol.t[end])))
        end

    end

    ylims!(ax,0.0,nothing)

    axislegend(ax)

    end # with_theme

    return fig

end


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