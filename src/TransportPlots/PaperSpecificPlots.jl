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