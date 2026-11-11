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

    """
        EquatorialPlaneMomentumDistributionPlot1D(type,PhaseSpace,sol,species;...)

    Returns a plot of the momentum distribution of a `species` along the equatorial plane (z=0 in cylindrical coordinates) as a function of time.
    """
    function EquatorialPlaneMomentumDistributionPlot1D(type::Static,PhaseSpace::PhaseSpaceStruct,sol,species::Vector{String};theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true)

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        CHAR_n = PhaseSpace.Characteristic.CHAR_number_density

        name_list = PhaseSpace.name_list
        Momentum = PhaseSpace.Momentum
        Space = PhaseSpace.Spacetime
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
        xr = Grids.xr
        yr = Grids.yr
        zr = Grids.zr

        if wide
            fig = Figure(size=(576,234)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
        if order == 1
            ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,\left[\text{m}^{-3}\right]\right)$"
        elseif order == 2
            ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,\left[\text{m}^{-3}\left(m_ec\right)\right]\right)$"
        elseif order == -2
            ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}V}\,\left[\text{m}^{-3}\left(m_ec\right)^{-3}\right]\right)$"
        elseif order == 100 
            ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,\left[\text{J}\,c^{-1}\,\text{m}^{-3}\right]\right)$"
        else
            ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,\left[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}\right]\right)$"
        end

        if order == 100 
            order_power = 2 
        else
            order_power = order
        end

        zlab = L"$\log_{10}\left(\rho\,[r_\text{s}]\right)$"

        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect(),xticks=-20:5:20,yticks=-20:5:20)
        ax.limits = plot_limits

        linestyles = [(:dash,:dense),:solid,(:dot),(:dashdot),(:dashdotdot)]
        max_total = -Inf
        p_min = Inf
        p_max = -Inf

        legend_elements = []
        line_labels = []

        min_dir = log10(1.0) #log10(yr[1] * zr[1])
        max_dir = log10(yr[end] * zr[end])

        for (species_idx, species_name) in enumerate(species) 

            species_index = findfirst(x->x==species_name,name_list)

            p_num = Momentum.px_num_list[species_index]  
            u_num = Momentum.py_num_list[species_index]
            h_num = Momentum.pz_num_list[species_index]
            dp = Grids.dpx_list[species_index]
            du = Grids.dpy_list[species_index]
            meanp = Grids.mpx_list[species_index]
            if PhaseSpace.Momentum.px_grid_list[species_index] == "l"
                mp_plot = log10.(meanp)
            elseif  PhaseSpace.Momentum.px_grid_list[species_index] == "u"
                mp_plot = meanp
            end
            meanu = Grids.mpy_list[species_index]
            meanh = Grids.mpz_list[species_index]
            p_r = Grids.pxr_list[species_index]
            u_r = Grids.pyr_list[species_index]
            h_r = Grids.pzr_list[species_index]
            mass = Grids.mass_list[species_index]

            f3D = zeros(Float32,p_num,u_num,h_num)
            f1D = zeros(Float32,p_num*u_num*h_num)

            pdNdp = zeros(Float64,p_num)

            p_min = min(p_min,p_r[1])
            p_max = max(p_max,p_r[end])

            @assert y_num == z_num "y_num must equal z_num for equatorial plane plot"

            for t_idx in [1,length(sol.t)] # only plot initial and final time steps

                if t_idx == 1 
                    linestyle_idx = 1
                elseif t_idx == length(sol.t)
                    linestyle_idx = 2
                end

            for space_idx in 1:y_num

                ρ = my[space_idx] * mz[space_idx] # radial coordinate

                f1D .= copy(LocationSpeciesToStateVector(sol.f[t_idx],PhaseSpace,species_index=species_index,x_idx=1,y_idx=space_idx,z_idx=space_idx))

                f3D .= reshape(f1D,(p_num,u_num,h_num))

                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order_power)) / dp[px] * CHAR_n
                end

                pdNdp = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                if order == 100 
                    pdNdp .*= PhaseSpace.Characteristic.CHAR_energy_density / DiplodocusTransport.CONST_c # convert to J c^-1 m^-3
                end

                color = theme.colormap[][(log10(ρ) - min_dir) / (max_dir - min_dir)]
                x = (ρ - min_dir) / (max_dir - min_dir)

                # plot individual lines for each position along the desired direction
                if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp)
                    lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = color,linestyle=linestyles[linestyle_idx])
                else
                    scatterlines!(ax,mp_plot,log10.(pdNdp),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[linestyle_idx])
                end  

                max_f = maximum(x for x in pdNdp if !isnan(x))
                max_total = max(max_f,max_total) 

            end

            end 

            push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
            if species_name == "Ele"
                push!(line_labels,"Electron")
            elseif species_name == "Pos"
                push!(line_labels,"Positron")
            elseif species_name == "Pho"
                push!(line_labels,"Photon")
            end

        end # species loop 

        t_unit_string = TimeUnits()

        if legend
            axislegend(ax,legend_elements,line_labels,position = :lt)
        end

        println("max_total = $max_total")
        if plot_limits == (nothing,nothing)
            xlims!(ax,(log10(p_min)-1.0,log10(p_max)))
            ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
        end
        #println("$((log10(p_min)-1.0,log10(p_max)+1.0))")
        #println("$((log10(max_total)-9.0,log10(max_total)+1.0))")

        Colorbar(fig[1,2],colormap = theme.colormap,limits=(min_dir,max_dir),label=zlab,vertical = true,flipaxis = true,tellheight=false)

        return fig

        end # with_theme

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


    ########### Angle Distribution ###########################
    ##########################################################

    function AngleDistributionWithPolarDistributionPlot(type::Static,PhaseSpace::PhaseSpaceStruct,sol,species::String,time_idx::Int64;theme=DiplodocusDark(),order::Int64=1,plot_limits=(nothing,nothing),wide=false,legend=true,u_indices::Union{Vector{Int64},Nothing}=nothing,x_idx::Int64=1,y_indices::Vector{Int64}=[1],z_idx::Int64=1,dis_range::Float64=10.0,total_line=true)

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        fig = Figure(size=(500,900)) # default single column 4:3 aspect ratio

        name_list = PhaseSpace.name_list
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        species_index = findfirst(x->x==species,name_list)
        CHAR_n = PhaseSpace.Characteristic.CHAR_number_density

        p_num = Momentum.px_num_list[species_index]  
        u_num = Momentum.py_num_list[species_index]
        h_num = Momentum.pz_num_list[species_index]
        dp = Grids.dpx_list[species_index]
        du = Grids.dpy_list[species_index]
        meanp = Grids.mpx_list[species_index]
        if PhaseSpace.Momentum.px_grid_list[species_index] == "l"
            mp_plot = log10.(meanp)
        elseif  PhaseSpace.Momentum.px_grid_list[species_index] == "u"
            mp_plot = meanp
        end
        meanu = Grids.mpy_list[species_index]
        meanh = Grids.mpz_list[species_index]
        p_r = Grids.pxr_list[species_index]
        u_r = Grids.pyr_list[species_index]
        h_r = Grids.pzr_list[species_index]
        mass = Grids.mass_list[species_index]

        xlab = L"$\log_{10}\left(p [m_ec]\right)$"
        if order == 1
            ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
        elseif order == 2
            ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
        elseif order == -2
            ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
        elseif order == 0
            ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
        else
            ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
        end
        ax1 = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,width=200,height=200#=,aspect=DataAspect()=#,xticks=-15:5:15,yticks=-15:5:15)
        ax1.limits = plot_limits
        hidexdecorations!(ax1,grid=false)
        ax2 = Axis(fig[2,1],xlabel=xlab,ylabel=ylab,width=200,height=200#=,aspect=DataAspect()=#,xticks=-15:5:15,yticks=-15:5:15)
        ax2.limits = plot_limits
        hidexdecorations!(ax2,grid=false)
        ax3 = Axis(fig[3,1],xlabel=xlab,ylabel=ylab,width=200,height=200#=,aspect=DataAspect()=#,xticks=-15:5:15,yticks=-15:5:15)
        ax3.limits = plot_limits

        linkxaxes!(ax1,ax2,ax3)

        axu1 = PolarAxis(fig[1,2],theta_0=-pi/2,direction=-1,radius_at_origin=log10(p_r[1])-1,rticks=0:5:5,width=165)
        thetalims!(axu1,0,pi)
        rlims!(axu1,log10(p_r[1]),log10(p_r[end]))

        axu2 = PolarAxis(fig[2,2],theta_0=-pi/2,direction=-1,radius_at_origin=log10(p_r[1])-1,rticks=0:5:5,width=165)
        thetalims!(axu2,0,pi)
        rlims!(axu2,log10(p_r[1]),log10(p_r[end]))

        axu3 = PolarAxis(fig[3,2],theta_0=-pi/2,direction=-1,radius_at_origin=log10(p_r[1])-1,rticks=0:5:5,width=165)
        thetalims!(axu3,0,pi)
        rlims!(axu3,log10(p_r[1]),log10(p_r[end]))

        u_as_theta_grid = zeros(Float64,length(u_r))
        u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
        u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
        u_as_theta_grid_tick_values_string[end] = "u=1.0" 
        u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
        @. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
        @. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

        linestyles = [(:dash,:dense),:solid,(:dot),(:dashdot),(:dashdotdot),(:dash,:dense),:solid,(:dot),(:dashdot),(:dashdotdot)]
        max_total = -Inf32
        p_min = Inf32
        p_max = -Inf32

        legend_elements = []
        line_labels = []
        legend_elements_angle = []
        line_labels_angle = []

        if !isnothing(u_indices)
            for (idx,u_idx) in enumerate(u_indices)
                u_val = PhaseSpace.Grids.mpy_list[1][u_idx]
                color = theme.colormap[][(u_val - u_r[1]) / (u_r[end] - u_r[1])]
                pushfirst!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
                pushfirst!(line_labels_angle,L"$u=%$(round(u_val,sigdigits=2))$")
            end
        end

        f3D = zeros(Float32,p_num,u_num,h_num)
        f1D = zeros(Float32,p_num*u_num*h_num)

        p_min = min(p_min,p_r[1])
        p_max = max(p_max,p_r[end])

        for (y_i,y_idx) in enumerate(y_indices) 

            if y_i == 1
                ax = ax3
                axu = axu3
                ax_col = fig[3,3]
            elseif y_i == 2
                ax = ax2
                axu = axu2
                ax_col = fig[2,3]
            elseif y_i == 3
                ax = ax1
                axu = axu1
                ax_col = fig[1,3]
            end

        f1D .= copy(LocationSpeciesToStateVector(sol.f[time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

        f3D .= reshape(f1D,(p_num,u_num,h_num))
        @. f3D = f3D*(f3D!=Inf)
        # scale by order
        # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
        for px in 1:p_num, py in 1:u_num, pz in 1:h_num
            f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px] * CHAR_n
        end

        for (i,u_idx) in enumerate(u_indices)

            u_val = meanu[u_idx]

            color = theme.colormap[][(u_val - u_r[1]) / (u_r[end] - u_r[1])]

            pdNdp = dropdims(sum(f3D, dims=(3)),dims=(3))[:,u_idx]

            if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp)
                lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = color,linestyle=:solid)
            else
                scatterlines!(ax,mp_plot,log10.(pdNdp),linewidth=2.0,color = color,markersize=0.0,linestyle=:solid)
            end
            max_f = maximum(x for x in pdNdp if !isnan(x))
            max_total = max(max_f,max_total)

        end # angle loop

        pt = 4/3

        if legend && y_i == 1
            axislegend(ax,legend_elements_angle,line_labels_angle,position = :rt,labelsize=7pt,rowgap=-10)
        end

        # angle distribution plot 
        dis = dropdims(sum(reshape(f3D,(p_num,u_num,h_num)), dims=(3)),dims=(3))
        for px in 1:p_num, py in 1:u_num
            dis[px,py] /= du[py]
        end
        replace!(dis,0.0 => NaN) # replace Inf with NaN for plotting

        max_dis = maximum(x for x in dis if !isnan(x))
        min_dis = minimum(x for x in dis if !isnan(x))
        col_range = (log10(max_dis)-dis_range+0.25,log10(max_dis)+0.25)

        

        hm = heatmap!(axu,u_as_theta_grid,log10.(p_r),log10.(dis'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-col_range[2]),lowclip=theme.colormap_var[][0.0])

        translate!(hm,0,0,-100)
        axu.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)

        if order == 1
            Colorbar(ax_col,hm,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}]\right)$",flipaxis=false#=,height=Relative(0.75)=#,tellheight=false)
        elseif order == -2
            Colorbar(ax_col,hm,label=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$",flipaxis=true#=,height=Relative(0.75)=#,tellheight=false)
        elseif order == 0
            Colorbar(ax_col,hm,label=L"$\log_{10}\left(\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
        else
            Colorbar(ax_col,hm,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=true,height=Relative(0.75))
        end
        
        text!(axu,L"$\log_{10}\left(p\,[m_ec]\right)$",position=(-3.05,log10(p_r[end])-1),rotation=pi/2,fontsize=9pt)

        if plot_limits == (nothing,nothing)
            xlims!(ax,(log10(p_min)-1.0,log10(p_max)+1.0))
            ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
        end

        end # y values

        colsize!(fig.layout,1,Aspect(1,1))
        colsize!(fig.layout,2,Fixed(100))
        #colsize!(fig.layout,3,Fixed(5))

        rowgap!(fig.layout,1,0.0)
        rowgap!(fig.layout,2,0.0)

        return fig

        end # with_theme

    end

##### MISC Undocumented Plots for papers ######
###############################################

    # ============== AM3 Test Plots ============== #

    function AM3_MomentumDistributionPlot(filePath,t_max,t_min,t_grid;plot_limits=(nothing,nothing),theme=DiplodocusDark(),wide=false)

        fileExist = isfile(filePath)

        if fileExist
            f = DC.jldopen(filePath,"r+");

            meanp_ele = f["meanp_ele"];
            f_ele = f["f_ele"];
            t_ele = f["t_ele"];

            meanp_pho = f["meanp_pho"];
            f_pho = f["f_pho"];
            t_pho = f["t_pho"];

            DC.close(f)  
        else
            error("no file at $filePath found")
        end

        # unit conversion 
        eV_to_mElec2 = 1.60217e-19 / 9.109e-31 / 2.9979e8^2
        cm3_to_m3 = 1e6 

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        fig = Figure()
        if wide
            fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        xlab = L"$\log_{10}\left(p [m_ec]\right)$"
        ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)]\right)$"
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
        ax.limits = plot_limits

        linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
        legend_elements = []
        line_labels = []
        
        for i in 1:length(t_pho)

            t = t_pho[i]
            print("t=$t \r")
            if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps
                if t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                elseif t_grid == "l"
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                end

                # sum along u and h directions
                pdNdp = f_pho[i,:][1]
                #println("$pdNdp")
                scatterlines!(ax,log10.(meanp_pho .* eV_to_mElec2),log10.(meanp_pho .* pdNdp .* cm3_to_m3 .* eV_to_mElec2),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
            end

        end

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
        push!(line_labels,"Pho")

        for i in 1:length(t_ele)

                t = t_ele[i]
                if log10(t) % 1 == 0.0 && t <= t_max # 10^n timesteps 
                    if t_grid == "u"
                        color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                    elseif t_grid == "l"
                        color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                    end

                    # sum along u and h directions
                    pdNdp = f_ele[i,:][1]
                    scatterlines!(ax,log10.(sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1)),log10.(sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1) .* pdNdp .* cm3_to_m3),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
                end

        end

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
        push!(line_labels,"Ele")


        if t_grid == "u"
            Colorbar(fig[1,2],colormap = theme.colormap,limits=(TimeUnits(t_min),TimeUnits(t_max)),label=L"$t$ $[\text{s} * \sigma_{T}c]$")
        elseif t_grid == "l"
            Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(t_min),log10(t_max)),label=L"$\log_{10}\left(t [\text{s}]\right)$")
        end

        axislegend(ax,legend_elements,line_labels,position = :lt)
    
        return fig

        end # with_theme

    end

    function AM3_DIP_Combo_MomentumDistributionPlot(filePath_AM3,sol_DIP,PhaseSpace_DIP,t_max,t_min,t_grid;plot_limits=(nothing,nothing),theme=DiplodocusDark(),ele_err=true,pos_err=true,yticks=1:2:15,pos=true)

        name_list = PhaseSpace_DIP.name_list
        Momentum = PhaseSpace_DIP.Momentum
        Grids = PhaseSpace_DIP.Grids
        Time = PhaseSpace_DIP.Time

        # load AM3 Data
        fileExist = isfile(filePath_AM3)
        if fileExist
            f = DC.jldopen(filePath_AM3,"r+");

            meanp_ele = f["meanp_ele"];
            f_ele = f["f_ele"];
            t_ele = f["t_ele"];

            if pos
                meanp_pos = f["meanp_pos"];
                f_pos = f["f_pos"];
                t_pos = f["t_pos"];
            end

            meanp_pho = f["meanp_pho"];
            f_pho = f["f_pho"];
            t_pho = f["t_pho"];

            DC.close(f)  
        else
            error("no file at $filePath_AM3 found")
        end

        # unit conversion 
        eV_to_mElec2 = 1.60217e-19 / 9.109e-31 / 2.9979e8^2
        cm3_to_m3 = 1e6 

        CairoMakie.activate!(inline=true) # plot in vs code window

        with_theme(theme) do

        fig = Figure(size=(576,460)) # double column 8:6 aspect ratio

        #g1 = fig[1,1] = GridLayout()
        g2 = fig[1,1]  = GridLayout()
        #g3 = fig[1,3] = GridLayout()

        #colsize!(g2.layout,0,Relative(0.05))
        #colsize!(g2.layout,1,Relative(0.9))
        #colsize!(g2.layout,2,Relative(0.05))

        xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
        ylab = L"$\log_{10}\left(p^2\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
        ylab_err = L"$\text{Frac. Diff.}$"
        
        ax_DIP = Axis(g2[1,1],aspect=DataAspect(),ylabel=ylab,yticks=yticks,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(2),xminorgridvisible=true,yminorticksvisible=true,yminorgridvisible=true)
        hidexdecorations!(ax_DIP,grid=false,minorgrid=false,minorticks=false)

        ax_AM3 = Axis(g2[2,1],aspect=DataAspect(),yticks=yticks,xminorticks=IntervalsBetween(5),yminorticks=IntervalsBetween(2),xminorgridvisible=true,yminorgridvisible=true,yminorticksvisible=true)
        hidexdecorations!(ax_AM3,grid=false,minorticks=false,minorgrid=false)

        ax_err = Axis(g2[3,1],xlabel=xlab,aspect=DataAspect(),ylabel=ylab_err,yticks=[-1.0,0.0,1.0],xminorticks=IntervalsBetween(5),xminorgridvisible=true,xminorticksvisible=true)

        #pt = 4/3
        #Label(g2[1:2,0],ylab,rotation = pi/2,fontsize=9pt)
        #Label(g2[3,0],ylab_err,rotation = pi/2,fontsize=9pt)

        linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
        legend_elements = []
        line_labels = []

        order = 2

        ######## Photon Lines ########
        ##############################

        species_name = "Pho"
        species_index = findfirst(x->x==species_name,name_list)
        p_num = Momentum.px_num_list[species_index]  
        u_num = Momentum.py_num_list[species_index]
        h_num = Momentum.pz_num_list[species_index]
        dp = Grids.dpx_list[species_index]
        du = Grids.dpy_list[species_index]
        meanp = Grids.mpx_list[species_index]
        meanu = Grids.mpy_list[species_index]
        meanh = Grids.mpz_list[species_index]
        p_r = Grids.pxr_list[species_index]
        u_r = Grids.pyr_list[species_index]
        h_r = Grids.pzr_list[species_index]
        mass = Grids.mass_list[species_index]

        f3D = zeros(Float32,p_num,u_num,h_num)
        f1D = zeros(Float32,p_num*u_num*h_num)
        
        for i in 1:length(t_pho)

            # AM3
            t = t_pho[i]
            #println("t=$t")
            if log10(t) % 1 == 0.0 && t <= t_max && t>=t_min# 10^n timesteps
                t_plot = t
                println("t=$t")

                if t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                elseif t_grid == "l"
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                end

                meanp_AM3 = meanp_pho .* eV_to_mElec2

                # AM3
                pdNdp_AM3 = meanp_AM3 .* f_pho[i,:][1] .* cm3_to_m3
                #println("$pdNdp")
                pdNdp_AM3_log=log10.(pdNdp_AM3)
                replace!(pdNdp_AM3_log,-Inf=>-45.0)
                scatterlines!(ax_AM3,log10.(meanp_AM3),pdNdp_AM3_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])

                # DIP
                t_idx = find_closest(sol_DIP.t,SIToCodeUnitsTime(t_plot))
                if isnothing(t_idx)
                    t_idx = 1 # works for both log and uniform diplodocus time stepping
                end          
                println(t_idx)
                f1D .= copy(LocationSpeciesToStateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
                f3D .= reshape(f1D,(p_num,u_num,h_num))
                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                end
                pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp)
                    lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[1])
                else
                    pdNdp_DIP_log=log10.(pdNdp_DIP)
                    replace!(pdNdp_DIP_log,-Inf=>-45.0)
                    scatterlines!(ax_DIP,log10.(meanp),pdNdp_DIP_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])
                end

                # error plot
                # regrid AM3 results (assuming it has finest grid) to grid size of DIP
                pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
                for p in eachindex(meanp)
                    p_idx = find_closest(meanp_AM3,meanp[p])
                    pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
                end
                err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
                replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
                replace!(err,Inf=>NaN) # remove values for which DIP does not have values

                scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[1])

            end

        end

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
        push!(line_labels,"Photon")

        ######## Electron Lines ########
        ################################

        species_name = "Ele"
        species_index = findfirst(x->x==species_name,name_list)
        p_num = Momentum.px_num_list[species_index]  
        u_num = Momentum.py_num_list[species_index]
        h_num = Momentum.pz_num_list[species_index]
        dp = Grids.dpx_list[species_index]
        du = Grids.dpy_list[species_index]
        meanp = Grids.mpx_list[species_index]
        meanu = Grids.mpy_list[species_index]
        meanh = Grids.mpz_list[species_index]
        p_r = Grids.pxr_list[species_index]
        u_r = Grids.pyr_list[species_index]
        h_r = Grids.pzr_list[species_index]
        mass = Grids.mass_list[species_index]

        f3D = zeros(Float32,p_num,u_num,h_num)
        f1D = zeros(Float32,p_num*u_num*h_num)

        for i in 1:length(t_ele)

            t = t_ele[i]
            if log10(t) % 1 == 0.0 && t <= t_max && t>=t_min# 10^n timesteps 
                t_plot = t
                if t_grid == "u"
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                elseif t_grid == "l"
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                end
                meanp_AM3 = sqrt.((meanp_ele .* eV_to_mElec2).^2 .-1)

                # AM3
                pdNdp_AM3 = meanp_AM3 .* f_ele[i,:][1] .* cm3_to_m3
                pdNdp_AM3_log=log10.(pdNdp_AM3)
                replace!(pdNdp_AM3_log,-Inf=>-45.0)
                scatterlines!(ax_AM3,log10.(meanp_AM3),pdNdp_AM3_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])

                # DIP
                t_idx = find_closest(sol_DIP.t,SIToCodeUnitsTime(t_plot))
                if isnothing(t_idx)
                    t_idx = 1 # works for both log and uniform diplodocus time stepping
                end 
                f1D .= copy(LocationSpeciesToStateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
                f3D .= reshape(f1D,(p_num,u_num,h_num))
                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                end
                pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                    idx = findfirst(!iszero,pdNdp)
                    lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[2])
                else
                    pdNdp_DIP_log=log10.(pdNdp_DIP)
                    replace!(pdNdp_DIP_log,-Inf=>-45.0)
                    scatterlines!(ax_DIP,log10.(meanp),pdNdp_DIP_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
                end

                if ele_err
                    # error plot
                    # regrid AM3 results (assuming it has finest grid) to grid size of DIP
                    pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
                    for p in eachindex(meanp)
                        p_idx = find_closest(meanp_AM3,meanp[p])
                        pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
                    end
                    err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
                    replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
                    replace!(err,Inf=>NaN) # remove values for which DIP does not have values

                    scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[2])
                end

            end

        end

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
        push!(line_labels,"Electron")

        ######## Positron Lines ########
        ################################

        species_name = "Pos"
        species_index = findfirst(x->x==species_name,name_list)
        if species_index != nothing
            p_num = Momentum.px_num_list[species_index]  
            u_num = Momentum.py_num_list[species_index]
            h_num = Momentum.pz_num_list[species_index]
            dp = Grids.dpx_list[species_index]
            du = Grids.dpy_list[species_index]
            meanp = Grids.mpx_list[species_index]
            meanu = Grids.mpy_list[species_index]
            meanh = Grids.mpz_list[species_index]
            p_r = Grids.pxr_list[species_index]
            u_r = Grids.pyr_list[species_index]
            h_r = Grids.pzr_list[species_index]
            mass = Grids.mass_list[species_index]

            f3D = zeros(Float32,p_num,u_num,h_num)
            f1D = zeros(Float32,p_num*u_num*h_num)

            for i in 1:length(t_pos)

                t = t_pos[i]
                if log10(t) % 1 == 0.0 && t <= t_max && t>=t_min# 10^n timesteps 
                    t_plot = t
                    if t_grid == "u"
                        color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                    elseif t_grid == "l"
                        color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                    end
                    meanp_AM3 = sqrt.((meanp_pos .* eV_to_mElec2).^2 .-1)

                    # AM3
                    pdNdp_AM3 = meanp_AM3 .* f_pos[i,:][1] .* cm3_to_m3
                    pdNdp_AM3_log=log10.(pdNdp_AM3)
                    replace!(pdNdp_AM3_log,-Inf=>-45.0)
                    scatterlines!(ax_AM3,log10.(meanp_AM3),pdNdp_AM3_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])

                    # DIP
                    t_idx = find_closest(sol_DIP.t,SIToCodeUnitsTime(t_plot))
                    if isnothing(t_idx)
                        t_idx = 1 # works for both log and uniform diplodocus time stepping
                    end 
                    f1D .= copy(LocationSpeciesToStateVector(sol_DIP.f[t_idx],PhaseSpace_DIP,species_index=species_index))
                    f3D .= reshape(f1D,(p_num,u_num,h_num))
                    @. f3D = f3D*(f3D!=Inf)
                    # scale by order
                    # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                    for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                        f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                    end
                    pdNdp_DIP = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                    if sum(@. !isnan(pdNdp_DIP) * !isinf(pdNdp_DIP) * !iszero(pdNdp_DIP)) == 1 # there is only one valid position so scatterlines doesn't work
                        idx = findfirst(!iszero,pdNdp)
                        lines!(ax,[log10(meanp[idx]), log10(meanp[idx])],[-20.0, log10(pdNdp_DIP[idx])],linewidth=2.0,color = color,linestyle=linestyles[3])
                    else
                        pdNdp_DIP_log=log10.(pdNdp_DIP)
                        replace!(pdNdp_DIP_log,-Inf=>-45.0)
                        scatterlines!(ax_DIP,log10.(meanp),pdNdp_DIP_log,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])
                    end

                    if pos_err
                    # error plot
                    # regrid AM3 results (assuming it has finest grid) to grid size of DIP
                    pdNdp_AM3_DIP_Grid = zeros(Float64,size(meanp))
                    for p in eachindex(meanp)
                        p_idx = find_closest(meanp_AM3,meanp[p])
                        pdNdp_AM3_DIP_Grid[p] = pdNdp_AM3[p_idx]
                    end
                    err = (pdNdp_AM3_DIP_Grid.-pdNdp_DIP)./pdNdp_DIP
                    replace!(err,-1.0=>NaN) # remove values for which AM3 does not have values
                    replace!(err,Inf=>NaN) # remove values for which DIP does not have values

                    scatterlines!(ax_err,log10.(meanp),err,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[3])
                end

                end
            end
            push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[3],linewidth = 2.0))
            push!(line_labels,"Positron")
        end

        if t_grid == "u"
            Colorbar(g2[1:2,2],colormap = theme.colormap,limits=(TimeUnits(t_min),TimeUnits(t_max)),label=L"$t\,$ $[\text{s} * \sigma_{T}c]$",height=Relative(0.66))
        elseif t_grid == "l"
            Colorbar(g2[1:2,2],colormap = theme.colormap,limits=(log10(t_min),log10(t_max)),label=L"$\log_{10}\left(t \,[\text{s}]\right)$",height=Relative(0.66))
        end


        axislegend(ax_DIP,legend_elements,line_labels,position = :lt)

        linkxaxes!(ax_DIP,ax_err)
        linkxaxes!(ax_DIP,ax_AM3)
        #ax_DIP.limits = plot_limits
        xlims!(ax_DIP,plot_limits[1][1],plot_limits[1][2]) 
        ylims!(ax_DIP,plot_limits[2][1],plot_limits[2][2]) 
        #ax_AM3.limits = plot_limits
        xlims!(ax_AM3,plot_limits[1][1],plot_limits[1][2])
        ylims!(ax_AM3,plot_limits[2][1],plot_limits[2][2]) 

        xlims!(ax_err,plot_limits[1][1],plot_limits[1][2])
        ylims!(ax_err,-1.0,1.0)
        
        rowsize!(g2,1,Relative(0.4))
        rowsize!(g2,2,Relative(0.4))
        rowsize!(g2,3,Relative(0.2))

        rowgap!(g2,1,0.0)
        rowgap!(g2,2,0.0)
    
        return fig

        end # with_theme

    end


    function TwoSolAngleDistributionPlot(twosol::Tuple{OutputStruct,OutputStruct},species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=false,legend=true,framerate=12,filename="TwoSolAngleDistribution.mp4",figure=nothing)

        CairoMakie.activate!(inline=true) # plot in vs code window
        with_theme(theme) do

        xlab = L"$\log_{10}\left(p [m_ec]\right)$"
        if order == 1
            ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
        elseif order == 2
            ylab = L"$\log_{10}\left(p^{2}\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)]\right)$"
        else 
            ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
        end

        if isnothing(figure)
            time_idx = Observable(1) # index of the current time step
            t = @lift(TimeUnits(twosol[1].t[$time_idx]))
            if wide
                fig = Figure(size=(500,216)) # double column 8:3 aspect ratio
            else
                fig = Figure() # default single column 4:3 aspect ratio
            end
            #fig = Figure(size =(3.25inch,3.25inch)) # 1:1 aspect ratio
            ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
        else
            fig, time_idx = figure # use the provided figure and time index
            ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
        end

        ax.limits = plot_limits

        t_unit_string = TimeUnits()
                
        text!(ax,@lift("t=$(round($(t), sigdigits = 3))  "),fontsize=12pt,align=(:right,:center),space=:relative,offset=(375.0,135.0))
        text!(ax,L"%$t_unit_string",fontsize=12pt,align=(:left,:center),space=:relative,offset=(375.0,135.0))

        linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
        line_labels = []
        legend_elements = []  
        legend_elements_angle = []
        line_labels_angle = []  

        name_list = PhaseSpace.name_list
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Time = PhaseSpace.Time

        counter = 1;

        for (species_idx, species_name) in enumerate(species) 

            species_index = findfirst(x->x==species_name,name_list)

            p_num = Momentum.px_num_list[species_index]  
            u_num = Momentum.py_num_list[species_index]
            h_num = Momentum.pz_num_list[species_index]
            dp = Grids.dpx_list[species_index]
            meanp = Grids.mpx_list[species_index]
            meanu = Grids.mpy_list[species_index]
            p_r = Grids.pxr_list[species_index]
            u_r = Grids.pyr_list[species_index]
            mass = Grids.mass_list[species_index]

            for i in 1:2

                color = theme.palette.color[][mod(2*i-1,7)+1]

                sol = twosol[i]

                pdNdp = @lift begin
                    f1D = zeros(Float32,p_num*u_num*h_num)
                    f1D .= copy(LocationSpeciesToStateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index))
                    f3D = zeros(Float32,p_num,u_num,h_num)
                    f3D .= reshape(f1D,(p_num,u_num,h_num))
                    @. f3D = f3D*(f3D!=Inf)
                    # scale by order
                    # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                    for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                        f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                    end
                    # sum along u and h directions
                    log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,2]
                end

                scatterlines!(ax,log10.(meanp),pdNdp,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

                if legend==true && counter == 1 
                    push!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
                    if i == 1
                        push!(line_labels_angle,L"\text{Iso}")
                    end
                    if i == 2 
                        push!(line_labels_angle,L"\text{Ani}")
                    end
                end

            end # sol loop 

            counter += 1

            push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
            push!(line_labels,species_name)

        end # species loop 

        if legend
            axislegend(ax,legend_elements,line_labels,position = :lt)
            axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
        end
        
        if !isnothing(filename)
            # recording the animation
            time_idxs = 1:length(twosol[1].t)
            record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
                print("generating frame: $frame \r")
                time_idx[] = frame
            end
        end

        end # with_theme

    end