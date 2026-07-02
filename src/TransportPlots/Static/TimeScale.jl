"""
    TimeScalePlot(PhaseSpace, scheme, state)

Function plots the time scales for interaction loss rates of a simulation given the system sulution `sol` at a time indices `t_idxs`.

Optional arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `p_timescale`: whether to plot the timescale for momentum magnitude loss or state vector loss, default is `false` i.e. plot state vector losses to aid in assessing time step limits for stability.
- `logt`: default is `false`. If a uniform time grid is used for the solution, this can be toggled to `true` to display the solution in log10 time steps rather than uniform time steps.
- `dt`: default is `false`. If `true` the time step `dt` corresponding to the current simulation time will also be plotted as a horizontal line.
"""
function TimeScalePlot(method::AbstractSteppingMethod,sol::OutputStruct,t_idxs::Vector{Int64},species::Vector{String};wide=false,paraperp::Bool=false,u_avg::Bool=true,logt::Bool=false,p_timescale=false,legend=true,plot_dt::Bool=false,plot_limits=(nothing,nothing),theme=DiplodocusDark(),TimeUnits::Function=CodeToCodeUnitsTime,direction::String="all")

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end

    t_unit_string = TimeUnits()

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if p_timescale
        ylab = L"$\log_{10}\left(\frac{p}{\mathrm{d}p/\mathrm{d}t}\,%$t_unit_string \right)$"
    else
        ylab = L"$\log_{10}\left(Timescale\,%$t_unit_string\right)$"
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
    legend_elements_species = []
    line_labels_species = []
    legend_elements_angle = []
    line_labels_angle = []

    PhaseSpace=method.PhaseSpace
    Momentum = PhaseSpace.Momentum
    Time = PhaseSpace.Time
    t_grid = PhaseSpace.Time.t_grid
    Grids = PhaseSpace.Grids
    tr = PhaseSpace.Grids.tr
    t_min = logt ? tr[2]/10 : tr[1]
    t_max = tr[end]
    p_num_list = Momentum.px_num_list 
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    dp_list = Grids.dpx_list
    mp_list = Grids.mpx_list
    name_list = PhaseSpace.name_list

    df = zeros(eltype(sol.f[1]),size(sol.f[1]))
    timescale = zeros(eltype(sol.f[1]),size(sol.f[1]))

    for (idx,t_idx) in enumerate(t_idxs)

        f = sol.f[t_idx]
        @. f = f*(f>=1f-20) # to align with solver cut
        t = sol.t[t_idx]
        
        if t_grid == "l" || logt
            color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
        elseif t_grid == "u"
            color = theme.colormap[][(t - t_min) / (t_max - t_min)]
        end

        t_idx_global = find_closest(tr,t)
        if isnothing(t_idx_global)
            t_idx_global = 1
        end
        if t_idx_global != 1
            t_idx_global -= 1 # -1 as sol is at t[t+1] except initial time step
        end 

        dt0 = tr[2] - tr[1]
        dt = PhaseSpace.Grids.dt[t_idx_global]

        method(df,f,dt0,dt,t)

        #= DIAG needs to be fixed as off diagonal terms can contribute to losses e.g. for compton scattering there could be a loss of electrons from photons. For a better estimate we should split method.temp into two parts: a positive part giving the rate of gain and a negative part giving the rate of loss.

        Could also split df into forces, emissive and absorption parts.
        =#

        println("Diag should be fixed!!!!!")

        if direction == "all"
            df .= DiplodocusTransport.diag(method.temp) .* f
            #df .= method.temp * f
        elseif direction == "I"
            df .= method.FluxM.Ap_Flux \ (DiplodocusTransport.diag(method.FluxM.I_Flux) .* (dt / dt0)) .* f
        end

        @. timescale =  dt * f / df
        timescale = timescale .* (timescale .<= 0.0) # only want loss rates 

        for (species_idx, species_name) in enumerate(species) 

            species_index = findfirst(x->x==species_name,name_list)

            name = name_list[species_index]
            p_num = p_num_list[species_index]
            u_num = u_num_list[species_index]
            h_num = h_num_list[species_index]
            dp = dp_list[species_index]
            mp = mp_list[species_index]
            if Momentum.px_grid_list[species_index] == "l"
                mp_plot = log10.(mp)
            elseif Momentum.px_grid_list[species_index] == "u"
                mp_plot = mp
            end
            mu = Grids.mpy_list[species_index]

            timescale1D = copy(LocationSpeciesToStateVector(timescale,PhaseSpace,species_index=species_index))
            timescale3D = reshape(timescale1D,(p_num,u_num,h_num))
            @. timescale3D = timescale3D*(timescale3D!=Inf)
            @. timescale3D = timescale3D*(timescale3D!=-Inf)
            #@. timescale3D = timescale3D*(timescale3D<=0.0)

            timescale2D = dropdims(sum(timescale3D, dims=(3)),dims=(3))
            # timescale is timescale for particle losses not particle energy losses, to get that we need to scale by dp/mp (i.e. dp/p)
            if p_timescale
                timescale2D = mp ./ dp .* timescale2D
            end

            if paraperp==true

                println("tscale=$(TimeUnits.(Float64.(abs.(timescale2D[10,ceil(Int64,u_num/2)]))))")

                scatterlines!(ax,mp_plot,log10.(TimeUnits.(Float64.(abs.(timescale2D[:,end])))),linewidth=2.0,color=theme.textcolor[],markersize=0.0,linestyle=linestyles[1])
                push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = linestyles[1],linewidth = 2.0))
                push!(line_labels_angle,L"\parallel")

                scatterlines!(ax,mp_plot,log10.(TimeUnits.(Float64.(abs.(timescale2D[:,ceil(Int64,u_num/2)])))),linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[2])
                push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = linestyles[2],linewidth = 2.0))
                push!(line_labels_angle,L"\perp")
                

            else
                u_avg_T = zeros(Float64,length(mp_plot))

                for u in 1:u_num

                    if u_avg
                        @. u_avg_T += timescale2D[:,u] / u_num
                    else
                        scatterlines!(ax,mp_plot,log10.(TimeUnits.(Float64.(abs.(timescale2D[:,u])))),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

                        #hlines!(ax,log10(1 / (1-mu[u]^2)),color = theme.palette.color[][mod(2*u-1,7)+1])

                        #=if species_index == 1
                            push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = :solid,linewidth = 2.0))
                            push!(line_labels_angle,L"%$(mu[u])")
                        end=#
                    end

                end 

                if u_avg
                    scatterlines!(ax,mp_plot,log10.(TimeUnits.(Float64.(abs.(u_avg_T)))),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

                end
            end

            #push!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = :solid,linewidth = 2.0))
            #push!(line_labels_angle,L"%$(mu[u])")

            if idx == 1
                push!(legend_elements_species,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
                if name == "Ele"
                    name = "Electron"
                elseif name == "Pho"
                    name = "Photon"
                elseif name == "Pos"
                    name = "Positron"
                end
                push!(line_labels_species,name)
                
            end

        end # species loop

        if plot_dt && idx == 1
            hlines!(ax,log10.(TimeUnits(dt)),color=theme.textcolor[],linewidth=2.0,linestyle = linestyles[length(species)+1])
            push!(legend_elements_species,LineElement(color = theme.textcolor[], linestyle = linestyles[length(species)+1],linewidth = 2.0))
            push!(line_labels_species,"dt")
        end

    end # t loop

        if paraperp
            axislegend(ax,legend_elements_angle,line_labels_angle,position = :lb)
        end
        if legend
            axislegend(ax,legend_elements_species,line_labels_species,position = :lb)
        end

    return fig

    end # with theme

end