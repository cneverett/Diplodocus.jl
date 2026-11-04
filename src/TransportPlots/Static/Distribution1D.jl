"""
    MomentumDistribution1DPlot(type,sol,species,PhaseSpace,time_idx;...)

Plots the angle averaged distribution function of a given vector of particle `species` as a function of distance along a single spatial dimension given by `direction` at a given time (`t_idx`) given by the `sol` based on the conditions held in `PhaseSpace`. 
    
The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
- `TimeUnits`: a function that converts the time given in code units to the desired units for plotting
- `plot_limits`: the limits of the x and y axes, default is `(nothing,nothing)` which sets the limits automatically based on the data.
- `wide`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), default is `false` (single column 4:3 aspect ratio).
- `legend`: if `true`, a legend is added to the plot, default is `true`.


Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumDistribution.mp4".
- `figure`: default is `nothing`, which creates a new figure. If a figure is provided, the plot is added to that figure instead of creating a new one, this should be of the form of a tuple of `figure` and `time_idx` from the main plot.
- `initial`: default is `false`, if `true` causes the initial distribution to remain on the plot
"""
function MomentumDistribution1DPlot(type::Static,sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,direction::String,t_idx::Int64;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true)

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
        fig = Figure(size=(576,432)) # double column 8:6 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
    elseif order == 2
        ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    else
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if direction == "x"
        zlab = L"$\left(x\,[\text{Code Units}]\right)$"
        n_dir = x_num
        grid_type = Space.x_grid
        if grid_type == "u"
            mdir_plot = mx
            min_dir = xr[1]
            max_dir = xr[end]
            zlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
            min_dir = log10(xr[1])
            max_dir = log10(xr[end])
            zlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
    elseif direction == "y"
        n_dir = y_num
        grid_type = Space.y_grid
        if grid_type isa StretchGrid
            mdir_plot = my
            min_dir = yr[1]
            max_dir = yr[end]
            zlab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(my)
            min_dir = log10(yr[1])
            max_dir = log10(yr[end])
            zlab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
    elseif direction == "z"
        n_dir = z_num
        grid_type = Space.z_grid
        if grid_type isa StretchGrid
            mdir_plot = mz
            min_dir = zr[1]
            max_dir = zr[end]
            zlab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
            min_dir = log10(zr[1])
            max_dir = log10(zr[end])
            zlab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
    max_total = -Inf
    p_min = Inf
    p_max = -Inf

    legend_elements = []
    line_labels = []

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

        pdNdp = zeros(Float64,p_num,n_dir)
        pdNdp_total = zeros(Float64,p_num)

        p_min = min(p_min,p_r[1])
        p_max = max(p_max,p_r[end])

        for x_idx in 1:x_num, y_idx in 1:y_num, z_idx in 1:z_num

            f1D .= copy(LocationSpeciesToStateVector(sol.f[t_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

            f3D .= reshape(f1D,(p_num,u_num,h_num))

            @. f3D = f3D*(f3D!=Inf)
            # scale by order
            # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
            for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px] * CHAR_n
            end

            if direction == "x"
                # sum along u and h directions
                @view(pdNdp[:,x_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            elseif direction == "y"
                # sum along u and h directions
                @view(pdNdp[:,y_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            elseif direction == "z"
                # sum along u and h directions
                @view(pdNdp[:,z_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3)) 
            end

            pdNdp_total .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3))

        end

        for i in 1:n_dir

            pdNdp_local = @view(pdNdp[:,i])

            pos = mdir_plot[i]

            color = theme.colormap[][(pos - min_dir) / (max_dir - min_dir)]

            # plot individual lines for each position along the desired direction
            if sum(@. !isnan(pdNdp_local) * !isinf(pdNdp_local) * !iszero(pdNdp_local)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp_local)
                lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp_local[idx])],linewidth=2.0,color = color,linestyle=linestyles[species_idx])
            else
                scatterlines!(ax,mp_plot,log10.(pdNdp_local),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])
            end  

        end

        max_f = maximum(x for x in pdNdp_total if !isnan(x))
        println("max_f for $species_name = $max_f")
        max_total = max(max_f,max_total)   
        
        # plot total line 
        if sum(@. !isnan(pdNdp_total) * !isinf(pdNdp_total) * !iszero(pdNdp_total)) == 1 # there is only one valid position so scatterlines doesn't work
            idx = findfirst(!iszero,pdNdp_total)
            lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp_total[idx])],linewidth=2.0,color = theme.textcolor[],linestyle=linestyles[species_idx])
        else
            scatterlines!(ax,mp_plot,log10.(pdNdp_total),linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[species_idx])
        end  

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
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
        xlims!(ax,(log10(p_min)-1.0,log10(p_max)+1.0))
        ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
    end
    #println("$((log10(p_min)-1.0,log10(p_max)+1.0))")
    #println("$((log10(max_total)-9.0,log10(max_total)+1.0))")

    Colorbar(fig[2,1],colormap = theme.colormap,limits=(min_dir,max_dir),label=zlab,vertical = false,flipaxis = false)

    return fig

    end # with_theme

end

function MomentumDistribution1DPlot(type::Animated,sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,direction::String;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,thermal=false,plot_limits=(nothing,nothing),wide=false,legend=true,framerate=12,filename="MomentumDistribution1D.mp4",figure=nothing)

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
    xr = Grids.xr
    yr = Grids.yr
    zr = Grids.zr

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}]\right)$"
    elseif order == 2
        ylab = L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    else
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if direction == "x"
        zlab = L"$\left(x\,[\text{Code Units}]\right)$"
        n_dir = x_num
        grid_type = Space.x_grid
        if grid_type == "u"
            mdir_plot = mx
            min_dir = xr[1]
            max_dir = xr[end]
            zlab = L"$\left(x\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
            min_dir = log10(xr[1])
            max_dir = log10(xr[end])
            zlab = L"$\log_{10}\left(x\,[\text{Code Units}]\right)$"
        end
    elseif direction == "y"
        n_dir = y_num
        grid_type = Space.y_grid
        if grid_type == "u"
            mdir_plot = my
            min_dir = yr[1]
            max_dir = yr[end]
            zlab = L"$\left(y\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(my)
            min_dir = log10(yr[1])
            max_dir = log10(yr[end])
            zlab = L"$\log_{10}\left(y\,[\text{Code Units}]\right)$"
        end
    elseif direction == "z"
        n_dir = z_num
        grid_type = Space.z_grid
        if grid_type == "u"
            mdir_plot = mz
            min_dir = zr[1]
            max_dir = zr[end]
            zlab = L"$\left(z\,[\text{Code Units}]\right)$"
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
            min_dir = log10(zr[1])
            max_dir = log10(zr[end])
            zlab = L"$\log_{10}\left(z\,[\text{Code Units}]\right)$"
        end
    end

    if isnothing(figure)
        time_idx = Observable(1)
        t = @lift(sol.t[$time_idx])
        if wide
            fig = Figure(size=(576,432)) # double column 8:6 aspect ratio
        else
            fig = Figure() # default single column 4:3 aspect ratio
        end
        ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    else
        fig, time_idx = figure
        ax = Axis(fig,xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    end
    
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]
    max_total = -Inf
    p_min = Inf
    p_max = -Inf

    legend_elements = []
    line_labels = []

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

        pdNdp = @lift begin
            f3D = zeros(Float32,p_num,u_num,h_num)
            f1D = zeros(Float32,p_num*u_num*h_num)
            pdNdp_t = zeros(Float64,p_num,n_dir)
            p_min = min(p_min,p_r[1])
            p_max = max(p_max,p_r[end])
            for x_idx in 1:x_num, y_idx in 1:y_num, z_idx in 1:z_num
                f1D .= copy(LocationSpeciesToStateVector(sol.f[$time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
                f3D .= reshape(f1D,(p_num,u_num,h_num))
                @. f3D = f3D*(f3D!=Inf)
                # scale by order
                # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
                for px in 1:p_num, py in 1:u_num, pz in 1:h_num
                    f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px]
                end
                if direction == "x"
                    # sum along u and h directions
                    @view(pdNdp_t[:,x_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                elseif direction == "y"
                    # sum along u and h directions
                    @view(pdNdp_t[:,y_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
                elseif direction == "z"
                    # sum along u and h directions
                    @view(pdNdp_t[:,z_idx]) .+= dropdims(sum(f3D, dims=(2,3)),dims=(2,3)) 
                end
            end
            log10.(pdNdp_t)
        end

        pdNdp_total = @lift log10.(dropdims(sum(10 .^ $pdNdp, dims=(2)),dims=(2)))

        for i in 1:n_dir

            pdNdp_local = @lift @view($pdNdp[:,i])

            pos = mdir_plot[i]

            color = theme.colormap[][(pos - min_dir) / (max_dir - min_dir)]

            scatterlines!(ax,mp_plot,pdNdp_local,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

        end
        
        # plot total line 
        scatterlines!(ax,mp_plot,pdNdp_total,linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[species_idx])

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
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


    Colorbar(fig[2,1],colormap = theme.colormap,limits=(min_dir,max_dir),label=zlab,vertical = false,flipaxis = false)

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
    MomentumAndPolarAngleDistribution1DPlot(type,sol,PhaseSpace,species,coordinates,time_idx;...)

Plots the distribution function given by `sol` of a given particle `species` as a function of momentum ``p`` and polar angle ``u`` at a given `time_idx` at three `coordinates` locations. 

The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.

Static arguments:
- `coordinates`: a NOT OPTIONAL `tuple` of three `tuples` of x,y,z coordinate indicies.
- `time_idx`: a NOT OPTIONAL index `Int64` of the time index of `sol` for which to plot the momentum distributions.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "MomentumAndPolarAngleDistribution.mp4".

"""
function MomentumAndPolarAngleDistribution1DPlot(type::Static,sol,PhaseSpace::PhaseSpaceStruct,species::String,coordinates::Tuple{Tuple{Int64,Int64,Int64},Tuple{Int64,Int64,Int64},Tuple{Int64,Int64,Int64}},time_idx::Int64;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime)

    CairoMakie.activate!(inline=true) # plot in vs code window

    with_theme(theme) do

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Time = PhaseSpace.Time

    species_index = findfirst(x->x==species,name_list)

    p_num = Momentum.px_num_list[species_index]  
    u_num = Momentum.py_num_list[species_index]
    h_num = Momentum.pz_num_list[species_index]
    dp = Grids.dpx_list[species_index]
    du = Grids.dpy_list[species_index]
    meanp = Grids.mpx_list[species_index]
    meanu = Grids.mpy_list[species_index]
    p_r = Grids.pxr_list[species_index]
    u_r = Grids.pyr_list[species_index]
    mass = Grids.mass_list[species_index]

    fig = Figure(size=(576,276)) # 8:3 aspect ratio

    (x_idx,y_idx,z_idx) = coordinates[1]
    f1 = copy(LocationSpeciesToStateVector(sol.f[time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
    (x_idx,y_idx,z_idx) = coordinates[2]
    f2 = copy(LocationSpeciesToStateVector(sol.f[time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
    (x_idx,y_idx,z_idx) = coordinates[3]
    f3 = copy(LocationSpeciesToStateVector(sol.f[time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

    #dis1 = dropdims(sum(reshape(sol.f[t_idx[1]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis2 = dropdims(sum(reshape(sol.f[t_idx[2]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)
    #dis3 = dropdims(sum(reshape(sol.f[t_idx[3]].x[species_index],(p_num,u_num,h_num)),dims=3),dims=3)

    dis1 = dropdims(sum(reshape(f1,(p_num,u_num,h_num)),dims=3),dims=3)
    dis2 = dropdims(sum(reshape(f2,(p_num,u_num,h_num)),dims=3),dims=3)
    dis3 = dropdims(sum(reshape(f3,(p_num,u_num,h_num)),dims=3),dims=3)

    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dpdu = f / dpdu and p^order * dN/dpdu = f * mp^order / dpdu
    for px in 1:p_num, py in 1:u_num
        dis1[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis2[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
        dis3[px,py] *= (meanp[px]^(order)) / dp[px] / du[py]
    end
    replace!(dis1,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis2,0.0 => NaN) # replace Inf with NaN for plotting
    replace!(dis3,0.0 => NaN) # replace Inf with NaN for plotting

    max_dis = maximum(x for x in [dis1; dis2; dis3] if !isnan(x))
    min_dis = minimum(x for x in [dis1; dis2; dis3] if !isnan(x))
    col_range = (log10(max_dis)-16.0,log10(max_dis))

    ax1 = PolarAxis(fig[1,1+1],theta_0=-pi/2,direction=-1,width=176)
    ax1.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax1,0,pi)

    ax2 = PolarAxis(fig[1,2+1],theta_0=-pi/2,direction=-1,width=176)
    ax2.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax2,0,pi)

    ax3 = PolarAxis(fig[1,3+1],theta_0=-pi/2,direction=-1,width=176)
    ax3.radius_at_origin = log10(p_r[1])-1.0
    thetalims!(ax3,0,pi)

    u_as_theta_grid = zeros(Float64,length(u_r))
    u_as_theta_grid_tick_values = Vector{Float64}(-1:0.5:1)
    u_as_theta_grid_tick_values_string = string.(u_as_theta_grid_tick_values)
    u_as_theta_grid_tick_values_string[end] = "u=1.0" 
    u_as_theta_grid_tick_locations = zeros(Float64,length(u_as_theta_grid_tick_values))
    @. u_as_theta_grid = pi - pi * (u_r+1)/2 # convert u grid to a set of theta values such that u can be plotted as polar angle
    @. u_as_theta_grid_tick_locations = pi - pi * (u_as_theta_grid_tick_values+1)/2 # convert u grid ticks to a set of theta values such that u can be plotted as polar angle

    hm1 = heatmap!(ax1,u_as_theta_grid,log10.(p_r),log10.(dis1'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
    hm2 = heatmap!(ax2,u_as_theta_grid,log10.(p_r),log10.(dis2'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))
    hm3 = heatmap!(ax3,u_as_theta_grid,log10.(p_r),log10.(dis3'),colormap=theme.colormap_var,colorrange=col_range,colorscale=x->asinh(x-log10(max_dis)))

    translate!(hm1,0,0,-100)
    translate!(hm2,0,0,-100)
    translate!(hm3,0,0,-100)

    rlims!(ax1,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax2,log10(p_r[1]),log10(p_r[end])+1.0)
    rlims!(ax3,log10(p_r[1]),log10(p_r[end])+1.0)
    ax1.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    ax2.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    ax3.thetaticks = (u_as_theta_grid_tick_locations,u_as_theta_grid_tick_values_string)
    #hidethetadecorations!(ax1, grid=false)
    #hidethetadecorations!(ax2, grid=false)
    #hidethetadecorations!(ax3, grid=false)

    if order == 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order == -2
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    elseif order != 1
        Colorbar(fig[1,1],hm1,label=L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$",flipaxis=false,height=Relative(0.75),tellheight=false)
    end

    t_unit_string = TimeUnits()

    pt = 4/3
    text!(ax1,L"$\log_{10}\left(p\,[m_ec]\right)$",position=(-3.05,log10(p_r[end])),rotation=pi/2,fontsize=9pt)
    ax1_label=Axis(fig[2,2])
    ax2_label=Axis(fig[2,3])
    ax3_label=Axis(fig[2,4])
    hidedecorations!(ax1_label)
    hidedecorations!(ax2_label)
    hidedecorations!(ax3_label)
    hidespines!(ax1_label)
    hidespines!(ax2_label)
    hidespines!(ax3_label)
    text!(ax1_label,L"$x,y,z=%$(coordinates[1])$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))
    text!(ax2_label,L"$x,y,z=%$(coordinates[2])$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))
    text!(ax3_label,L"$x,y,z=%$(coordinates[3])$",space=:relative,position=(0.5,0.5),fontsize=10pt,align=(:center,:center))

    colsize!(fig.layout,1,Relative(0.1))
    colsize!(fig.layout,2,Relative(0.3))
    colsize!(fig.layout,3,Relative(0.3))
    colsize!(fig.layout,4,Relative(0.3))
    rowsize!(fig.layout,2,Relative(0.06))
    rowgap!(fig.layout,1,0.0)
    
    return fig

    end # with_theme

end




################# Plots that only work with Paraboloidal coordinates ############### 
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