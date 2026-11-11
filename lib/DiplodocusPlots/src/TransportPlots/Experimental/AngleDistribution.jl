"""
    AngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::PlotType;angle_step=1,order=1,plot_limits=(nothing,nothing),theme=DiplodocusDark())

Plots the momentum averaged distribution function of of a given vector of particle `species` as a function of time given by the `sol` based on the conditions held in `PhaseSpace`. 
    
The plot can be either static, animated or interactive depending on the `type` argument. `Static` and `Animated` plots are generated using CairoMakie and are best for publications and presentations, while `Interactive` plots are generated using GLMakie and allow for user interaction with the plot.

Common arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `order`: the order of p in p^order * dN/dp dV, default is 1, i.e. number density spectrum. 2 is "energy" density spectrum.
- `TimeUnits`: a function that converts the time given in code units to the desired units for plotting
- `plot_limits`: the limits of the x and y axes, default is `(nothing,nothing)` which sets the limits automatically based on the data.
- `wide`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), default is `false` (single column 4:3 aspect ratio).
- `legend`: if `true`, a legend is added to the plot, default is `true`.
- `angle_step`: the angular grid step to reduce plotting lines.

Static arguments:
- `time_idx`: the time index to plot at.

Animated arguments:
- `framerate`: the frame rate of the animation, default is 12 fps.
- `filename`: the name of the file to save the animation to, default is "AngleDistribution.mp4".
- `figure`: default is `nothing`, which creates a new figure. If a figure is provided, the plot is added to that figure instead of creating a new one, this should be of the form of a tuple of `figure` and `time_idx` from the main plot.
"""
function AngleDistributionPlot(type::Static,PhaseSpace::PhaseSpaceStruct,sol,species::Vector{String},time_idx::Int64;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Tuple{Float64,String}=(1.0,"\\text{ units}"),plot_limits=(nothing,nothing),wide=false,legend=true,angle_step::Int64=1,x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1,total_line=false)

    CairoMakie.activate!(inline=true) # plot in vs code window

    CHAR_n = PhaseSpace.Characteristic.CHAR_number_density

    with_theme(theme) do

    if wide
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
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
    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab,aspect=DataAspect())
    ax.limits = plot_limits

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot,:dense),(:dashdotdot,:dense)]
    max_total = -Inf32
    p_min = Inf32
    p_max = -Inf32

    legend_elements = []
    line_labels = []
    legend_elements_angle = []
    line_labels_angle = []

    counter = 1

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

    p_min = min(p_min,p_r[1])
    p_max = max(p_max,p_r[end])

    f1D .= copy(LocationSpeciesToStateVector(sol.f[time_idx],PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

    f3D .= reshape(f1D,(p_num,u_num,h_num))
    @. f3D = f3D*(f3D!=Inf)
    # scale by order
    # f = dN/dpdudh * dpdudh therefore dN/dp = f / dp and p^order * dN/dp = f * mp^order / dp
    for px in 1:p_num, py in 1:u_num, pz in 1:h_num
        f3D[px,py,pz] = f3D[px,py,pz] * (meanp[px]^(order)) / dp[px] * CHAR_n
    end

    for i in 1:angle_step:u_num

        u_val = meanu[i]

        color = theme.colormap[][(u_val - u_r[1]) / (u_r[end] - u_r[1])]

        pdNdp = dropdims(sum(f3D, dims=(3)),dims=(3))[:,i]

        if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
            idx = findfirst(!iszero,pdNdp)
            lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = color,linestyle=linestyles[species_idx])
        else
            scatterlines!(ax,mp_plot,log10.(pdNdp),linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])
        end
        max_f = maximum(x for x in pdNdp if !isnan(x))
        max_total = max(max_f,max_total)

        if legend==true && counter == 1 
            u_val = meanu[i]
            color = theme.colormap[][(u_val - u_r[1]) / (u_r[end] - u_r[1])]
            pushfirst!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
            pushfirst!(line_labels_angle,L"$u=%$(round(u_val,sigdigits=2))$")
        end

    end # angle loop

        total = 1
        if total_line == true

            pdNdp = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))
            if sum(@. !isnan(pdNdp) * !isinf(pdNdp) * !iszero(pdNdp)) == 1 # there is only one valid position so scatterlines doesn't work
                idx = findfirst(!iszero,pdNdp)
                lines!(ax,[mp_plot[idx], mp_plot[idx]],[-20.0, log10(pdNdp[idx])],linewidth=2.0,color = theme.textcolor[],linestyle=linestyles[species_idx])
            else
                scatterlines!(ax,mp_plot,log10.(pdNdp),linewidth=2.0,color = theme.textcolor[],markersize=0.0,linestyle=linestyles[species_idx+1])
            end

            if total ==1 
                pushfirst!(legend_elements_angle,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx+1],linewidth = 2.0))
                pushfirst!(line_labels_angle,L"$\text{Total}$")
            end

            total +=1
        end

        counter += 1

        push!(legend_elements,LineElement(color = theme.textcolor[], linestyle = linestyles[species_idx],linewidth = 2.0))
        push!(line_labels,species_name)

    end # species loop 


    if legend
        axislegend(ax,legend_elements,line_labels,position = :rt)
        Legend(fig[1,2],legend_elements_angle,line_labels_angle)
    end

    if plot_limits == (nothing,nothing)
        xlims!(ax,(log10(p_min)-1.0,log10(p_max)+1.0))
        ylims!(ax,(log10(max_total)-9.0,log10(max_total)+1.0)) 
    end

    return fig

    end # with_theme

end

function AngleDistributionPlot(sol,species::Vector{String},PhaseSpace::PhaseSpaceStruct,type::Animated;theme=DiplodocusDark(),order::Int64=1,TimeUnits::Function=CodeToCodeUnitsTime,plot_limits=(nothing,nothing),wide=false,legend=true,angle_step::Int64=1,framerate=12,filename="AngleDistribution.mp4",figure=nothing)

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    xlab = L"$\log_{10}\left(p [m_ec]\right)$"
    if order == 1
        ylab = L"$\log_{10}\left(p\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}]\right)$"
    elseif order == -2
        ylab=L"$\log_{10}\left(\frac{\mathrm{d}N}{p^2\mathrm{d}p\mathrm{d}u\mathrm{d}V}\,[\text{m}^{-3}\left(m_ec\right)^{-3}]\right)$"
    elseif order != 1
        ylab = L"$\log_{10}\left(p^{%$(order)}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V} [\text{m}^{-3}\left(m_ec\right)^{%$(order-1)}]\right)$"
    end

    if isnothing(figure)
        time_idx = Observable(1) # index of the current time step
        t = @lift(TimeUnits(sol.t[$time_idx]))
        if wide
            fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
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

    for i in ceil(Int64,u_num/2):angle_step:u_num

        u_val = meanu[i]
        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]

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
        log10.(dropdims(sum(f3D, dims=(3)),dims=(3)))[:,i]
        end

        scatterlines!(ax,log10.(meanp),pdNdp,linewidth=2.0,color = color,markersize=0.0,linestyle=linestyles[species_idx])

        if legend==true && counter == 1 
        u_val = meanu[i]
        color = theme.colormap[][(u_val - u_r[ceil(Int64,u_num/2)]) / (u_r[end] - u_r[ceil(Int64,u_num/2)])]
        push!(legend_elements_angle,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        push!(line_labels_angle,L"$\theta=%$(round(acos(u_val)/pi,sigdigits=2)) \pi$")
        end

    end # angle loop 

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
        time_idxs = 1:length(sol.t)
        record(fig,filename,time_idxs,framerate=framerate,backend=CairoMakie) do frame
            print("generating frame: $frame \r")
            time_idx[] = frame
        end
    end

    end # with_theme

end