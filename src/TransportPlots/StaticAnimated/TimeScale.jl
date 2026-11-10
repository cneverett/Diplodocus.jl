"""
    ImprovedTimeScalePlot(PhaseSpace, sol,species,t_idx,EmiInteractions,BinInteractions,Forces;advection=true,dt=1.0)

Function plots the time scales for gain and loss rates of a simulation given the system solution `sol` at a time index `t_idx` for a set of `species`, and the desired sets of emissive interactions `EmiInteractions`, binary interactions `BinInteractions`, and forces `Forces`, with an to also display the advection timescales and time step.

Optional arguments:
- `theme`: the colour theme to use for the plot, default is `DiplodocusDark()`.
- `dt`: default is `nothing`. !nothing then the time step `dt` (in code units) will also be plotted as a horizontal line.
- `x_idx`, `y_idx`, `z_idx`: default is `1`. The spatial indices of the cell to plot.
- `u_indices`: default is `nothing`. If nothing then only the minimum angular timescale will be plotted, otherwise a vector of angle indices can be provided to plot the timescales for specific angles.
- `plot_limits`: default is `(nothing,nothing)`. The limits of the plot axes.
- `rate`: default is `false`. The default of `false` plots the gain and loss timescales of the specified interactions (τ), if set to `true` the rate will be plotted instead (1/τ).
- `p_timescale`: default is `false`. If `true`, the plot will show the timescale for momentum loss rather than that of number loss.
"""
function ImprovedTimeScalePlot(PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,t_idx::Int64,species::Vector{String},DataDirectory::String;forces::Union{Nothing,Vector{AbstractForce}}=nothing,emissions::Union{Nothing,Vector{EmissiveInteraction}}=nothing,binaries::Union{Nothing,Vector{BinaryInteraction}}=nothing,advection=true,dt::Union{Nothing,Float64}=nothing,theme=DiplodocusDark(),x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1,u_indices::Union{Nothing,Vector{Int64}}=nothing,y_limits_time=(nothing,nothing),y_limits_rate=(nothing,nothing),y_limits_dis=(nothing,nothing),p_timescale::Bool=false,net_rate::Bool=false,rate_order::Int64=2,figsize::Tuple{Int64,Int64}=(576,1080))

    CairoMakie.activate!(inline=true) # plot in vs code window
    with_theme(theme) do

    fig = Figure(size=figsize) 

    Momentum = PhaseSpace.Momentum
    Spacetime = PhaseSpace.Spacetime
    Grids = PhaseSpace.Grids
    dt0 = Spacetime.dt0
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    n_space = x_num * y_num * z_num
    name_list = PhaseSpace.name_list
    p_num_list = Momentum.px_num_list
    u_num_list = Momentum.py_num_list
    h_num_list = Momentum.pz_num_list
    n_momentum = sum(p_num_list.*u_num_list.*h_num_list)
    dp_list = Grids.dpx_list
    mp_list = Grids.mpx_list
    pr_list = Grids.pxr_list
    mu_list = Grids.mpy_list
    ur_list = Grids.pyr_list
    CHAR_n = PhaseSpace.Characteristic.CHAR_number_density
    CHAR_e = PhaseSpace.Characteristic.CHAR_energy_density
    CONST_c = DiplodocusTransport.CONST_c

    if rate_order != 2 
        @warn "Rate order is not 2, this may lead to incorrect units for the rate plot axes."
    end

    xlab = L"$\log_{10}\left(p\,[m_ec]\right)$" 
    ylab_rate_gain = L"\log_{10}\left(R_\text{gain}\right)"#L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N_\text{gain}}{\mathrm{d}p\mathrm{d}V\mathrm{d}t}\,\left[\text{m}^{-3}\,\left(m_e c\right)\,\left(r_\text{s}/c\right)^{-1}\right]\right)$"
    ylab_rate_loss = L"\log_{10}\left(R_\text{loss}\right)"#L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N_\text{loss}}{\mathrm{d}p\mathrm{d}V\mathrm{d}t}\,\left[\text{m}^{-3}\,\left(m_e c\right)\,\left(r_\text{s}/c\right)^{-1}\right]\right)$"
    ylab_time_gain = L"$\log_{10}\left(\tau_\text{gain}\right)$"# L"$\log_{10}\left(\tau_\text{gain} \,\left[\left(r_\text{s}/c\right)^{-1}\right]\right)$"
    ylab_time_loss = L"$\log_{10}\left(\tau_\text{loss}\right)$"#L"$\log_{10}\left(\tau_\text{loss} \,\left[\left(r_\text{s}/c\right)^{-1}\right]\right)$"

    ylab_dis = L"\log_{10}\left(\text{Mom. dis.}\right)"#L"$\log_{10}\left(p^{2}\,\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}V}\,\left[\,\text{m}^{-3}\,\left(m_e c\right)\right]\right)$"

    xlim1 = log10.((pr_list[1][1],pr_list[1][end]))
    xlim2 = log10.((pr_list[2][1],pr_list[2][end]))

    #=ax1_dis = Axis(fig[1,1],xlabel=xlab,ylabel=ylab_dis,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax1_dis.limits = (xlim1,y_limits_dis)
    hidexdecorations!(ax1_dis,grid=false)
    ax2_dis = Axis(fig[1,2],xlabel=xlab,ylabel=ylab_dis,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    hideydecorations!(ax2_dis,grid=false)
    hidexdecorations!(ax2_dis,grid=false)
    ax2_dis.limits = (xlim2,y_limits_dis)=#

    ax1_rate_gain = Axis(fig[1,1],xlabel=xlab,ylabel=ylab_rate_gain,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax1_rate_gain.limits = (xlim1,y_limits_rate)
    hidexdecorations!(ax1_rate_gain,grid=false)
    ax2_rate_gain = Axis(fig[1,2],xlabel=xlab,ylabel=ylab_rate_gain,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    hideydecorations!(ax2_rate_gain,grid=false)
    hidexdecorations!(ax2_rate_gain,grid=false)
    ax2_rate_gain.limits = (xlim2,y_limits_rate)

    ax1_rate_loss = Axis(fig[2,1],xlabel=xlab,ylabel=ylab_rate_loss,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax1_rate_loss.limits = (xlim1,y_limits_rate)
    hidexdecorations!(ax1_rate_loss,grid=false)
    ax2_rate_loss = Axis(fig[2,2],xlabel=xlab,ylabel=ylab_rate_loss,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax2_rate_loss.limits = (xlim2,y_limits_rate)
    hideydecorations!(ax2_rate_loss,grid=false)
    hidexdecorations!(ax2_rate_loss,grid=false)

    ax1_time_gain = Axis(fig[3,1],xlabel=xlab,ylabel=ylab_time_gain,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax1_time_gain.limits = (xlim1,y_limits_time)
    hidexdecorations!(ax1_time_gain,grid=false)
    ax2_time_gain = Axis(fig[3,2],xlabel=xlab,ylabel=ylab_time_gain,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    hideydecorations!(ax2_time_gain,grid=false)
    hidexdecorations!(ax2_time_gain,grid=false)
    ax2_time_gain.limits = (xlim2,y_limits_time)

    ax1_time_loss = Axis(fig[4,1],xlabel=xlab,ylabel=ylab_time_loss,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    ax1_time_loss.limits = (xlim1,y_limits_time)
    ax2_time_loss = Axis(fig[4,2],xlabel=xlab,ylabel=ylab_time_loss,aspect=DataAspect(),xticks=-10:5:5,yticks=-20:5:20)
    hideydecorations!(ax2_time_loss,grid=false)
    ax2_time_loss.limits = (xlim2,y_limits_time)

    f = copy(sol.f[t_idx])
    # remove small values 
    #@. f = ifelse(f < 1e-10, zero(eltype(f)), f) 
    df_gain = zeros(eltype(f),size(f))
    df_loss = zeros(eltype(f),size(f))

    # we work with rates as these add linearly 
    rate_gain = zeros(eltype(f),size(f))
    rate_loss = zeros(eltype(f),size(f))

    # calculate number of interactions for colors 
    # colors are then ordered as advection, forces, emissions, binaries
    num_interactions = 0 
    if !isnothing(advection)
        num_interactions += 1
        num_advection = 1
    else
        num_advection = 0
    end
    if !isnothing(forces)
        num_interactions += length(forces)
        num_force = length(forces)
    else
        num_force = 0
    end
    if !isnothing(emissions)
        num_interactions += length(emissions)
        num_emission = length(emissions)
    else
        num_emission = 0
    end
    if !isnothing(binaries)
        num_interactions += length(binaries)
        num_binary = length(binaries)
    else
        num_binary = 0
    end

    # linestyles for angles 
    linestyles = [(:dash,:dense),:solid,(:dot),(:dashdot),(:dashdotdot),(:dash,:dense),:solid,(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []
    legend_elements_angle = []
    line_labels_angle = []

    if !isnothing(u_indices)
        for (idx,u_idx) in enumerate(u_indices)
            u_val = PhaseSpace.Grids.mpy_list[1][u_idx]
            color = theme.textcolor[]
            linestyle = linestyles[idx]
            pushfirst!(legend_elements_angle,LineElement(color = color, linestyle = linestyle,linewidth = 2.0))
            pushfirst!(line_labels_angle,L"$u=%$(round(u_val,sigdigits=2))$")
        end
    end

    if !isnothing(dt)

        hlines!(ax1_time_gain,log10(dt),color=theme.textcolor[],linewidth=1.0,linestyle=:dashdot)
        hlines!(ax2_time_gain,log10(dt),color=theme.textcolor[],linewidth=1.0,linestyle=:dashdot)
        hlines!(ax1_time_loss,log10(dt),color=theme.textcolor[],linewidth=1.0,linestyle=:dashdot)
        hlines!(ax2_time_loss,log10(dt),color=theme.textcolor[],linewidth=1.0,linestyle=:dashdot)

        pushfirst!(legend_elements,LineElement(color = theme.textcolor[], linestyle = :dashdot,linewidth = 2.0))
        pushfirst!(line_labels,L"$\text{Time Step}$")

    end

    # distribution plots 
    #=for (species_idx, species_name) in enumerate(species)

        color = theme.textcolor[]

        ax_dis = species_idx == 1 ? ax1_dis : ax2_dis

        species_index = findfirst(x->x==species_name,name_list)
        name = name_list[species_index]
        p_num = p_num_list[species_index]
        u_num = u_num_list[species_index]
        h_num = h_num_list[species_index]
        dp = dp_list[species_index]
        mp = mp_list[species_index]
        mu = mu_list[species_index]
        ur = ur_list[species_index]
        if Momentum.px_grid_list[species_index] == "l"
            mp_plot = log10.(mp)
        elseif Momentum.px_grid_list[species_index] == "u"
            mp_plot = mp
        end

        f1D = copy(LocationSpeciesToStateVector(f,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
        f3D = reshape(f1D,(p_num,u_num,h_num))

        # remove infinite values 
        replace!(f1D,0.0=>NaN)
        replace!(f1D,-0.0=>NaN)

        # angle plots 
        if isnothing(u_indices)
            #= minimum angular timescale =#
            pdNdp = dropdims(sum(f3D, dims=(2,3)),dims=(2,3))

            # timescale is timescale for particle losses not particle energy losses, to get that we need to scale by dp/mp (i.e. dp/p for τ and p/dp for rate)
            if p_timescale
                @. pdNdp *= mp^2 / dp * CHAR_n # convert to (m_e c)^{order-1} m^-3 / time units 
            end

            # plot lines

            scatterlines!(ax_dis,mp_plot,log10.(abs.(pdNdp)),linewidth=2.0,color=color,markersize=0.0)

        else

            pdNdp = dropdims(sum(f3D, dims=(3)),dims=(3))

            if p_timescale
                @. pdNdp *= mp^2 / dp * CHAR_n
            end

            # plot lines for each angle in u_indices

            for (idx,u_idx) in enumerate(u_indices)

                scatterlines!(ax_dis,mp_plot,log10.(pdNdp[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

            end
        end

    end # species loop =#

    if advection

        color = theme.colormap[][(1-1) / (num_interactions-1)] # first color for advection
        println("Advection color = $(color)")

        pushfirst!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
        pushfirst!(line_labels,L"$\text{Advection}$")

        #= 

            A * df = - Flux * f 
            df/dt = - Flux * f / A / dt0

            rate = (df/dt) = (- Flux * f / A / dt0) = - (Flux * f) / (A * dt0) 
            time = f / (df/dt)

        =#

        Flux = DiplodocusTransport.BuildFluxMatrices(PhaseSpace,AbstractForce[]) # only generates the advection flux matrices

        LossFlux = spdiagm(0 => diag(Flux.X_Flux))      # loss terms are on the diagonal
        GainFlux = Flux.X_Flux - LossFlux  # gain terms are off diagonal 

        mul!(df_gain, GainFlux, f)
        mul!(df_loss, LossFlux, f)

        @. rate_gain = - df_gain / (Flux.Ap_Flux * dt)
        @. rate_loss = - df_loss / (Flux.Ap_Flux * dt0)

        if net_rate 
            rate_net = rate_gain + rate_loss
            @. rate_gain = ifelse(rate_net > 0.0, rate_net,0.0)
            @. rate_loss = ifelse(rate_net < 0.0, -rate_net,0.0)
        end

        # plotting specific cell for a specific species
        for (species_idx, species_name) in enumerate(species)

            ax_rate_gain = species_idx == 1 ? ax1_rate_gain : ax2_rate_gain
            ax_rate_loss = species_idx == 1 ? ax1_rate_loss : ax2_rate_loss
            ax_time_gain = species_idx == 1 ? ax1_time_gain : ax2_time_gain 
            ax_time_loss = species_idx == 1 ? ax1_time_loss : ax2_time_loss

            species_index = findfirst(x->x==species_name,name_list)
            name = name_list[species_index]
            p_num = p_num_list[species_index]
            u_num = u_num_list[species_index]
            h_num = h_num_list[species_index]
            dp = dp_list[species_index]
            mp = mp_list[species_index]
            mu = mu_list[species_index]
            ur = ur_list[species_index]
            if Momentum.px_grid_list[species_index] == "l"
                mp_plot = log10.(mp)
            elseif Momentum.px_grid_list[species_index] == "u"
                mp_plot = mp
            end

            rate_gain_location1D = LocationSpeciesToStateVector(rate_gain,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            rate_gain_location3D = reshape(rate_gain_location1D,(p_num,u_num,h_num))
            rate_loss_location1D = LocationSpeciesToStateVector(rate_loss,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
            rate_loss_location3D = reshape(rate_loss_location1D,(p_num,u_num,h_num))

            f_location1D = copy(LocationSpeciesToStateVector(f,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
            f_location3D = reshape(f_location1D,(p_num,u_num,h_num))

            # remove infinite values 
            replace!(rate_gain_location3D,0.0=>NaN)
            replace!(rate_loss_location3D,0.0=>NaN)
            replace!(rate_gain_location3D,-0.0=>NaN)
            replace!(rate_loss_location3D,-0.0=>NaN)

            # angle plots 
            if isnothing(u_indices)
                #= minimum angular timescale =#
                rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(2,3)),dims=(2,3))
                rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(2,3)),dims=(2,3))

                f_location = dropdims(sum(f_location3D, dims=(2,3)),dims=(2,3))

                τ_gain_location = f_location ./ rate_gain_location 
                τ_loss_location = f_location ./ rate_loss_location

                # timescale is timescale for particle losses not particle energy losses, to get that we need to scale by dp/mp (i.e. dp/p for τ and p/dp for rate)
                if p_timescale
                    @. rate_gain_location *= mp^rate_order / dp * CHAR_n # convert to J c^-1 m^-3 / time units 
                    @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                end

                # plot lines

                scatterlines!(ax_rate_gain,mp_plot,log10.(abs.(rate_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                scatterlines!(ax_time_gain,mp_plot,log10.(abs.(τ_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location)),linewidth=2.0,color=color,markersize=0.0)

            else

                rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(3)),dims=(3))
                rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(3)),dims=(3))

                f_location = dropdims(sum(f_location3D, dims=(3)),dims=(3))

                τ_gain_location = f_location ./ rate_gain_location
                τ_loss_location = f_location ./ rate_loss_location

                if p_timescale
                    @. rate_gain_location *= mp^rate_order / dp * CHAR_n
                    @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                end

                # plot lines for each angle in u_indices

                for (idx,u_idx) in enumerate(u_indices)

                    scatterlines!(ax_rate_gain,mp_plot,log10.(rate_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                    scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                    scatterlines!(ax_time_gain,mp_plot,log10.(τ_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                    scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                end
            end

        end # species loop 

    end # advection 

    if !isnothing(forces)

        #= 

            A * df = - Flux * f 
            df/dt = - Flux * f / A / dt0

            rate = (df/dt) = - Flux * f / A / dt0 = - (Flux * f) / (A * dt0 ) 
            time = f / (df/dt)

        =#

        for (force_idx,force) in enumerate(forces) 

            color = theme.colormap[][(force_idx+num_advection-1) / (num_interactions-1)]

            pushfirst!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
            if force isa CoordinateForce
                pushfirst!(line_labels,L"$\text{Ricci{-}Rot. Coeffs.}$")
            elseif force isa FirstOrderGuidingCentre
                pushfirst!(line_labels,L"$\text{1$^{\text{st}}$-Order Lorentz}$")
            elseif force isa SyncRadReact
                pushfirst!(line_labels,L"$\text{Synchrotron}$")
            end

            Flux = DiplodocusTransport.BuildFluxMatrices(PhaseSpace,AbstractForce[force]) # only generates the force flux

            LossFlux = spdiagm(0 => diag(Flux.P_Flux))      # loss terms are on the diagonal
            GainFlux = Flux.P_Flux - LossFlux  # gain terms are off diagonal 

            println("sum(GainFlux) = $(sum(GainFlux)), sum(LossFlux) = $(sum(LossFlux))")            

            mul!(df_gain, GainFlux, f)
            mul!(df_loss, LossFlux, f)

            @. rate_gain = - df_gain / (Flux.Ap_Flux * dt0 #=* f=#)
            @. rate_loss = - df_loss / (Flux.Ap_Flux * dt0 #=* f=#)

            if net_rate 
                rate_net = rate_gain + rate_loss
                @. rate_gain = ifelse(rate_net > 0.0, rate_net,0.0)
                @. rate_loss = ifelse(rate_net < 0.0, -rate_net,0.0)
            end

            # plotting specific cell for a specific species
            for (species_idx, species_name) in enumerate(species)

                ax_rate_gain = species_idx == 1 ? ax1_rate_gain : ax2_rate_gain
                ax_rate_loss = species_idx == 1 ? ax1_rate_loss : ax2_rate_loss
                ax_time_gain = species_idx == 1 ? ax1_time_gain : ax2_time_gain 
                ax_time_loss = species_idx == 1 ? ax1_time_loss : ax2_time_loss

                species_index = findfirst(x->x==species_name,name_list)
                name = name_list[species_index]
                p_num = p_num_list[species_index]
                u_num = u_num_list[species_index]
                h_num = h_num_list[species_index]
                dp = dp_list[species_index]
                mp = mp_list[species_index]
                mu = mu_list[species_index]
                ur = ur_list[species_index]
                if Momentum.px_grid_list[species_index] == "l"
                    mp_plot = log10.(mp)
                elseif Momentum.px_grid_list[species_index] == "u"
                    mp_plot = mp
                end

                rate_gain_location1D = LocationSpeciesToStateVector(rate_gain,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_gain_location3D = reshape(rate_gain_location1D,(p_num,u_num,h_num))
                rate_loss_location1D = LocationSpeciesToStateVector(rate_loss,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_loss_location3D = reshape(rate_loss_location1D,(p_num,u_num,h_num))

                f_location1D = copy(LocationSpeciesToStateVector(f,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
                f_location3D = reshape(f_location1D,(p_num,u_num,h_num))

                # remove infinite values 
                replace!(rate_gain_location3D,-0.0=>NaN)
                replace!(rate_loss_location3D,-0.0=>NaN)
                replace!(rate_gain_location3D,0.0=>NaN)
                replace!(rate_loss_location3D,0.0=>NaN)

                # angle plots 
                if isnothing(u_indices)
                    #= minimum angular timescale =#
                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(2,3)),dims=(2,3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(2,3)),dims=(2,3))

                    f_location = dropdims(sum(f_location3D, dims=(2,3)),dims=(2,3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^rate_order / dp * CHAR_n 
                        @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                    end

                    # plot lines
                    scatterlines!(ax_rate_gain,mp_plot,log10.(abs.(rate_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                    scatterlines!(ax_time_gain,mp_plot,log10.(abs.(τ_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                else

                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(3)),dims=(3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(3)),dims=(3))

                    f_location = dropdims(sum(f_location3D, dims=(3)),dims=(3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^rate_order / dp * CHAR_n
                        @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                    end
                    # plot lines for each angle in u_indices

                    for (idx,u_idx) in enumerate(u_indices)

                        scatterlines!(ax_rate_gain,mp_plot,log10.(rate_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                        scatterlines!(ax_time_gain,mp_plot,log10.(τ_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx]) 

                    end
                end

            end # species loop 

        end

    end # forces

    if !isnothing(emissions)

        #= 

            A * df = MEmi * f 
            df/dt = MEmi * f / A / dt0

            rate = (df/dt) = MEmi * f / A / dt0 = (MEmi * f) / (A * dt0) 
            time = f / (df/dt)

        =#

        for (emission_idx,emission) in enumerate(emissions) 

            color = theme.colormap[][(emission_idx+num_force+num_advection-1) / (num_interactions-1)]

            pushfirst!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
            if emission.EmiName == "Sync"
                pushfirst!(line_labels,L"$\text{Synchrotron}$")                
            end

            Bin = BuildBinaryMatrices(PhaseSpace,BinaryInteraction[],Int[],DataDirectory) # empty binary arrays but needed for domain
            Emi = BuildEmissionMatrices(PhaseSpace,[emission],DataDirectory,Bin;Emi_corrected=true)

            Flux = DiplodocusTransport.BuildFluxMatrices(PhaseSpace,AbstractForce[])

            # Build df_gain and df_loss based on spatially local emission matrix
            off_space = (x_idx-1)*y_num*z_num+(y_idx-1)*z_num+z_idx-1
            MEmi = Emi.M_Emi[off_space+1]
            LossEmi = spdiagm(0 => diag(MEmi))      # loss terms are on the diagonal
            GainEmi = MEmi - LossEmi  # gain terms are off diagonal
            
            start_idx = n_momentum*off_space+1
            end_idx = n_momentum*(off_space+1)
            f_local = @view(f[start_idx:end_idx])
            df_gain_local = @view(df_gain[start_idx:end_idx])
            df_loss_local = @view(df_loss[start_idx:end_idx])

            mul!(df_gain_local, GainEmi, f_local)
            mul!(df_loss_local, LossEmi, f_local)

            # build global rate_gain and rate_loss but only local cells will be non-zero
            fill!(rate_gain,0.0)
            fill!(rate_loss,0.0)
            @. rate_gain = df_gain / (Flux.Ap_Flux * dt0 #=* f=#)
            @. rate_loss = df_loss / (Flux.Ap_Flux * dt0 #=* f=#)

            if net_rate 
                rate_net = rate_gain + rate_loss
                @. rate_gain = ifelse(rate_net > 0.0, rate_net,0.0)
                @. rate_loss = ifelse(rate_net < 0.0, -rate_net,0.0)
            end

            # plotting specific cell for a specific species
            for (species_idx, species_name) in enumerate(species)

                ax_rate_gain = species_idx == 1 ? ax1_rate_gain : ax2_rate_gain
                ax_rate_loss = species_idx == 1 ? ax1_rate_loss : ax2_rate_loss
                ax_time_gain = species_idx == 1 ? ax1_time_gain : ax2_time_gain 
                ax_time_loss = species_idx == 1 ? ax1_time_loss : ax2_time_loss

                species_index = findfirst(x->x==species_name,name_list)
                name = name_list[species_index]
                p_num = p_num_list[species_index]
                u_num = u_num_list[species_index]
                h_num = h_num_list[species_index]
                dp = dp_list[species_index]
                mp = mp_list[species_index]
                mu = mu_list[species_index]
                ur = ur_list[species_index]
                if Momentum.px_grid_list[species_index] == "l"
                    mp_plot = log10.(mp)
                elseif Momentum.px_grid_list[species_index] == "u"
                    mp_plot = mp
                end

                rate_gain_location1D = LocationSpeciesToStateVector(rate_gain,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_gain_location3D = reshape(rate_gain_location1D,(p_num,u_num,h_num))
                rate_loss_location1D = LocationSpeciesToStateVector(rate_loss,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_loss_location3D = reshape(rate_loss_location1D,(p_num,u_num,h_num))

                f_location1D = copy(LocationSpeciesToStateVector(f,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))
                f_location3D = reshape(f_location1D,(p_num,u_num,h_num))

                # remove infinite values 
                replace!(rate_gain_location3D,-0.0=>NaN)
                replace!(rate_loss_location3D,-0.0=>NaN)
                replace!(rate_gain_location3D,0.0=>NaN)
                replace!(rate_loss_location3D,0.0=>NaN)

                # angle plots 
                if isnothing(u_indices)
                    #= minimum angular timescale =#
                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(2,3)),dims=(2,3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(2,3)),dims=(2,3))

                    f_location = dropdims(sum(f_location3D, dims=(2,3)),dims=(2,3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^rate_order / dp * CHAR_n
                        @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                    end

                    # plot lines
                    scatterlines!(ax_rate_gain,mp_plot,log10.(abs.(rate_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                    scatterlines!(ax_time_gain,mp_plot,log10.(abs.(τ_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                else

                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(3)),dims=(3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(3)),dims=(3))

                    f_location = dropdims(sum(f_location3D, dims=(3)),dims=(3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^rate_order / dp * CHAR_n 
                        @. rate_loss_location *= mp^rate_order / dp * CHAR_n
                    end

                    # plot lines for each angle in u_indices
                    for (idx,u_idx) in enumerate(u_indices)

                        scatterlines!(ax_rate_gain,mp_plot,log10.(rate_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                        scatterlines!(ax_time_gain,mp_plot,log10.(τ_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                    end
                end

            end # species loop 

        end

    end # emissions

    if !isnothing(binaries)

        #= 

            A * df = MBin * f * f * Vol
            df/dt = MBin * f * f * Vol / A / dt0

            rate = (df/dt) = MBin * f * f * Vol / A / dt0 = (MBin * f * f * Vol) / (A * dt0) 
            time = f / (df/dt)

        =#

        for (binary_idx,binary) in enumerate(binaries) 

            fill!(df_gain,0.0)
            fill!(df_loss,0.0)

            color = theme.colormap[][(binary_idx+num_emission+num_force+num_advection-1) / (num_interactions-1)]

            pushfirst!(legend_elements,LineElement(color = color, linestyle = :solid,linewidth = 2.0))
            if binary.name1 == "Ele" && binary.name2 == "Pho" && binary.name3 == "Ele" && binary.name4 == "Pho"
                pushfirst!(line_labels,L"$\text{Compton}$")                
            elseif binary.name1 == "Pho" && binary.name2 == "Pho" && binary.name3 == "Ele" && binary.name4 == "Pos"
                pushfirst!(line_labels,L"$\text{Pair Prod.}$")
            end

            Bin_gain = BuildBinaryMatrices(PhaseSpace,BinaryInteraction[binary],Int[],DataDirectory;Bin_corrected=true,Bin_sparse=true,symmetric=true,gain_loss="gain") 
            GainBin = Bin_gain.M_Bin
            Bin_loss = BuildBinaryMatrices(PhaseSpace,BinaryInteraction[binary],Int[],DataDirectory;Bin_corrected=true,Bin_sparse=true,symmetric=true,gain_loss="loss") 
            LossBin = Bin_loss.M_Bin
            Flux = DiplodocusTransport.BuildFluxMatrices(PhaseSpace,AbstractForce[]) # needed for volume element

            # Build df_gain and df_loss based on spatially local emission matrix
            off_space = (x_idx-1)*y_num*z_num+(y_idx-1)*z_num+z_idx-1
            Vol = Flux.Vol[off_space+1]
            M_Bin_Mul_Step = zeros(n_momentum,n_momentum)
            M_Bin_Mul_Step_reshape = reshape(M_Bin_Mul_Step,n_momentum^2)
            
            start_idx = n_momentum*off_space+1
            end_idx = n_momentum*(off_space+1)
            f_local = @view(f[start_idx:end_idx])
            df_gain_local = @view(df_gain[start_idx:end_idx])
            df_loss_local = @view(df_loss[start_idx:end_idx])

            println("maximum f_local = $(maximum(f_local)), minimum f_local = $(minimum(f_local))")

            mul!(M_Bin_Mul_Step_reshape, GainBin, f_local)
            mul!(df_gain_local, M_Bin_Mul_Step, f_local)

            println("maximum df_gain_local = $(maximum(df_gain_local)), minimum df_gain_local = $(minimum(df_gain_local))")

            mul!(M_Bin_Mul_Step_reshape, LossBin, f_local)
            mul!(df_loss_local, M_Bin_Mul_Step, f_local)

            println("maximum df_loss_local = $(maximum(df_loss_local)), minimum df_loss_local = $(minimum(df_loss_local))")

            # build global rate_gain and rate_loss but only local cells will be non-zero
            fill!(rate_gain,0.0)
            fill!(rate_loss,0.0)
            @. rate_gain = df_gain * Vol / (Flux.Ap_Flux * dt0 #=* f=#)
            @. rate_loss = df_loss * Vol / (Flux.Ap_Flux * dt0 #=* f=#)

            println("sum(rate_gain) = $(sum(rate_gain)), sum(rate_loss) = $(sum(rate_loss))")

            if net_rate 
                rate_net = rate_gain + rate_loss
                @. rate_gain = ifelse(rate_net > 0.0, rate_net,0.0)
                @. rate_loss = ifelse(rate_net < 0.0, -rate_net,0.0)
            end

            println("sum(rate_gain) = $(sum(rate_gain)), sum(rate_loss) = $(sum(rate_loss))")

            # plotting specific cell for a specific species
            for (species_idx, species_name) in enumerate(species)

                ax_rate_gain = species_idx == 1 ? ax1_rate_gain : ax2_rate_gain
                ax_rate_loss = species_idx == 1 ? ax1_rate_loss : ax2_rate_loss
                ax_time_gain = species_idx == 1 ? ax1_time_gain : ax2_time_gain 
                ax_time_loss = species_idx == 1 ? ax1_time_loss : ax2_time_loss

                species_index = findfirst(x->x==species_name,name_list)
                name = name_list[species_index]
                p_num = p_num_list[species_index]
                u_num = u_num_list[species_index]
                h_num = h_num_list[species_index]
                dp = dp_list[species_index]
                mp = mp_list[species_index]
                mu = mu_list[species_index]
                ur = ur_list[species_index]
                if Momentum.px_grid_list[species_index] == "l"
                    mp_plot = log10.(mp)
                elseif Momentum.px_grid_list[species_index] == "u"
                    mp_plot = mp
                end

                rate_gain_location1D = LocationSpeciesToStateVector(rate_gain,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_gain_location3D = reshape(rate_gain_location1D,(p_num,u_num,h_num))
                rate_loss_location1D = LocationSpeciesToStateVector(rate_loss,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                rate_loss_location3D = reshape(rate_loss_location1D,(p_num,u_num,h_num))

                f_location1D = LocationSpeciesToStateVector(f,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx)
                f_location3D = reshape(f_location1D,(p_num,u_num,h_num))

                # remove infinite values 
                replace!(rate_gain_location3D,-0.0=>NaN)
                replace!(rate_loss_location3D,-0.0=>NaN)
                replace!(rate_gain_location3D,0.0=>NaN)
                replace!(rate_loss_location3D,0.0=>NaN)

                # angle plots 
                if isnothing(u_indices)
                    #= minimum angular timescale =#
                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(2,3)),dims=(2,3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(2,3)),dims=(2,3))

                    f_location = dropdims(sum(f_location3D, dims=(2,3)),dims=(2,3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^2 / dp * CHAR_n
                        @. rate_loss_location *= mp^2 / dp * CHAR_n
                    end

                    # plot lines
                    scatterlines!(ax_rate_gain,mp_plot,log10.(abs.(rate_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                    scatterlines!(ax_time_gain,mp_plot,log10.(abs.(τ_gain_location)),linewidth=2.0,color=color,markersize=0.0)
                    scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location)),linewidth=2.0,color=color,markersize=0.0)

                else

                    rate_gain_location = dropdims(sum(rate_gain_location3D, dims=(3)),dims=(3))
                    rate_loss_location = dropdims(sum(rate_loss_location3D, dims=(3)),dims=(3))

                    f_location = dropdims(sum(f_location3D, dims=(3)),dims=(3))

                    τ_gain_location = f_location ./ rate_gain_location
                    τ_loss_location = f_location ./ rate_loss_location

                    if p_timescale
                        @. rate_gain_location *= mp^2 / dp * CHAR_n
                        @. rate_loss_location *= mp^2 / dp * CHAR_n
                    end

                    # plot lines for each angle in u_indices

                    for (idx,u_idx) in enumerate(u_indices)

                        scatterlines!(ax_rate_gain,mp_plot,log10.(rate_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_rate_loss,mp_plot,log10.(abs.(rate_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                        scatterlines!(ax_time_gain,mp_plot,log10.(τ_gain_location[:,u_idx]),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])
                        scatterlines!(ax_time_loss,mp_plot,log10.(abs.(τ_loss_location[:,u_idx])),linewidth=2.0,color=color,markersize=0.0,linestyle=linestyles[idx])

                    end
                end

            end # species loop 

        end

    end # emissions

    # legends 
    if !isnothing(u_indices)
        Legend(fig[5,1],legend_elements_angle,line_labels_angle#=,tellwidth=false,tellheight=false=#,rowgap=-10)
    end
    Legend(fig[5,2],legend_elements,line_labels,nbanks = 2#=,tellwidth=false,tellheight=false=#,rowgap=-10)

    colgap!(fig.layout,1,0.0)
    rowgap!(fig.layout,1,0.0)
    rowgap!(fig.layout,2,0.0)
    rowgap!(fig.layout,3,0.0)
    rowgap!(fig.layout,4,0.0)

    colsize!(fig.layout,1,Aspect(1,1))
    colsize!(fig.layout,2,Aspect(1,19/10))

    #rowsize!(fig.layout,1,Fixed(0))
    rowsize!(fig.layout,5,Fixed(90))

    linkyaxes!(ax1_time_gain, ax2_time_gain)
    linkyaxes!(ax1_time_loss, ax2_time_loss)

    linkyaxes!(ax1_rate_gain, ax2_rate_gain)
    linkyaxes!(ax1_rate_loss, ax2_rate_loss)
    #linkyaxes!(ax1_dis, ax2_dis)

    linkxaxes!(ax1_time_gain, ax1_time_loss,ax1_rate_gain, ax1_rate_loss#=,ax1_dis=#)
    linkxaxes!(ax2_time_gain, ax2_time_loss,ax2_rate_gain, ax2_rate_loss#=,ax2_dis=#)

    resize_to_layout!(fig)

    return fig

    end # theme 

end