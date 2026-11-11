"""
    NumberDensityPlot1D(type,PhaseSpace,sol,species,x_dir,y_dir,z_dir;...)

Returns a plot of the Number Density of `species` as a function of time along a single coordinate direction. The direction is dictated by the values of `x_dir`,`y_dir`, and `z_dir`. The coordinate direction to plot is set by specifying one of these three to be the string `dir`. The other two may then either be set to `sum` to sum average along these directions or a single coordinate index defining the surface along which the quantity is measured. 

Arguments:
- `type::PlotType` determines the type of plot to generate, it can be either `Static` or `Animated`.
- `PhaseSpace::PhaseSpaceStruct` is the structure containing the phase space information of the simulation.
- `sol::OutputStruct` is the solution object containing the distribution function of all particles and time stepping information of the simulation. 
- `species::Vector{String}` is a vector of the abbreviated three letter names of the particle species to plot.
- `x_dir::Union{String,Int},y_dir::Union{String,Int},z_dir::Union{String,Int}` if `dir` then the plot will be along this coordinate direction, if `sum` then the plot will average along this direction, if `Int` then this coordinate will be kept constant at this coordinate index.

Common keyword arguments: 
- `theme=DiplodocusDark()`: the colour theme to use for the plot, can be either `DiplodocusDark()` for dark mode or `DiplodocusLight()` for light mode.
- `TimeUnits::Tuple{Float64,String}=(1.0, "\\text{Code Units}")`: a tuple that converts the time given in code units to the desired units for plotting. The first entry is the conversion factor and the second is a string that will be converted into a LaTeX string for the time label.
- `plot_limits=(nothing,nothing)`: the limits of the x and y axes, default sets the limits automatically based on the data.
- `wide=false`: if `true`, the plot is generated in a wide format (double column 8:3 aspect ratio), if `false` (single column 4:3 aspect ratio).
- `logt=false`: If `true` time will be converted to `log10` format for plotting.
"""
function NumberDensityPlot1D(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,species::Vector{String},x_dir::Union{String,Int},y_dir::Union{String,Int},z_dir::Union{String,Int};fig=nothing,theme=DiplodocusDark(),title=nothing,TimeUnits::Function=CodeToCodeUnitsTime,step::Int=1,wide=false,logt::Bool=false,plot_limits=(nothing,nothing))

    @assert (x_dir != "dir" && y_dir != "dir") || (y_dir != "dir" && z_dir != "dir") || (x_dir != "dir" && z_dir != "dir") "More than one coordinate direction cannot be the plot direction"

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
        fig = Figure(size=(576,216)) # double column 8:3 aspect ratio
    else
        fig = Figure() # default single column 4:3 aspect ratio
    end
    ylab = L"$n\, [\text{Code Units}]$"

    if x_dir == "dir"
        xlab = L"$x\, [\text{Code Units}]$"
        n_dir = x_num
        grid_type = Spacetime.x_grid
        if grid_type == "u"
            mdir_plot = mx
        elseif grid_type == "l"
            mdir_plot = log10.(mx)
        end
    elseif y_dir == "dir"
        xlab = L"$y\, [\text{Code Units}]$"
        n_dir = y_num
        grid_type = Spacetime.y_grid
        if grid_type == "u"
            mdir_plot = my
        elseif grid_type == "l"
            mdir_plot = log10.(my)
        end
    elseif z_dir == "dir"
        xlab = L"$z\, [\text{Code Units}]$"
        n_dir = z_num
        grid_type = Spacetime.z_grid
        if grid_type == "u"
            mdir_plot = mz
        elseif grid_type == "l"
            mdir_plot = log10.(mz)
        end
    end

    ax = Axis(fig[1,1],xlabel=xlab,ylabel=ylab)
    ax.limits = plot_limits

    linestyles = [:solid,(:dash,:dense),(:dot),(:dashdot),(:dashdotdot)]

    legend_elements = []
    line_labels = []

    t_save = length(sol.t)
    t_plot = ceil(Int64,t_save/step)

    values = (1:t_save)*step .+ 2 # add 2 to skip initial and kernel steps

    num = zeros(Float64,n_dir)
    num_total = zeros(Float64,n_dir)

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

        t_min = logt ? sol.t[2]/10 : sol.t[1]
        t_max = sol.t[end]

        for i in 1:t_save

            if (i in values || i == 1 || i == 2) # plot first step for initial conds, second for kernel 

                t = sol.t[i]
                #println("t=$(CodeToSIUnitsTime(t))")
                if logt
                    color = theme.colormap[][(log10(t) - log10(t_min)) / (log10(t_max) - log10(t_min))]
                else
                    color = theme.colormap[][(t - t_min) / (t_max - t_min)]
                end

                # sum over other 2 spatial directions to get average number density along desired direction
                fill!(num,0.0)

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    Nᵃ = DiplodocusTransport.FourFlow(sol.f[i],PhaseSpace,species_index;x_idx=x,y_idx=y,z_idx=z)
                    #Ua = HydroFourVelocity(Na)
                    Uₐ = [-1.0,0.0,0.0,0.0] # static observer
                    n = DiplodocusTransport.ScalarNumberDensity(Nᵃ,Uₐ) * CHAR_n

                    #println("n = $n at (x,y,z) = ($x,$y,$z)")

                    if x_dir == "dir" 
                        if y_dir isa Int && z_dir isa Int
                            if y != y_dir || z != z_dir 
                                continue 
                            else
                                num[x] += n
                            end
                        elseif y_dir isa Int && !(z_dir isa Int)
                            if y != y_dir 
                                continue 
                            else
                                num[x] += n * dz[z] / sum(dz) 
                            end
                        elseif !(y_dir isa Int) && (z_dir isa Int)
                            if z != z_dir 
                                continue 
                            else
                                num[x] += n * dy[y] / sum(dy) 
                            end
                        else
                            num[x] += n * dy[y] * dz[z] / (sum(dy)*sum(dz))
                        end
                    elseif y_dir == "dir" 
                        if x_dir isa Int && z_dir isa Int
                            if x != x_dir || z != z_dir 
                                continue 
                            else
                                num[y] += n 
                            end
                        elseif x_dir isa Int && !(z_dir isa Int)
                            if x != x_dir 
                                continue 
                            else
                                num[y] += n * dz[z] / sum(dz) 
                            end
                        elseif !(x_dir isa Int) && (z_dir isa Int)
                            if z != z_dir 
                                continue 
                            else
                                num[y] += n * dx[x] / sum(dx) 
                            end
                        else
                            num[y] += n * dx[x] * dz[z] / (sum(dx)*sum(dz))
                        end
                    elseif z_dir == "dir" 
                        if x_dir isa Int && y_dir isa Int
                            if x != x_dir || y != y_dir 
                                continue 
                            else
                                num[z] += n 
                            end
                        elseif x_dir isa Int && !(y_dir isa Int)
                            if x != x_dir 
                                continue 
                            else
                                num[z] += n * dy[y] / sum(dy) 
                            end
                        elseif !(x_dir isa Int) && (y_dir isa Int)
                            if y != y_dir 
                                continue 
                            else
                                num[z] += n * dx[x] / sum(dx) 
                            end
                        else
                            num[z] += n * dx[x] * dy[y] / (sum(dx)*sum(yz))
                        end
                    end

                end

                replace!(num,0.0=>NaN)

                scatterlines!(ax,mdir_plot,num,color = color,markersize=5.0)

            end # if times should be plotted
        
        end # time loop
       
    end # species loop 

    t_unit_string = TimeUnits[2]
    t_unit_scale = TimeUnits[1]

    if logt
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(log10(round(t_units_scale*sol.t[1],sigdigits=5)),log10(round(t_units_scale*sol.t[end],sigdigits=5))),label=L"$\log_{10}\left(t\,[%$t_unit_string] \right)$")
    else
        Colorbar(fig[1,2],colormap = theme.colormap,limits=(sol.t[1] * t_unit_scale,sol.t[end] * t_unit_scale),label=L"$t\, [%$t_unit_string]$")
    end
    
    #axislegend(ax,position = :rb)

    end # with_theme

    return fig

end


