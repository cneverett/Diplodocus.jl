function find_closest(A::AbstractArray{T}, b::T) where {T<:Real}
    if length(A) <= 1
        return firstindex(A)
    end

    i = searchsortedfirst(A, b)

    if i == firstindex(A)
        return i
    elseif i > lastindex(A)
        return lastindex(A)
    else
        prev_dist = b - A[i-1]
        next_dist = A[i] - b

        if prev_dist < next_dist
            return i - 1
        else
            return i
        end
    end
end

function SimToPlotGrid(PhaseSpace::PhaseSpaceStruct,surface::Tuple{String,Int64},plot_coordinates::AbstractCoordinates,sim_coordinates::AbstractCoordinates;x_units::String="\\text{Code Units}",y_units::String="\\text{Code Units}")

    surface_string = surface[1]
    surface_index = surface[2]
    Spacetime = PhaseSpace.Spacetime
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    if plot_coordinates isa Cartesian

        x_string = "x"
        y_string = "y"
        z_string = "z"
        
    elseif plot_coordinates isa Cylindrical

        x_string = "ρ"
        y_string = "φ"
        z_string = "z"

        if sim_coordinates isa Paraboloidal
            if surface_string == "xz" # ρ-z plane maps to u-v plane
                nx_dir = y_num
                ny_dir = z_num
            end
        end

     elseif plot_coordinates isa Paraboloidal

        x_string = "ϕ"
        y_string = "u"
        z_string = "v"

        if sim_coordinates isa Paraboloidal
            if surface_string == "xy" # ϕ-u plane maps to x-y plane
                nx_dir = x_num
                ny_dir = y_num
            elseif surface_string == "xz" # ϕ-v plane maps to x-z plane
                nx_dir = x_num
                ny_dir = z_num
            elseif surface_string == "yz" # u-v plane maps to y-z plane
                nx_dir = y_num
                ny_dir = z_num
            end
        end

    else

        error("Invalid plot coordinates specified. Must be one of Cartesian, Cylindrical, or Paraboloidal.")

    end
   

    if surface_string == "xy"
        xlab = L"$%$(x_string)\,[%$(x_units)]$"
        ylab = L"$%$(y_string)\,[%$(y_units)]$"
    elseif surface_string == "yz"
        xlab = L"$%$(y_string)\,[%$(x_units)]$"
        ylab = L"$%$(z_string)\,[%$(y_units)]$"
    elseif surface_string == "xz"
        xlab = L"$%$(x_string)\,[%$(x_units)]$"
        ylab = L"$%$(z_string)\,[%$(y_units)]$"
    else
        error("Invalid surface specified. Must be one of 'xy', 'yz', or 'xz'.")
    end

    return xlab, ylab, nx_dir, ny_dir

end

function PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,::Paraboloidal,::Cylindrical)

    surface_string = surface[1]
    surface_index = surface[2]

    if surface_string == "xz" # ρ-z plane maps to u-v plane
        x = surface_index
        y = x_plot_idx
        z = y_plot_idx
    else
        error("Invalid surface specified for plotting a Paraboloidal coordinates simulation in Cylindrical coordinates. Only 'xz' is supported.")
    end

    return x,y,z

end

function PlotIndexTransform2D!(x_plot_idx,y_plot_idx,surface,::Paraboloidal,::Paraboloidal)

    surface_string = surface[1]
    surface_index = surface[2]

    if surface_string == "xz" # ϕ-v plane maps to ϕ-v plane
        y = surface_index
        x = x_plot_idx
        z = y_plot_idx
    elseif surface_string == "yz" # u-v plane maps to u-v plane
        x = surface_index
        y = x_plot_idx
        z = y_plot_idx
    elseif surface_string == "xy" # ϕ-u plane maps to ϕ-u plane
        z = surface_index
        x = x_plot_idx
        y = y_plot_idx
    else
        error("Invalid surface specified for plotting a Paraboloidal coordinates simulation in Cylindrical coordinates. Only 'xz', 'yz', and 'xy' are supported.")
    end

    return x,y,z

end