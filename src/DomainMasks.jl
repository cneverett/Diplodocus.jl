"""
    InclusiveDomainMask(PhaseSpace;x_dom,y_dom,z_dom)

Returns a vector of offset indices (`off_space`) of spatial grid points that are included within (inclusive) the specified spatial domain defined by x_dom, y_dom, and z_dom. x_dom,y_dom,z_dom are tuples of the form (lower_bound, upper_bound) defining the spatial domain in each direction.
"""
function InclusiveDomainMask(PhaseSpace::PhaseSpaceStruct;x_dom=(-Inf, Inf),y_dom=(-Inf, Inf),z_dom=(-Inf, Inf))

    x_dom_up = x_dom[2]
    x_dom_low = x_dom[1]
    y_dom_up = y_dom[2]
    y_dom_low = y_dom[1]
    z_dom_up = z_dom[2]
    z_dom_low = z_dom[1]
    
    Spacetime = PhaseSpace.Spacetime
    Grids = PhaseSpace.Grids

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    mx = Grids.mx
    my = Grids.my
    mz = Grids.mz

    offset_space_indices::Vector{Int64} = Int64[]

    for x in eachindex(mx), y in eachindex(my), z in eachindex(mz)
        if (x_dom_low <= mx[x] <= x_dom_up) && (y_dom_low <= my[y] <= y_dom_up) && (z_dom_low <= mz[z] <= z_dom_up)

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
            push!(offset_space_indices,off_space)

        end
    end

    return offset_space_indices

end

"""
    ExclusiveDomainMask(PhaseSpace;x_dom,y_dom,z_dom)

Returns a vector of offset indices (`off_space`) of spatial grid points that are outside of (exclusive) the specified spatial domain defined by x_dom, y_dom, and z_dom. x_dom,y_dom,z_dom are tuples of the form (lower_bound, upper_bound) defining the spatial domain in each direction.
"""
function ExclusiveDomainMask(PhaseSpace::PhaseSpaceStruct;x_dom=(-Inf, Inf),y_dom=(-Inf, Inf),z_dom=(-Inf, Inf))

    x_dom_up = x_dom[2]
    x_dom_low = x_dom[1]
    y_dom_up = y_dom[2]
    y_dom_low = y_dom[1]
    z_dom_up = z_dom[2]
    z_dom_low = z_dom[1]
    
    Spacetime = PhaseSpace.Spacetime
    Grids = PhaseSpace.Grids

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    mx = Grids.mx
    my = Grids.my
    mz = Grids.mz

    offset_space_indices::Vector{Int64} = Int64[]

    for x in eachindex(mx), y in eachindex(my), z in eachindex(mz)
        if !(x_dom_low <= mx[x] <= x_dom_up) && !(y_dom_low <= my[y] <= y_dom_up) && !(z_dom_low <= mz[z] <= z_dom_up)

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
            push!(offset_space_indices,off_space)

        end
    end

    return offset_space_indices

end