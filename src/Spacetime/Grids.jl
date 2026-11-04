abstract type AbstractSpacetimeGrid end

"""
    UniformGrid(up,low,num)

Creates a uniformly spaced grid from `low` to `up` with `num` divisions.

    | a | a | a | a | a | 
    low                  up
    a = (up - low)/num

"""
struct UniformGrid <: AbstractSpacetimeGrid
    up::Float64
    low::Float64
    num::Int64
end

""" 
    Log10Grid(up,low,num)

Creates a log10 spaced grid from `10^low` to `10^up` with `num` divisions.
    
    | 1 | 10 | 100 | 1000 | 10000 | 
    10^low                        10^up

"""
struct Log10Grid <: AbstractSpacetimeGrid
    up::Float64
    low::Float64
    num::Int64
end

"""
    StretchGrid(low,up,sub,num)

Creates a stretched grid that starts at `low` and passes through `up` with `sub` grid points between `low` and `up`, the total number of grid cells is `num`. If `num` = `sub+2` then the grid ends at `up`, otherwise the grid continues beyond `up`. The spacing between grid points is determined by a geometric progressions with `ratio = 2^(1/(sub+1))`
"""
struct StretchGrid <: AbstractSpacetimeGrid
    low::Float64
    up::Float64 
    sub::Int64
    num::Int64
end

"""
    SpacetimeGrid!(grid::AbstractSpacetimeGrid)

Returns a `num+1` long `Vector{Float}` of grid bounds for a given spacetime coordinate based on the type of `grid` specified.
"""
SpacetimeGrid!(grid::AbstractSpacetimeGrid) = error("Spacetime grid function not defined for this grid $(typeof(grid)).")


function SpacetimeGrid(grid::UniformGrid) 
    #= uniform spacing
        | a | a | a | a | a | 
        low                  up
        a = (up - low)/num
    =#
    return [range(grid.low,grid.up,grid.num+1);]
end

function SpacetimeGrid(grid::Log10Grid) 
    #= log10 spacing
        | 1 | 10 | 100 | 1000 | 10000 | 
        10^low                        10^up
    =#
    return 10 .^[range(grid.low,grid.up,grid.num+1);]
end

function SpacetimeGrid(grid::StretchGrid) 
    #= stretch spacing 
        grid starts at `low` and passes through `up` with `sub` grid points between `low` and `up`, the total number of grid cells is `num`. If `num` = `sub+2` then the grid ends at `up`, otherwise the grid continues beyond `up`. The spacing between grid points is determined by a geometric progressions with `ratio = 2^(1/(sub+1))`
    =#
    r = 2^(1/(grid.sub+1))
    vec = [grid.low + (grid.up - grid.low) * (r^i-1) / (r^(grid.sub+1)-1) for i in 0:(grid.num)]
    return vec
end


grid = StretchGrid(0.0,1.0,1,8)
SpacetimeGrid(grid)