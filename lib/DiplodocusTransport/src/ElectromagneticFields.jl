"""
    ElectromagneticFieldGrid(Spacetime,tetrad,Grids)

Returns the average values of the magnetic and electric fields within a given spatial grid cell (x_index,y_index,z_index) as measured by an observer associated with the `tetrad`, which also defines the type of field. Returning two arrays for the magnetic field B(x_index,y_index,z_index) and electric field E(x_index,y_index,z_index) respectively. Note: this is currently only defined for "force-free" tetrads the **magnetic is measured with respect to the orthonormal observer frame of the tetrad** and is taken to be along the local momentum z-direction and the electric field is zero. Tetrads that allow an E parallel to B will be implemented in the future.
"""
ElectromagneticFieldGrid(Spacetime::SpacetimeStruct,tetrad::ForceFreeTetrad,Grids::GridsStruct) = error("Electromagnetic field not defined for this tetrad $(typeof(tetrad)).")


function ElectromagneticFieldGrid(Spacetime::SpacetimeStruct,tetrad::UniformElectromagneticFieldTetrad,Grids::GridsStruct) 

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    B0 = tetrad.B0 # measured with respect to stationary observer NOT tetrad frame
    E0 = tetrad.E0 # measured with respect to stationary observer NOT tetrad frame

    B = sqrt(B0^2-E0^2)

    B_field = zeros(Float64,x_num,y_num,z_num)
    E_field = zeros(Float64,x_num,y_num,z_num)

    #=
     B field is ALWAYS in local momentum z-direction and E field in local momentum y-direction using local orthonormal basis (T^α,ϵ^αβγδn_βE_γB_δ,E^α,B^α), in this frame E is zero.

     B field is taken to be the average B field in a coordinate grid cell i.e. ∫B(x,y,z)χ dxdydz / ∫χ dxdydz where χ is the volume element of the grid cell. For a uniform B field and E field this is just B and E.
    =#

    for ix in 1:x_num, iy in 1:y_num, iz in 1:z_num
        B_field[ix,iy,iz] = B
        E_field[ix,iy,iz] = 0.0
    end

    return B_field, E_field

end

function ElectromagneticFieldGrid(Spacetime::SpacetimeStruct,tetrad::ParabolicForceFreeFieldTetrad,Grids::GridsStruct)
    
    metric = Spacetime.metric
    coordinates = Spacetime.coordinates

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    B0 = tetrad.B0
    @inline Bfunction(pos) = tetrad.Bfunction(pos,coordinates) *  VolumeElement(pos,metric,coordinates)
    @inline Vfunction(pos) = VolumeElement(pos,metric,coordinates)

    B_field = zeros(Float64,x_num,y_num,z_num)
    E_field = zeros(Float64,x_num,y_num,z_num)

    #=
     B field is ALWAYS in local momentum z-direction and E field in local momentum y-direction using local orthonormal basis (T^α,ϵ^αβγδn_βE_γB_δ,E^α,B^α), in this frame E is zero.

     B field is taken to be the average B field in a coordinate grid cell i.e. ∫B(x,y,z)χ dtdxdydz / ∫χ dtdxdydz where χ is the volume element of the grid cell. For a uniform B field and E field this is just B and E.
    =#

    t0 = Grids.tr[1]
    t1 = Grids.tr[2]
    n::SVector{4,Int64} = [2,2,16,16] 

    for ix in 1:x_num, iy in 1:y_num, iz in 1:z_num
        x0 = Grids.xr[ix]
        x1 = Grids.xr[ix+1]
        y0 = Grids.yr[iy]
        y1 = Grids.yr[iy+1]
        z0 = Grids.zr[iz]
        z1 = Grids.zr[iz+1]
        a::SVector{4,Float64} = [t0,x0,y0,z0]
        b::SVector{4,Float64} = [t1,x1,y1,z1]
        B_field[ix,iy,iz] = B0 * Simpson4D(Bfunction,a,b,n) / Simpson4D(Vfunction,a,b,n)
        E_field[ix,iy,iz] = 0.0
    end

    return B_field, E_field

end

#="""
    ElectroMagneticFieldFunction_InvZDecay(Space::SpaceStruct,Momentum::MomentumStruct,Characteristic::CharacteristicStruct,Grids::GridsStruct,parameters::Vector{Float64})

Returns the values of the magnetic and electric fields at a given spatial grid point (x_index,y_index,z_index) for a magnetic field that has the form B=B0*L0/z (Note this is only for testing as grad(B)!=0), with the parameters `B0` and `L0` defined by the `ElectroMagneticField_InvZDecay` struct. Note: L0 and B0 are defined in terms of the characteristic length and magnetic field scales respectively.
"""
function ElectroMagneticFieldFunction_InvZDecay(Space::SpaceStruct,Momentum::MomentumStruct,Characteristic::CharacteristicStruct,Grids::GridsStruct,parameters::Vector{Float64}) 

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    Characteristic = Characteristic
    Grids = Grids

    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num

    L = Characteristic.CHAR_length
    B0 = parameters[1]
    L0 = parameters[2]

    mz = Grids.mz

    B_field = zeros(Float64,x_num,y_num,z_num)
    E_field = zeros(Float64,x_num,y_num,z_num)

    # B field is ALWAYS in local momentum z-direction and E field in local momentum y-direction using local orthonormal basis (n^α,ϵ^αβγδn_βE_γB_δ,E^α,B^α)

    if space_coords isa Cartesian # B along z
        for ix in 1:x_num, iy in 1:y_num, iz in 1:z_num
            B_field[ix,iy,iz] = B0*L0/mz[iz]
            E_field[ix,iy,iz] = 0e0
        end
    elseif space_coords isa Cylindrical # B along z
        for ix in 1:x_num, iy in 1:y_num, iz in 1:z_num
            B_field[ix,iy,iz] = B0*L0/mz[iz]
            E_field[ix,iy,iz] = 0e0
        end
    end

    return B_field, E_field

end=#