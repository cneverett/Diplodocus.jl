#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT then:
        - for A fluxes the dimensions of `flux` are L^3
        - for B, C, D fluxes the dimensions of `flux` are L^2*cT
        - for I, J, and K fluxes the dimensions of `flux` are cT*L^2
        - for volume elements the dimensions are L^3*cT

=#

# =============== ABCD Vol  ===================== #

    function AFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        if plus_minus == "plus"
            t = Grids.tr[t_idx+1]
        elseif plus_minus == "minus"
            t = Grids.tr[t_idx]
        else
            error("plus_minus string not recognised.")
        end
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the A surface i.e. surface of constant t

        flux::Float64 = 1.0

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= (-x0 + x1) * (-y0 + y1) * (-z0 + z1) 
            flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            flux *= (-x0^2 + x1^2) / 2 * (-y0 + y1) * (-z0 + z1) 
            flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            flux *= -(1/3) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))
            flux *= (p0 - p1) * (u0 - u1) * (h0 - h1)

        else
            error("A Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function BFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        if plus_minus == "plus"
            x = Grids.xr[x_idx+1]
        elseif plus_minus == "minus"
            x = Grids.xr[x_idx]
        else
            error("plus_minus string not recognised.")
        end
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the B surface i.e. surface of constant x coordinate: cartesian x, spherical r, cylindrical r

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi
        
            flux *=  (t0 - t1) * (y0 - y1) * (z0 - z1)
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(h0) - sinpi(h1)) / 2

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * x * (y0 - y1) * (z0 - z1)
            flux *= (-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * (sinpi(γ - h0) - sinpi(γ - h1))

        elseif space_coords isa Spherical && momentum_coords isa Spherical 

            h1 = h1/pi
            h0 = h0/pi
            y1 = y1/pi
            y0 = y0/pi

            if space_coords.xp_BC isa Escape
                # TODO: This was a quick fix, check later
                flux *= (1/3) * (x^3) * (z0 - z1) * (cospi(y0) - cospi(y1))
                flux *= (p0 - p1) * (u0 - u1) * (h0 - h1)
            else
                flux *= (t0 - t1) * x^2 * (z0 - z1) * (cospi(y0) - cospi(y1)) / 2
                flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * pi*(h0 - h1)
            end

        else
            error("B Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function CFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        if plus_minus == "plus"
            y = Grids.yr[y_idx+1]
        elseif plus_minus == "minus"
            y = Grids.yr[y_idx]
        else
            error("plus_minus string not recognised.")
        end
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the C surface i.e. surface of constant y coordinate: cartesian y, spherical θ, cylindrical ϕ

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *=  (t0 - t1) * (x0 - x1) * (z0 - z1) 
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (cospi(h0) - cospi(h1)) / 2 

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (-t0 + t1) * (-x0 + x1) * (-z0 + z1) 
            flux *= ((-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * cospi(β) * (cospi(γ - h0) - cospi(γ - h1)) + u0^2 * h0 * sinpi(β) - u1^2 * h0 * sinpi(β) - u0^2 * h1 * sinpi(β) + u1^2 * h1 * sinpi(β))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * sin(y) / 2
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(h0) - sinpi(h1)) / 2

        else 
            error("C Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function DFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        if plus_minus == "plus"
            z = Grids.zr[z_idx+1]
        elseif plus_minus == "minus"
            z = Grids.zr[z_idx]
        else
            error("plus_minus string not recognised.")
        end
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the D surface i.e. surface of constant z coordinate: cartesian z, spherical ϕ, cylindrical z

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) 
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1) / 2

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            flux *= (1/4) * (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) 
            flux *= ((-u0^2 + u1^2) * (h0 - h1) * cospi(β) + (-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * (cospi(γ - h0) - cospi(γ - h1)) * sinpi(β))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) / 2
            flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0) - cospi(h1)) / 2

        else
            error("D Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function VolFunction(PhaseSpace::PhaseSpaceStruct,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64)

        Grids = PhaseSpace.Grids
        metric = PhaseSpace.Spacetime.metric
        coordinates = PhaseSpace.Spacetime.coordinates

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        # Evaluate the spacetime volume element

        vol::Float64 = 1.0

        if metric isa Minkowski
            if coordinates isa Cartesian

                vol *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)

            elseif coordinates isa Cylindrical

                vol *= (-t0 + t1) * (-x0^2 + x1^2) * (-y0 + y1) * (-z0 + z1) / 2

            elseif coordinates isa Spherical

                y1 = y1/pi
                y0 = y0/pi

                vol *= -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0) - cospi(y1))

            elseif coordinates isa Paraboloidal

                vol *= 1/8 * (t0 - t1) * (x0 - x1) * (y0^2 - y1^2) * (z0^2 - z1^2) * (y0^2 + y1^2 + z0^2 + z1^2)

            else
                error("Volume function not implemented metric $(typeof(metric)) with coordinates $(typeof(coordinates)).")
            end
        else
            error("Volume function not implemented metric $(typeof(metric)).")
        end

        return vol

    end

# ======================================================== #

"""
    IFluxFunction(force,PhaseSpace,species_idx,plus_minus,t_idx,x_idx,y_idx,z_idx,px_idx,py_idx,pz_idx)

Returns the flux through the I surface (surface of constant p coordinate) for a given analytic force. plus_minus is a string that can be either "plus" or "minus" to indicate which p coordinate surface the flux is being evaluated at for a given px_idx. The indices t_idx, x_idx, y_idx, z_idx, px_idx, py_idx, and pz_idx specify the location in the phase space grid where the flux is being evaluated.   
"""
IFluxFunction(force::AnalyticForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64) = error("I Flux Function not defined for this analytic force type $(typeof(force))).")

"""
    JFluxFunction(force,PhaseSpace,species_idx,plus_minus,t_idx,x_idx,y_idx,z_idx,px_idx,py_idx,pz_idx)

Returns the flux through the J surface (surface of constant u coordinate) for a given analytic force. plus_minus is a string that can be either "plus" or "minus" to indicate which u coordinate surface the flux is being evaluate at for a given py_idx. The indices t_idx, x_idx, y_idx, z_idx, px_idx, py_idx, and pz_idx specify the location in the phase space grid where the flux is being evaluated.   
"""
JFluxFunction(force::AnalyticForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64) = error("J Flux Function not defined for this analytic force type $(typeof(force))).")

"""
    KFluxFunction(force,PhaseSpace,species_idx,plus_minus,t_idx,x_idx,y_idx,z_idx,px_idx,py_idx,pz_idx)

Returns the flux through the K surface (surface of constant phi coordinate) for a given analytic force. plus_minus is a string that can be either "plus" or "minus" to indicate which phi coordinate surface the flux is being evaluate at for a given pz_idx. The indices t_idx, x_idx, y_idx, z_idx, px_idx, py_idx, and pz_idx specify the location in the phase space grid where the flux is being evaluated.   
"""
KFluxFunction(force::AnalyticForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64) = error("K Flux Function not defined for this analytic force type $(typeof(force))).")

# =============== IJK Cartesian Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        if plus_minus == "plus"
            p = Grids.pxr_list[species_idx][px_idx+1]
        elseif plus_minus == "minus"
            p = Grids.pxr_list[species_idx][px_idx]
        else
            error("plus_minus string not recognised.")
        end
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            flux *= 0e0

        else
            error("IFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function JFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        if plus_minus == "plus"
            u = Grids.pyr_list[species_idx][py_idx+1]
        elseif plus_minus == "minus"
            u = Grids.pyr_list[species_idx][py_idx]
        else
            error("plus_minus string not recognised.")
        end
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1) 
            flux *= cospi(γ - h0/2 - h1/2) * sinpi(β) * sinpi((h0 - h1)/2) * (-4u*sqrt(1 - u^2)*sinpi(β) -4(-1 + u^2)*cospi(β)*cospi((h0 - h1)/2)*sinpi(γ - h0/2 - h1/2))


        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi
            y1 = y1/pi
            y0 = y0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * (cospi(y0) - cospi(y1)) / 2
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1+u^2) * pi*(h0 - h1)

        else
            error("JFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function KFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        if plus_minus == "plus"
            h = Grids.pzr_list[species_idx][pz_idx+1]
        elseif plus_minus == "minus"
            h = Grids.pzr_list[species_idx][pz_idx]
        else
            error("plus_minus string not recognised.")
        end

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("KFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h = h/pi

            flux *= -(1/16) * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            flux *= (6u0^2*sinpi(2*β) - 6u1^2*sinpi(2*β) - 8*asin(u0)*sinpi(γ - h) + 8*asin(u1)*sinpi(γ - h) - u0^2*sinpi(2*(β + γ - h)) + u1^2*sinpi(2*(β + γ - h)) - 4*u0*sqrt(1 - u0^2)*sinpi(2*β + γ - h) + 4*u1*sqrt(1 - u1^2)*sinpi(2*β + γ - h) - u0^2*sinpi(2*(β - γ + h)) + u1^2*sinpi(2*(β - γ + h)) + 4*u0*sqrt(1 - u0^2)*sinpi(2*β - γ + h) - 4*u1*sqrt(1 - u1^2)*sinpi(2*β - γ + h))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h = h/pi
            y1 = y1/pi
            y0 = y0/pi

            flux *=  (t0 - t1) * (x0^2 - x1^2) * (sinpi(y0) - sinpi(y1)) * (z0 - z1) / 2
            flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * ((u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1)) * sinpi(h)) / 2

        else
            error("KFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

# ======================================================== #

function acot_mod(x)

    # ArcCot[(1+x)/Sqrt[1-x^2]]
    # this is actually just acos(x)/2

    #=if x == -1
        return pi/2
    elseif x == 1
        return 0.0
    else
        return acot((1 + x)/sqrt(1 - x^2))
    end=#

    return acos(x)/2

end

# ============ CoordinateForce and CoordinateFlux Momentum Integrals =========== #

function CoordinateFluxMomentumIntegral!(FluxMomentumArray::MVector{4,Float64},PhaseSpace::PhaseSpaceStruct,species_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

    Grids = PhaseSpace.Grids
    m = Grids.mass_list[species_idx]
    p0 = Grids.pxr_list[species_idx][px_idx]
    p1 = Grids.pxr_list[species_idx][px_idx+1]
    u0 = Grids.pyr_list[species_idx][py_idx]
    u1 = Grids.pyr_list[species_idx][py_idx+1]
    h0 = Grids.pzr_list[species_idx][pz_idx]
    h1 = Grids.pzr_list[species_idx][pz_idx+1]

    h0divpi = h0/pi
    h1divpi = h1/pi
    sh0,ch0 = sincospi(h0divpi)
    sh1,ch1 = sincospi(h1divpi)

    FluxMomentumArray[1] = (-h0 + h1) * (-p0 + p1) * (-u0 + u1)
    FluxMomentumArray[2] = 0.5*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - acos(u0) + acos(u1)) * (sh0 - sh1)
    FluxMomentumArray[3] = 0.5*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - acos(u0) + acos(u1)) * (ch0 - ch1)
    FluxMomentumArray[4] = 0.5*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1)

end

function CoordinateForceMomentumIntegral!(ForceMomentumArray::MArray{Tuple{4,4,4,3,2},Float64,5,384},PhaseSpace::PhaseSpaceStruct,species_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

    Grids = PhaseSpace.Grids
    m = Grids.mass_list[species_idx]
    p0 = Grids.pxr_list[species_idx][px_idx]
    p1 = Grids.pxr_list[species_idx][px_idx+1]
    u0 = Grids.pyr_list[species_idx][py_idx]
    u1 = Grids.pyr_list[species_idx][py_idx+1]
    h0 = Grids.pzr_list[species_idx][pz_idx]
    h1 = Grids.pzr_list[species_idx][pz_idx+1]

    fill!(ForceMomentumArray,zero(Float64))

    CoordinateForceIFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,1,1),m,p1,u0,u1,h0,h1)
    CoordinateForceIFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,1,2),m,p0,u0,u1,h0,h1)

    CoordinateForceJFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,2,1),m,p0,p1,u1,h0,h1)
    CoordinateForceJFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,2,2),m,p0,p1,u0,h0,h1)

    CoordinateForceKFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,3,1),m,p0,p1,u0,u1,h1)
    CoordinateForceKFluxMomentumIntegral!(view(ForceMomentumArray,:,:,:,3,2),m,p0,p1,u0,u1,h0)

end

function CoordinateForceIFluxMomentumIntegral!(ForceMomentumArrayView::AbstractArray{T,3},m::T,p::T,u0::T,u1::T,h0::T,h1::T) where T

    h0divpi = h0/pi
    h1divpi = h1/pi
    sh0,ch0 = sincospi(h0divpi)
    sh1,ch1 = sincospi(h1divpi)

    M121::T = -(1/4) * sqrt(m^2 + p^2) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - acos(u0) + acos(u1)) * (sh0 - sh1)
    M131::T = 1/4 * sqrt(m^2 + p^2) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - acos(u0) + acos(u1)) * (ch0 - ch1)
    M141::T = -(1/4) * (h0 - h1) * sqrt(m^2 + p^2) * (u0^2 - u1^2)

    M122::T = 1/12 * p * (-3u0 + u0^3 + 3u1 - u1^3) * (h0 - h1 + ch0 * sh0 - ch1 * sh1)
    M132::T = -(1/12) * p * (-3u0 + u0^3 + 3u1 - u1^3) * (ch0^2 - ch1^2)
    M142::T = 1/6 * p * ((1 - u0^2)^(3/2) - (1 - u1^2)^(3/2)) * (sh0 - sh1)

    M123::T = M132
    M133::T = 1/12 * p * (-3u0 + u0^3 + 3u1 - u1^3) * (h0 - h1 - ch0 * sh0 + ch1 * sh1)
    M143::T = 1/6 * p * (-(1 - u0^2)^(3/2) + (1 - u1^2)^(3/2)) * (ch0 - ch1)

    M124::T = M142
    M134::T = M143
    M144::T = 1/6 * (h0 - h1) * p * (-u0^3 + u1^3)

    ForceMomentumArrayView[1,2,1] = M121
    ForceMomentumArrayView[1,3,1] = M131
    ForceMomentumArrayView[1,4,1] = M141
    ForceMomentumArrayView[1,2,2] = M122
    ForceMomentumArrayView[1,3,2] = M132
    ForceMomentumArrayView[1,4,2] = M142
    ForceMomentumArrayView[1,2,3] = M123
    ForceMomentumArrayView[1,3,3] = M133
    ForceMomentumArrayView[1,4,3] = M143
    ForceMomentumArrayView[1,2,4] = M124
    ForceMomentumArrayView[1,3,4] = M134
    ForceMomentumArrayView[1,4,4] = M144

    ForceMomentumArrayView[2,1,1] = -M121
    ForceMomentumArrayView[3,1,1] = -M131
    ForceMomentumArrayView[4,1,1] = -M141
    ForceMomentumArrayView[2,1,2] = -M122
    ForceMomentumArrayView[3,1,2] = -M132
    ForceMomentumArrayView[4,1,2] = -M142
    ForceMomentumArrayView[2,1,3] = -M123
    ForceMomentumArrayView[3,1,3] = -M133
    ForceMomentumArrayView[4,1,3] = -M143
    ForceMomentumArrayView[2,1,4] = -M124
    ForceMomentumArrayView[3,1,4] = -M134
    ForceMomentumArrayView[4,1,4] = -M144
    
end

function CoordinateForceJFluxMomentumIntegral!(ForceMomentumArrayView::AbstractArray{T,3},m::T,p0::T,p1::T,u::T,h0::T,h1::T) where T

    h0divpi = h0/pi
    h1divpi = h1/pi
    sh0,ch0 = sincospi(h0divpi)
    sh1,ch1 = sincospi(h1divpi)
    sqrt_mp0 = sqrt(m^2 + p0^2)
    sqrt_mp1 = sqrt(m^2 + p1^2)

    M121::T = 0.5 * u * sqrt(1 - u^2) * (sqrt_mp0 - sqrt_mp1 - m*asinh(m/p0) + m*asinh(m/p1)) * (sh0 - sh1)
    M131::T = -0.5 * u * sqrt(1 - u^2) * (sqrt_mp0 - sqrt_mp1 - m*asinh(m/p0) + m*asinh(m/p1)) * (ch0 - ch1)
    M141::T = -0.5 * (h0 - h1) * (1 - u^2) * (sqrt_mp0 - sqrt_mp1 - m*asinh(m/p0) + m*asinh(m/p1))
    M241::T = -0.5 * (p0 - p1) * sqrt(1 - u^2) * (sh0 - sh1)
    M341::T = 0.5 * (p0 - p1) * sqrt(1 - u^2) * (ch0 - ch1)

    M122::T = 0.25 * (p0 - p1) * u * (1 - u^2) * (h0 - h1 + ch0 * sh0 - ch1 * sh1)
    M132::T = -0.25 * (p0 - p1) * u * (1 - u^2) * (ch0^2 - ch1^2)
    M142::T = -0.5 * (p0 - p1) * (1 - u^2)^(3/2) * (sh0 - sh1)
    M242::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (1 - u^2) * (h0 - h1 + ch0 * sh0 - ch1 * sh1)
    M342::T = 0.25 * (sqrt_mp0 - sqrt_mp1) * (1 - u^2) * (ch0^2 - ch1^2)

    M123::T = -0.25 * (p0 - p1) * u * (1 - u^2) * (ch0^2 - ch1^2)
    M133::T = 0.25 * (p0 - p1) * u * (1 - u^2) * (h0 - h1 - ch0 * sh0 + ch1 * sh1)
    M143::T = 0.5 * (p0 - p1) * (1 - u^2)^(3/2) * (ch0 - ch1)
    M243::T = 0.25 * (sqrt_mp0 - sqrt_mp1) * (1 - u^2) * (ch0^2 - ch1^2)
    M343::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (1 - u^2) * (h0 - h1 - ch0 * sh0 + ch1 * sh1)

    M124::T = 0.5 * (p0 - p1) * u^2 * sqrt(1 - u^2) * (sh0 - sh1)
    M134::T = -0.5 * (p0 - p1) * u^2 * sqrt(1 - u^2) * (ch0 - ch1)
    M144::T = -0.5 * (h0 - h1) * (p0 - p1) * u * (1 - u^2)
    M244::T = -0.5 * (sqrt_mp0 - sqrt_mp1) * u * sqrt(1 - u^2) * (sh0 - sh1)
    M344::T = 0.5 * (sqrt_mp0 - sqrt_mp1) * u * sqrt(1 - u^2) * (ch0 - ch1)

    ForceMomentumArrayView[1,2,1] = M121
    ForceMomentumArrayView[1,3,1] = M131
    ForceMomentumArrayView[1,4,1] = M141
    ForceMomentumArrayView[1,2,2] = M122
    ForceMomentumArrayView[1,3,2] = M132
    ForceMomentumArrayView[1,4,2] = M142
    ForceMomentumArrayView[1,2,3] = M123
    ForceMomentumArrayView[1,3,3] = M133
    ForceMomentumArrayView[1,4,3] = M143
    ForceMomentumArrayView[1,2,4] = M124
    ForceMomentumArrayView[1,3,4] = M134
    ForceMomentumArrayView[1,4,4] = M144

    ForceMomentumArrayView[2,1,1] = -M121
    ForceMomentumArrayView[3,1,1] = -M131
    ForceMomentumArrayView[4,1,1] = -M141
    ForceMomentumArrayView[2,1,2] = -M122
    ForceMomentumArrayView[3,1,2] = -M132
    ForceMomentumArrayView[4,1,2] = -M142
    ForceMomentumArrayView[2,1,3] = -M123
    ForceMomentumArrayView[3,1,3] = -M133
    ForceMomentumArrayView[4,1,3] = -M143
    ForceMomentumArrayView[2,1,4] = -M124
    ForceMomentumArrayView[3,1,4] = -M134
    ForceMomentumArrayView[4,1,4] = -M144

    ForceMomentumArrayView[2,4,1] = M241
    ForceMomentumArrayView[3,4,1] = M341
    ForceMomentumArrayView[2,4,2] = M242
    ForceMomentumArrayView[3,4,2] = M342
    ForceMomentumArrayView[2,4,3] = M243
    ForceMomentumArrayView[3,4,3] = M343
    ForceMomentumArrayView[2,4,4] = M244
    ForceMomentumArrayView[3,4,4] = M344

    ForceMomentumArrayView[4,2,1] = -M241
    ForceMomentumArrayView[4,3,1] = -M341
    ForceMomentumArrayView[4,2,2] = -M242
    ForceMomentumArrayView[4,3,2] = -M342
    ForceMomentumArrayView[4,2,3] = -M243
    ForceMomentumArrayView[4,3,3] = -M343
    ForceMomentumArrayView[4,2,4] = -M244
    ForceMomentumArrayView[4,3,4] = -M344
end

function CoordinateForceKFluxMomentumIntegral!(ForceMomentumArrayView::AbstractArray{T,3},m::T,p0::T,p1::T,u0::T,u1::T,h::T) where T

    sqrt_mp0 = sqrt(m^2 + p0^2)
    sqrt_mp1 = sqrt(m^2 + p1^2)
    hdivpi = h / pi
    sh,ch = sincospi(hdivpi) # high accuracy than sincos 

    M121::T = -0.5 * (acos(u0) - acos(u1)) * (sqrt_mp0 - sqrt_mp1 - m * asinh(m/p0) + m * asinh(m/p1)) * sh
    M131::T = 0.5 * (acos(u0) - acos(u1)) * (sqrt_mp0 - sqrt_mp1 - m * asinh(m/p0) + m * asinh(m/p1)) * ch

    M231::T = -0.5 * (p0 - p1) * (u0 - u1)
    M241::T = 0.5 * (p0 - p1) * (sqrt(1 - u0^2) - sqrt(1 - u1^2)) * sh
    M341::T = -0.5 * (p0 - p1) * (sqrt(1 - u0^2) - sqrt(1 - u1^2)) * ch

    M122::T = 0.5 * (p0 - p1) * (u0 - u1) * sh * ch
    M132::T = -0.5 * (p0 - p1) * (u0 - u1) * ch^2

    M232::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * ch
    M242::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (u0^2 - u1^2) * sh * ch
    M342::T = 0.25 * (sqrt_mp0 - sqrt_mp1) * (u0^2 - u1^2) * ch^2

    M123::T = 0.5 * (p0 - p1) * (u0 - u1) * sh^2
    M133::T = -M122

    M233::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * sh
    M243::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (u0^2 - u1^2) * sh^2
    M343::T = -M242

    M124::T = -M241
    M134::T = -M341

    M234::T = -0.25 * (sqrt_mp0 - sqrt_mp1) * (u0^2 - u1^2)
    M244::T = -M233
    M344::T = M232

    # TODO: Assign correct values!!!!
    ForceMomentumArrayView[1,2,1] = M121
    ForceMomentumArrayView[1,3,1] = M131
    ForceMomentumArrayView[2,3,1] = M231
    ForceMomentumArrayView[2,4,1] = M241
    ForceMomentumArrayView[3,4,1] = M341

    ForceMomentumArrayView[1,2,2] = M122
    ForceMomentumArrayView[1,3,2] = M132
    ForceMomentumArrayView[2,3,2] = M232
    ForceMomentumArrayView[2,4,2] = M242
    ForceMomentumArrayView[3,4,2] = M342

    ForceMomentumArrayView[1,2,3] = M123
    ForceMomentumArrayView[1,3,3] = M133
    ForceMomentumArrayView[2,3,3] = M233
    ForceMomentumArrayView[2,4,3] = M243
    ForceMomentumArrayView[3,4,3] = M343

    ForceMomentumArrayView[1,2,4] = M124
    ForceMomentumArrayView[1,3,4] = M134
    ForceMomentumArrayView[2,3,4] = M234
    ForceMomentumArrayView[2,4,4] = M244
    ForceMomentumArrayView[3,4,4] = M344

    # antisymmetric parts
    ForceMomentumArrayView[2,1,1] = -M121
    ForceMomentumArrayView[3,1,1] = -M131
    ForceMomentumArrayView[3,2,1] = -M231
    ForceMomentumArrayView[4,2,1] = -M241
    ForceMomentumArrayView[4,3,1] = -M341

    ForceMomentumArrayView[2,1,2] = -M122
    ForceMomentumArrayView[3,1,2] = -M132
    ForceMomentumArrayView[3,2,2] = -M232
    ForceMomentumArrayView[4,2,2] = -M242
    ForceMomentumArrayView[4,3,2] = -M342

    ForceMomentumArrayView[2,1,3] = -M123
    ForceMomentumArrayView[3,1,3] = -M133
    ForceMomentumArrayView[3,2,3] = -M233
    ForceMomentumArrayView[4,2,3] = -M243
    ForceMomentumArrayView[4,3,3] = -M343

    ForceMomentumArrayView[2,1,4] = -M124
    ForceMomentumArrayView[3,1,4] = -M134
    ForceMomentumArrayView[3,2,4] = -M234
    ForceMomentumArrayView[4,2,4] = -M244
    ForceMomentumArrayView[4,3,4] = -M344

end



# OUTDATED BELOW HERE
#= # ======= IJK Cylindrical Ricci aligned to B field ======= #

    #=
    B field is formed of a axial field and azimuthal field, such that B^α = (0,0,b1/r^2,b2), local coordinates are aligned to this field such that B^a = (0,0,0,B) where B= sqrt(b1^2/r^2+b2^2)
    =#

    function IFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,h0,h1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        b1::Float64 = space_coords.b1
        b2::Float64 = space_coords.b2

        B::Float64 = b1/b2

        if b1 == 0 # no azimuthal field
            flux = 0e0
        elseif b2 == 0 # no axial field
            flux = (1/(48*sqrt(m^2 + p^2)))*p^2*(t0 - t1)*(x0 - x1)*(y0 - y1)*(z0 - z1)*(sinpi(h0/pi) - sinpi(h1/pi))*(-2* u0*sqrt(1 - u0^2)*(-11 + x0^2 + x0*x1 + x1^2) + 4*u0^3*sqrt(1 - u0^2)*(-7 + x0^2 + x0*x1 + x1^2) + 2* u1*sqrt(1 - u1^2)*(-11 + x0^2 + x1*(x0 + x1) - 2*u1^2*(-7 + x0^2 + x0*x1 + x1^2)) + (u0*sqrt(1 - u0^2)*(-5 + 2*u0^2) + u1*(5 - 2*u1^2)*sqrt(1 - u1^2))*(cos(2*h0) + cos(2*h1) - 2*sinpi(h0/pi)*sinpi(h1/pi)) - 2*acot_mod(u0)*(-3*cos(2*h0) - 3*cos(2*h1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(h0/pi)*sinpi(h1/pi))) + 2*acot_mod(u1)*(-3*cos(2*h0) - 3*cos(2*h1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(h0/pi)*sinpi(h1/pi))))
        else
            #flux = 
        end

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,h0,h1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        return 0f0

    end

    function KFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)
        phipi = Float64(phi/pi) # using sinpi rather than sin to avoid float issues 

        m = eval(Symbol("CONST_mu"*name))

        flux = (1/2)*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*sinpi(phipi)
        flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)

        return flux

    end

# ======================================================== # =#