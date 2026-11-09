#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT.
    All force fluxes i.e. `flux` calculated in the functions below have dimensions of: either c*T*L^2 or L^3T/T_c depending on if the flux depends on external parameters that denote a characteristic time (e.g. magnetic field in synchrotron radiation reaction). 

=#

# ======= IJK Sync Rad Reaction Cylindrical ============== #

    function IFluxFunction(force::SyncRadReact,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Spacetime.coordinates
        momentum_coords = PhaseSpace.Momentum.coordinates

        @assert momentum_coords isa ModifiedSpherical "SyncRadReact force only implemented for modified spherical momentum coordinates"

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

        if m == 0 || Z == 0
            return 0.0
        end

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

        # Evaluate the flux through the I surface i.e. surface of constant p

        mode = force.mode
        if isnothing(force.B) && isnothing(force.B_sampled)
            B = Grids.B_field[x_idx,y_idx,z_idx]
        elseif isnothing(force.B) && !isnothing(force.B_sampled)
            B_true = Grids.B_field[x_idx,y_idx,z_idx]
            B_idx = findmin(abs.(force.B_sampled .- B_true))[2]
            B = force.B_sampled[B_idx]
        else
            B = force.B
        end

        # fluxScale must have no dimensions
        invTc = (Z^4*B^2)/(CONST_μ0*m^3*CONST_mEle*CONST_c^2) * CONST_σT * CONST_c # old normalisation / (CONST_σT*CONST_c) 
        T = Characteristic.CHAR_time
        
        flux::Float64 = T*invTc

        # momentum part (independent of space)
        if typeof(mode) == Ani
            flux *= 1/(3*m)
            flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
        elseif typeof(mode) == Axi # average over azimuthal angles
            flux *= 1/(3*m)
            flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
        elseif typeof(mode) == Iso # averaged force over all angles
            flux *= -2/(3*m)
            flux *= p*sqrt(m^2 + p^2) * (h0 - h1) * (u0 - u1)
        else
            error("Synchrotron mode not recognised.")
        end

        # space part (independent of momentum)
        if space_coords isa Cartesian
            flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1) 
        elseif space_coords isa Cylindrical
            flux *= (1/2) * (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
        elseif space_coords isa Spherical
            flux *= -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
        elseif space_coords isa Paraboloidal
            flux *= 1/8 * (t0 - t1) * (x0 - x1) * (y0^2 - y1^2) * (z0^2 - z1^2) * (y0^2 + y1^2 + z0^2 + z1^2)
        else
            error("Space co-ordinate system not recognised.")
        end

        return flux

    end

    function JFluxFunction(force::SyncRadReact,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Spacetime.coordinates
        momentum_coords = PhaseSpace.Momentum.coordinates

        @assert momentum_coords isa ModifiedSpherical "SyncRadReact force only implemented for modified spherical momentum coordinates"

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

        if m == 0 || Z == 0
            return 0.0
        end

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

        # Evaluate the flux through the J surface i.e. surface of constant u

        mode = force.mode
        if isnothing(force.B)
            B = Grids.B_field[x_idx,y_idx,z_idx]
        else
            B = force.B
        end

        invTc = (Z^4*B^2)/(CONST_μ0*m^3*CONST_mEle*CONST_c^2) * CONST_σT * CONST_c # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        # momentum part (independent of space)
        if typeof(mode) == Ani
            flux *= (1/2) * m * u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
        elseif typeof(mode) == Axi
            flux *= (1/2) * m * u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
        elseif typeof(mode) == Iso
            flux *= 0e0
        end

        # space part (independent of momentum)
        if space_coords isa Cartesian
            flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
        elseif space_coords isa Cylindrical
            flux *= (1/2) * (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
        elseif space_coords isa Spherical
            flux *= -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
        elseif space_coords isa Paraboloidal
            flux *= 1/8 * (t0 - t1) * (x0 - x1) * (y0^2 - y1^2) * (z0^2 - z1^2) * (y0^2 + y1^2 + z0^2 + z1^2)
        else
            error("Space co-ordinate system not recognised.")
        end

        return flux

    end

    function KFluxFunction(force::SyncRadReact,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Spacetime.coordinates
        momentum_coords = PhaseSpace.Momentum.coordinates

        @assert momentum_coords isa ModifiedSpherical "SyncRadReact force only implemented for modified spherical momentum coordinates"

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

        if m == 0 || Z == 0
            return 0.0
        end

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

        # Evaluate the flux through the K surface i.e. surface of constant phi

        mode = force.mode
        if isnothing(force.B)
            B = Grids.B_field[x_idx,y_idx,z_idx]
        else
            B = force.B
        end

        invTc = (Z^4*B^2)/(CONST_μ0*m^3*CONST_mEle*CONST_c^2) * CONST_σT * CONST_c # old normalisation / (CONST_σT*CONST_c)  
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        # momentum part (independent of space)
        if typeof(mode) == Ani
            flux *= 0e0
        elseif typeof(mode) == Axi
            flux *= 0e0
        elseif typeof(mode) == Iso
            flux *= 0e0
        end

        # space part (independent of momentum)
        if space_coords isa Cartesian
            flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
        elseif space_coords isa Cylindrical
            flux *= (1/2) * (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
        elseif space_coords isa Spherical
            flux *= -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
        elseif space_coords isa Paraboloidal
            flux *= 1/8 * (t0 - t1) * (x0 - x1) * (y0^2 - y1^2) * (z0^2 - z1^2) * (y0^2 + y1^2 + z0^2 + z1^2)
        else
            error("Space co-ordinate system not recognised.")
        end

        return flux

    end
