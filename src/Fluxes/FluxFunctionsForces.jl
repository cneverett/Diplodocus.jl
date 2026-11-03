#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT.
    All force fluxes i.e. `flux` calculated in the functions below have dimensions of: either c*T*L^2 or L^3T/T_c depending on if the flux depends on external parameters that denote a characteristic time (e.g. magnetic field in synchrotron radiation reaction). 

=#

"""
    SpaceScalarForceFunction!(::SpaceScalarForce,PhaseSpace,x_idx,y_idx,z_idx)

Returns the spatial component of the force flux through the I,J,K surface for a given space scalar force. 
"""
SpaceScalarForceFunction!(force,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    MomentumScalarForceFunction!(::SpaceScalarForce,MScaForce::MMatrix{3,2,T},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 2D array of momentum components of the force flux through the I,J,K surface for a given space scalar force. The last dimension of the array corresponds to the plus and minus fluxes through the surface.  
"""
MomentumScalarForceFunction!(force,MScaForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    SpaceVectorForceFunction!(::SpaceVectorForce,SVecForce::MVector{4,T},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 1D array of spatial components of the force flux through the I,J,K surface for a given space vector force. 
"""
SpaceVectorForceFunction!(force,SVecForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    MomentumScalarForceFunction!(::SpaceVectorForce,MVecForce::MArray{Tuple{3,4,2},T,3,24},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 3D array of momentum components of the force flux through the I,J,K surface for a given space scalar force. The last dimension of the array corresponds to the plus and minus fluxes through the surface.   
"""
MomentumVectorForceFunction!(force,MVecForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    SpaceMatrixForceFunction!(::SpaceMatrixForce,SMatForce::MMatrix{4,4,T},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 2D array of spatial components of the force flux through the I,J,K surface for a given space vector force. 
"""
SpaceMatrixForceFunction!(force,SMatForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    MomentumMatrixForceFunction!(::SpaceMatrixForce,MMatForce::MArray{Tuple{3,4,4,2},T,4,96},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 4D array of momentum components of the force flux through the I,J,K surface for a given space scalar force. The last dimension of the array corresponds to the plus and minus fluxes through the surface.   
"""
MomentumMatrixForceFunction!(force,MMatForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    Space3TensorForceFunction!(::Space3TensorForce,S3TenForce::MArray{Tuple{4,4,4},T,3,64},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 3D array of spatial components of the force flux through the I,J,K surface for a given space vector force. 
"""
Space3TensorForceFunction!(force,S3TenForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")

"""
    Momentum3TensorForceFunction!(::Space3TensorForce,M3TenForce::MArray{Tuple{3,4,4,4,2},T,5,384},PhaseSpace,x_idx,y_idx,z_idx)

Modifies the 5D array of momentum components of the force flux through the I,J,K surface for a given space scalar force. The last dimension of the array corresponds to the plus and minus fluxes through the surface.
"""
Momentum3TensorForceFunction!(force,M3TenForce,PhaseSpace,x_idx,y_idx,z_idx) = error("Function not defined for this force type $(typeof(force))).")


# ======= Local Lorentz Force Zeroth Order ============ # 

    function SpaceMatrixForceFunction!(force::LocalLorentzForceZerothOrder,SMatForce::MMatrix{4,4,T},PhaseSpace::PhaseSpaceStruct,x_idx::Int64,y_idx::Int64,z_idx::Int64) where T

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad
        coordinates = PhaseSpace.Spacetime.coordinates

        t0 = Grids.tr[1]
        t1 = Grids.tr[2]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        @inline B_func(pos_local) = tetrad.Bfunction(pos_local,coordinates)
        @inline SpaceIntegrand!(pos_local,SpaceMatrix_local) = LocalLorentzForceZerothOrderSpaceIntegrand!(pos_local,SpaceMatrix_local,metric,coordinates)

        # flux must have no dimensions
        flux_non_dimensional::T = CONST_c * Characteristic.CHAR_time / Characteristic.CHAR_length

        a::SVector{4,T} = [t0,x0,y0,z0]
        b::SVector{4,T} = [t1,x1,y1,z1]
        n::SVector{4,Int64} = [2,16,16,16] # number of points for integration, can be changed to increase accuracy

        # Integrate space part of force
        Simpson4D!(SpaceIntegrand!,SMatForce,a,b,n)

        # put a tolerance on integration to avoid small numerical errors causing a non-zero flux when the flux should be zero
        tolerance::T = 100*eps(T)
        SMatForce[abs.(SMatForce) .< tolerance * max(abs.(SMatForce))] .= zero(T)

        SMatForce .*= flux_non_dimensional

    end

    @inline function LocalLorentzForceZerothOrderSpaceIntegrand!(pos::MVector{4,T},SMatForce::MMatrix{4,4,T},metric::AbstractMetric,coordinates::AbstractCoordinates) where T

        #= zeroth order Lorentz force

            F^a = q/m (F^{ab}p_b)
            so spatial integral is 
            ∫ χ F^{ab} dtdxdydz 
            and momentum integral is 
            ∫ D^i_j R^j_a q/m p_b * m/p^0

            F^{ab} = (0,0,0,0 ; 0,0,-B,0; 0,B,0,0; 0,0,0,0) for local observer
        =#

        fill!(SMatForce, zero(T))
        χ::T = VolumeElement(pos,metric,coordinates)

        SMatForce[1,2] = -χ * B_func(pos)
        SMatForce[2,1] = χ * B_func(pos)
        
        return nothing

    end

    function MomentumMatrixForceFunction!(force::LocalLorentzForceZerothOrder,MMatForce::MArray{Tuple{3,4,4,2},T,3,96},PhaseSpace::PhaseSpaceStruct,species_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64) where T

        #= zeroth order Lorentz force

            F^a = q/m (F^{ab}p_b)
            so spatial integral is 
            ∫ χ F^{ab} dtdxdydz 
            and momentum integral is 
            ∫ D^i_j R^j_a q/m p_b * m/p^0

            F^{ab} = (0,0,0,0 ; 0,0,-B,0; 0,B,0,0; 0,0,0,0) for local observer
        =#

        fill!(MMatForce, zero(T))

        Grids = PhaseSpace.Grids

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        h0divpi = h0/pi
        h1divpi = h1/pi
        sh0, ch0 = sincospi(h0divpi)
        sh1, ch1 = sincospi(h1divpi)

        if m != 0 && Z != 0

            # I_plus
            p = p1
                # due to symmetries these will all cancel out with Space Matrix when summed

            # I_minus
            p = p0
                # due to symmetries these will all cancel out with Space Matrix when summed 

            # J_plus 
            u = u1
                # due to symmetries these will all cancel out with Space Matrix when summed 
  
            # J_minus
            u = u0
                # due to symmetries these will all cancel out with Space Matrix when summed 

            # K_plus
            h = h1
            hdivpi = h1divpi
                # Only elements 2,3 and 3,2 are needed as F^{ab} is only non-zero for 2,3 and 3,2 
            MMatForce[3,2,3,1] = (-u0 + u1)*(asinh(p0/m) - asinh(p1/m)) * sinpi(h1divpi)^2
            MMatForce[3,3,2,1] = (-u0 + u1)*(-asinh(p0/m) + asinh(p1/m)) * cospi(h1divpi)^2

            # K_minus
            h = h0
            hdivpi = h0divpi
                # Only elements 2,3 and 3,2 are needed as F^{ab} is only non-zero for 2,3 and 3,2 
            MMatForce[3,2,3,1] = (-u0 + u1)*(asinh(p0/m) - asinh(p1/m)) * sinpi(h0divpi)^2
            MMatForce[3,3,2,1] = (-u0 + u1)*(-asinh(p0/m) + asinh(p1/m)) * cospi(h0divpi)^2

        end

        # scale by q and mEle as q/m is in dimensional units but the non-dimensional part of m cancels out with m/p^0 in momentum integral
        MMatForce .*= Z * CONST_q / CONST_mEle

        return nothing

    end

# ======= First Order Guiding Centre ============== #

    function SpaceVectorForceFunction!(force::FirstOrderGuidingCentre,SVecForce::MVector{4,T},PhaseSpace::PhaseSpaceStruct,x_idx::Int64,y_idx::Int64,z_idx::Int64) where T

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad
        coordinates = PhaseSpace.Spacetime.coordinates

        t0 = Grids.tr[1]
        t1 = Grids.tr[2]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        pos = MVector{4,T}(zeros(T,4))
        e::MMatrix{4,4,T,16} = zeros(T,4,4)
        ∂lnB = MVector{4,T}(zeros(T,4))

        @inline lnB_func(pos_local) = log(tetrad.Bfunction(pos_local,coordinates))
        cfg::ForwardDiff.GradientConfig = ForwardDiff.GradientConfig(lnB_func,pos)
        @inline SpaceIntegrand!(pos_local,IFluxSpaceVector_local) = GuidingCentreForceSpaceIntegrand!(pos_local,IFluxSpaceVector_local,metric,coordinates,tetrad,e,∂lnB,lnB_func,cfg)

        # flux must have no dimensions
        flux_non_dimensional::T = CONST_c * Characteristic.CHAR_time / Characteristic.CHAR_length

        a::SVector{4,T} = [t0,x0,y0,z0]
        b::SVector{4,T} = [t1,x1,y1,z1]
        n::SVector{4,Int64} = [2,16,16,16] # number of points for integration, can be changed to increase accuracy

        # Integrate space part of force
        Simpson4D!(SpaceIntegrand!,SVecForce,a,b,n)

        # put a tolerance on integration to avoid small numerical errors causing a non-zero flux when the flux should be zero
        tolerance::T = 100*eps(T)
        SVecForce[abs.(SVecForce) .< tolerance * max(abs.(SVecForce))] .= zero(T)

        SVecForce .*= flux_non_dimensional

        return nothing

    end

    @inline function GuidingCentreForceSpaceIntegrand!(pos::MVector{4,T},SVecForce::MVector{4,T},metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad,e,∂lnB,func,cfg::ForwardDiff.GradientConfig) where T

        fill!(SVecForce, zero(T))
        TetradComponents!(pos,e,metric,coordinates,tetrad)
        ForwardDiff.gradient!(∂lnB, func, pos, cfg)

        χ::T = VolumeElement(pos,metric,coordinates)

        @inbounds for a in 1:4, b in 1:4
            SVecForce[a] -= χ == zero(T) ? zero(T) : e[a,b]*∂lnB[b]*χ/2
        end

        return nothing

    end

    function MomentumVectorForceFunction!(force::FirstOrderGuidingCentre,MVecForce::MArray{Tuple{3,4,2},T,3,24},PhaseSpace::PhaseSpaceStruct,species_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64) where T

        fill!(MVecForce, zero(T))

        Grids = PhaseSpace.Grids

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        h0divpi = h0/pi
        h1divpi = h1/pi
        sh0, ch0 = sincospi(h0divpi)
        sh1, ch1 = sincospi(h1divpi)

        if m != 0 && Z != 0

            # I_plus
            MVecForce[1,2,1] = -p1^2*(-5u0*sqrt(1 - u0^2) + 2u0^3*sqrt(1 - u0^2) + 5u1*sqrt(1 - u1^2) - 2u1^3*sqrt(1 - u1^2) + 3acos(u0) - 3acos(u1))*(sh0 - sh1)/(8*sqrt(m^2 + p1^2))
            MVecForce[1,3,1] = p1^2*(-5u0*sqrt(1 - u0^2) + 2u0^3*sqrt(1 - u0^2) + 5u1*sqrt(1 - u1^2) - 2u1^3*sqrt(1 - u1^2) + 3acos(u0) - 3acos(u1))*(ch0 - ch1)/(8*sqrt(m^2 + p1^2))
            MVecForce[1,4,1] = -(h0 - h1)*p1^2*(-2u0^2 + u0^4 + 2u1^2 - u1^4)/(4*sqrt(m^2 + p1^2))
            # I_minus
            MVecForce[1,2,2] = -p0^2*(-5u0*sqrt(1 - u0^2) + 2u0^3*sqrt(1 - u0^2) + 5u1*sqrt(1 - u1^2) - 2u1^3*sqrt(1 - u1^2) + 3acos(u0) - 3acos(u1))*(sh0 - sh1)/(8*sqrt(m^2 + p0^2))
            MVecForce[1,3,2] = p0^2*(-5u0*sqrt(1 - u0^2) + 2u0^3*sqrt(1 - u0^2) + 5u1*sqrt(1 - u1^2) - 2u1^3*sqrt(1 - u1^2) + 3acos(u0) - 3acos(u1))*(ch0 - ch1)/(8*sqrt(m^2 + p0^2))
            MVecForce[1,4,2] = -(h0 - h1)*p0^2*(-2u0^2 + u0^4 + 2u1^2 - u1^4)/(4*sqrt(m^2 + p0^2))
            # J_plus 
            MVecForce[2,2,1] = -(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u1*(1 - u1^2)^(3/2)*(sh0 - sh1)
            MVecForce[2,3,1] = (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u1*(1 - u1^2)^(3/2)*(ch0 - ch1)
            MVecForce[2,4,1] = (h0 - h1)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(-1 + u1^2)^2
            # J_minus
            MVecForce[2,2,2] = -(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u0*(1 - u0^2)^(3/2)*(sh0 - sh1)
            MVecForce[2,3,2] = (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u0*(1 - u0^2)^(3/2)*(ch0 - ch1)
            MVecForce[2,4,2] = (h0 - h1)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(-1 + u0^2)^2
            # K_plus
            MVecForce[3,2,1] = -0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * sh1
            MVecForce[3,3,1] = 0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * ch1
            MVecForce[3,4,1] = 0.0
            # K_minus
            MVecForce[3,2,1] = -0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * sh0
            MVecForce[3,3,1] = 0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * ch0
            MVecForce[3,4,1] = 0.0

        end

        return nothing

    end

    function JFluxFunction(force::FirstOrderGuidingCentre,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)
        
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad
        coordinates = PhaseSpace.Spacetime.coordinates
        momentum_coords = PhaseSpace.Momentum.coordinates

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        pos = MVector{4,Float64}(zeros(Float64,4))
        e::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
        ∂lnB = MVector{4,Float64}(zeros(Float64,4))
        JFluxSpaceVector = MVector{3,Float64}(zeros(Float64,3))

        @inline lnB_func(pos_local) = log(tetrad.Bfunction(pos_local,coordinates))
        cfg::ForwardDiff.GradientConfig = ForwardDiff.GradientConfig(lnB_func,pos)
        @inline ForceSpaceIntegrand!(pos_local,JFluxSpaceVector_local) = GuidingCentreSpaceFluxVector(pos_local,JFluxSpaceVector_local,metric,coordinates,tetrad,e,∂lnB,lnB_func,cfg)

        flux::Float64 = c * T / L # L0 is the characteristic length in terms of L 
        #flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            JFluxMomentumVector = [
            -(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u*(1 - u^2)^(3/2)*(sin(h0) - sin(h1)),
            (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*u*(1 - u^2)^(3/2)*(cos(h0) - cos(h1)),
            (h0 - h1)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(-1 + u^2)^2
            ]
            a::SVector{4,Float64} = [t0,x0,y0,z0]
            b::SVector{4,Float64} = [t1,x1,y1,z1]
            n::SVector{4,Int64} = [2,16,16,16] # number of points for integration, can be changed to increase accuracy

            # Integrate space part of force
            #Simpson4D!(ForceSpaceIntegrand!,JFluxSpaceVector,a,b,n)
            JFluxSpaceVector[1] = rand(Float64) # placeholder for testing

            @inbounds for i in 1:3
                flux *= JFluxMomentumVector[i] * JFluxSpaceVector[i]
            end

        end

        return flux

    end

    function KFluxFunction(force::FirstOrderGuidingCentre,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad
        coordinates = PhaseSpace.Spacetime.coordinates
        momentum_coords = PhaseSpace.Momentum.coordinates

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        pos = MVector{4,Float64}(zeros(Float64,4))
        e::MMatrix{4,4,Float64,16} = zeros(Float64,4,4)
        ∂lnB = MVector{4,Float64}(zeros(Float64,4))
        KFluxSpaceVector = MVector{3,Float64}(zeros(Float64,3))

        @inline lnB_func(pos_local) = tetrad.Bfunction(pos_local,coordinates)
        cfg::ForwardDiff.GradientConfig = ForwardDiff.GradientConfig(lnB_func,pos)
        @inline ForceSpaceIntegrand!(pos_local,KFluxSpaceVector_local) = GuidingCentreSpaceFluxVector(pos_local,KFluxSpaceVector_local,metric,coordinates,tetrad,e,∂lnB,lnB_func,cfg)

        #flux::Float64 = T*invTc
        flux::Float64 = c * T / L # L0 is the characteristic length in terms of L 

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            KFluxMomentumVector = [
            -0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * sin(h),
            0.5*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * cos(h),
            0.0
            ]
            a::SVector{4,Float64} = [t0,x0,y0,z0]
            b::SVector{4,Float64} = [t1,x1,y1,z1]
            n::SVector{4,Int64} = [2,8,8,8] # number of points for integration, can be changed to increase accuracy

            # Integrate space part of force
            #Simpson4D!(ForceSpaceIntegrand!,KFluxSpaceVector,a,b,n)
            KFluxSpaceVector[1] = rand(Float64) # placeholder for testing

            @inbounds for i in 1:3
                flux *= KFluxMomentumVector[i] * KFluxSpaceVector[i]
            end

        end

        return flux

    end

# ======= IJK GradB Inv Z Decay Cylindrical ============== #

    function IFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        #flux::Float64 = T*invTc
        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            flux *= -(1/(32(m^2 + p^2)^(1/2))) * p^2 * (2(u0 - u1)*(u0 + u1)*(-2 + u0^2 + u1^2)*(h0 - h1)*cospi(β) + 2(u0*sqrt(1 - u0^2)*(-5 + 2u0^2) + u1*(5 - 2u1^2)*sqrt(1 - u1^2) - 3asin(u0) + 3asin(u1))*sinpi(β)*sinpi((h0 - h1)/(2pi))*sinpi(γ - h0/(2pi) - h1/(2pi))) 

            # space part (independent of momentum)
            flux *= (log(z0) - log(z1))
            flux *= (t0 - t1) * (x0 - x1) * (x0 + x1) * (y0 - y1)

        end

        return flux

    end

    function JFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)
        
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 
        #flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            # momentum part (independent of space)
            #flux *= m^2 * ((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1 + u^2) * ((-1 + u^2) * (h0 - h1) * cospi(β) + u * sqrt(1 - u^2) * (cospi(γ - h0/pi) - cospi(γ - h1/pi)) * sinpi(β)))/(2 * sqrt((m^2 + p0^2) * (m^2 + p1^2))) 

            flux *= ((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1 + u^2) * ((-1 + u^2) * (h0 - h1) * cospi(β) + u * sqrt(1 - u^2) * (cospi(γ - h0/pi) - cospi(γ - h1/pi)) * sinpi(β)))/(4) 


            # space part (independent of momentum)
            flux *= (log(z0) - log(z1)) 
            flux *= (t0 - t1) * (x0 - x1) * (x0 + x1) * (y0 - y1)

        end

        return flux

    end

    function KFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]        

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        #flux::Float64 = T*invTc
        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            # momentum part (independent of space)
            #flux *= -(m^2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))  * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * cospi(γ - h/pi) * sinpi(β))/(4 * sqrt((m^2 + p0^2) * (m^2 + p1^2)))

            flux *= -((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))  * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * cospi(γ - h/pi) * sinpi(β))/(8)

            # space part (independent of momentum)
            flux *= (log(z0) - log(z1))
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1)
        end

        return flux

    end


# ================= IJK ExB Cartesian ==================== #

    function IFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Roation = PhaseSpace.GlobalToLocalRotation
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

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else
            # momentum part (independent of space)
            flux *= -(1/2) * E0* (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0/pi) - cospi(h1/pi))

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

    function JFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time    

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else

            # momentum part (independent of space)
            flux *= E0 * log(p1/p0) * u*sqrt(1 - u^2) * (cospi(h0/pi) - cospi(h1/pi))

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

    function KFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time    

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else

            # momentum part (independent of space)
            flux *= -E0 * (-2acot_mod(u0)+2acot_mod(u1)) * cospi(phi/pi) * log(p1/p0) + B0 * (u0 - u1) * (-asinh(p0/m) + asinh(p1/m)) # sign of E0 flipped to make sure direction is correct ?????

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

# ======================================================== #


# TODO: Update this function for aximuthal B and radial E to meet standards of functions above
# =============== IJK ExB Cylindrical ================== #

    function IFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c) 
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= 1/2 * E0 * (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1) * xi
            flux *= (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1))
            flux *= fluxScale
        end

        return flux

    end

    function JFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time
        
        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) * xi
            flux *= E0 * u * sqrt(1 - u^2) * log(p1/p0)*(-sinpi(h0/pi) + sinpi(h1/pi))
            flux *= fluxScale
        end

        return flux

    end

    function KFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) * xi
            flux *= (1/2 * B0 * (-u0 + u1) * (-2asinh(p0/m) + 2asinh(p1/m)) + 2*E0 * acot_mod(u0) * log(p1/p0) * sinpi(phi/pi) - 2*E0 * acot_mod(u1) * log(p1/p0) * sinpi(phi/pi))
            flux *= fluxScale
        end

        return flux

    end

# ======================================================== 