"""
    CharacteristicStruct()

A struct for storing the characteristic scales (in SI units) for the simulation, from which all other values are normalised to.
"""
@kwdef struct CharacteristicStruct

    # fundamental scales (SI units)
    CHAR_mass::Float64 = CONST_mEle
    CHAR_charge::Float64 = CONST_q
    CHAR_speed::Float64 = CONST_c
    CHAR_magnetic_field::Float64 = 1.0 # 1 Tesla

    # user defined scales (SI units)
    CHAR_length::Float64 = (isdefined(Main,Symbol("CHAR_length")) && isdefined(Main,Symbol("CHAR_time"))) ? error("Both CHAR_time and CHAR_length are defined. Only one of these should be defined as the other is derived from it using CHAR_speed.") : isdefined(Main,Symbol("CHAR_length")) ? getfield(Main,Symbol("CHAR_length")) : isdefined(Main,Symbol("CHAR_time")) ? getfield(Main,Symbol("CHAR_time")) * CHAR_speed : CONST_pcs
    #    catch
    #        if isdefined(Main,Symbol("CHAR_time"))
    #            getfield(Main,Symbol("CHAR_time")) * CHAR_speed
    #        else
    #            CONST_pcs
    #        end
    #    end

    CHAR_time::Float64 = isdefined(Main,Symbol("CHAR_time")) ? getfield(Main,Symbol("CHAR_time")) : CHAR_length / CHAR_speed 

    CHAR_number_density::Float64 = isdefined(Main,Symbol("CHAR_number_density")) ? getfield(Main,Symbol("CHAR_number_density")) : 1.0 # 1 particle per m^3    

    # derived scales (SI units)
    CHAR_momentum::Float64 = CHAR_mass*CHAR_speed
    CHAR_energy::Float64 = CHAR_mass*CHAR_speed^2
    CHAR_mass_density::Float64 = CHAR_mass*CHAR_number_density
    CHAR_energy_density::Float64 = CHAR_mass_density*CHAR_speed^2
    CHAR_pressure::Float64 = CHAR_energy_density
    CHAR_electric_field::Float64 = CHAR_magnetic_field*CONST_c

    # scales used in DiplodocusCollisions and for normalisation of Big and Flux Matrices
    Bin_Norm::Float64 = CONST_σT * CONST_c * CHAR_time * CHAR_number_density
    Emi_Norm::Float64 = CONST_σT * CONST_c * CHAR_time
    Flux_Norm::Float64 = CONST_c * CHAR_time / CHAR_length
    
end


"""
    TimeStruct()

A struct for storing the time domain of the simulation.
"""
@kwdef struct TimeStruct

    dt0::Float64 = 1.0 # dt used for building flux terms

    #t_low::Float64   = getfield(Main,Symbol("t_low"))
    #t_num::Int64     = getfield(Main,Symbol("t_num"))
    #t_grid::String   = getfield(Main,Symbol("t_grid"))

end

"""
    SpacetimeStruct()

A struct for storing the space domain of the simulation.
"""
@kwdef struct SpacetimeStruct

    metric::AbstractMetric = getfield(Main,Symbol("metric"))
    coordinates::AbstractCoordinates = getfield(Main,Symbol("coordinates"))
    tetrad::AbstractTetrad = getfield(Main,Symbol("tetrad"))

    scheme::String = try 
            getfield(Main,Symbol("space_scheme"))
        catch
            "upwind"
        end

    dt0::Float64 = 1.0 # dt used for building flux terms

    x_grid::AbstractSpacetimeGrid  = getfield(Main,Symbol("x_grid"))
    y_grid::AbstractSpacetimeGrid  = getfield(Main,Symbol("y_grid"))
    z_grid::AbstractSpacetimeGrid  = getfield(Main,Symbol("z_grid"))

    x_num::Int64 = x_grid.num
    y_num::Int64 = y_grid.num
    z_num::Int64 = z_grid.num

end

"""
    MomentumStruct()

A struct for storing the momentum domain of the simulation.
"""
@kwdef struct MomentumStruct 

    coordinates::AbstractCoordinates = try 
            getfield(Main,Symbol("momentum_coords"))
        catch
            ModifiedSpherical()
        end
    scheme::String = try 
            getfield(Main,Symbol("momentum_scheme"))
        catch
            "upwind"
        end

    px_up_list::Vector{Float64}     = getfield(Main,Symbol("px_up_list"))
    px_low_list::Vector{Float64}    = getfield(Main,Symbol("px_low_list"))
    px_grid_list::Vector{String}    = getfield(Main,Symbol("px_grid_list"))
    px_num_list::Vector{Int64}      = getfield(Main,Symbol("px_num_list"))

    py_up_list::Vector{Float64}     = getfield(Main,Symbol("py_up_list"))
    py_low_list::Vector{Float64}    = getfield(Main,Symbol("py_low_list"))
    py_grid_list::Vector{String}    = getfield(Main,Symbol("py_grid_list"))
    py_num_list::Vector{Int64}      = getfield(Main,Symbol("py_num_list"))

    pz_up_list::Vector{Float64}     = getfield(Main,Symbol("pz_up_list"))
    pz_low_list::Vector{Float64}    = getfield(Main,Symbol("pz_low_list"))
    pz_grid_list::Vector{String}    = getfield(Main,Symbol("pz_grid_list"))
    pz_num_list::Vector{Int64}      = getfield(Main,Symbol("pz_num_list"))

end

"""
    GridsStruct()

A struct for storing the grid values for each particle in the simulation.
"""
@kwdef mutable struct GridsStruct <: Function

    n_space::Int64 
    n_momentum::Int64

    mass_list::Vector{Float64}
    charge_list::Vector{Float64}

    tr::Vector{Float64}

    xr::Vector{Float64}
    yr::Vector{Float64}
    zr::Vector{Float64}

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    mx::Vector{Float64}
    my::Vector{Float64}
    mz::Vector{Float64}

    pxr_list::Vector{Vector{Float64}}
    pyr_list::Vector{Vector{Float64}}
    pzr_list::Vector{Vector{Float64}}

    dpx_list::Vector{Vector{Float64}}
    dpy_list::Vector{Vector{Float64}}
    dpz_list::Vector{Vector{Float64}}

    mpx_list::Vector{Vector{Float64}}
    mpy_list::Vector{Vector{Float64}}
    mpz_list::Vector{Vector{Float64}}

    dE_list::Vector{Vector{Float64}}

    momentum_species_offset::Vector{Int64}

    B_field::Array{Float64,3} # magnetic field values on the x,y,z grid
    E_field::Array{Float64,3} # electric field values on the x,y,z grid

    function GridsStruct(name_list,spacetime::SpacetimeStruct,momentum::MomentumStruct)

        self = new()

        # number of grid points in space and momentum domains

            self.n_space = spacetime.x_num * spacetime.y_num * spacetime.z_num
            self.n_momentum = sum(momentum.px_num_list .* momentum.py_num_list .* momentum.pz_num_list)

        # spacetime domain

            self.tr = [0.0,spacetime.dt0] # only first set

            x_grid  = spacetime.x_grid
            y_grid  = spacetime.y_grid
            z_grid  = spacetime.z_grid

            self.xr = SpacetimeGrid(x_grid)
            self.yr = SpacetimeGrid(y_grid)
            self.zr = SpacetimeGrid(z_grid)

            self.dx = deltaVector(self.xr)
            self.dy = deltaVector(self.yr)
            self.dz = deltaVector(self.zr)

            self.mx = meanVector(self.xr)
            self.my = meanVector(self.yr)
            self.mz = meanVector(self.zr)

        # momentum domain grids

            px_up_list      = momentum.px_up_list
            px_low_list     = momentum.px_low_list
            px_grid_list    = momentum.px_grid_list
            px_num_list     = momentum.px_num_list

            py_up_list      = momentum.py_up_list
            py_low_list     = momentum.py_low_list
            py_grid_list    = momentum.py_grid_list
            py_num_list     = momentum.py_num_list

            pz_up_list      = momentum.pz_up_list
            pz_low_list     = momentum.pz_low_list
            pz_grid_list    = momentum.pz_grid_list
            pz_num_list     = momentum.pz_num_list

            num_species     = length(name_list);

            self.mass_list  = Vector{Float64}(undef,num_species)
            self.charge_list=Vector{Float64}(undef,num_species)

            self.pxr_list   = Vector{Vector{Float64}}(undef,num_species);
            self.pyr_list   = Vector{Vector{Float64}}(undef,num_species);
            self.pzr_list   = Vector{Vector{Float64}}(undef,num_species);

            self.dpx_list   = Vector{Vector{Float64}}(undef,num_species);
            self.dpy_list   = Vector{Vector{Float64}}(undef,num_species);
            self.dpz_list   = Vector{Vector{Float64}}(undef,num_species);

            self.mpx_list   = Vector{Vector{Float64}}(undef,num_species);
            self.mpy_list   = Vector{Vector{Float64}}(undef,num_species);
            self.mpz_list   = Vector{Vector{Float64}}(undef,num_species);

            self.dE_list    = Vector{Vector{Float64}}(undef,num_species);

            self.momentum_species_offset = Vector{Int64}(undef,num_species);

            for i in eachindex(name_list)

                self.mass_list[i]   = eval(Symbol("CONST_mu"*name_list[i]));
                self.charge_list[i] = eval(Symbol("CONST_z"*name_list[i]));

                self.pxr_list[i]    = bounds(px_low_list[i],px_up_list[i],px_num_list[i],px_grid_list[i]);
                self.pyr_list[i]    = bounds(py_low_list[i],py_up_list[i],py_num_list[i],py_grid_list[i]);
                self.pzr_list[i]    = bounds(pz_low_list[i],pz_up_list[i],pz_num_list[i],pz_grid_list[i]);

                self.dpx_list[i]    = deltaVector(self.pxr_list[i]);
                self.dpy_list[i]    = deltaVector(self.pyr_list[i]);
                self.dpz_list[i]    = deltaVector(self.pzr_list[i]);
                self.mpx_list[i]    = meanVector(self.pxr_list[i]);
                self.mpy_list[i]    = meanVector(self.pyr_list[i]);
                self.mpz_list[i]    = meanVector(self.pzr_list[i]);

                self.dE_list[i]     = deltaEVector(self.pxr_list[i],self.mass_list[i]) ./ self.dpx_list[i];

                if i == 1
                    self.momentum_species_offset[i] = 0
                else
                    self.momentum_species_offset[i] = self.momentum_species_offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
                end

            end

            if spacetime.tetrad isa ElectromagneticTetrad
                println("Building electromagnetic field grids...")
                # build electromagnetic field grids
                self.B_field, self.E_field = ElectromagneticFieldGrid(spacetime,spacetime.tetrad,self)
            else
                self.B_field = zeros(Float64,0,0,0)
                self.E_field = zeros(Float64,0,0,0)
            end

        return self

    end

end

"""
    PhaseSpaceStruct(name_list,time,space,momentum,Binary_list,Emi_list,forces)

A struct for storing the phase space of the simulation.
"""
@kwdef struct PhaseSpaceStruct

    # list of particle names
    name_list::Vector{String}           = getfield(Main,Symbol("name_list"))
    Characteristic::CharacteristicStruct= CharacteristicStruct()
    Spacetime::SpacetimeStruct          = SpacetimeStruct()    
    Momentum::MomentumStruct            = MomentumStruct()

    # grids
    Grids::GridsStruct  = GridsStruct(name_list,Spacetime,Momentum)

end
