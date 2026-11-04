mutable struct OutputStruct{T<:AbstractFloat,V<:AbstractVector{T}}
    
    f::Vector{V}
    t::Vector{T}

    function OutputStruct(f0::V,n_save::Int64) where {T<:Union{Float32,Float64},V<:AbstractVector{T}}

        self = new{T,V}()
        self.f = [similar(f0) for _ in 1:n_save]
        self.t = Vector{T}(undef,n_save)

        return self
    
    end

end


function Solve(method::AbstractSteppingMethod,dt_initial::T,t_save::Vector{T};progress::Bool=false,fileName::String=nothing,fileLocation::String=pwd(),Verbose::Int64=1) where T

    if isdir(fileLocation) == false
        mkpath(fileLocation)
    end

    PhaseSpace = method.PhaseSpace

    f=zeros(eltype(method.f_init),size(method.f_init))
    copyto!(f,method.f_init)
    n_save = length(t_save)
    output = OutputStruct(f,n_save)

    # save initial state (step 1)
    copyto!(output.f[1], f)
    output.t[1] = t_save[1]

    # progress bar
    if progress
        p = Progress(n_save)
    end

    copyto!(method.f,method.f_init) # reset f to initial condition at start of each solve (important for multiple solves in same session)
    #println(sum(method.f))
    method.step = 0
    
    dt::T = dt_initial

    if method isa ExpRBKIOPSStruct 
        fill!(method.dt_guess,ldexp(dt_initial, -5)) # initial guess for dt is 1/32 of the initial dt
    elseif method isa ERBEKrylovStruct
        fill!(method.dt_guess,ldexp(dt_initial, -10)) # initial guess for dt is 1/32 of the initial dt
    elseif method isa ExponentialRosenbrockEulerKrylovStruct
        fill!(method.dt_guess,ldexp(dt_initial, -6)) # initial guess for dt is 1/32 of the initial dt
    elseif method isa ExponentialRosenbrockEulerKrylovMixedStruct
        @. method.dt_guess = ldexp(dt_initial * method.momentum_step, -6) # initial guess for dt is 1/32 of the initial dt multiplied by the size of the local step allowed by X_flux
        # check dt <= dt_min from X_Flux 
        if dt > method.dt_min
            @warn "Chosen dt may lead to unstable space updates"
        end
        if !ispow2(dt)
            @warn "Chosen dt is not of the form 2.0^n, may lead to floating point issues in timestep"
        end
        if !all(iszero,mod.(t_save,dt))
            @warn "Chosen t_save are not all integer multiples of dt"
        end
        # largest momentum step update 
        if method.LocalMomentumUpdate
            max_dt = method.dt_min * maximum(method.momentum_step)
            # check end step is a multiple of max_dt 
            if mod(t_save[end],max_dt) != 0
                @warn "simulation end time is not a multiple of the largest advection timescale, may result in inaccuracy"
            end
            if dt != method.dt_min 
                @assert "For this method dt_Initial: $(dt) must equal dt_min: $(method.dt_min)"
            end
        end
    end

    for i in 2:n_save # start at 2 since initial state already saved

        t_start = t_save[i-1]
        t_stop = t_save[i]

        save = false
        # perform timestep
        while !save
            dt,save = method(t_start,t_stop,dt,Verbose)
            t_start += dt
        end

        # update t and dt and save state after timestep
        copyto!(output.f[i],method.f)
        output.t[i] = t_stop
        
        if progress
            next!(p)
        end
    end
    if progress
        finish!(p)
    end

    if !isnothing(fileName)
        # save output to file
        println("Saving Simulation output")

        SolutionFileSave(output,PhaseSpace,fileLocation,fileName);

    end


    return output

end

