"""
    ForwardEuler(t_start,t_stop,dt,Verbose)

Forward (Explicit) Euler time-stepping method for the transport equation.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # update momentum and space space using f at time t

    mul!(method.df_Flux,method.F_Flux,method.f)
    @. method.df = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # create df_Emi due to emission terms
    if method.Emission_Interactions
        mul!(method.df_Emi,method.M_Emi,method.f)
        @. method.df += method.df_Emi
    end
        
    # create df_Bin due to binary interactions
    if method.Binary_Interactions
        update_Big_Bin!(method,method.f)
        @. method.df += method.df_Bin
    end

    #println("df: $(method.df)")

    @. method.df *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
    # add injection term
    @. method.df += method.df_Inj * dt_scale 

    if !isnothing(method.df_mask)
        @. method.df *= method.df_mask
    end 

    if !isfinite(sum(method.df))
        println("non-finite value in df calculation")
    end

    # remove df that are small compared to f
    #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)

    # Cr (CFL) condition check
    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 3

        Cr_Bin = 0.0
        Cr_Emi = 0.0
        Cr_Flux = 0.0

        sum_f = sum(method.f)
        if sum_f != 0.0

            # Binary CFL
            if method.Binary_Interactions
                @. method.df_tmp = method.df_Bin / method.f 
                @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
                Cr_Bin = -minimum(method.df_tmp)
            end

            # Emission CFL
            if method.Emission_Interactions
                @. method.df_tmp = method.df_Emi / method.f 
                @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
                Cr_Emi = -minimum(method.df_tmp)
            end

            # Flux CFL
            @. method.df_tmp = -method.df_Flux / method.f 
            @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
            Cr_Flux = -minimum(method.df_tmp)

        end   

    end

    # Adaptive timestep control based on CFL (Courant) number `Cr`.
    # Strategy: be conservative when increasing dt (small, gradual growth)
    # but allow immediate reduction when Cr rises. This avoids oscillation
    # caused by aggressive grow/shrink symmetry.
    dt_old = dt
    target_Cr = 1f0        # desired Courant number
    safety = 0.8f0         # safety factor
    max_growth = 1.1f0     # hard cap on growth per step
    max_shrink = 0.67f0    # hard cap on shrink per step
    gentle_growth = 1.1f0  # soft cap for routine increases
    grow_trigger_ratio = 1.1e0  # require Cr_prev/Cr >= this to allow growth

    if !(isfinite(Cr) && Cr > 0.0)
        adaptive_factor = 1.0f0
    else
        raw_factor = safety * (target_Cr / max(Cr, eps(eltype(Cr))))
        clamped = clamp(raw_factor, max_shrink, max_growth)

        if raw_factor > 1.0f0 
            # growth: only if Cr has dropped sufficiently compared to previous Cr
            if method.step == 1
                # no history: do not grow on first step
                adaptive_factor = 1.0f0
            else
                # reduce growth factor smoothly as Cr approaches target: when Cr ≈ safety*target_Cr, factor → 1.0
                target_threshold = safety * target_Cr
                reduction_factor = min(1.0f0, Cr / max(target_threshold, eps(eltype(Cr))))
                clamped = 1.0f0 + (clamped - 1.0f0) * (1.0f0 - reduction_factor)
                adaptive_factor = clamped
                println("1: ", raw_factor, " ", adaptive_factor)
            end
        else
            # shrink: apply immediately (within hard cap)
            adaptive_factor = (Cr_prev+safety*target_Cr)/2Cr#raw_factor
            println("2: ", raw_factor, " ", adaptive_factor)
        end
    end

    if method.step == 1
        adaptive_factor = 1.0
    end

    #=if method.step > 2
    if (Cr < 1e-3)
        adaptive_factor = min(1.1,(Cr_prev/Cr))
        println("5: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    elseif Cr > Cr_prev
        adaptive_factor = max(0.67,0.5Cr_prev/Cr)
        println("7: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    else
        adaptive_factor = 1.0
        println("7: ", Cr_prev, " ", Cr, " ", adaptive_factor)
    end
    end=#

    adaptive_factor = 1.0  # TODO: this prevents adaption

    # apply adaptive factor and update scaling
    dt = dt * adaptive_factor
    dt_scale = dt / dt0

    if method.step == 1
        method.Cr = Cr
    else
        method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
    end

    # will we reached the next t_save?

    t_next = t_start + dt
    if abs(t_next - t_stop) <= eps(t_stop) * 10 #eps(max(abs(t_next), abs(t_stop)))
        save = true
    elseif t_next >= t_stop
        adaptive_factor *= (t_stop - t_start) / dt # adjust adaptive factor for final time step to ensure we end exactly at t_stop 
        dt = t_stop - t_start
        save = true
    else
        save = false
    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)), Cr_adapted = $(round(Cr*adaptive_factor,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # remove df that are small compared to f again after scaling time
    @. method.df *= adaptive_factor
    #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
    #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

    #@. method.df_tmp = method.df / method.f
    #@. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
    #println(-minimum(method.df_tmp) )

    # update state vector f 
    @. method.f += method.df  

    # remove masked off domain regions
    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    # removing negative values (values less than 1f-28 for better stability)
    @. method.f = method.f * (method.f>method.n_cut) * sign(method.f)
    # hacky fix for inf values
    @. method.f = method.f * (method.f!=Inf)

    return dt,save

end

function update_Big_Bin!(method::ExplicitSteppingMethod,f::AbstractVector)
    
    n_momentum=size(method.M_Bin,2)

    for off_space in method.Bin_Domain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fView = @view(f[start_idx:end_idx])

        df_BinView = @view method.df_Bin[start_idx:end_idx]
        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

        @inbounds vol = method.Vol[off_space+1]
        mul!(df_BinView,method.M_Bin_Mul_Step,fView,vol,zero(eltype(method.df_Bin)))

    end

    return nothing

end

"""
    ForwardSymplecticEuler(t_start,t_stop,dt,Verbose)

Forward Symplectic (Semi-Implicit) Euler time-stepping method for the transport equation. Symplectic integrator updates momentum space first then physical space, i.e. a full step in momentum followed by a full step in physical space using the updated distribution function.

# Explicit method:
Evaluates `dg`, ``dg = g^{t+1}-g^{t}`` from the following expression:
```math
dg = \\left[-\\left(\\mathcal{A}^{+}+\\mathcal{A}^{-}+\\mathcal{B}+\\mathcal{C}+\\mathcal{D}+\\mathcal{I}+\\mathcal{J}+\\mathcal{K}+\\right)g^{t}+\\left(M_\\text{Emi}g^{t}+M_\\text{Bin}g^{t}g^{t}\\right)\\mathcal{V}\\right]/ \\mathcal{A}^+
````
"""
function (method::ForwardSymplecticEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    dt0 = method.dt0

    # will we reached the next t_save? TODO: adjust this for adaptive time stepping, currently assumes constant time stepping

        if isapprox(t_start + dt, t_stop)
            save = true
        elseif t_start + dt >= t_stop
            dt = t_stop - t_start
            save = true
        else
            save = false
        end

    # scaling of time stepping

        dt_scale = dt / dt0

    # update momentum space using f at time t

        # create df_PFlux due to momentum flux terms
        mul!(method.df_PFlux,method.P_Flux,method.f)
        @. method.df_Momentum = -method.df_PFlux # minus sign as flux terms are on RHS of transport equation, also resets df_Momentum

        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f)
            @. method.df_Momentum += method.df_Emi
        end
            
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method)
            @. method.df_Momentum += method.df_Bin
        end

        @. method.df_Momentum *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector

        # update f_momentum (f after momentum update)
        @. method.f_Momentum = method.f + method.df_Momentum
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Momentum = method.f_Momentum*(method.f_Momentum>=1f-10)

    # update physical space using f_Momentum (f after momentum update)

        # create df_XFlux due to space flux terms
        mul!(method.df_XFlux,method.X_Flux,method.f_Momentum)
        @. method.df_Space = -method.df_XFlux * method.invAp_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # update f_Space (f after momentum and space updates)
        @. method.f_Space = method.f_Momentum + method.df_Space
        # removing negative values (values less than 0f0 for better stability)
        @. method.f_Space = method.f_Space*(method.f_Space>=1f-20)

    # CFL condition check TODO: add adaptive time stepping based on CFL condition 

        if Verbose == 1 || Verbose == 2 || Verbose == 3

            Cr = 0.0
            Cr_Bin = 0.0
            Cr_Emi = 0.0
            Cr_PFlux = 0.0
            Cr_Momentum = 0.0
            Cr_Space = 0.0

            # Cr (CFL) condition check
            if sum(method.f) != 0.0

                if Verbose == 3
            
                    # Binary CFL
                    if method.Binary_Interactions
                        @. method.df_tmp = method.df_Bin * method.invAp_Flux * dt_scale / method.f
                        Cr_Bin = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    end

                    # Emission CFL
                    if method.Emission_Interactions
                        @. method.df_tmp = method.df_Emi * method.invAp_Flux * dt_scale / method.f
                        Cr_Emi = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    end

                    # P Flux CFL
                    @. method.df_tmp = -method.df_PFlux * method.invAp_Flux * dt_scale / method.f
                    Cr_PFlux = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback

                    # Momentum CFL
                    @. method.df_tmp = method.df_Momentum /  method.f
                    Cr_Momentum = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                    # Space CFL
                    @. method.df_tmp = method.df_Space / method.f_Momentum
                    Cr_Space = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) # non-allocating and GPU compatible without CPU fallback
                end

                # Cr is calculated for the entire time step (momentum and space updates)
                # f_tmp - f = df of the momentum step
                @. method.df_tmp = (method.df_Space + method.df_Momentum) / method.f 
                Cr = -mapreduce(x -> isnan(x) ? Inf : x, min, method.df_tmp;init=0.0) 

            end   

            if Verbose == 1 && Cr > 1.0
                println("Cr = $(round(Cr, sigdigits=3)), t=$t_start, dt=$dt, system may be unstable")
            elseif Verbose == 2
                println("Cr = $(round(Cr, sigdigits=3)), t=$t_start, dt=$dt")
            elseif Verbose == 3
                println("Cr = $(round(Cr, sigdigits=3)), Cr_Bin = $(round(Cr_Bin, sigdigits=3)), Cr_Emi = $(round(Cr_Emi, sigdigits=3)), Cr_PFlux = $(round(Cr_PFlux, sigdigits=3)), Cr_Momentum = $(round(Cr_Momentum, sigdigits=3)), Cr_Space = $(round(Cr_Space, sigdigits=3)),  t=$t_start, t_save=$t_stop , dt=$dt")
            end
        end

    # update state vector f with momentum, space and injection updates
    
        @. method.f = method.f_Space + method.df_Inj * dt_scale
        # removing negative values (values less than p_cut for better stability) and ensure positivity for CFL calculations
        @. method.f = method.f*(method.f>=method.p_cut)*sign(method.f)
        # hacky fix for inf values
        @. method.f = method.f*(method.f!=Inf)

    # return dt and save

    return dt,save

end

"""
    ExponentialRosenbrockEulerKrylov(dg,g,t,dt)

2nd order Exponential integration time-stepping method for the transport equation.

An operator splitting method is used whereby the spatial flux terms are updated explicitly for a half-step, then the momentum distribution is updated using an exponential Rosenbrock method, and finally the spatial flux terms are updated again using a half-step.

This method uses a host+device architecture for momentum updates where binary interactions maybe evaluated on a device (e.g. GPU) while multiple spatial cells may have their momentum distributions updated on the host in parallel. Spatial updates are performed on the device. 

If `LocalMomentumUpdate` is set to `true` then local time stepping is used for momentum updates with dt based on the minimum advection time for each spatial cell. 

"""
function (method::ExponentialRosenbrockEulerKrylovStruct)(t_start::T,t_stop::T,dt::T,Verbose::Int64) where T

    DevicePrecision = method.DevicePrecision
    HostPrecision = method.HostPrecision

    @assert T == HostPrecision "Host precision and timing precision do not match"

    method.step += 1

    dt0::HostPrecision = method.dt0

    # will we reached the next t_save?

    t_next::HostPrecision = t_start + dt
    if abs(t_next - t_stop) <= sqrt(eps(t_stop)) #eps(max(abs(t_next), abs(t_stop)))
        dt = t_stop - t_start # make it exact to avoid round off error in time stepping
        save = true
    elseif t_next >= t_stop
        @warn "t_next ($t_next) >= t_stop ($t_stop), with dt=$dt and t_start=$t_start"
        #adaptive_factor *= (t_stop - t_start) / dt # adjust adaptive factor for final time step to ensure we end exactly at t_stop 
        #dt = t_stop - t_start
        save = true
    else
        save = false
    end

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale_d = method.DevicePrecision(dt / dt0)
    dt_scale = method.HostPrecision(dt / dt0)

    f_d = method.f 
    fstep_d = method.fstep_d
    #ftmp_d = method.ftmp_d
    df_d = method.df_d
    df_tmp_d = method.df_tmp_d
    E_long_d = method.E_long_d

    # set fstep to intial value

        copyto!(fstep_d,f_d)

    # half space update

        mul!(df_d,method.XFlux_d,fstep_d)
        @. df_d *= -method.invAFlux_d * dt_scale_d / 2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask_d)
            @. df_d *= method.df_mask_d
        end 

        @. fstep_d += df_d

        @. fstep_d = ifelse(fstep_d < method.n_cut_d && fstep_d .* E_long_d < method.n_cut_d,zero(DevicePrecision),fstep_d)

    # momentum update

        update_momentum!(method,dt_scale;every_step=true) # this MUST updates method.fstep_d in place

    # half space update

        mul!(df_d,method.XFlux_d,fstep_d)
        @. df_d *= -method.invAFlux_d * dt_scale_d / 2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask_d)
            @. df_d *= method.df_mask_d
        end 

        @. fstep_d += df_d

        @. fstep_d = ifelse(fstep_d < method.n_cut_d && fstep_d .* E_long_d < method.n_cut_d,zero(DevicePrecision),fstep_d)

    # momentum update (only if LocalMomentumUpdate)

        if method.LocalMomentumUpdate
            update_momentum!(method,dt_scale;every_step=false) # this MUST updates method.fstep_d in place
        end

    # Cr calculation

    Cr = 0.0
    sum_f = sum(f_d) # sum of previous state
    if sum_f != zero(DevicePrecision)
        @. df_d = fstep_d - f_d
        @. df_tmp_d = df_d / f_d
        @. df_tmp_d = ifelse(isnan(df_tmp_d), DevicePrecision(Inf), df_tmp_d)
        Cr = -minimum(df_tmp_d) 
    end

    if Verbose == 1 && Cr > 1.0
        println("\r step=$(method.step), t_start=$(round(t_start,sigdigits=16)), t_end=$(round(t_start+dt,sigdigits=16)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 2
        println("\r step=$(method.step), t_start=$(round(t_start,sigdigits=16)), t_end=$(round(t_start+dt,sigdigits=16)), Cr = $(round(Cr,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    # update state vector
    @. f_d = fstep_d

    # remove masked off domain regions

    if !isnothing(method.f_mask_d)
        @. f_d *= method.f_mask_d
    end

    return dt, save

end

function update_momentum!(method::ExponentialRosenbrockEulerKrylovStruct,dt::T;every_step::Bool=true) where T

    HostPrecision = method.HostPrecision

    @assert T == HostPrecision "Host precision and timing precision do not match"

    WorkerPool = method.WorkerPool
    done = Channel{Nothing}(length(method.ActiveDomain))
    tasks_this_step::Int64 = 0

    for off_space in method.ActiveDomain

        #= Control when to update momenum  
            1. if LocalMomentumUpdate == false or every_step == true, then this updates momentum in between half space updates for every time step.
            2. if LocalMomentumUpdate == true and every_step == false, then this only updates momentum after 2 half space updates but halfway in terms of global steps towards momentum_step 
            
            e.g. for case 1:
            step:                   1         2         3         4         ...
            space half steps:       | 1/2 1/2 | 1/2 1/2 | 1/2 1/2 | 1/2 1/2 ...
            momentum updates:       |    1    |    1    |    1    |    1    ...

            e.g. for case 2 with momentum step = 2:
            step:                   1         2         3         4         ...
            space half steps:       | 1/2 1/2 | 1/2 1/2 | 1/2 1/2 | 1/2 1/2 ...
            momentum updates:       |         2         |         2         ...
        =#

        momentum_step = method.momentum_step[off_space+1]

        if method.LocalMomentumUpdate
            if every_step && momentum_step == 1 # evaluate momentum on every step if momentum_step = 1
                dt_cell = HostPrecision(momentum_step * method.dt_initial)
                put!(WorkerPool.jobs,(off_space,dt_cell,done))
                tasks_this_step += 1
            elseif !every_step && momentum_step > 1 && mod(method.step+Int(momentum_step / 2),momentum_step) == 0 # or halfway to every momentum_step if momentum_step > 1
                # i.e. if momentum_step = 64, momentum should first be evaluated at step 32 for dt=64*dt and then on every 64 steps after that, i.e. 32, 96, 160 etc.
                dt_cell = HostPrecision(momentum_step * method.dt_initial)
                put!(WorkerPool.jobs,(off_space,dt_cell,done))
                tasks_this_step += 1
            end
        elseif every_step
            put!(WorkerPool.jobs,(off_space,dt,done))
            tasks_this_step += 1
        end
        
    end

    for task in WorkerPool.tasks
    if istaskfailed(task)
        try
            fetch(task)
        catch err
            showerror(stdout, err, catch_backtrace())
            println()
        end
    end
    end

    if tasks_this_step > 0 # if there are tasks to do
        for i in 1:tasks_this_step

            #println("Waiting for worker to finish off_space $i")
            #flush(stdout)
            take!(done) # when worker is done it will put a value in the done channel, so we can wait for all workers to finish. take waits for a value to be available in the done channel and then removes if from the channel. This is a blocking operation, so it will wait until a worker is done before continuing.
        end
    end

    return nothing

end

function worker!(method::ExponentialRosenbrockEulerKrylovStruct,worker::Int,jobs::Channel{Tuple{Int,T,Channel{Nothing}}}) where T

    println("Worker $worker started")

    wait(jobs) # wait for jobs to be available in the channel before assigning the variables below to prevent them from intialiseing when `WorkerPool` is created and before `method` is fully constructed.

    HostPrecision = method.HostPrecision
    DevicePrecision = method.DevicePrecision

    @assert T == HostPrecision "Host precision and timing precision do not match"

    n_momentum = method.PhaseSpace.Grids.n_momentum
    n_space = method.PhaseSpace.Grids.n_space

    J = method.J[worker]            # on host
    J_d = method.J_d[worker]        # on device
    Jsparse = method.Jsparse[worker]# on host
    F = method.F[worker]            # on host
    F_d = method.F_d[worker]        # on device
    JF = method.JF[worker]          # on host

    MBin = method.MBin
    MBinMulStepReshape = method.MBinMulStepReshape[worker]
    MBinMulStep = method.MBinMulStep[worker]

    fold = method.fold[worker]      # on host
    fold_d = method.fold_d[worker]  # on device
    fout = method.fout[worker]      # on host

    mB = method.mB 
    mL = method.mL
    KsB = method.KsB[worker]        # on host
    KsL = method.KsL[worker]        # on host
    ϕcacheB = method.ϕcacheB[worker]# on host
    ϕcacheL = method.ϕcacheL[worker]# on host

    ϕ = method.ϕ[worker]            # on host
    δ = method.δ[worker]            # on host

    E = method.E                    # on host
    E_d = method.E_d                # on device
    Einv = method.Einv              # on host
    Einv_d = method.Einv_d          # on device

    S = method.S[worker]            # on host 
    Sinv = method.Sinv[worker]      # on host

    # arnoldi settings and tolerances
    iop = 0
    reorthogonalize = true
    arnoldi_tol = 1e-8
    ηKrylov = method.ηKrylov
    correct = true

    for (off_space, dt, done) in jobs # each job is a space index to be processed

        idx_range = n_momentum*off_space+1:n_momentum*(off_space+1)

        @inbounds HasInjection = method.HasInjection[off_space+1]

        fstep_d = @view(method.fstep_d[idx_range]) # on device
        df_Inj = @view(method.df_Inj[idx_range]) # on host
        df_Inj_d = @view(method.df_Inj_d[idx_range]) # on device

        if !HasInjection && sum(fstep_d) == zero(DevicePrecision) 
            println("space: ", off_space," no change"," dt: ", dt)
            put!(done, nothing) # signal that this job is done
            continue
        end

        @inbounds Vol = method.Vol[off_space+1]
        Vol_d = DevicePrecision(Vol)
        @inbounds invA = method.invA[off_space+1]
        invA_d = DevicePrecision(invA)

        @inbounds dt_guess = method.dt_guess[off_space+1]

        @inbounds MEmiPFlux = method.MEmiPFlux[off_space+1]

        if method.BinaryInteractions && off_space in method.BinaryDomain # mix of host and device

            if method.J_d isa CuArray
            @assert MEmiPFlux isa CuArray "MEmiPFlux must be a dense for binary interactions"
            end

            copyto!(fold_d,fstep_d)
            copyto!(fold,fstep_d) # copy to host for Arnoldi

            t = 0.0
            dt_local = dt_guess
            dt_old = dt_local

            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = HostPrecision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local
                dtscale_d = DevicePrecision(dtscale)
                if dtscale < HostPrecision(1e-8) 
                    @warn "Time step small $dt_local in cell $off_space, injection = $HasInjection, type: binary"
                end

                kϕold = kϕ
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                mul!(MBinMulStepReshape,MBin,fold_d,Vol_d,zero(DevicePrecision))
                # Form J
                @. J_d = MBinMulStep * 2
                @. J_d += MEmiPFlux
                @. MBinMulStep += MEmiPFlux

                @. J_d *= dtscale_d * invA_d

                # Form F
                mul!(F_d,MBinMulStep,fold_d,invA_d * dtscale_d,zero(DevicePrecision))

                @. F_d += df_Inj_d * dtscale_d 
                if norm(F_d) < DevicePrecision(1e-30) # no change in distribution or just small numerical round off
                    break
                end

                # Diagonal scaling of J and F to improve conditioning for Krylov subspace approximation
                J_d .*= (Einv_d') # right mul
                J_d .*= (E_d) # left mul
                F_d .*= (E_d) # left mul

                copyto!(J,J_d) # copy to host for Arnoldi and scaling
                copyto!(F,F_d) # copy to host for Arnoldi and scaling

                tau = 1e-8
                mul!(JF,J,F)
                @. S = max(abs(F),k*abs(JF)/2)
                floor = tau * maximum(S)
                @. S = max(S, floor)
                @. Sinv = 1.0 / S

                J .*= S' # right mul
                J .*= Sinv # left mul
                F .*= Sinv # left mul

                arnoldi!(KsB, J, F;m=mB,reorthogonalize=reorthogonalize,remove_drift=false,tol=arnoldi_tol,iop=iop)

                V = @view(ExponentialUtilities.getV(KsB)[:,1:end-1])

                if iop == 0 && norm(V' * V - I) > 1e-5 
                    dt_local = ldexp(dt_local, -1) # reduce time step by factor of 2
                    @warn "Krylov subspace has poor orthogonalisation $(norm(V' * V - I)), dimension $(KsB.m), max dimension $(KsB.maxiter), reducing time step to $dt_local"
                    continue
                end

                # Compute φ functions of H
                phiv!(ϕ,k,KsB,2;cache=ϕcacheB,correct=correct,errest=false) # TODO: This allocates                  
                @. δ = S * Einv * @view(ϕ[:,2]) * k
                @. fout = fold + δ
                
                @. fout = ifelse(fout < zero(HostPrecision), zero(HostPrecision),fout)

                dt_old = dt_local

                kϕ = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                errest = Inf
                while errest > ηKrylov # k from ϕv errestimate, can be more strict than energy error estimate
                    _, errest = phiv!(ϕ,kϕ,KsB,2;cache=ϕcacheB,correct=correct,errest=true) # TODO: This allocates
                    if errest > ηKrylov
                        kϕ *= 0.77 * 0.5(tanh(log10(errest/ηKrylov)+1)+1)
                    end
                end
                k = min(prevpow2_float(kϕ),2.0) # maximum increase of 2 and set to only 2^n values

                dt_local = k*dt_old
                if t + dt_local/dt > 1.0
                    dt_local = (1.0 - t) * dt
                    println("space: ", off_space," region: Binary", " t_local/dt: ", t, " dt: ", dt, " dt_local: ", dt_local)
                    t = 1.0
                else
                    println("space: ", off_space," region: Binary", " t_local/dt: ", t, " dt: ", dt, " dt_local: ", dt_local)
                    t += dt_local/dt
                end
                k = dt_local / dt_old

                if k != 1.0

                    phiv!(ϕ,k,KsB,2;cache=ϕcacheB,correct=correct,errest=false) # TODO: This allocates
                    @. δ = S * Einv * @view(ϕ[:,2]) * k
                    @. fout = fold + δ

                end

                # number and energy cut 
                @. fout = ifelse(fout < method.n_cut && fout * E < method.n_cut, zero(HostPrecision),fout)

                copyto!(fold, fout)
                copyto!(fold_d, fout) # needed to update J and F on device

            end # while not accurate

            copyto!(fstep_d, fold_d)

            @inbounds method.dt_guess[off_space+1] = dt_old # update dt_guess with last plus one dt (avoid last as this could be limited by (1.0-t))

        else # if no binary interaction this is done purely on host

            @assert MEmiPFlux isa SparseMatrixCSC "MEmiPFlux must be a SparseMatrixCSC for linear terms"

            copyto!(fold,fstep_d)

            t = 0.0

            dt_local = dt_guess #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness
            dt_old = dt_local

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = HostPrecision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                if dtscale < HostPrecision(1e-8) 
                    @warn "Time step small $dt_local in cell $off_space, injection = $HasInjection, type: linear"
                end

                kϕold = kϕ
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                # Form J
                copyto!(Jsparse,MEmiPFlux)
                @. Jsparse.nzval *= dtscale * invA

                if HasInjection
                    # Form F
                    mul!(F,Jsparse,fold)
                    @. F += df_Inj * dtscale
                    if norm(F) < HostPrecision(1e-30) # no change in distribution
                        break
                    end
                else 
                    # no injection, just linear Jacobian so use exp over phi
                    @. F = fold
                end

                Jsparse .*= Einv' # right mul
                Jsparse .*= E # left mul
                F .*= E # left mul

                tau = 1e-8
                mul!(JF,Jsparse,F)
                @. S = max(abs(F),k*abs(JF)/2)
                floor = tau * maximum(S)
                @. S = max(S, floor)
                @. Sinv = 1.0 / S

                Jsparse .*= S' # right mul
                Jsparse .*= Sinv # left mul
                F .*= Sinv # left mul
            
                arnoldi!(KsL,Jsparse,F;m=mL,reorthogonalize=reorthogonalize,tol=arnoldi_tol,iop=iop)

                V = @view(ExponentialUtilities.getV(KsL)[:,1:end-1])

                norm_unscaled = norm(V' * V - I)
                if isnan(norm_unscaled) # no J or F means no Krylov subspace generated, so just accept step as is and move on,
                    break
                end

                if iop == 0 && norm_unscaled > 1e-5 
                    dt_local = ldexp(dt_local, -1) # reduce time step by factor of 2
                    @warn "Krylov subspace has poor orthogonalisation $norm_unscaled), reducing time step to $dt_local"
                    continue
                end

                # Compute φ functions of H
                if HasInjection 
                    if KsL.m != 1 # only do ϕv error estimate if injection is present and J has entries
                        phiv!(ϕ,k,KsL,2;cache=ϕcacheL,correct=correct,errest=false) # TODO: This allocates
                        @. δ = S * Einv * @view(ϕ[:,2]) * k
                        @. fout = fold + δ
                    else
                        @. fout = fold + df_Inj * k * dtscale
                    end
                else # no injection, just linear Jacobian so use exp over phi
                    phiv!(ϕ,k,KsL,2;cache=ϕcacheL,correct=correct,errest=false) # TODO: This allocates
                    @. fout = S * Einv * @view(ϕ[:,1])
                end

                # adaptive time stepping

                @. fout = ifelse(fout < zero(HostPrecision), zero(HostPrecision),fout)

                dt_old = dt_local

                kϕ = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                errest = Inf
                if HasInjection && KsL.m != 1 # only do ϕv error estimate if injection is present and J has entries
                    kϕ = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                    errest = Inf
                    while errest > ηKrylov # k from ϕv errestimate, can be more strict than energy error estimate
                        _, errest = phiv!(ϕ,kϕ,KsL,2;cache=ϕcacheL,correct=correct,errest=true) # TODO: This allocates
                        if errest > ηKrylov
                            kϕ *= 0.77 * 0.5(tanh(log10(errest/ηKrylov)+1)+1)
                        end
                    end
                else
                    while errest > ηKrylov # k from ϕv errestimate, can be more strict than energy error estimate
                        _, errest = phiv!(ϕ,kϕ,KsL,2;cache=ϕcacheL,correct=correct,errest=true) # TODO: This allocates
                        if errest > ηKrylov
                            kϕ *= 0.77 * 0.5(tanh(log10(errest/ηKrylov)+1)+1)
                        end
                    end
                end
                k = min(prevpow2_float(kϕ),2.0) # maximum increase of 2 and set to only 2^n values

                dt_local = k*dt_old
                if t + dt_local/dt > 1.0
                    dt_local = (1.0 - t) * dt
                    println("space: ", off_space," region: Linear", " t_local/dt: ", t, " dt: ", dt, " dt_local: ", dt_local)
                    t = 1.0
                else
                    println("space: ", off_space," region: Linear", " t_local/dt: ", t, " dt: ", dt, " dt_local: ", dt_local)
                    t += dt_local/dt
                end
                k = dt_local / dt_old

                if k != 1.0
                    if HasInjection 
                        if KsL.m != 1 # only do ϕv error estimate if injection is present and J has entries
                            phiv!(ϕ,k,KsL,2;cache=ϕcacheL,correct=correct,errest=false) # TODO: This allocates
                            @. δ = S * Einv * @view(ϕ[:,2]) * k
                            @. fout = fold + δ
                        else
                            @. fout = fold + df_Inj * k * dtscale
                        end
                    else # no injection, just linear Jacobian so use exp over phi
                        phiv!(ϕ,k,KsL,2;cache=ϕcacheL,correct=correct,errest=false) # TODO: This allocates
                        @. fout = S * Einv * @view(ϕ[:,1])
                    end

                end

                @. fout = ifelse(fout < method.n_cut && fout * E < method.n_cut, zero(HostPrecision),fout)

                @. fold = fout

            end # while not accurate

            copyto!(fstep_d, fold)

            @inbounds method.dt_guess[off_space+1] = dt_old # update dt_guess with last plus one dt (avoid last as this could be limited by (1.0-t))

        end

        put!(done, nothing) # signal that this job is done

    end

    return nothing

end
