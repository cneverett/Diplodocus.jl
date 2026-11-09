
"""
    Heun(t_start,t_stop,dt,Verbose)

Heun's time-stepping method for the transport equation.
"""
function (method::HeunStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # first step

        mul!(method.df_Flux,method.F_Flux,method.f)
        @. method.df_step = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux
        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f)
            @. method.df_step += method.df_Emi
        end
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f)
            @. method.df_step += method.df_Bin
        end
        @. method.df_step *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
        # add injection term
        @. method.df_step += method.df_Inj * dt_scale 

        if !isfinite(sum(method.df_step))
            println("non-finite value in df_step calculation")
        end
        # update state vector f_step 
        @. method.f_step += method.df_step  

        # removing negative values

        @. method.f_step = method.f_step * (method.f_step>method.n_cut) * sign(method.f_step)

    # second step

        mul!(method.df_Flux,method.F_Flux,method.f_step)
        @. method.df = -method.df_Flux # minus sign as flux terms are on RHS of transport equation, also resets df_Flux
        # create df_Emi due to emission terms
        if method.Emission_Interactions
            mul!(method.df_Emi,method.M_Emi,method.f_step)
            @. method.df += method.df_Emi
        end
        # create df_Bin due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f_step)
            @. method.df += method.df_Bin
        end
        @. method.df *= method.invAp_Flux * dt_scale # Assumes Ap_flux is diagonal and stored as a vector
        # add injection term
        @. method.df += method.df_Inj * dt_scale 

        if !isfinite(sum(method.df))
            println("non-finite value in df calculation")
        end

    # update state vector df for both steps
    
        @. method.df += method.df_step

        # mask of df regions
        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        if !isfinite(sum(method.df))
            println("non-finite value in df calculation")
        end

    # Cr (CFL) condition check

        Cr = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=method.n_cut, Inf,method.df / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
            Cr = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        if method.step == 1
            method.Cr = Cr
        else
            method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        end

        # remove df that are small compared to f again after scaling time
        @. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

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

        if dt < 0.0
            error("Negative time step calculated, something went wrong with the CFL condition calculation")
        end

    # output verbose information 

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

    # update state vector f 

        two = eltype(method.f)(2)
        @. method.f += method.df / two

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing negative values (values less than 1f-28 for better stability)

        @. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end


"""
    MPE(t_start,t_stop,dt,Verbose)

Modified Patankar Euler time-stepping method for the transport equation.
"""
function (method::MPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = dt / dt0

    # Create Q matrix for Patankar update

        println("1")
        Q = copy(method.Q) # method.Q is (F+MEmi) with structural zeros for binary interactions 
        println("max diag: ",maximum(diag(Q)))
        println("min diag: ",minimum(diag(Q)))
        println(length(Q.nzval))
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))
        println("2")
        # add M_Bin_Mul_Step to Q due to binary interactions
        if method.Binary_Interactions
            update_Big_Bin!(method,method.f,Q)
        end # now Q = (F+MEmi+MBin)
        println(length(Q.nzval))
        println(maximum(abs.(Q)))
        println("max diag: ",maximum(diag(Q)))
        println("min diag: ",minimum(diag(Q)))
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))

    # set up Ay=b system (where A=Q, y=f(t+1), b=A^-f(t)+df_Inj)

        println("4")
        Q = method.invA_Flux * Q # Assumes Ap_flux is diagonal and stored as a vector, also scales Q by dt_scale
        println(length(Q.nzval)) 
        
        println("5")
        @. Q.nzval *= dt_scale

        println(maximum(diag(Q)))

        println("6")
        Q = I - Q 
        println(length(Q.nzval))
        
        println("maxnon diag: ",maximum(Q - spdiagm(diag(Q))))
        println("minnon diag: ",minimum(Q - spdiagm(diag(Q))))

        println("7")
        method.f_tmp .= method.f + method.df_Inj .* dt_scale #+ method.invA_Flux * (method.MEmi * method.f) .* dt_scale

    # do inversion, gives f(t+1) in f_tmp using GMRES with ILU preconditioner

        println(maximum(diag(Q)))
        println("8")
        Qpre = ilu(Q)
        (x,stats) = gmres(Q,method.f_tmp, M=Qpre,ldiv=true,rtol=eps(eltype(Q)),atol=eps(eltype(Q)),restart=true)
        display(stats)
        method.f_tmp .= x

        #Qfact = lu(Q)
        #ldiv!(Qfact,method.f_tmp)

    # mask of df regions

        println("9")
        @. method.df = method.f_tmp - method.f 

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check

        Cr = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=method.n_cut, Inf,method.df / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
            Cr = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        if method.step == 1
            method.Cr = Cr
        else
            method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        end

        # remove df that are small compared to f again after scaling time
        @. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

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

        if dt < 0.0
            error("Negative time step calculated, something went wrong with the CFL condition calculation")
        end

    # output verbose information 

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

    # update state vector f 

        if !isnothing(method.df_mask)
            @. method.f = (method.f_tmp * (method.df_mask)) + method.f * (1-method.df_mask)
        else
            @. method.f = method.f_tmp
        end
        #@. method.f += method.df

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        #@. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

function updateMomentumPatankarEuler!(method::SymplecticMPEStruct,f::AbstractVector{T},dt_scale::T) where T

    Spacetime = method.PhaseSpace.Spacetime
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    n_space = x_num*y_num*z_num
    
    n_momentum=size(method.Q,2)

    I_matrix = Diagonal(ones(method.Precision, n_momentum))
    Dl = Diagonal(ones(Float64, n_momentum))
    Dr = Diagonal(ones(Float64, n_momentum))
    Qs = zeros(Float64, size(method.Q))
    Qtmp = zeros(Float64, size(method.Q))
    Qflux = zeros(Float64, size(method.Q))
    Qbin = zeros(Float64, size(method.Q))
    fstep = zeros(Float64, n_momentum)
    df_Inj = zeros(Float64, n_momentum)
    ftmp = zeros(Float64, n_momentum)
    fs = zeros(Float64,size(fstep))
    sigma = Diagonal(zeros(Float64, n_momentum))
    fdivsigma = Diagonal(zeros(Float64, n_momentum))
    ftildedivsigma = Diagonal(zeros(Float64, n_momentum))
    tmpdiag = Diagonal(zeros(Float64, n_momentum))
    onesvec = ones(Float64, n_momentum)
    source = Diagonal(zeros(Float64, n_momentum))

    for off_space in 0:n_space-1 

        fill!(method.Q,zero(method.Precision)) # reset Q
        fill!(tmpdiag.diag,zero(eltype(tmpdiag.diag))) # reset tmpdiag
        fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        fill!(Qflux,zero(eltype(Qflux))) # reset Qflux

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)
        println("off space: $off_space, start idx: $start_idx, end idx: $end_idx")

        P_Flux_View = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        invA_Flux_View = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        if method.Emission_Interactions
            M_Emi_View = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        end
        df_Inj .= @view(method.df_Inj[start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx]) 
        @. fdivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        #fstep .+= (df_Inj + invA_Flux_View * M_Emi_View * fstep) * dt_scale
        #fstep[fstep .< 1e16eps(0.0)] .= 0.0 # remove negative values for stability, also sets small values to zero which allows for more sparsity and faster solves in the Patankar update
        @. sigma.diag = max(fstep,eps(Float64))
        @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0) #fstep ./ sigma.diag # this is just fstep but with a cutoff to avoid division by zero, also allows for more sparsity and faster solves in the Patankar update
        fill!(Qbin,zero(eltype(Qbin))) # reset Qbin NEEDED otherwise binary interactions from previous off_space will carry over and contaminate current off_space update
        @. sigma.diag = ifelse(fstep > 0.0, fstep,eps(Float64))
        Qflux .= P_Flux_View

        # first do advection update as binary interaction can be be large in comparision.

        #=    mul!(tmpdiag.diag,transpose(Qflux),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
            Qtmp .= invA_Flux_View * (Qflux * ftildedivsigma - tmpdiag) * dt_scale

            #Qtmp .-= invA_Flux_View * M_Emi_View * ftildedivsigma * dt_scale # minus as M_Emi is on LHS of transport equation and Q is on RHS

            #mul!(source.diag, M_Emi_View, fstep)
            #@. source.diag /= sigma.diag
            #Qbin .-= invA_Flux_View * source * dt_scale # minus as M_Bin is on LHS of transport equation and Q is on RHS

            #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            Qtmp .+= I_matrix # add identity for Patankar update 

            ftmp .= fstep

            QLU = lu(Qtmp)
            ldiv!(QLU,ftmp)

            fstep .= ftmp
            =#
            # TODO: if dt too big this leads to poor energy conservation due to excess gain of photons and fast cooling of electrons 

        # now do binary interaction update
        
        if method.Binary_Interactions && off_space in method.Bin_Domain

            @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        
            @inbounds vol = method.Vol[off_space+1]

            #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            # loss
            #mul!(ftmp,method.Lij,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape
            #tmpdiag.diag .= ftmp .* ftildedivsigma.diag
            #Qbin .-= tmpdiag # minus as M_Bin is on LHS of transport equation and Q is on RHS
            mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstep,vol,zero(method.Precision)) 
            mul!(tmpdiag.diag,transpose(method.M_Bin_Mul_Step),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag

            Qbin .= tmpdiag # minus as M_Bin is on LHS of transport equation and Q is on RHS

            # gain 
            #mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            Qbin .-= method.M_Bin_Mul_Step * ftildedivsigma # minus as M_Bin is on LHS of transport equation and Q is on RHS


            cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            colsum = sum(Qbin,dims=1)
            rowsum = sum(Qbin,dims=2)
            println("Qbin - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))

            Qtmp .= invA_Flux_View * (Qbin) * dt_scale

            Qtmp .+= I_matrix # add identity for Patankar update

            norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
            Qs .= Dl * Qtmp * Dr

            ftmp .= fstep

            fs .= Dl * ftmp

            QLU = lu(Qs)
            ldiv!(QLU,fs)

            s = svdvals(Qs)
            println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

            ftmp .= Dr * fs

            #fstep .= patankar_energy_fix!(method,ftmp,fstep)

            #fstep .= patankar_energy_fix_quadratic!(method,ftmp,fstep)

            #fstep .= patankar_energy_fix_linear!(method,ftmp,fstep)

            fstep .= patankar_energy_fix_linear_filtered!(method,ftmp,fstep)

            #fstep .= ftmp

            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")
            
        end

       #= if method.Emission_Interactions
            mul!(source.diag, M_Emi_View, fstep)
            @. source.diag /= sigma.diag # divide by fstep for source term to be in correct units for Patankar update, also sets source to zero where fstep is zero which allows for more sparsity and faster solves in the Patankar update
            println("max source diag: ", findmax(source.diag), " min source diag: ", findmin(source.diag), " sum source diag: ", sum(source.diag))
        else
            fill!(source.diag, zero(eltype(source.diag)))
        end=#
 



        #=println("dt scale: $dt_scale")
        Qtmp .+= invA_Flux_View * (Qbin) * dt_scale
        colsum = sum(Qtmp,dims=1)
        rowsum = sum(Qtmp,dims=2)
        println("Qtmp - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
        display(Qtmp[115:125,120]')

        colsum = sum(Qtmp,dims=1)
        rowsum = sum(Qtmp,dims=2)
        println("Qtmp - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
        display(Qtmp[115:125,120]')

        cutoff_offdiagonal_to_diag_relative!(Qtmp, 1e-6)
        display(Qtmp[115:125,120]')

        println("sum diag: ", sum(diag(Qtmp)), " sum off diag: ", sum(Qtmp) - sum(diag(Qtmp)))
        display(sparse(Qtmp))

        #fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        Qtmp .+= I_matrix # add identity for Patankar update

        kwo = cond(Qtmp)
        println("off space: $off_space, condition without equilibration: $kwo")

        norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=3)
        Qs .= Dl * Qtmp * Dr

        kwo = cond(Qs)
        println("off space: $off_space, condition with equilibration: $kwo")

        #=if kwo > sqrt(1/eps(eltype(Qtmp)))
            k = kwo
            step = zero(eltype(method.Precision))
            fstep .= fstep
            dstep = one(eltype(method.Precision))
            a = one(eltype(method.Precision))
            while step < one(eltype(method.Precision))
                k = method.Precision(Inf)
                while k > sqrt(1/eps(eltype(Qtmp)))

                    s = svdvals(Qtmp)
                    println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

                    if method.Binary_Interactions && off_space in method.Bin_Domain   
            
                        @inbounds vol = method.Vol[off_space+1]

                        mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fstep,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                        @. Qbin = -method.M_Bin_Mul_Step # minus as M_Bin is on LHS of transport equation and Q is on RHS

                        Qtmp .= invA_Flux_View * (Qbin + Qflux) * dt_scale * a + I_matrix # add identity for Patankar update
                
                    else

                        Qtmp .= invA_Flux_View * Qflux * dt_scale * a + I_matrix # add identity for Patankar update

                    end

                    norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=3)
                    Qs .= Dl * Qtmp * Dr

                    k = cond(Qs)

                    println("off space: $off_space, cond without: $kwo, condition with right and left: $k, a = $a, step = $step, dstep = $dstep")

                    if k > sqrt(1/eps(eltype(Qtmp)))
                        a *= method.Precision(1/2)
                        dstep *= method.Precision(1/2)
                    end
                    
                end

                if method.Emission_Interactions
                    ftmp .= fstep + df_Inj .* dt_scale * a + invA_Flux_View * M_Emi_View * fstep .* dt_scale * a # add M_Emi contribution to ftmp for Patankar update
                else
                    ftmp .= fstep + df_Inj .* dt_scale * a
                end

                s = svdvals(Qs)
                println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")

                fs .= Dl * ftmp
                    println("debug ftmp: sum=", sum(ftmp), " max=", maximum(ftmp), " sum(Dl*ftmp)=", sum(Dl * ftmp))

                QLU = lu(Qs)
                ldiv!(QLU,fs)
                fstep .= Dr * fs
                # diagnostics: check mass conservation and residual for this block
                mass_before = sum(@view(f[start_idx:end_idx]))
                mass_after = sum(fstep)
                res = norm(Qs * fs - Dl * ftmp) / max(one(eltype(fs)), norm(Dl * ftmp))
                println("off_space=$off_space mass_before=", mass_before, " mass_after=", mass_after, " Δ=", mass_after-mass_before, " res=", res)
                @. sigma.diag = ifelse(fstep != 0.0, 1.0,0.0)

                #gmres!(method.GMRESWorkspace,Qs,fs;rtol = 1f-8,atol = 0f0)
                #fstep .= Dr * method.GMRESWorkspace.x

                step += dstep

            end

        else=#

            if method.Emission_Interactions
                ftmp .= fstep #+ df_Inj .* dt_scale + invA_Flux_View * (M_Emi_View * fstep) .* dt_scale # add M_Emi contribution to ftmp for Patankar update
                #println("min ftmp: ", findmin(ftmp), " max ftmp: ", findmax(ftmp), " sum ftmp: ", sum(ftmp))
                #ftmp[ftmp .< 1e16eps(0.0)] .= 0.0
            else
                ftmp .= fstep #+ df_Inj .* dt_scale
            end

            #norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
            #Qs .= Dl * Qtmp * Dr

            #k = cond(Qs)

            #println("off space: $off_space, cond without: $kwo, condition with right and left: $k")


             
            fs .= Dl * ftmp

            QLU = lu(Qs)
            ldiv!(QLU,fs)
            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")

            #QLU = lu(Qtmp)
            #ldiv!(QLU,ftmp)

            fstep .= Dr * fs
            #fstep .= ftmp

            #gmres!(method.GMRESWorkspace,Qtmp,ftmp;rtol = 1f-10,atol = 0f0)
            #fstep .= method.GMRESWorkspace.x

        #end
        =#

        @view(method.f_step[start_idx:end_idx]) .= fstep

    end

    return nothing

end

"""
    SymplecticMPE(t_start,t_stop,dt,Verbose)

Symplectic Modified Patankar Euler time-stepping method for the transport equation. Momentum transport is updated first using the implicit MPE scheme then space transport is performed using simple Explicit Forward Euler. 
"""
function (method::SymplecticMPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    # Create Q matrix for Patankar update and solve for momentum transport (includes injection)

        updateMomentumPatankarEuler!(method,method.f,dt_scale)
        if !isnothing(method.df_mask)
            @. method.f_step = (method.f_step * method.df_mask) + method.f * (1-method.df_mask)
        end 
        if !isnothing(method.f_mask)
            @. method.f_step *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f_step = ifelse(method.f_step<=method.n_cut,zero(eltype(method.f_step)),method.f_step)

    # Cr (CFL) condition check for momentum update

        Cr_momentum = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=0.0, Inf,(method.f_step - method.f) / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_momentum = -minimum(method.df_tmp)
        end

    # space transport

        #@. method.f_step = method.f + method.df_Inj * dt_scale # add injection term to f_step for space transport
        mul!(method.df_Flux,method.X_Flux,method.f_step)
        method.df .= -method.invA_Flux * method.df_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check space update

        Cr_space = 0.0
        dt_old = dt
        sum_f = sum(method.f_step)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f_step<=0.0, Inf,(method.df / method.f_step))
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_space = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        #if method.step == 1
        #    method.Cr = Cr
        #else
        #    method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        #end

        # remove df that are small compared to f again after scaling time
        #@. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

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

        if dt < 0.0
            error("Negative time step calculated, something went wrong with the CFL condition calculation")
        end

    # output verbose information 

        if Verbose == 1 && Cr > 1.0
            println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
        elseif Verbose == 2
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr_momentum = $(round(Cr_momentum,sigdigits=3)), Cr_space = $(round(Cr_space,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        elseif Verbose == 3
            println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        end
        if Verbose > 0
            flush(stdout)
        end

    # update state vector f 

        #println(sum(method.f[1:288]), " ", sum(method.df[1:288]), " ", sum(method.f_step[1:288]))
        #println(sum(method.f[289:888]), " ", sum(method.df[289:888]), " ", sum(method.f_step[289:888]))
        println(sum(method.f), " ", sum(method.df), " ", sum(method.f_step))

        @. method.f = method.df + method.f_step

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        #@. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

function (method::SymplecticSymmetricMPEStruct)(t_start,t_stop,dt,Verbose::Int64)

    @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f)

    method.step += 1

    dt0 = method.dt0
    Cr_prev = method.Cr

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    # Create Q matrix for Patankar update and solve for momentum transport (includes injection)

        updateMomentumPatankarEulerSymmetric!(method,method.f,dt_scale)
        if !isnothing(method.df_mask)
            @. method.f_step = (method.f_step * method.df_mask) + method.f * (1-method.df_mask)
        end 
        if !isnothing(method.f_mask)
            @. method.f_step *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f_step = ifelse(method.f_step<=method.n_cut,zero(eltype(method.f_step)),method.f_step)

    # Cr (CFL) condition check for momentum update

        Cr_momentum = 0.0
        dt_old = dt
        sum_f = sum(method.f)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f<=0.0, Inf,(method.f_step - method.f) / method.f)
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_momentum = -minimum(method.df_tmp)
        end

    # space transport

        #@. method.f_step = method.f + method.df_Inj * dt_scale # add injection term to f_step for space transport
        mul!(method.df_Flux,method.X_Flux,method.f_step)
        method.df .= -method.invA_Flux * method.df_Flux * dt_scale # minus sign as flux terms are on RHS of transport equation, also resets df_Flux

    # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

    # Cr (CFL) condition check space update

        Cr_space = 0.0
        dt_old = dt
        sum_f = sum(method.f_step)
        if sum_f != 0.0
            @. method.df_tmp = ifelse(method.f_step<=0.0, Inf,(method.df / method.f_step))
            test = findmin(method.df_tmp)
            println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]]+method.f_step[test[2]] - method.f[test[2]], " min f: ", method.f[test[2]])
            Cr_space = -minimum(method.df_tmp)
        end

        adaptive_factor = 1.0

    # apply adaptive factor and update scaling

        dt = dt * adaptive_factor
        dt_scale = dt / dt0

        #if method.step == 1
        #    method.Cr = Cr
        #else
        #    method.Cr = Cr_prev + (Cr-Cr_prev) * adaptive_factor
        #end

        # remove df that are small compared to f again after scaling time
        #@. method.df *= adaptive_factor
        #@. method.df = ifelse(abs(method.df) < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)

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

        if dt < 0.0
            error("Negative time step calculated, something went wrong with the CFL condition calculation")
        end

    # output verbose information 

        if Verbose == 1 && Cr > 1.0
            println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
        elseif Verbose == 2
            println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr_momentum = $(round(Cr_momentum,sigdigits=3)), Cr_space = $(round(Cr_space,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        elseif Verbose == 3
            println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
        end
        if Verbose > 0
            flush(stdout)
        end

    # update state vector f 

        #println(sum(method.f[1:288]), " ", sum(method.df[1:288]), " ", sum(method.f_step[1:288]))
        #println(sum(method.f[289:888]), " ", sum(method.df[289:888]), " ", sum(method.f_step[289:888]))
        println(sum(method.f), " ", sum(method.df), " ", sum(method.f_step))

        @. method.f = method.df + method.f_step

    # remove masked off domain regions

        if !isnothing(method.f_mask)
            @. method.f *= method.f_mask
        end

    # removing n_cut values for saving

        @. method.f = ifelse(method.f<=method.n_cut,zero(eltype(method.f)),method.f) #method.f * (method.f>method.n_cut) * sign(method.f)
        #@. method.f = ifelse(method.f<method.n_cut,method.n_cut,method.f)
        # hacky fix for inf values
        #@. method.f = method.f * (method.f!=Inf)

    return dt,save

end

function updateMomentumPatankarEulerSymmetric!(method::SymplecticSymmetricMPEStruct,f::AbstractVector{T},dt_scale::T) where T

    Spacetime = method.PhaseSpace.Spacetime
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    n_space = x_num*y_num*z_num
    
    n_momentum=size(method.Q,2)

    I_matrix = Diagonal(ones(method.Precision, n_momentum))
    Dl = Diagonal(ones(Float64, n_momentum))
    Dr = Diagonal(ones(Float64, n_momentum))
    Qs = zeros(Float64, size(method.Q))
    Qtmp = zeros(Float64, size(method.Q))
    Qflux = zeros(Float64, size(method.Q))
    Qbin = zeros(Float64, size(method.Q))
    fstep = zeros(Float64, n_momentum)
    df_Inj = zeros(Float64, n_momentum)
    ftmp = zeros(Float64, n_momentum)
    fs = zeros(Float64,size(fstep))
    sigma = Diagonal(zeros(Float64, n_momentum))
    fdivsigma = Diagonal(zeros(Float64, n_momentum))
    ftildedivsigma = Diagonal(zeros(Float64, n_momentum))
    tmpdiag = Diagonal(zeros(Float64, n_momentum))
    onesvec = ones(Float64, n_momentum)
    source = Diagonal(zeros(Float64, n_momentum))

    for off_space in 0:n_space-1 

        fill!(method.Q,zero(method.Precision)) # reset Q
        fill!(tmpdiag.diag,zero(eltype(tmpdiag.diag))) # reset tmpdiag
        fill!(Qtmp,zero(eltype(Qtmp))) # reset Qtmp
        fill!(Qflux,zero(eltype(Qflux))) # reset Qflux

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)
        println("off space: $off_space, start idx: $start_idx, end idx: $end_idx")

        P_Flux_View = @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        invA_Flux_View = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        if method.Emission_Interactions
            M_Emi_View = @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        end
        df_Inj .= @view(method.df_Inj[start_idx:end_idx])
        fstep .= @view(f[start_idx:end_idx]) 
        @. fdivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0)
        #fstep .+= (df_Inj + invA_Flux_View * M_Emi_View * fstep) * dt_scale
        #fstep[fstep .< 1e16eps(0.0)] .= 0.0 # remove negative values for stability, also sets small values to zero which allows for more sparsity and faster solves in the Patankar update
        @. sigma.diag = max(fstep,eps(Float64))
        @. ftildedivsigma.diag = ifelse(fstep > 0.0, 1.0,0.0) #fstep ./ sigma.diag # this is just fstep but with a cutoff to avoid division by zero, also allows for more sparsity and faster solves in the Patankar update
        fill!(Qbin,zero(eltype(Qbin))) # reset Qbin NEEDED otherwise binary interactions from previous off_space will carry over and contaminate current off_space update
        @. sigma.diag = ifelse(fstep > 0.0, fstep,eps(Float64))
        Qflux .= P_Flux_View

        # first do advection update as binary interaction can be be large in comparision.

            #=mul!(tmpdiag.diag,transpose(Qflux),onesvec)
            tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
            Qtmp .= invA_Flux_View * (Qflux * ftildedivsigma - tmpdiag) * dt_scale

            #Qtmp .-= invA_Flux_View * M_Emi_View * ftildedivsigma * dt_scale # minus as M_Emi is on LHS of transport equation and Q is on RHS

            #mul!(source.diag, M_Emi_View, fstep)
            #@. source.diag /= sigma.diag
            #Qbin .-= invA_Flux_View * source * dt_scale # minus as M_Bin is on LHS of transport equation and Q is on RHS

            #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e-6)

            Qtmp .+= I_matrix # add identity for Patankar update 

            ftmp .= fstep

            QLU = lu(Qtmp)
            ldiv!(QLU,ftmp)

            fstep .= ftmp=#
            
            # TODO: if dt too big this leads to poor energy conservation due to excess gain of photons and fast cooling of electrons 

        # now do binary interaction update
        
        if method.Binary_Interactions && off_space in method.Bin_Domain

            fstepE = fstep #.* method.E

            @. ftildedivsigma.diag = ifelse(fstepE > 0.0, 1.0,0.0)
            @. sigma.diag = ifelse(fstepE > 0.0, 1/fstepE,1/eps(Float64))
        
            @inbounds vol = method.Vol[off_space+1]

            #mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fView,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

            # loss ij
                fill!(Qbin,zero(eltype(Qbin))) # reset Qbin for loss calculation

                mul!(method.M_Bin_Mul_Step_reshape,method.Lijk,fstepE,vol,zero(method.Precision))
                mul!(tmpdiag.diag,method.M_Bin_Mul_Step,fstepE) 

                println(isdiag(method.M_Bin_Mul_Step))

                Qbin .-= tmpdiag * sigma #./2 #* sigma #=+ method.M_Bin_Mul_Step * ftildedivsigma)/2=# # minus as M_Bin is on LHS of transport equation and Q is on RHS

                #mul!(method.M_Bin_Mul_Step_reshape,method.Lijk,fstepE,vol,zero(method.Precision))
                #mul!(tmpdiag.diag,method.M_Bin_Mul_Step,fstepE) 

                #Qbin .-= sigma * tmpdiag ./2 #* sigma #=+ method.M_Bin_Mul_Step * ftildedivsigma)/2=# # minus as M_Bin is on LHS of transport equation and Q is on RHS

                # gain 
                mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,fstepE,vol,zero(method.Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                Qbin .-= (method.M_Bin_Mul_Step * ftildedivsigma) #./2 #+ ftildedivsigma * method.M_Bin_Mul_Step) / 2 # minus as M_Bin is on LHS of transport equation and Q is on RHS

                #mul!(method.M_Bin_Mul_Step_reshape,method.Gijk,ftildedivsigma.diag,vol,zero(method.Precision))
                #Qbin .-= (method.M_Bin_Mul_Step * diagm(fstepE)) ./2
                # together
                #mul!(method.M_Bin_Mul_Step_reshape,method.Aijk,fstepE,vol,zero(method.Precision))
                #Qbin .= method.M_Bin_Mul_Step
                #Qbin .= -Qbin * ftildedivsigma 

                #cutoff_offdiagonal_to_diag_relative!(Qbin, 1e0)

                #mul!(tmpdiag.diag,transpose(Qflux),onesvec)
                #tmpdiag.diag .= tmpdiag.diag .* ftildedivsigma.diag
                #Qbin .+=  (Qflux * ftildedivsigma - tmpdiag)

                display(sparse(Qbin))

                println("has_positive_offdiag(Qbin): $(has_positive_offdiag(Qbin))")

                Qtmp .= invA_Flux_View * (Qbin) * dt_scale

                colsum = sum(Qtmp,dims=1)
                rowsum = sum(Qtmp,dims=2)
                println("Qbin - max col sum: ", findmax(colsum), " max row sum: ", findmax(rowsum), " min col sum: ", findmin(colsum), " min row sum: ", findmin(rowsum))
                println("Qbin - max: ", findmax(Qtmp[:,findmax(colsum)[2][2]]), " min: ", findmin(Qtmp[:,findmax(colsum)[2][2]]))
                #display(reshape(Qtmp[:,findmax(colsum)[2][2]],37,24))
                display(reshape(Qtmp[:,150],37,24))
                println("sum: $(sum(Qtmp[:,150])), Esum: $(sum(method.E .* Qtmp[:,150]))")

                println("sumQf ", sum(Qtmp * fstep), " sumEQf ", sum(method.E .*(Qtmp * fstep)), " maxQf ", maximum(Qtmp * fstep), " minQf ", minimum(Qtmp * fstep))
                println("sumQf ele ", sum(Qtmp[1:288,:] * fstep), " sumEQf ele ", sum(method.E[1:288] .*(Qtmp[1:288,:] * fstep)), "sumQf pho ", sum(Qtmp[289:end,:] * fstep), "  sumEQf pho ", sum(method.E[289:end] .*(Qtmp[289:end,:] * fstep)))

                println("max timescale: ", maximum(replace!(abs.(max.(fstep,eps(Float64)) ./ (Qtmp * fstep)),NaN=>0.0,0.0=>Inf)), " min timescale: ", minimum(replace(abs.(fstep ./ (Qtmp * fstep)),NaN=>Inf,0.0=>Inf)))

                Qtmp .+= I_matrix # add identity for Patankar update

                norm_equilibration_matrices!(Dl,Dr,Qtmp; iters=0)
                Qs .= Dl * Qtmp * Dr

                println("max 1TQ ", maximum(transpose(onesvec) * Qtmp - transpose(onesvec)), " min 1TQ ", minimum(transpose(onesvec) * Qtmp - transpose(onesvec)))
                println("max ETQ ", maximum(transpose(method.E) * Qtmp - transpose(method.E))," min ETQ ", minimum(transpose(method.E) * Qtmp - transpose(method.E)))

                println("has_positive_offdiag(Qs): $(has_positive_offdiag(Qs))")

                ftmp .= fstepE

                fs .= Dl * ftmp

                QLU = lu(Qs)
                ldiv!(QLU,fs)

                s = svdvals(Qs)
                println("σmax = $(maximum(s)), σmin = $(minimum(s)), cond = $(maximum(s)/minimum(s))")
                println("cond = $(cond(Qs))")

                ftmp .= (Dr * fs) #./ method.E

            println("num before $(sum(fstep)), num after $(sum(ftmp)), Δ = $(sum(ftmp)-sum(fstep))")
            println("eng before $(sum(method.E .* fstep)), eng after $(sum(method.E .* ftmp)), Δ = $(sum(method.E .* ftmp) - sum(method.E .* fstep))")

            fstep .= ftmp

            println("any negative fstep: ", any(fstep .< 0.0))
            println("sum of negative fstep: ", sum((fstep .< 0.0) .* fstep))
            println("sum of positive fstep: ", sum((fstep .> 0.0) .* fstep))

            println("off space: $off_space, Cr_momentum = $(-minimum(@. ifelse(fstep<=0.0, Inf,(Dr * fs - fstep) / fstep)))")
            
        end


        @view(method.f_step[start_idx:end_idx]) .= fstep

    end

    return nothing

end

"""
    BackwardsEuler(dg,g,t,dt)

Backwards (Implicit) Euler time-stepping method for the transport equation. 

"""
function (method::BackwardEulerStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)
        
    # create df_Bin due to binary interactions (jacobians are added in update_Big_Bin! if implicit)
    if method.Binary_Interactions
        update_Big_Bin!(method,dt_scale)
    end


    
    # removing negative values
    @. method.fstep = method.fstep * (method.fstep>=method.n_cut)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    return dt,save

end

function update_Big_Bin!(method::BackwardEulerStruct,dt_scale::T) where T
    
    n_momentum=size(method.M_Bin,2)
    
    f_iter = zeros(T,n_momentum)
    f_guess = zeros(T,n_momentum)
    f_tmp = zeros(T,n_momentum)
    df_tmp = zeros(T,n_momentum)
    J = method.J
    F = method.F

    tol = T(1e-6)
    maxiter = 10
    damping = false

    for off_space in method.Bin_Domain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        f = @view(method.f[start_idx:end_idx])
        fstep = @view(method.fstep[start_idx:end_idx])
        @inbounds vol = method.Vol[off_space+1] * dt_scale
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        
        f_guess .= f

        for λ in [T(1e-12), T(1e-10), T(1e-8), T(1e-6), T(1e-4), T(1e-2), T(1.0)]

            f_iter .= f_guess

            # Newton iteration
            for it in 1:maxiter

                println("Newton iteration: $it, λ: $λ")

                build_residual_jacobian!(method, f_iter, f, invA_Flux, λ * vol)

                if norm(F, Inf) < tol
                    f_guess .= f_iter
                    break
                end

                println("Residual norm: $(norm(F, Inf))")

                println("cond J: $(cond(J))")

                δ = J \ F

                δ .= δ.* (abs.(δ) .> 0.0) # remove small values
                
                display(reshape(δ, 37,24))
                display(reshape(δ ./ f_iter, 37,24))

                f_iter .-= δ

                f_iter .= f_iter .* (f_iter .> 0.0) # enforce positivity

                # Simple backtracking line search on residual norm
                #=α = one(T)
                Fn = norm(F, 2)

                while α > T(1e-12)
                    trial = f_guess .+ α .* δ
                    build_residual_jacobian!(method, trial, f, invA_Flux, vol)

                    if norm(F, 2) <= (one(T) - T(1e-4) * α) * Fn
                        f_iter .= trial
                        break
                    end

                    α *= T(0.5)
                end

                if α <= T(1e-12)
                    error("Newton line search failed to reduce the residual.")
                end=#

            end

            #error("Newton did not converge in $maxiter iterations.")

            f_guess .= f_iter

        end

        fstep .= f_iter

    end

    return nothing

end

function build_residual_jacobian!(method::BackwardEulerStruct, f_iter, f, invA_Flux, vol)

    J = method.J
    F = method.F

    I_matrix = Diagonal(ones(method.Precision, length(F)))

    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,f_iter,vol,zero(eltype(f_iter))) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

    mul!(J,invA_Flux,method.M_Bin_Mul_Step,-2.0,0.0)

    @. J += I_matrix

    mul!(F,method.M_Bin_Mul_Step,f_iter)

    tmp = f_iter .- f .- (invA_Flux * F) 

    F .= tmp

end

"""
    ES2(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. 

"""
function (method::ES2Struct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = zeros(eltype(method.f),length(method.f))

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # half momentum update

        #=mul!(ftmp,method.invImMP,method.fstep+method.df_Inj * dt_scale/2) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end =#
        
    # binary update

        if method.Binary_Interactions
            update_momentum!(method,dt_scale)
        end

    # half momentum update

        #=mul!(ftmp,method.invImMP,method.fstep+method.df_Inj * dt_scale/2) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end=#

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end

function update_momentum!(method::ES2Struct,dt_scale::T) where T

    Precision = method.Precision
    
    n_momentum=method.PhaseSpace.Grids.n_momentum
    n_space=method.PhaseSpace.Grids.n_space
    
    fold = zeros(Precision,n_momentum)
    fxdiv2= zeros(Precision,n_momentum)
    fx = zeros(Precision,n_momentum)
    fz = zeros(Precision,n_momentum)
    fnonzero = Diagonal(zeros(Precision,n_momentum))
    J = method.J

    n_species = length(method.PhaseSpace.name_list)

    M = zeros(Precision,(n_species+1),n_momentum)

    #LinearFlux = zeros(Precision,n_momentum,n_momentum)

    @view(M[1:n_species,:]) .= method.N 
    @view(M[n_species+1,:]) .= method.E

    #A = [ zeros(Precision,n_momentum,n_momentum) zeros(Precision,n_momentum) ; zeros(Precision,n_momentum)' zero(Precision) ]

    consv = zeros(n_species+1)
    consv[end] = one(Precision)

    for off_space in 0:n_space-1

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        fin = zeros(Precision,n_momentum)
        fout = zeros(Precision,n_momentum)
        df_Inj = @view(method.df_Inj[start_idx:end_idx])
        @inbounds vol = method.Vol[off_space+1]
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        A_Flux = @view(method.A_Flux[start_idx:end_idx,start_idx:end_idx])

        @. fnonzero.diag = ifelse(fold > zero(Precision), one(Precision),zero(Precision))

        #LinearFlux .= @view(method.P_Flux[start_idx:end_idx,start_idx:end_idx])
        #if method.Emission_Interactions
        #    LinearFlux .-= @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])
        #end

        if method.Binary_Interactions && off_space in method.Bin_Domain

            fold .= @view(method.f[start_idx:end_idx]) 

            #if sum(fold) == Precision(0)
            #    continue
            #end

            α = Precision(1.0)
            accurate = false
            ΔE = zero(Precision)
            step = zero(Float64)
            dstep = one(Float64)
            ϵ = eps(Precision)

            while !accurate || step < 1.0

                #fold .+= df_Inj * α * dt_scale

            # linear flux update (should conserve energy directly) so do single exponential step
            #LinearFlux .= invA_Flux * LinearFlux * fnonzero

            #display(reshape(real.(eigvals(-LinearFlux .* dt_scale)) .* fnonzero,37,24))
            #Println("max 1/τ: ", maximum(-real.(eigvals(-LinearFlux .* dt_scale) .* fnonzero) ))
            #display(fnonzero)

            #fold .= exp(-LinearFlux .* dt_scale) * fold

            #@. fnonzero.diag = ifelse(fold > 0.0, 1.0,0.0)

            # Second order in time ES2 method
                # x/2 step

                    @time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fold

                    #s = zeros(Precision,n_momentum)
                    #D = Diagonal(s)
                    #s .= max.(abs.(fin), 0.5*dt_scale*α * abs.(method.M_Bin_Mul_Step * fin), ϵ)
                    #D = Diagonal(100 .*s)

                    #println("max D: ", maximum(D), " min D : ", minimum(D))
                    #println("max D: ", maximum(s), " min s : ", minimum(s))

                    #display(D)

                    #println("max J: ", maximum(method.M_Bin_Mul_Step), " min J : ", minimum(method.M_Bin_Mul_Step))

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])

                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    A = [J (D \ df_Inj) ; zeros(Precision,n_momentum+1)']

                    #println("max J: ", maximum(J), " min J : ", minimum(J))

                    #display(sparse(A))

                    fin .= D \ fin 

                    f = [fin ; one(Precision)]

                    #fxdiv2 .= exp(method.M_Bin_Mul_Step .* dt_scale) * fin
                    #fxdiv2 .= expv(0.5*dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)
                    f .= expv(dt_scale*α,A,f;tol=Precision(1e-12),m=500)

                    fxdiv2 .= D * @view(f[1:n_momentum])

                    #display(reshape(fxdiv2,37,24))
                    #display(reshape(df_Inj*α*dt_scale,37,24))

                    #F = schur(method.M_Bin_Mul_Step)
                    #w = F.Z' * fin
                    #fxdiv2 .= F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    #println(maximum(real.(F.values)))



                    #F = schur(J)
                    #w = F.Z' * fin
                    #fxdiv2 .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    #display(F)
                    #display(sum(F.values))
                    #println(maximum(real.(F.values)))

                    fxdiv2 .= invA_Flux * fxdiv2

                    @. fxdiv2 = ifelse(fxdiv2<=method.n_cut,zero(eltype(fxdiv2)),fxdiv2)
                
                # z step 

                    #=@time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fxdiv2,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fold 

                    s = max.(abs.(fin), dt_scale*α * abs.(method.M_Bin_Mul_Step*fin), ϵ)
                    D = Diagonal(s)

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])

                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    #fz .= expv(dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)

                    fin .= D \ fin

                    fz .= D * expv(dt_scale*α,J,fin;tol=Precision(1e-8),m=100)

                    #F = schur(J)
                    #w = F.Z' * fin
                    #fz .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    @. fz = ifelse(fz<=method.n_cut,zero(eltype(fz)),fz)

                    fz .= invA_Flux * fz=#

                # x step 

                    #=@time mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fz,vol,zero(Precision)) # temp is linked to M_Bin_Mul_Step so it gets edited while maintaining is 2D shape

                    fin .= A_Flux * fxdiv2

                    s = max.(abs.(fin), 0.5*dt_scale*α * abs.(method.M_Bin_Mul_Step*fin), ϵ)
                    D = Diagonal(s)

                    method.M_Bin_Mul_Step += @view(method.M_Emi[start_idx:end_idx,start_idx:end_idx])  
                    
                    D = Diagonal(1 ./ method.E)

                    J .= D \ method.M_Bin_Mul_Step * D

                    #fx .= expv(0.5*dt_scale*α,method.M_Bin_Mul_Step,fin;tol=Precision(1e-8),m=100)

                    fin .= D \ fin

                    fx .= D * expv(0.5*dt_scale*α,J,fin;tol=Precision(1e-8),m=100)

                    #F = schur(J)
                    #w = F.Z' * fin
                    #fx .= D * F.Z * expv(0.5*dt_scale*α,F.T,w;tol=Precision(1e-8),m=100)

                    @. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)

                    fx .= invA_Flux * fx=#

                # combine 
               
                #@. fx = (fx + fz) / 2

                #@. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)

                @. fx = ifelse(fxdiv2<=method.n_cut,zero(eltype(fxdiv2)),fxdiv2)

                #println(sum(fx), " ", sum(fold), " ", sum(fx-fold), " ", sum(df_Inj * α * dt_scale) )

                #println(sum(fx[1:288]), " ", sum(fold[1:288]), " ", sum((fx-fold)[1:288]), " ", sum((df_Inj * α * dt_scale)[1:288]))

                # energy correction  
                ΔE = dot(method.E, fx) - dot(method.E, fold .+ df_Inj * α *dt_scale)

                Eold = dot(method.E, fold .+ df_Inj * α *dt_scale)
                Eerr = abs(dot(method.E, fx) / Eold - one(Precision))

                #=if Eerr < 1e-3
                    accurate = true
                    step += dstep * α
                    # error correction step to better improve accuracy``
                    #ϵ = eps(Float64(method.n_cut))
                    #Winv = method.E .* max.(fx, ϵ)
                    #q = Winv .* method.E
                    #q ./= dot(method.E, q)
                    #@. fx -= ΔE * q
  
                    #=ϵ = eps(Float64(method.n_cut))
                    Winv = method.E .* max.(fx, ϵ)
                    B = M * Diagonal(Winv) * M'
                    # precondition B with diagonal of B
                    D = Diagonal(1.0 ./ sqrt.(diag(B)))
                    Bs = D * B * D
                    q = Diagonal(Winv) * (M' * D * (Bs \ (D * consv)))
                    @. fx -= ΔE*q
                    @. fx = ifelse(fx<=method.n_cut,zero(eltype(fx)),fx)=#
                    fold .= fx
                    println("step: ", step, " α: ", α, " Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fx))
                    if Eerr < 1e-4 
                        if step + dstep * 4.0 * α < 1.0 
                            # can fit two 2x steps
                            α *= 2.0
                        end
                    end 
                else
                    accurate = false
                    println("Energy error: ", Eerr, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fx), " reducing α: ", α)
                    α *= 0.5
                end=#

                    ϵ = eps(Precision(method.n_cut))

                    Winv = #=method.E .*=# max.(fx, ϵ)

                    MW = M .* Winv'

                    L = MW * M'

                    display(M)
                    display(M * fx)
                    display([dot(method.N[1,:],fx) ; dot(method.N[2,:],fx) ; Eold])

                    r = M * fx - [dot(method.N[1,:],fx) ; dot(method.N[2,:],fx) ; Eold]
                    r = [0f0 ; 0f0 ; ΔE]

                    println("r: ", r, " ΔE: ", ΔE)

                    λ = L \ r

                    δu = Winv .* (M' * λ)

                    fx .-= δu

                    accurate = true
                    step += dstep * α

            end # while not accurate

            @. fstep = fx

        else #

            # linear flux update (should conserve energy directly) so do single exponential step
            #@view(A[1:n_momentum,1:n_momentum]) .= -invA_Flux * LinearFlux

            #fin = [fold ; one(Precision)]
            #fout = similar(fin)

            #fout .= expv(dt_scale,A,fin;tol=Precision(1e-6),m=50)

            #fstep .= @view(fout[1:n_momentum])

        end

    end

    return nothing

end

"""
    ExpRBKIOPS(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. Using the KIOPS algorithm for the exponential action approximation. 

"""
function (method::ExpRBKIOPSStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = similar(method.f)

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end 

    # binary update

        update_momentum!(method,dt_scale)

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end

function update_momentum!(method::ExpRBKIOPSStruct,dt::T) where T

    Precision = method.Precision
    
    n_momentum = method.PhaseSpace.Grids.n_momentum
    #n_space = method.PhaseSpace.Grids.n_space
    momentum_species_offset = method.PhaseSpace.Grids.momentum_species_offset
    name_list = method.PhaseSpace.name_list
    num_species = length(name_list)
    
    fold = method.fold
    fout = method.fout
    J = method.J
    F = method.F
    g = method.g
    fscale = method.fscale #::Vector{Precision} = zeros(Precision,n_momentum)
    onevec = copy(F)
    fill!(onevec,one(Precision))
    D = method.D #Diagonal(ones(Precision,n_momentum)) # 1/E
    Dinv = method.Dinv #Diagonal(ones(Precision,n_momentum)) # E
    E = method.E

    δ = method.δ #::Vector{Precision} = zeros(Precision,n_momentum)

    EmiTrue::Bool = true

    KIOPS_workspace = method.KIOPS_workspace
    mmax = KIOPS_workspace.mmax
    mmin = KIOPS_workspace.mmin

    for off_space in method.ActiveDomain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        df_Inj = @view(method.df_Inj[start_idx:end_idx])

        has_injection = sum(df_Inj) > zero(Precision)

        if !has_injection && sum(fstep) == zero(Precision) 
            continue
        end

        @inbounds vol = method.Vol[off_space+1]
        @inbounds invA = method.invA[off_space+1]
        @inbounds dt_guess = method.dt_guess[off_space+1]

        if method.Binary_Interactions && off_space in method.Bin_Domain

            if isassigned(method.M_Emi, off_space+1)
                EmiTrue = true
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt_guess # initial guess for local time step, can adjust based on desired accuracy and problem stiffness
            # TODO: have each cell its own dt initial guess that can then be updated each completed timestep.
            dt_next = dt_local 

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                #dt_local = dt_next

                dtscale = Precision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive
                dt_old = dt_local

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision))
                    # Form J
                    @. J = Precision(2) * method.M_Bin_Mul_Step
                    if EmiTrue
                        @. J += M_Emi
                        @. method.M_Bin_Mul_Step += M_Emi
                    end

                    println("max M_Emi: $(maximum(M_Emi)), max M_Bin: $(maximum(method.M_Bin_Mul_Step)), off_space: $off_space, dt_local: $dt_local, dtscale: $dtscale, vol: $vol, invA: $invA")

                    J .*= dtscale * invA
                    # Form F
                    mul!(F,method.M_Bin_Mul_Step,fold)
                    F *= invA

                    # Energy error estimate from F
                    ηEest = abs(dot(E,F)) * dtscale / dot(E,fold)
                    println("Energy error estimate from F: ", ηEest)

                    @. F += #=A *=# df_Inj #* dtldt0
                    #res = sum(abs, F)
                    #display(F)
                    if all(iszero, F) # no change in distribution
                        break
                    end
                    lmul!(dtscale,F)

                    
                    # get residuals 
                    mul!(g,J,fold)
                    println("norm g: $(norm(g)), norm J: $(norm(J))")
                    LinearAlgebra.axpby!(Precision(1),F,Precision(-1),g)


                    # Diagonal scaling of J and g to improve conditioning for Krylov subspace approximation
                    J .*= (D)' # right mul
                    J .*= (Dinv) # left mul
                    F .*= (Dinv) # left mul
                    g .*= (Dinv) # left mul

                    # form U
                    fscale .= (Dinv) .* fold
                    # KIOPS 

                        tau = one(Precision)
                        #_, stats = kiops_roe!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)
                        _, stats = kiops_roe_noview!(KIOPS_workspace,tau,J,fscale,F;tol=1e-16)
                        #_, stats = kiops_roe_panel!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)
                        #_, stats = kiops_roe_panel_new!(KIOPS_workspace,tau,J,fscale,F;tol=1e-7)
                        
                        if stats.m_final >= mmax+1
                            # reject step and reduce time step
                            println("KIOPS subspace dimension is large $(stats.m_final), reducing time step by 0.5")
                            dt_local = ldexp(dt_local, -1)
                            continue
                        end

                        println("KIOPS stats: ", stats)

                        println("sum w: ", sum(D .* KIOPS_workspace.w), " norm w: ", norm(D .* KIOPS_workspace.w), " max w: ", maximum(D .* KIOPS_workspace.w), " min w: ", minimum(D .* KIOPS_workspace.w), "dtscale:", dtscale)

                        fout .= fold .+ (D) .* KIOPS_workspace.w

                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                    # direct matrix exponential for testing
                        #=Atilde = KIOPS_workspace.Atilde
                        # form augmented matrix Atilde = [A u_phi; 0 0] for Rosenbrock-Euler
                        @views copyto!(Atilde[1:n_momentum, 1:n_momentum], J)
                        @views copyto!(Atilde[1:n_momentum, end], g)
                        expm_kiops_pade!(KIOPS_workspace.F, Atilde, KIOPS_workspace)
                        mul!(fout,@view(KIOPS_workspace.F[1:n_momentum, 1:n_momentum]),fscale)
                        fout .+= @view(KIOPS_workspace.F[1:n_momentum,end])
                        fout .= (D) .* fout # right mul
                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)=#

                    # time step adaption based on energy error 
                        # energy error
                        ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * dtscale)
                        Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout), " inj E: ", dot(E, df_Inj * dtscale))
                        
                        if isnan(ηE) || ηE > 1e-4
                            # reject step and reduce time step
                            println("Rejecting first step due to large energy error: ", ηE, " reducing time step by 0.5")
                            dt_local = ldexp(dt_local, -1)
                            continue
                        end

                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    @. fold = fout

                    if ηE < 1e-5 #= (ηE < 5e-5 && stats.m_final <= mmax / 2) || stats.m_final == 10=#
                        #println("KIOPS good $(stats.m_final), increasing time step by 2")
                        dt_next = ldexp(dt_old, 1)
                        dt_local = ldexp(dt_local, 1)
                    else
                        dt_next = dt_old
                    end

            end # while not accurate

            @. fstep = fold

            @inbounds method.dt_guess[off_space+1] = dt_next

        else #

            if isassigned(method.M_Emi, off_space+1)
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt_guess #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffnes
            dt_next = dt_local

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = Precision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive
                dt_old = dt_local

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    if EmiTrue
                        copyto!(Jsparse,M_Emi)
                        Jsparse *= dtscale * invA
                    else
                        fill!(Jsparse,zero(Precision))
                        dropzeros!(Jsparse)
                    end

                    # Form F
                    mul!(F,Jsparse,fold)
                    @. F += #=A *=# df_Inj * dtscale
                    if all(iszero, F) # no change in distribution
                        break
                    end

                    # get residuals 
                    mul!(g,J,fold)
                    LinearAlgebra.axpby!(Precision(1),F,Precision(-1),g)

                    Jsparse .*= D' # right mul
                    Jsparse .*= Dinv # left mul
                    g .*= Dinv # left mul

                    # form U
                    fscale .= (Dinv) .* fold

                    tau = one(Precision)
                    #_, stats = kiops_roe!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)
                    #_, stats = kiops_roe_noview!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)
                    #_, stats = kiops_roe_panel!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)
                    _, stats = kiops_roe_panel_new!(KIOPS_workspace,tau,J,fscale,g;tol=1e-7)

                    if stats.m_final >= 17 #mmax+1
                        # reject step and reduce time step
                        println("KIOPS subspace dimension is large $(stats.m_final), reducing time step by 0.5")
                        dt_local = ldexp(dt_local, -1)
                        continue
                    end

                    println("KIOPS stats: ", stats)

                    fout .= (D) .* KIOPS_workspace.w

                    # adaptive time stepping

                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * dtscale)
                        Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout), " inj E: ", dot(E, df_Inj * dtscale))
                        
                        if isnan(ηE) || ηE > 1e-4
                            # reject step and reduce time step
                            println("Rejecting first step due to large energy error: ", ηE, " reducing time step by 0.5, $kE")
                            dt_local = ldexp(dt_local, -1)
                            continue
                        end

                       if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

                    if ηE < 1e-5 #= (ηE < 5e-5 && stats.m_final <= mmax / 2) || stats.m_final == 10=#
                        #println("KIOPS good $(stats.m_final), increasing time step by 2")
                        dt_next = ldexp(dt_old, 1)
                    else
                        dt_next = dt_old
                    end

            end # while not accurate

            @. fstep = fold

            @inbounds method.dt_guess[off_space+1] = dt_next

        end

    end

    return nothing

end

function (method::ExponentialRosenbrockEulerLejaStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = zeros(eltype(method.f),length(method.f))

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end 

    # binary update

        update_momentum!(method,dt_scale)

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        test = findmin(method.df_tmp)
        println("loc: ", test[2], " Minimum df_tmp: ", test[1], " min df: ", method.df[test[2]], " min f: ", method.f[test[2]], " min df_tmp old: ", (@. ifelse(isnan(method.df / method.f), Inf, method.df / method.f))[test[2]])
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end

function update_momentum!(method::ExponentialRosenbrockEulerLejaStruct,dt::T) where T

    Precision = method.Precision
    
    n_momentum = method.PhaseSpace.Grids.n_momentum
    #n_space = method.PhaseSpace.Grids.n_space
    momentum_species_offset = method.PhaseSpace.Grids.momentum_species_offset
    name_list = method.PhaseSpace.name_list
    num_species = length(name_list)
    
    fold = method.fold
    fout = method.fout
    Estate::Vector{Precision} = zeros(Precision,n_momentum)
    J = method.J
    #Jtmp = copy(J)
    F = method.F
    fscale = method.fscale #::Vector{Precision} = zeros(Precision,n_momentum)
    J64 = zeros(Float64,n_momentum,n_momentum)
    F64 = zeros(Float64,n_momentum)
    #Ftmp = copy(F)
    ϕ = method.ϕ #zeros(Precision,n_momentum,2)  
    D = method.D #Diagonal(ones(Precision,n_momentum)) # 1/E
    Dinv = method.Dinv #Diagonal(ones(Precision,n_momentum)) # E
    E = method.E

    #D.diag .= one(Precision)
    #Dinv.diag .= one(Precision)

    δ = method.δ #::Vector{Precision} = zeros(Precision,n_momentum)

    m = method.m # dimension of Krylov subspace for exponential action approximation, can adjust based on problem size and desired accuracy

    maxdeg = method.m # maximum degree of polynomial interpolation for exponential action approximation, can adjust based on problem size and desired accuracy
    Leja_nodes = method.Leja_nodes # precomputed Leja nodes for polynomial interpolation, can adjust based on problem size and desired accuracy

    Ks = method.Ks
    ϕcache = method.ϕcache

    Ks64 = KrylovSubspace{Float64}(n_momentum, m)

    for off_space in method.ActiveDomain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        df_Inj = @view(method.df_Inj[start_idx:end_idx])

        has_injection = sum(df_Inj) > zero(Precision)

        if !has_injection && sum(fstep) == zero(Precision) 
            continue
        end

        @inbounds vol = method.Vol[off_space+1]
        invA_Flux = @view(method.invA_Flux[start_idx:end_idx,start_idx:end_idx])
        A_Flux = @view(method.A_Flux[start_idx:end_idx,start_idx:end_idx])
        invA = invA_Flux[1,1] # assumes the diagonal is just constant, should be true for a single spacial cell
        A = A_Flux[1,1]

        if method.Binary_Interactions && off_space in method.Bin_Domain

            if isassigned(method.M_Emi, off_space+1)
                EmiTrue = true
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision))
                    # Form J
                    @. J = Precision(2) * method.M_Bin_Mul_Step 
                    @. J += M_Emi
                    lmul!(dtscale,J)
                    @. J *= invA # as we are evolving Af
                    # Form F
                    @. method.M_Bin_Mul_Step += M_Emi
                    mul!(F,method.M_Bin_Mul_Step,fold)
                    @. F *= invA

                    # Energy error estimate from F
                    ηEest = abs(dot(E,F)) * dtscale / dot(E,fold)
                    println("Energy error estimate from F: ", ηEest)

                    @. F += df_Inj 
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end
                    lmul!(dtscale,F)

                    # Diagonal scaling of J and F to improve conditioning for Krylov subspace approximation
                    rmul!(J,D)
                    lmul!(Dinv,J)
                    lmul!(Dinv,F)

                    ϕ1Av, info = phi1_leja_dense(J,F,Leja_nodes;dt = k,gamma = 1.5,maxdeg = maxdeg,tol = 1e-6,use_directional_mu = true,use_dense_mu2 = false)

                    println("ϕ1Av info: ", info)

                    # Compute φ functions of H
                    #phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true,errest=false) # TODO: This allocates
                    mul!(δ,D,ϕ1Av,k#=*invA=#,zero(Precision))
                    #@time δ .= invA_Flux * δ # TODO: this allocates
                    @. fout = fold + δ
                    
                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * k*dtscale)
                        Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout), " inj E: ", dot(E, df_Inj * k * dtscale))

                        if !isfinite(ηE)
                            println("ηE: $ηE, sum(δ): $(sum(δ))")
                            @warn "Energy error is Inf or NaN, may be unstable, consider reducing time step or adjusting ηtarget"
                            dt_local *= convert(typeof(dt_local), 1/64)
                            continue
                        end

                        if info.err_est > 1e-6 || ηE > 1e-5
                            @warn "ϕ1Av error estimate $(info.err_est) is large, may be unstable, consider reducing time step or adjusting ηtarget"
                            dt_local = ldexp(dt_local, -1)
                            continue
                        elseif info.err_est < 1e-7 && ηE < 1e-6
                            @warn "ϕ1Av error estimate $(info.err_est) is small, may be inefficient, consider increasing #time step or adjusting ηtarget"
                            dt_local = ldexp(dt_local, 1)
                            continue
                        end

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        ηtarget = 1e-16

                        #kE = (1e-4/(ηE+eps(1e-4)))^(1.0/(order+1)) # k from energy error estimate 
                        #kϕmax = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                        #kϕ = min(kE,kϕmax) # max k is 2.0
                        #errest = Inf
                        #while errest > ηtarget # k from ϕv errestimate, can be more strict than energy error estimate
                        #    _, errest = phiv!(ϕ,kϕ,Ks64,1;cache=ϕcache,correct=true,errest=true) # TODO: This allocates
                        #    println("ϕv error estimate during: ", errest)
                        #    if errest > ηtarget
                        #        kϕ *= 0.77 * 0.5(tanh(log10(errest)-log10(ηtarget)+1)+1)
                        #    end
                        #end
                        #println("ϕv error estimate after: ", errest)
                        #k = min(kE,2.0)
                        #println("kE: ", kE, " kϕ: ", kϕ, " errest: ", errest)

                        #k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        dt_local = k*dt_old
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old


                    #=if k != 1.0

                        println("here")

                        ϕ1Av, info = phi1_leja_dense(J,F,Leja_nodes;dt = k,gamma = 1.5,maxdeg = maxdeg,tol = 1e-6,use_directional_mu = false,use_dense_mu2 = false)
                        mul!(δ,D,ϕ1Av,k#=*invA=#,zero(Precision))

                        #phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true,errest=false) # TODO: This allocates
                        #mul!(δ,D,@view(ϕ[:,2]),k#=*invA=#,zero(Precision))
                        @. fout = fold + δ

                        # energy after correction  
                        ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * k * dtscale)
                        Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error after adaption: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))
                    end=#

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)


                    #=if ηE > 2e-4
                        @warn "Energy error in EXPRB1 step is large $ηE, may be unstable, consider reducing time step or adjusting ηtarget"
                    end=#

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        else #

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    @. J = M_Emi
                    lmul!(dtscale*invA,J) # scale as evolving Af
                    # Form F
                    mul!(F,M_Emi,fold)
                    @. F *= invA 
                    @. F += df_Inj 
                    lmul!(dtscale,F)
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end

                    rmul!(J,D)
                    lmul!(Dinv,J)
                    lmul!(Dinv,F)

                    @. J64 = Float64(J)

                    if has_injection 
                        @. F64 = Float64(F)
                        arnoldi!(Ks64,J64,F64;m=m,reorthogonalize=true)
                        #arnoldi!(Ks,J,F;m=m,reorthogonalize=true)
                    else # no injection, just linear Jacobian so use exp over phi
                        f64 = Float64.(Dinv * fold)
                        arnoldi!(Ks64,J64,f64;m=m,reorthogonalize=true)
                        #mul!(fscale,Dinv,fold)
                        #arnoldi!(Ks,J,fscale;m=m,reorthogonalize=true)
                    end

                    #V = ExponentialUtilities.getV(Ks)[:,1:end-1]
                    #H = ExponentialUtilities.getH(Ks)[1:end-1,1:end]
                    
                    V = ExponentialUtilities.getV(Ks64)[:,1:end-1]
                    H = ExponentialUtilities.getH(Ks64)[1:end-1,1:end]

                    norm_unscaled = norm(V' * V - I)
                    if isnan(norm_unscaled) # no J or F means no Krylov subspace generated, so just accept step as is and move on,
                        break
                    end

                    if #=cond((I - dt_local*H)) > 1f3 ||=#  norm(V' * V - I) > 1f-5 
                        dt_local = dt_local * 0.5
                        @warn "Krylov subspace has poor orthogonalisation $(norm(V' * V - I)), reducing time step to $dt_local"
                        continue
                    end

                    # Compute φ functions of H
                    if has_injection 
                        phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true,errest=false) # TODO: This allocates
                        mul!(δ,D,@view(ϕ[:,2]),k#=*invA=#,zero(Precision))
                        @. fout = fold + δ
                    else # no injection, just linear Jacobian so use exp over phi
                        expv!(δ,k#=*invA=#,Ks64)
                        @. fout = D.diag * δ
                    end


                    # adaptive time stepping

                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        if has_injection
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * dtscale)
                            Eold = dot(E, fold .+ df_Inj * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        #println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        if isinf(ηE)
                            @warn "Energy error is Inf or NaN, may be unstable, consider reducing time step or adjusting ηtarget"
                            dt_local *= convert(typeof(dt_local), 1/64)
                            continue
                        end

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        ηtarget = 1e-5

                        kE = (1e-4/(ηE+eps(1e-4)))^(1.0/(order+1)) # k from energy error estimate
                        if has_injection
                            kϕmax = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                            kϕ = min(kE,kϕmax) # max k is 2.0
                            errest = Inf
                            while errest > ηtarget # k from ϕv errestimate, can be more strict than energy error estimate
                                _, errest = phiv!(ϕ,kϕ,Ks64,1;cache=ϕcache,correct=true,errest=true) # TODO: This allocates
                                println("ϕv error estimate during: ", errest)
                                if errest > ηtarget
                                    kϕ *= 0.77 * 0.5(tanh(log10(errest)-log10(ηtarget)+1)+1)
                                end
                            end
                            println("ϕv error estimate after: ", errest)
                            k = min(kE,kϕ,2.0)
                            println("kE: ", kE, " kϕ: ", kϕ, " errest: ", errest)
                        else
                            k = min(kE,2.0)
                        end
                        #k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        dt_local = k*dt_old
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                    if k != 1.0
                        if has_injection 
                            phiv!(ϕ,k,Ks64,1;cache=ϕcache,correct=true)
                            mul!(δ,D,@view(ϕ[:,2]),k#=*invA=#,zero(Precision))
                            @. fout = fold + δ
                        else # no injection, just linear Jacobian so use exp over phi
                            expv!(δ,k#=*invA=#,Ks64)
                            @. fout = D.diag * δ
                        end

                        # energy after correction  
                        #=if has_injection 
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * k * dtscale)
                            Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))=#
                        #println("Energy error: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        #=if ηE > 2e-4
                            @warn "Energy error in EXPRB1 step is large $ηE, may be unstable, consider reducing time step or adjusting ηtarget"
                        end=#
                        
                    end

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        end

    end

    return nothing

end


has_positive_offdiag(M) = any((i != j && M[i, j] > 0) for i in axes(M, 1), j in axes(M, 2))

function patankar_energy_fix!(method::SymplecticMPEStruct,fstep::AbstractVector{T},fprev::AbstractVector{T}) where T

    E = method.E
    df = similar(fstep)
    fill!(df,zero(eltype(df)))
    ftmp = copy(fstep)

    w = zeros(eltype(df), length(df))
    df_active = zeros(eltype(df), length(df))
    df_inactive = zeros(eltype(df), length(df))

    active = trues(length(fstep))
    nonpositive = true

    while nonpositive

        fill!(w, zero(eltype(w)))
        fill!(df_active, zero(eltype(df_active)))
        fill!(df_inactive, zero(eltype(df_inactive)))

    # fixed inactive correction
    df_inactive .= -ftmp .* (.!active)

    # remaining defects caused by clipping
    ΔM = sum(ftmp[.!active])

    ΔE = sum(E .* fprev) - sum(E .* (ftmp .+ df_inactive))

    # weights on active set
    eta = E
    w[active] .= ftmp[active] #.* eta

    S0 = sum(w)
    Ebar = sum(E .* w) / S0

    Ẽ = E .- Ebar

    S2 = sum((Ẽ.^2) .* w)

    λ = -2ΔM / S0
    μ = -2*(ΔE - Ebar*ΔM) / S2
    # active correction ONLY
    df_active[active] .= -0.5 .* w[active] .* (λ .+ μ .* Ẽ[active])

    # total correction
    df = df_inactive + df_active

    println(sum(df_active))
    println("mass defect = ", ΔM)
    println("energy defect = ", ΔE)
    println("constraint mass = ", sum(df_active))
    println("constraint energy = ", sum(E .* df_active))
    println("minimum corrected state = ", minimum(ftmp .+ df))

    ucorr = ftmp .+ df

    new_active = ucorr .> 1e-10

    active .&= new_active

    nonpositive = any(ucorr .< 0.0)

    println("sum(df_inactive) = ", sum(df_inactive))
    println("sum(df_active)   = ", sum(df_active))
    println("total correction = ", sum(df))

    end

    println("mass defect = ", sum(ftmp .+ df) - sum(fprev))
    println("energy defect = ", sum(E .* (ftmp .+ df)) - sum(E .* fprev))


    return ftmp .+ df

end

function patankar_energy_fix_linear!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,) where T

    E = method.E

    #
    # Uncorrected Patankar state
    #
    ftmp = copy(fstep)

    N = length(ftmp)

    #
    # Corrections
    #
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    #
    # Active-set mask
    #
    active = ftmp .> tol

    #
    # Weights
    #
    w = zeros(T, N)

    for iter in 1:maxiter

        old_active = copy(active)

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))

        #
        # Inactive states are clipped to zero
        #
        df_inactive .= -ftmp .* (.!active)

        #
        # Remaining conservation defects
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Active weights
        #
        w[active] .= ftmp[active]

        #
        # Weighted energy centroid
        #
        W0 = sum(w)

        if W0 <= tol
            error("No active states remain in energy correction.")
        end

        Ebar = sum(E .* w) / W0

        #
        # Centered energies for conditioning
        #
        Ẽ = E .- Ebar

        #
        # Weighted moments
        #
        S0 = sum(w)
        S2 = sum((Ẽ.^2) .* w)

        #
        # Solve:
        #
        # δu_i = w_i (a + b Ẽ_i)
        #
        # Constraints:
        #
        # Σ δu_i      = ΔM
        # Σ E_i δu_i = ΔE
        #
        # Because Ẽ is centered:
        #
        # Σ w_i Ẽ_i = 0
        #
        # giving a diagonalized system.
        #

        a = ΔM / S0

        #
        # Energy defect relative to centroid
        #
        ΔẼ = ΔE - Ebar * ΔM

        b = ΔẼ / S2

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ)

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        #
        # Corrected state
        #
        ucorr = ftmp .+ df

        #
        # Update active set
        #
        active .&= (ucorr .> 0.0)

        #
        # Convergence test
        #
        if all(active .== old_active)
            break
        end

        if iter == maxiter
            @warn "Linear active-set energy projection reached maxiter."
        end
    end

    fcorrected = ftmp .+ df

    #
    # Diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
end

function patankar_energy_fix_linear_filtered!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,
    energy_window_orders = 4,) where T

    E = method.E

    #
    # Uncorrected Patankar state
    #
    ftmp = copy(fstep)

    N = length(ftmp)

    #
    # Correction vectors
    #
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    #
    # Working arrays
    #
    w      = zeros(T, N)
    filter = zeros(T, N)

    #
    # Permanent inactive mask
    #
    inactive = falses(N)

    #
    # Initial positivity check
    #
    inactive .|= (ftmp .<= tol)

    for iter in 1:maxiter

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))
        fill!(filter,      zero(T))

        #
        # Active set is complement of inactive set
        #
        active = .!inactive

        #
        # Permanently clipped states
        #
        df_inactive .= -ftmp .* inactive

        #
        # Remaining defects AFTER clipping
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Determine occupied energy range
        #
        occupied =
            (ftmp .+ df_inactive) .> tol

        if !any(occupied)
            error("No occupied states remain.")
        end

        Emax_occ = maximum(E[occupied])

        #
        # Smooth logarithmic filter
        #
        for i in eachindex(E)

            if E[i] <= 0 || Emax_occ <= 0
                filter[i] = zero(T)
            else

                x =
                    log10(E[i] / Emax_occ) +
                    energy_window_orders

                filter[i] =
                    clamp(x, zero(T), one(T))
            end
        end

        #
        # Active correction weights
        #
        w[active] .=
            ftmp[active] .* filter[active]

        #
        # Prevent singular systems
        #
        S0 = sum(w)

        if S0 <= tol
            error("Correction weights vanished.")
        end

        #
        # Center energies
        #
        Ebar = sum(E .* w) / S0

        Ẽ = E .- Ebar

        #
        # Weighted variance
        #
        S2 = sum((Ẽ.^2) .* w)

        if S2 <= tol
            error("Energy variance too small.")
        end

        #
        # Linear correction:
        #
        # δu_i = w_i (a + b Ẽ_i)
        #
        a = ΔM / S0

        ΔẼ = ΔE - Ebar * ΔM

        b = ΔẼ / S2

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ)

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        #
        # Corrected state
        #
        ucorr = ftmp .+ df

        #
        # Detect NEW negative states
        #
        newly_inactive =
            (ucorr .<= tol) .& (.!inactive)

        #
        # Diagnostics
        #
        println("iter = ", iter)

        println("ΔM target         = ", ΔM)
        println("ΔE target         = ", ΔE)

        println("constraint mass   = ", sum(df))
        println("constraint energy = ", sum(E .* df))

        println("minimum corrected = ", minimum(ucorr))

        #
        # Converged?
        #
        if !any(newly_inactive)
            break
        end

        #
        # Permanently accumulate clipped states
        #
        inactive .|= newly_inactive

        if iter == maxiter
            @warn "Energy projection reached maxiter."
        end
    end

    #
    # Final corrected state
    #
    fcorrected = ftmp .+ df

    #
    # Final diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
end

function patankar_energy_fix_quadratic!(
    method,
    fstep::AbstractVector{T},
    fprev::AbstractVector{T};
    tol = 1e-14,
    maxiter = 50,) where T

    E = method.E

    # unconstrained Patankar result
    ftmp = copy(fstep)

    N = length(fstep)

    # correction vectors
    df          = zeros(T, N)
    df_active   = zeros(T, N)
    df_inactive = zeros(T, N)

    # active set
    active = ftmp .> tol

    # weights
    w = zeros(T, N)

    for iter in 1:maxiter

        old_active = copy(active)

        fill!(df,          zero(T))
        fill!(df_active,   zero(T))
        fill!(df_inactive, zero(T))
        fill!(w,           zero(T))

        #
        # Inactive states are clipped to zero
        #
        df_inactive .= -ftmp .* (.!active)

        #
        # Remaining conservation defects
        #
        ΔM =
            sum(fprev) -
            sum(ftmp .+ df_inactive)

        ΔE =
            sum(E .* fprev) -
            sum(E .* (ftmp .+ df_inactive))

        #
        # Active-set weights
        #
        w[active] .= ftmp[active]

        #
        # Center energies for conditioning
        #
        W0 = sum(w)

        if W0 <= tol
            error("No active states remain in correction.")
        end

        Ebar = sum(E .* w) / W0

        Ẽ = E .- Ebar

        #
        # Build weighted moments
        #
        M0 = sum(w)
        M1 = sum(Ẽ .* w)
        M2 = sum((Ẽ.^2) .* w)
        M3 = sum((Ẽ.^3) .* w)
        M4 = sum((Ẽ.^4) .* w)

        #
        # We seek:
        #
        # δu_i = w_i (a + b Ẽ_i + c Ẽ_i^2)
        #
        # subject to:
        #
        # Σ δu_i = ΔM
        # Σ Ẽ_i δu_i = ΔẼ
        #
        # and minimizing weighted norm.
        #
        # This gives the normal equations:
        #
        # [M0 M1 M2] [a]   [ΔM ]
        # [M1 M2 M3] [b] = [ΔẼ]
        # [M2 M3 M4] [c]   [0   ]
        #
        # Last row enforces minimum-norm quadratic correction.
        #

        ΔẼ = ΔE - Ebar * ΔM

        A = [
            M0 M1 M2
            M1 M2 M3
            M2 M3 M4
        ]

        rhs = T[ΔM, ΔẼ, 0]

        coeffs = A \ rhs

        a, b, c = coeffs

        #
        # Active correction
        #
        df_active .=
            w .* (a .+ b .* Ẽ .+ c .* (Ẽ.^2))

        #
        # Total correction
        #
        df .= df_inactive .+ df_active

        ucorr = ftmp .+ df

        #
        # Update active set
        #
        active .&= (ucorr .> 0.0)

        #
        # Converged?
        #
        if all(active .== old_active)
            break
        end

        if iter == maxiter
            @warn "Active-set projection reached maxiter."
        end
    end

    #
    # Final corrected state
    #
    fcorrected = ftmp .+ df

    #
    # Diagnostics
    #
    println("Final mass defect   = ",
        sum(fcorrected) - sum(fprev))

    println("Final energy defect = ",
        sum(E .* fcorrected) - sum(E .* fprev))

    println("Minimum state       = ",
        minimum(fcorrected))

    return fcorrected
end

function pow2_equilibration_matrices(A::AbstractMatrix{T}; iters=3) where {T<:AbstractFloat}
    n, m = size(A)

    dl = ones(T, n)
    dr = ones(T, m)

    As = copy(A)

    for k in 1:iters
        # Left scaling: rows
        for i in 1:n
            s = maximum(abs, @view As[i, :])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[i, :] .*= α
                dl[i] *= α
            end
        end

        # Right scaling: columns
        for j in 1:m
            s = maximum(abs, @view As[:, j])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[:, j] .*= α
                dr[j] *= α
            end
        end
    end

    return Diagonal(dl), Diagonal(dr)
end

function pow2_equilibration_matrices!(Dl::Diagonal{T,Vector{T}}, Dr::Diagonal{T,Vector{T}}, A::AbstractMatrix{T}; iters=3) where {T<:AbstractFloat}

    n, m = size(A)

    Dl.diag .= ones(T)
    Dr.diag .= ones(T)

    As = copy(A)

    for k in 1:iters
        # Left scaling: rows
        for i in 1:n
            s = maximum(abs, @view As[i, :])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[i, :] .*= α
                Dl.diag[i] *= α
            end
        end

        # Right scaling: columns
        for j in 1:m
            s = maximum(abs, @view As[:, j])
            if s > zero(T)
                α = T(2.0)^(-round(Int, log2(s)))
                As[:, j] .*= α
                Dr.diag[j] *= α
            end
        end
    end

    return nothing
end

function norm_equilibration_matrices(A::AbstractMatrix{T}; iters=5) where {T<:AbstractFloat}
    n, m = size(A)

    dl = ones(T, n)
    dr = ones(T, m)

    As = copy(A)

    for k in 1:iters
        # Row norm scaling
        for i in 1:n
            s = norm(@view As[i, :])
            if s > zero(T)
                α = one(T) / T(s)
                As[i, :] .*= α
                dl[i] *= α
            end
        end

        # Column norm scaling
        for j in 1:m
            s = norm(@view As[:, j])
            if s > zero(T)
                α = one(T) / T(s)
                As[:, j] .*= α
                dr[j] *= α
            end
        end
    end

    return Diagonal(dl), Diagonal(dr)
end

function norm_equilibration_matrices!(Dl::Diagonal{T,Vector{T}}, Dr::Diagonal{T,Vector{T}}, A::AbstractMatrix{T}; iters=5) where {T<:AbstractFloat}
    n, m = size(A)

    Dl.diag .= ones(T)
    Dr.diag .= ones(T)
    As = copy(A)

    for k in 1:iters
        # Row norm scaling
        for i in 1:n
            s = norm(@view As[i, :])
            if s > zero(T)
                α = one(T) / T(s)
                As[i, :] .*= α
                Dl.diag[i] *= α
            end
        end

        # Column norm scaling
        for j in 1:m
            s = norm(@view As[:, j])
            if s > zero(T)
                α = one(T) / T(s)
                As[:, j] .*= α
                Dr.diag[j] *= α
            end
        end
    end

    return nothing

end

function cutoff_offdiagonal_to_diag_relative!(A, scale)
    n, m = size(A)
    @assert n == m "A must be square"

    for j in 1:n
        colmax = maximum(A[i, j] for i in 1:n)
        cutoff = scale *eps(colmax)

        for i in 1:n
            if i != j && abs(A[i, j]) < cutoff
                removed = A[i, j]
                A[i, j] = zero(eltype(A))
                A[j, j] -= removed  # Add the value to diagonal
            end
        end
    end

    return nothing
end

"""
    threshold_Jrows!(J; rtol=1e-14)

Removes weak coupling in the rows of the jacobian matrix J, i.e. elements of the equation for each state variable that are small compared to the largest element of that equation.
"""
function threshold_Jrows!(J; rtol=1e-14)
    m, n = size(J)

    for i in 1:m
        rowmax = maximum(abs, @view J[i,:])

        if rowmax == 0
            continue
        end

        thresh = rtol * rowmax

        for j in 1:n
            if i!=j && sign(J[i,j]) == 1 && abs(J[i,j]) < thresh # doesn't remove diagonal entries 
                J[i,j] = zero(eltype(J))
            end
        end
    end

    return nothing
end


"""
    threshold_Jcols!(J; rtol=1e-14)

Removes weak coupling in the columns of the jacobian matrix J, i.e. elements of the equation for each state variable that are small compared to the largest element of that equation.
"""
function threshold_Jcols!(J; rtol=1e-14)

    m, n = size(J)

    for j in 1:n

        colmax = maximum(abs, @view J[:,j])

        if colmax == 0
            continue
        end

        thresh = rtol * colmax

        for i in 1:m
            if i != j #=&& sign(J[i,j]) == 1=# && abs(J[i,j]) < thresh # postive and small
                J[i,j] = zero(eltype(J))
            end
        end
    end

    return J
end

"""
    threshold_F!(F; rtol=1e-14)

Removes weak transition rates in F, i.e. terms df_i/dt = F_i that are small and by small I mean f^{t+1}=f^{t}+ΔtF so if ΔtF < eps(f^{t})*rtol.
"""
function threshold_F!(F,f,dt; rtol=1e-2)
    
    @. F = ifelse(dt*abs(F) < f*rtol, zero(eltype(F)), F)
    #@. F = ifelse(abs(F) < rtol, zero(eltype(F)), F)

    return nothing
end

"""
    remove_small_J_contributions!(J,f; rtol=1e-14)

Removes weak coupling in the rows of the Jacobian matrix J when scaled by the state vector f
"""
function remove_small_J_contributions!(J,Jtmp,f; tol=1e-14)

    m, n = size(J)

    Jtmp .= abs.(J) .* reshape(abs.(f), 1, :)   # each entry's contribution magnitude
    rowmax  = vec(maximum(Jtmp; dims = 2))      # max contribution per row
    mask    = Jtmp .>= tol .* rowmax

    @. J *= mask

    return J
end

"""
    balance_similarity(J; maxiter=50, tol=1e-3)

Compute a diagonal scaling D such that

    Js = D \\ J * D

has approximately balanced row and column 2-norms.

Returns (D, Js).
"""
function balance_similarity(J; maxiter=50, tol=1e-3)

    n = size(J,1)
    @assert size(J,2) == n

    d = ones(eltype(J), n)

    Js = copy(J)

    for k = 1:maxiter

        maxerr = 0.0

        for i = 1:n

            r = norm(view(Js, i, :))
            c = norm(view(Js, :, i))

            if r == 0 || c == 0
                continue
            end

            ratio = c / r

            ratio = clamp(ratio, 1e-16, 1e16) # prevent extreme scaling

            # desired scaling
            s = sqrt(ratio)

            d[i] *= s
            d[i] = clamp(d[i], 1e-8, 1e8) # prevent extreme scaling

            # similarity update
            Js[:,i] .*= s
            Js[i,:] ./= s

            maxerr = max(maxerr, abs(log(ratio)))
        end

        maxerr < tol && break
    end

    return Diagonal(d), Js
end


"""
    ERBEKrylov(dg,g,t,dt)

2nd order Expontial integration time-stepping method for the transport equation. 

"""
function (method::ERBEKrylovStruct)(t_start,t_stop,dt,Verbose::Int64)

    method.step += 1

    dt0 = method.dt0

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

    if dt < 0.0
        error("Negative time step calculated, something went wrong with the CFL condition calculation")
    end

    # scaling of time stepping

    dt_scale = method.Precision(dt / dt0)

    ftmp = similar(method.f)

    # set fstep to intial value

        @. method.fstep = method.f

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        # mask of df regions

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df

    # Half injection 

        @. method.fstep += method.df_Inj * dt_scale / 2

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end 

    # binary update

        update_momentum!(method,dt_scale)

    # half momentum update

        mul!(ftmp,method.invImMP,method.fstep) 

        if !isnothing(method.df_mask)
            @. method.fstep = ftmp * method.df_mask + method.fstep * (1-method.df_mask)
        else
            @. method.fstep = ftmp
        end

    # half injection 

        @. method.fstep += method.df_Inj * dt_scale / 2

    # half space update

        mul!(method.df,method.X_Flux,method.fstep)
        method.df .= - method.invA_Flux * method.df .* dt_scale/2  # minus sign as flux terms are on RHS of transport equation, also resets df_Space

        if !isnothing(method.df_mask)
            @. method.df *= method.df_mask
        end 

        @. method.fstep += method.df
    
    # removing negative values

        @. method.fstep = ifelse(method.fstep<=method.n_cut,zero(eltype(method.fstep)),method.fstep)


    Cr = 0.0
    sum_f = sum(method.f)
    if sum_f != 0.0

        #@. method.df_tmp = ifelse(method.f + method.df < method.n_cut, zero(eltype(method.df)), method.df)
        #@. method.df_tmp = method.df_tmp / method.f
        @. method.df = method.fstep-method.f
        @. method.df_tmp = method.df / method.f
        @. method.df_tmp = ifelse(isnan(method.df_tmp), Inf, method.df_tmp)
        Cr = -minimum(method.df_tmp) 

    end

    if Verbose == 1 && Cr > 1.0
        println("step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3)), dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3)) system may be unstable")
    elseif Verbose == 2
        println("\r step=$(method.step), t=$(round(t_start,sigdigits=4)), Cr = $(round(Cr,sigdigits=3))")
    elseif Verbose == 3
        println("step=$(method.step), Cr = $(round(Cr,sigdigits=3)),Cr_Bin = $(round(Cr_Bin,sigdigits=3)), Cr_Emi = $(round(Cr_Emi,sigdigits=3)), Cr_Flux = $(round(Cr_Flux,sigdigits=3)), t=$t_start, t_save =$t_stop, dt_attempted=$(round(dt_old,sigdigits=3)), dt_adapted = $(round(dt,sigdigits=3))")
    end
    if Verbose > 0
        flush(stdout)
    end

    method.f .= method.fstep

    # remove masked off domain regions

    if !isnothing(method.f_mask)
        @. method.f *= method.f_mask
    end

    return dt,save

end


function update_momentum!(method::ERBEKrylovStruct,dt::T) where T

    Precision = method.Precision
    
    n_momentum = method.PhaseSpace.Grids.n_momentum
    #n_space = method.PhaseSpace.Grids.n_space
    momentum_species_offset = method.PhaseSpace.Grids.momentum_species_offset
    name_list = method.PhaseSpace.name_list
    num_species = length(name_list)
    
    fold = method.fold
    fout = method.fout
    J = method.J
    Jsparse = sparse(copy(J))
    #Jtmp = copy(J)
    F = method.F
    fscale = method.fscale #::Vector{Precision} = zeros(Precision,n_momentum)

    D = method.D #Diagonal(ones(Precision,n_momentum)) # 1/E
    Dinv = method.Dinv #Diagonal(ones(Precision,n_momentum)) # E
    E = method.E

    δ = method.δ #::Vector{Precision} = zeros(Precision,n_momentum)

    m = method.m # dimension of Krylov subspace for exponential action approximation, can adjust based on problem size and desired accuracy


    EmiTrue::Bool = true

    for off_space in method.ActiveDomain

        start_idx = n_momentum*off_space+1
        end_idx = n_momentum*(off_space+1)

        fstep = @view(method.fstep[start_idx:end_idx])
        df_Inj = @view(method.df_Inj[start_idx:end_idx])

        has_injection = sum(df_Inj) > zero(Precision)

        if !has_injection && sum(fstep) == zero(Precision) 
            continue
        end

        @inbounds vol = method.Vol[off_space+1]
        @inbounds invA = method.invA[off_space+1]
        @inbounds dt_guess = method.dt_guess[off_space+1]

        if method.Binary_Interactions && off_space in method.Bin_Domain

            if isassigned(method.M_Emi, off_space+1)
                EmiTrue = true
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0 # t is ratio of current t to dt. i.e. t=1.0 is end of timestep

            dt_local = dt_guess # initial guess for local time step, based on previous time steps
            dt_old = dt_local
            dt_next = dt_local

            # k's scale dt_local
            kE = 1.0
            kH = 1.0
            k = 1.0

            while t < 1.0

                dtscale = dt_local / method.dt0 # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    mul!(method.M_Bin_Mul_Step_reshape,method.M_Bin,fold,vol,zero(Precision))
                    # Form J
                    @. J = Precision(2) * method.M_Bin_Mul_Step
                    #fill!(J,zero(Precision))
                    #fill!(method.M_Bin_Mul_Step,zero(Precision)) 
                    if EmiTrue
                        @. J += M_Emi
                        @. method.M_Bin_Mul_Step += M_Emi
                    end

                    @. J *= dtscale * invA

                    println("max M_Emi: ", maximum(M_Emi), " min M_Emi : ", minimum(M_Emi))
                    println("max J: ", maximum(J), " min J : ", minimum(J))
                    println("max fold: ", maximum(fold), " min fold : ", minimum(fold))
                    println("max df_Inj: ", maximum(df_Inj), " min df_Inj : ", minimum(df_Inj))

                    # Form F
                    mul!(F,method.M_Bin_Mul_Step,fold)
                    @. F *= invA

                    # Energy error estimate using F 
                    ηEest = abs(dot(E,F) * dtscale) / dot(E,fold)
                    println("Energy error estimate: ", ηEest)

                    #@. F += #=A *=# df_Inj
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end
                    lmul!(dtscale,F)

                    # Diagonal scaling of J and F to improve conditioning for Krylov subspace approximation
                    J .*= (D)' # right mul
                    J .*= (Dinv) # left mul
                    F .*= (Dinv) # left mul

                    # TODO: define J and F as views of Aaug to avoid copying data
                    #Aaug = method.Aaug
                    #fill!(Aaug,zero(Precision))
                    #@view(Aaug[1:n_momentum,1:n_momentum]) .= J
                    #@view(Aaug[1:n_momentum,n_momentum+1]) .= F

                    # scalar
                        _, hmp1m, β = arnoldi_scalar!(J,F,method.V,method.Hs,method.w;reorth=false)
                    # block
                        #Rlast = arnoldi_block4!(J,F,method.V,method.Hb,method.W,method.R0;reorth=false)

                    order = 1.0 # energy error is order 1
                    ηHtarget = 1e-5
                    ηEtarget = 1e-5
                    k = min(1.0,(1.0-t)*dt/dt_local) # initial guess for k, limited to not overshoot final time 
                    k_old = 1.0
                    k_next = 1.0
                    q_init = m/4
                    q = q_init
                    ηH_old = 1.0
                    γ = 0.6
                    substeps = 0
                    ηE = Inf 
                    ηH = Inf
                    non_finite = false

                    while (!isfinite(ηE) || !isfinite(ηH)) || (ηE > ηEtarget || ηH > 1.4ηHtarget) 

                        # scalar
                            expm_pade!(k, method.Hsexp, method.Hs, method.PadeExpWorkspace_s)
                            #display(method.Hs)
                            #display(method.Hsexp)
                            #display(exp(method.Hs))
                            ferr = method.Hsexp[m,m+1]
                            #println("dtscale: ", dtscale, " k: ", k, " ferr: ", ferr)
                            ηH = abs(k * dtscale * β * hmp1m * ferr) # error estimate from Krylov subspace approximation of exponential action
                            @views mul!(δ, method.V[1:n_momentum, 1:m], method.Hsexp[1:m, m+1],β*k,zero(Precision)) # δ is energy change
                            δ .*= D # scale δ back to original space
                        # block
                            #expm_pade!(k, method.Hbexp, method.Hb, method.PadeExpWorkspace_s)
                            #ferr = @view(method.Hbexp[m-3:m,m+1])
                            #ηH = k * norm(Rlast * ferr)
                            #@views mul!(δ, method.V[1:n_momentum, 1:m], method.Hbexp[1:m, m+1],one(Precision),zero(Precision)) # δ is energy change
                            #δ .*= D # scale δ back to original space

                            #display(method.Hb)
                            #display(method.Hbexp)

                        @. fout = fold + δ
                        #@. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        Eold = dot(E, fold #=.+ df_Inj * k * dtscale=#)
                        ηE = abs(dot(E, fout) / Eold - one(Precision))

                        # KIOPS adaptive time 
                        #=if substeps == 0 || substeps == 1
                            q = q_init
                            q_old = q 
                        else
                            q_old = q
                            q = log(k/k_old) / log(ηH/ηH_old)
                        end=#

                        kE = k*(ηEtarget/(ηE+eps(ηEtarget)))^(1.0/(order+1)) # k from energy error estimate 
                        kH = k*(k*γ/((ηH+eps(ηHtarget))/ηHtarget))^(1/q) #(0.77 * 0.5(tanh(log10(ηH/ηHtarget)+1)+1)) # TODO: update this with a more specific method

                        k_old = k

                        if !isfinite(ηE) || !isfinite(ηH)
                            @warn "Rejecting step due to non-finite error estimate: ηE: $ηE, ηH: $ηH reducing time step by 0.5"
                            k *= 0.5
                            non_finite = true
                            break
                        elseif ηE < ηEtarget && ηH < 1.4ηHtarget # good step
                            # define how large to guess the next step can be
                            k_next = min(sqrt(kE),sqrt(kH),1.5) # max k is 2.0 or 1.0-t to avoid overshooting final time
                            println("Accepting step with k: ", k, " k_next: ", k_next, " t: ", t, " dt_local: ", dt_local," ηE: ", ηE, " ηH: ", ηH, " dtscale: ", dtscale)
                            break
                        else # 
                            k = min(0.9*kE,kH)
                            k = max(k,0.5*k_old)
                            println("Rejecting step with k: ", k, " t: ", t, " dt_local: ", dt_local," ηE: ", ηE, " ηH: ", ηH)
                        end

                        substeps += 1

                    end

                    if non_finite
                        dt_local = 0.5*dt_local
                        continue # skip to next spatial point if non-finite error estimate
                    end

                    t += k*dt_local / dt
                    println("space: ", off_space," region: Binary", " t: ", t, " dt: ", dt, " dt_local: ", k*dt_local)
                    dt_old = k*dt_local # used for updating dt_guess
                    dt_local = k_next*dt_local

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

            @inbounds method.dt_guess[off_space+1] = dt_old # update dt_guess with last plus one dt (avoid last as this could be limited by (1.0-t))

        else #

            if isassigned(method.M_Emi, off_space+1)
                @inbounds M_Emi = method.M_Emi[off_space+1]
            else
                EmiTrue = false
            end

            fold .= fstep 

            #if sum(fold) == Precision(0)
            #    continue
            #end
            t = 0.0

            dt_local = dt #/ 2 # initial guess for local time step, can adjust based on desired accuracy and problem stiffness

            kE = 1.0
            kϕ = 1.0
            k = 1.0

            while t < 1.0

                dtscale = Precision(dt_local / method.dt0) # scale for Jacobian as `vol` is calculated using `dt0` then the time step dt is just k as k*dt_local

                kold = k 
                kϕold = kϕ
                kEold = kE
                k = 1.0 # time step as a ratio of dt_local to dt_local before adaptive

                # EXPRB First Order Exponential Rosenbrock method with adaptive timestepping

                    # Form J
                    if EmiTrue
                        copyto!(Jsparse,M_Emi)
                        Jsparse *= dtscale * invA
                    else
                        fill!(Jsparse,zero(Precision))
                        dropzeros!(Jsparse)
                    end

                    # Form F
                    mul!(F,Jsparse,fold)
                    @. F += #=A *=# df_Inj * dtscale
                    if norm(F) == zero(Precision) # no change in distribution
                        break
                    end

                    Jsparse .*= D' # right mul
                    Jsparse .*= Dinv # left mul
                    F .*= Dinv # left mul

                    #@. J64 = Float64(J)

                    if has_injection 
                        #@. F64 = Float64(F)
                        #arnoldi!(Ks64,J64,F64;m=m,reorthogonalize=true)
                        arnoldi!(Ks,Jsparse,F;m=m,reorthogonalize=true)
                    else # no injection, just linear Jacobian so use exp over phi
                        #f64 = Float64.(Dinv * fold)
                        #arnoldi!(Ks64,J64,f64;m=m,reorthogonalize=true)
                        mul!(fscale,Dinv,fold)
                        arnoldi!(Ks,Jsparse,fscale;m=m,reorthogonalize=true)
                    end

                    V = ExponentialUtilities.getV(Ks)[:,1:end-1]
                    H = ExponentialUtilities.getH(Ks)[1:end-1,1:end]
                    #V = ExponentialUtilities.getV(Ks64)[:,1:end-1]
                    #H = ExponentialUtilities.getH(Ks64)[1:end-1,1:end]

                    norm_unscaled = norm(V' * V - I)
                    if isnan(norm_unscaled) # no J or F means no Krylov subspace generated, so just accept step as is and move on,
                        break
                    end

                    if #=cond((I - dt_local*H)) > 1f3 ||=#  norm(V' * V - I) > 1f-5 
                        dt_local = dt_local * 0.5
                        @warn "Krylov subspace has poor orthogonalisation $(norm(V' * V - I)), reducing time step to $dt_local"
                        continue
                    end

                    # Compute φ functions of H
                    if has_injection 
                        phiv!(ϕ,k,Ks,1;cache=ϕcache,correct=true,errest=false) # TODO: This allocates
                        @. δ = D * @view(ϕ[:,2]) * k
                        @. fout = fold + δ
                    else # no injection, just linear Jacobian so use exp over phi
                        expv!(δ,k#=*invA=#,Ks)
                        @. fout = D * δ
                    end


                    # adaptive time stepping

                        @. fout = ifelse(fout <= zero(Precision), zero(Precision),fout)

                        # energy error
                        if has_injection
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * dtscale)
                            Eold = dot(E, fold .+ df_Inj * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        #println("Energy error before adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        if isinf(ηE)
                            @warn "Energy error is Inf or NaN, may be unstable, consider reducing time step or adjusting ηtarget"
                            dt_local *= convert(typeof(dt_local), 1/64)
                            continue
                        end

                        dt_old = dt_local
                        order = 1.0 # energy error is order 1

                        ηtarget = 1e-16

                        kE = (1e-4/(ηE+eps(1e-4)))^(1.0/(order+1)) # k from energy error estimate
                        if has_injection
                            kϕmax = kϕold < 1.0 ? 1.0 + kϕold : 2.0
                            kϕ = min(kE,kϕmax) # max k is 2.0
                            errest = Inf
                            while errest > ηtarget # k from ϕv errestimate, can be more strict than energy error estimate
                                _, errest = phiv!(ϕ,kϕ,Ks,1;cache=ϕcache,correct=true,errest=true) # TODO: This allocates
                                println("ϕv error estimate during: ", errest)
                                if errest > ηtarget
                                    kϕ *= 0.77 * 0.5(tanh(log10(errest)-log10(ηtarget)+1)+1)
                                end
                            end
                            println("ϕv error estimate after: ", errest)
                            k = min(kE,kϕ,2.0)
                            println("kE: ", kE, " kϕ: ", kϕ, " errest: ", errest)
                        else
                            k = min(kE,2.0)
                        end
                        #k = 1.0
                        #println("dt_local: ", dt_local, " k: ", k, " new: ", dt_old*k, " old: ", dt_old)
                        dt_local = k*dt_old
                        #dt_local = max(dt_local, dt_old/1.01)
                        #println(dt)
                        if t + dt_local/dt > 1.0
                            dt_local = (1.0 - t) * dt
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t = 1.0
                        else
                            println("space: ", off_space," region: Linear", " t: ", t, " dt: ", dt, " dt_local: ", dt_local)
                            t += dt_local/dt
                        end
                        k = dt_local / dt_old

                    if k != 1.0
                        if has_injection 
                            phiv!(ϕ,k,Ks,1;cache=ϕcache,correct=true)
                            mul!(δ,D,@view(ϕ[:,2]),k#=*invA=#,zero(Precision))
                            @. fout = fold + δ
                        else # no injection, just linear Jacobian so use exp over phi
                            expv!(δ,k#=*invA=#,Ks)
                            @. fout = D.diag * δ
                        end

                        # energy after correction  
                        if has_injection 
                            ΔE = dot(E, fout) - dot(E, fold .+ df_Inj * k * dtscale)
                            Eold = dot(E, fold .+ df_Inj * k * dtscale)
                        else
                            ΔE = dot(E, fout) - dot(E, fold)
                            Eold = dot(E, fold)
                        end
                        ηE = abs(dot(E, fout) / Eold - one(Precision))
                        println("Energy error after adaptive: ", ηE, " ΔE: ", ΔE, " Eold: ", Eold, " Enew: ", dot(method.E, fout))

                        #=if ηE > 2e-4
                            @warn "Energy error in EXPRB1 step is large $ηE, may be unstable, consider reducing time step or adjusting ηtarget"
                        end=#
                        
                    end

                    @. fout = ifelse(fout <= method.n_cut, zero(Precision),fout)

                    #@. fold = ifelse(fout <= method.n_cut, zero(Precision),fout)
                    #@. fold = max(fout, zero(Precision))

                    #println(minimum(fold), " ", maximum(fold))

                    @. fold = fout

            end # while not accurate

            @. fstep = fold

        end

    end

    return nothing

end