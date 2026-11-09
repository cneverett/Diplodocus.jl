function EmissionCorrection!(PhaseSpace::PhaseSpaceStruct,GainMatrix3::AbstractArray{Float64,6},Parameters)

    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,Ext) = Parameters

    Characteristic = PhaseSpace.Characteristic

    if type == "Sync"

        force = SyncRadReact(mode=Ani(),B=Ext)
        Momentum = PhaseSpace.Momentum
        scheme = Momentum.scheme
        name_list = PhaseSpace.name_list

        Grids = PhaseSpace.Grids
        dE_list = Grids.dE_list

        name1_loc = findfirst(==(name1),name_list)
        name3_loc = findfirst(==(name3),name_list)

        dE1 = dE_list[name1_loc]
        dE3 = dE_list[name3_loc]

        p3r = Grids.pxr_list[name3_loc]
        p1m = Grids.mpx_list[name1_loc]      

        for px in axes(GainMatrix3, 4), py in axes(GainMatrix3,5), pz in axes(GainMatrix3,6) # loop over p1 states

            # critical frequency
            ω0 = abs((z1*1.6e-19*Ext))/(p1m[px]*9.11e-31)
            pc = 1.054e-34*ω0/(9.11e-31*3e8^2)*(p1m[px])^3
            pmin = p3r[1]
        
            GainSumE3 = zero(Float64)
            LossSumE1 = zero(Float64)

            px_num = px1_num
            pxp = RightBound(px,px_num,Closed())
            pxm = LeftBound(px,px_num,Closed())

            I_plus = zero(Float64)
            I_minus = zero(Float64)

            I_plus += IFluxFunction(force,PhaseSpace,name1_loc,"plus",1,1,1,1,px,py,pz)
            I_minus -= IFluxFunction(force,PhaseSpace,name1_loc,"minus",1,1,1,1,px,py,pz)

            # scheme
            (i_plus_right, i_plus_left, i_minus_right, i_minus_left) = SchemeCoefficients(scheme,I_plus,I_minus)

            #=
            ________________________________
            a |  I_m   | I_m+I_p | I_p    |
            __|________|_________|________|_
                    b-1       b        b+1  
            =#

            Mom_norm = MomentumSpaceNorm(Grids,name1_loc,px,py,pz)

            # normalised fluxes
            if px != pxp
                #LossSumE1 += convert(Float64,(I_plus * i_plus_right) / ((pxr[pxp+1]-pxr[pxp])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) * dE1[pxp]
                LossSumE1 += convert(Float64,(I_plus * i_plus_left) / Mom_norm)  * (dE1[px]-dE1[px+1])
            end
            if px != pxm
                LossSumE1 += convert(Float64,(I_minus * i_minus_right) / Mom_norm) * (dE1[px]-dE1[px-1])
                #LossSumE1 += convert(Float64,(I_minus * i_minus_left) / ((pxr[pxm+1]-pxr[pxm])*(pyr[py+1]-pyr[py])*(pzr[pz+1]-pzr[pz]))) * dE1[pxm]
            end

            # calculate total rate of energy gain from p1
            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                GainSumE3 += GainMatrix3[p3,u3,h3,px,py,pz] * dE3[p3]
            end

            vol = VolFunction(PhaseSpace,1,1,1,1)

            GainSumE3 *= vol

            
            if px != pxp && px != pxm
            #println(i_plus_right,i_plus_left,i_minus_right,i_minus_left,"I_plus = $I_plus, I_minus = $I_minus"," Ep1 = $(dE1[px+1]), E1 = $(dE1[px]), Em1 = $(dE1[px-1]), Gain: $GainSumE3, Loss: $LossSumE1, Correction = $(LossSumE1/GainSumE3), pc = $pc, pmin = $pmin, px = $px, py = $py, pz = $pz, Norm: $Mom_norm, Norm2: $(MomentumSpaceNorm(Grids,name3_loc,px+1,py,pz))")
            end

            if GainSumE3 != 0e0
                Correction = LossSumE1/GainSumE3
                if Correction < 0.0 
                    println("Negative correction factor, check flux calculations, setting correction to zero")
                    Correction = 0.0
                    @view(GainMatrix3[:,:,:,px,py,pz]) .= Correction * @view(GainMatrix3[:,:,:,px,py,pz])
                elseif pc < pmin # Peak of spectrum is below minimum photon momentum. To conserve energy we correct only the lowest photon momentum bin and assume isotropic emission as electron momentum is low. 
                    # calculate total rate of energy gain from lowest energy state
                    GainSumE3Min = zero(Float64)
                    if sum(@view(GainMatrix3[1,:,:,px,py,pz])) == 0.0 
                        @view(GainMatrix3[1,:,:,px,py,pz]) .= 1e-10 # add a random value that will be scaled
                        GainSumE3 += 1e-10 * dE3[1] * vol * py3_num * pz3_num # correct for new bins
                        Correction = LossSumE1/GainSumE3 
                        for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                            GainSumE3Min += GainMatrix3[1,u3,h3,px,py,pz] * dE3[1] * vol
                        end
                        @view(GainMatrix3[1,:,:,px,py,pz]) .= @view(GainMatrix3[1,:,:,px,py,pz]) * (Correction-1.0) * GainSumE3 / GainSumE3Min
                    else
                        for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                            GainSumE3Min += GainMatrix3[1,u3,h3,px,py,pz] * dE3[1] * vol
                        end
                        @view(GainMatrix3[1,:,:,px,py,pz]) .= @view(GainMatrix3[1,:,:,px,py,pz]) * (Correction-1.0) * GainSumE3 / GainSumE3Min
                    end
                else
                    if Correction > 1e2 || Correction < 1e-2 # if outside this range kernel is inaccurate or sync critical frequency out of range.
                        @warn "Correction factor large, but pc>pmin, check sync kernel calculations Ext = $Ext, p1x = $px, p1y = $py, p1z = $pz, pc = $pc, pmin = $pmin, GainE = $GainSumE3, LossE = $LossSumE1, Correction = $Correction, vol = $vol"
                    end
                    @view(GainMatrix3[:,:,:,px,py,pz]) .= Correction * @view(GainMatrix3[:,:,:,px,py,pz])
                end
            end

            # filter values below threshold after correction 
            # maximum gain value 
            GainMax = maximum(@view(GainMatrix3[:,:,:,px,py,pz]))
            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                if GainMatrix3[p3,u3,h3,px,py,pz] < eps(Float64) * GainMax
                    GainMatrix3[p3,u3,h3,px,py,pz] = 0.0
                end
                if isnan(GainMatrix3[p3,u3,h3,px,py,pz])
                    println("NaN value in GainMatrix3 after correction, setting to zero, p3 = $p3, u3 = $u3, h3 = $h3, px = $px, py = $py, pz = $pz, $GainMax, $GainSumE3")
                    GainMatrix3[p3,u3,h3,px,py,pz] = 0.0
                end
            end

        end # loop over p1 states

    else
        return nothing

    end

end