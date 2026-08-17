#function JetIntensityImage(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,observerDistance::Float64,observerInclination::Float64,constantCoordinateIndex::Int64)


using DiplodocusPlots
using LinearAlgebra
using StaticArrays
using CairoMakie
using GeometryBasics
using GeometryOps
using ColorSchemes


fileName = "M87_test_10_Spin.jld2"
fileName = "M87_t10_10x10_med_a5.jld2"
fileName = "M87_t50_10x10_med_a5.jld2"
#fileName = "M87_t100_10x10_med_a5.jld2"
fileName = "M87_t100_10x10_med_a9.jld2"
fileName = "M87_t200_10x10_med_a1.jld2"
fileName = "M87_t300_10x10_med_a1.jld2"
fileName = "M87_t300_10x10_med_a5.jld2"
fileName = "M87_t300_10x10_med_a9.jld2"

fileLocation = joinpath(pwd(),"..","DiplodocusTransport","TestingGround","Data");
(PhaseSpace,sol) = DiplodocusPlots.SolutionFileLoad(fileLocation,fileName);

constantCoordinateIndex_list = [2,3,4] # index of the constant coordinate (z=0 plane)
observerDistance = 16.8e6 * DiplodocusPlots.DiplodocusTransport.CONST_pcs / (6.5e9*DiplodocusPlots.DiplodocusTransport.CONST_rs_Sun) # distance from observer to origin in code units (16.8 Mpc for M87)
observerInclination = (90-17) * π / 180 # inclination angle of the observer in radians (0 = looking down the z-axis, π/2 = looking along the x-y plane)
screen_rotate = pi/4+(0.0)
screen_distance = 1.0 # distance from observer to image plane
ϕsamples = 41 
u_sub_samples = 8
two_jet = true # if true plots both approaching and receding jets

DiplodocusPlots.DiplodocusTransport.CONST_pcs
6.5e9*DiplodocusPlots.DiplodocusTransport.CONST_rs_Sun

################# define frequencies for images ######################

# 86 GHz = 6.96e-10 m_e, roughly EHT
p0_EHT = 6.96e-10 
# 10^15 Hz = 8.0e-6 m_e, roughly optical HST 
p0_HST = 8.0e-6
# 10^18 Hz = 8.0e-3 m_e, roughly X-ray Chandra
p0_Chandra = 8.0e-3
# 10^24 Hz = 8.0e3 m_e, roughly gamma-ray Fermi
p0_Fermi = 8.0e3
# 10^26 Hz = 8.0e5 m_e, roughly gamma-ray HESS, Veritas, MAGIC
p0_HESS = 8.0e5

white_purple_black = ColorScheme([ColorSchemes.Colors.RGB(1.0,1.0,1.0), ColorSchemes.Colors.RGB(223/255, 0, 255/255), ColorSchemes.Colors.RGB(0.0,0.0,0.0)]);

white_orange_black = ColorScheme([ColorSchemes.Colors.RGB(1.0,1.0,1.0), ColorSchemes.Colors.RGB(255/255, 127/255, 0.0), ColorSchemes.Colors.RGB(0.0,0.0,0.0)]);

white_green_black = ColorScheme([ColorSchemes.Colors.RGB(1.0,1.0,1.0), ColorSchemes.Colors.RGB(0.0, 255/255, 31/255), ColorSchemes.Colors.RGB(0.0,0.0,0.0)]);

p0_list = [p0_EHT,p0_HST,p0_Chandra,p0_Fermi,p0_HESS] # sample momenta for the intensity calculation
p0_color_list = [:gist_heat,:ice,ColorSchemes.ColorScheme(get(white_purple_black,range(1.0,0.0,length=256))),ColorSchemes.ColorScheme(get(white_orange_black,range(1.0,0.0,length=256))),ColorSchemes.ColorScheme(get(white_green_black,range(1.0,0.0,length=256)))] # colors for the intensity plots
p0_list_names = ["EHT","HST","Chandra","Fermi","HESS"]
p0_SED_color_list = [:red,:blue,:purple,:orange,:green]
p0_ColorLimits_list = [(5,10),(5,10),(5,10),(0,5),(-10,0)]


####### Figure setup ####################### 

fig = Figure(size=(720,1080),backgroundcolor=:black);
ax_EHT = Axis(fig[1,1];title="$(p0_list_names[1])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_HST = Axis(fig[2,1];title="$(p0_list_names[2])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_Chandra = Axis(fig[3,1];title="$(p0_list_names[3])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_Fermi = Axis(fig[4,1];title="$(p0_list_names[4])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_HESS = Axis(fig[5,1];title="$(p0_list_names[5])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_list = [ax_EHT, ax_HST, ax_Chandra, ax_Fermi, ax_HESS]
ax_SED = Axis(fig[6,1];title="SED",aspect=DataAspect())
hidedecorations!.(ax_list)

####### Plots on figure ####################

for freq_idx in eachindex(p0_list)

    p0 = p0_list[freq_idx]
    p0_color = p0_color_list[freq_idx]
    p0_name = p0_list_names[freq_idx]
    p0_SED_color = p0_SED_color_list[freq_idx]
    #p0_ColorLimits = p0_ColorLimits_list[freq_idx]
    ax = ax_list[freq_idx]

    points_I = Vector{Float64}(undef, 0) # Intensity points on the image plane to be plotted
    points_poly = Vector{Polygon}(undef, 0) # polygons for the surface to be plotted

    num_constantcoords = length(constantCoordinateIndex_list)
    α = 1.0/num_constantcoords

for coord_idx in eachindex(constantCoordinateIndex_list)

    constantCoordinateIndex = constantCoordinateIndex_list[coord_idx]

    ###########  Define Coordinates and Transforms ############

        sim_coordinates = PhaseSpace.Spacetime.coordinates
        coordinates = DiplodocusPlots.DiplodocusTransport.Cartesian()

    ########### Define Image Surface and Resolution ###########

        Grids = PhaseSpace.Grids
        xr_sim = Grids.xr
        yr_sim = Grids.yr
        zr_sim = Grids.zr
        mx_sim = Grids.mx
        my_sim = Grids.my
        mz_sim = Grids.mz

        # ϕ grid for image surface
        xr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.UniformGrid(2pi,0.0,ϕsamples))
        mx_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(xr_points_sim)

        # u grid for image surface
        a = 1/sqrt(2)
        ratio = sqrt(2)
        nsub = u_sub_samples
        cells = 10
        yr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.StretchGrid(0.0,a*ratio^((1-2^nsub)/2^nsub),ratio^(1/2^nsub),10*(2^nsub)))
        my_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(yr_points_sim)

        # v grid for image surface
        mz_point_sim = mz_sim[constantCoordinateIndex]

        #### Observer Position and Orientation ####

        pos_obs = zeros(Float64,3)
        pos_obs[1] = observerDistance*sin(observerInclination)
        pos_obs[2] = 0.0
        pos_obs[3] = observerDistance*cos(observerInclination)

        # form x-y basis of image plane perpendicular to the observer's line of sight
        screen_normvec = pos_obs/norm(pos_obs) # normal vector of the image plane
        screen_xvec = [0,0,1] - dot([0,0,1],screen_normvec)*screen_normvec # x basis vector of the image plane
        screen_xvec = screen_xvec/norm(screen_xvec) # normalize x basis vector
        screen_yvec = cross(screen_normvec,screen_xvec) # y basis vector of the image plane  
        screen_yvec = screen_yvec/norm(screen_yvec) # normalize y basis vector    

        @inline function normalToSurface!(pos_plot_normal::Vector{Float64},pos_sim::Vector{Float64})
            ϕ = pos_sim[1]
            u = pos_sim[2]
            v = pos_sim[3]
            norm = sqrt(u^2+v^2)
            z = (u^2 - v^2)/2 
            if z > 0
                pos_plot_normal[1] = u*cos(ϕ)/norm 
                pos_plot_normal[2] = u*sin(ϕ)/norm
                pos_plot_normal[3] = -v/norm
            else
                pos_plot_normal[1] = v*cos(ϕ)/norm 
                pos_plot_normal[2] = v*sin(ϕ)/norm
                pos_plot_normal[3] = -u/norm
            end
            return nothing
        end

        # calculate plotting points excluding ones that cannot be seen by the observer
        
        points_sim = Vector{Vector{Float64}}(undef, 0) # vector from sim origin to surface point
        points_plot = Vector{Vector{Float64}}(undef, 0)# vector from plot origin to surface point
        points_vertices = Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}(undef, 0) # vertices of the surface to be plotted
        points_surface_to_obs = Vector{Vector{Float64}}(undef, 0) # vector from surface point to observer
        points_vertices_surface_to_obs = Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}(undef, 0) # vectors from vertices of the surface to observer
        indices_sim = Vector{Vector{Int64}}(undef, 0)

        pos_plot = zeros(Float64,4)
        pos_plot_normal = zeros(Float64,3)
        pos_surface_to_obs = zeros(Float64,3)

        for j in 1:length(my_points_sim), i in 1:length(mx_points_sim)

            #= steps: 
                - load point in simulation coordinates
                - convert plot point to Cartesian coordinates 
                - evaluate vector from surface point to observer
                - work out normal vector to surface at that point
                if visible 
                    - push point to list of points to be plotted 
                    - push index of point to list of indices to be plotted
                    - push vector from surface point to observer to list of vectors to be plotted
                    - convert vertices of plot surface to Cartesian coordinates
                    - push vertices of plot surface to list of vertices to be plotted
                    - evaluate vector from vertices of plot surface to observer
            =# 

            x_sim = mx_points_sim[i]
            y_sim = my_points_sim[j]
            z_sim = mz_point_sim
            
            # convert sim coordinates to cartesian coordinates
            DiplodocusPlots.DiplodocusTransport.CoordinateTransform!([1,x_sim,y_sim,z_sim],pos_plot,sim_coordinates,coordinates)

            x_plot = pos_plot[2]
            y_plot = pos_plot[3]
            z_plot = pos_plot[4]

            normalToSurface!(pos_plot_normal,[x_sim,y_sim,z_sim])

            pos_surface_to_obs = pos_obs - [x_plot,y_plot,z_plot]
            visible = dot(pos_surface_to_obs,pos_plot_normal) > 0.0 # check if the surface point is visible to the observer

            
            if visible && z_plot > 0.0
                push!(points_sim, [x_sim,y_sim,z_sim])
                push!(indices_sim, [1, ceil(Int64,j/(length(my_points_sim)/10)), constantCoordinateIndex])
                push!(points_plot, [x_plot,y_plot,z_plot])
                push!(points_surface_to_obs, pos_surface_to_obs)

                # convert vertices of plot surface to Cartesian coordinates
                v1_sim = [1,xr_points_sim[i],yr_points_sim[j],mz_point_sim]
                v2_sim = [1,xr_points_sim[i+1],yr_points_sim[j],mz_point_sim]
                v3_sim = [1,xr_points_sim[i+1],yr_points_sim[j+1],mz_point_sim]
                v4_sim = [1,xr_points_sim[i],yr_points_sim[j+1],mz_point_sim]

                # convert sim coordinates to cartesian coordinates
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v1_sim,pos_plot,sim_coordinates,coordinates)
                x_v1 = pos_plot[2]
                y_v1 = pos_plot[3]
                z_v1 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v2_sim,pos_plot,sim_coordinates,coordinates)
                x_v2 = pos_plot[2]
                y_v2 = pos_plot[3]
                z_v2 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v3_sim,pos_plot,sim_coordinates,coordinates)
                x_v3 = pos_plot[2]
                y_v3 = pos_plot[3]
                z_v3 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v4_sim,pos_plot,sim_coordinates,coordinates)
                x_v4 = pos_plot[2]
                y_v4 = pos_plot[3]
                z_v4 = pos_plot[4]
                push!(points_vertices, ([x_v1,y_v1,z_v1],[x_v2,y_v2,z_v2],[x_v3,y_v3,z_v3],[x_v4,y_v4,z_v4]))

                pos_v1_to_obs = pos_obs - [x_v1,y_v1,z_v1]
                pos_v2_to_obs = pos_obs - [x_v2,y_v2,z_v2]
                pos_v3_to_obs = pos_obs - [x_v3,y_v3,z_v3]
                pos_v4_to_obs = pos_obs - [x_v4,y_v4,z_v4]
                push!(points_vertices_surface_to_obs, (pos_v1_to_obs,pos_v2_to_obs,pos_v3_to_obs,pos_v4_to_obs))
            end
        end

        if two_jet
        # repeat above with u and v i.e. sim y and z swapped to plot the receding jet    
            # v grid for image surface
            a = 1/sqrt(2)
            ratio = sqrt(2)
            nsub = u_sub_samples
            cells = 10
            zr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.StretchGrid(0.0,a*ratio^((1-2^nsub)/2^nsub),ratio^(1/2^nsub),10*(2^nsub)))
            mz_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(zr_points_sim)
            # u grid for image surface
            my_point_sim = my_sim[constantCoordinateIndex]

            for k in 1:length(mz_points_sim), i in 1:length(mx_points_sim)

            #= steps: 
                - load point in simulation coordinates
                - convert plot point to Cartesian coordinates 
                - evaluate vector from surface point to observer
                - work out normal vector to surface at that point
                if visible 
                    - push point to list of points to be plotted 
                    - push index of point to list of indices to be plotted
                    - push vector from surface point to observer to list of vectors to be plotted
                    - convert vertices of plot surface to Cartesian coordinates
                    - push vertices of plot surface to list of vertices to be plotted
                    - evaluate vector from vertices of plot surface to observer
            =# 

            x_sim = mx_points_sim[i]
            y_sim = my_point_sim
            z_sim = mz_points_sim[k]
            
            # convert sim coordinates to cartesian coordinates
            DiplodocusPlots.DiplodocusTransport.CoordinateTransform!([1,x_sim,y_sim,z_sim],pos_plot,sim_coordinates,coordinates)

            x_plot = pos_plot[2]
            y_plot = pos_plot[3]
            z_plot = pos_plot[4]

            normalToSurface!(pos_plot_normal,[x_sim,y_sim,z_sim])

            pos_surface_to_obs = pos_obs - [x_plot,y_plot,z_plot]
            visible = dot(pos_surface_to_obs,pos_plot_normal) > 0.0 # check if the surface point is visible to the observer

            
            if visible && z_plot < 0.0
                push!(points_sim, [x_sim,y_sim,z_sim])
                push!(indices_sim, [1, constantCoordinateIndex,ceil(Int64,k/(length(mz_points_sim)/10))])
                push!(points_plot, [x_plot,y_plot,z_plot])
                push!(points_surface_to_obs, pos_surface_to_obs)

                # convert vertices of plot surface to Cartesian coordinates
                v1_sim = [1,xr_points_sim[i],my_point_sim,zr_points_sim[k]]
                v2_sim = [1,xr_points_sim[i+1],my_point_sim,zr_points_sim[k]]
                v3_sim = [1,xr_points_sim[i+1],my_point_sim,zr_points_sim[k+1]]
                v4_sim = [1,xr_points_sim[i],my_point_sim,zr_points_sim[k+1]]

                # convert sim coordinates to cartesian coordinates
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v1_sim,pos_plot,sim_coordinates,coordinates)
                x_v1 = pos_plot[2]
                y_v1 = pos_plot[3]
                z_v1 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v2_sim,pos_plot,sim_coordinates,coordinates)
                x_v2 = pos_plot[2]
                y_v2 = pos_plot[3]
                z_v2 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v3_sim,pos_plot,sim_coordinates,coordinates)
                x_v3 = pos_plot[2]
                y_v3 = pos_plot[3]
                z_v3 = pos_plot[4]
                DiplodocusPlots.DiplodocusTransport.CoordinateTransform!(v4_sim,pos_plot,sim_coordinates,coordinates)
                x_v4 = pos_plot[2]
                y_v4 = pos_plot[3]
                z_v4 = pos_plot[4]
                push!(points_vertices, ([x_v1,y_v1,z_v1],[x_v2,y_v2,z_v2],[x_v3,y_v3,z_v3],[x_v4,y_v4,z_v4]))

                pos_v1_to_obs = pos_obs - [x_v1,y_v1,z_v1]
                pos_v2_to_obs = pos_obs - [x_v2,y_v2,z_v2]
                pos_v3_to_obs = pos_obs - [x_v3,y_v3,z_v3]
                pos_v4_to_obs = pos_obs - [x_v4,y_v4,z_v4]
                push!(points_vertices_surface_to_obs, (pos_v1_to_obs,pos_v2_to_obs,pos_v3_to_obs,pos_v4_to_obs))
            end
            end
        end

    ########### Sample Image Surface and Compute Intensity ####

        points_x = Vector{Float64}(undef, 0) # x points on the image plane to be plotted
        points_y = Vector{Float64}(undef, 0) # y points on the image plane to be plotted
        points_x_vertices = Vector{Tuple{Float64,Float64,Float64,Float64}}(undef, 0) # x vertices of the surface to be plotted
        points_y_vertices = Vector{Tuple{Float64,Float64,Float64,Float64}}(undef, 0) # y vertices of the surface to be plotted


        pvec_plot = zeros(Float64,4) # momentum vector of photon emitted to observer in plot frame
        pvec_sim = zeros(Float64,4) # momentum vector of photon emitted to observer in the sim frame
        pvec_tmp = zeros(Float64,4) # momentum vector of photon emitted to observer in the sim frame

        inve_sim = @MMatrix zeros(Float64,4,4) # inverse tetrad matrix at the surface point in the sim frame
        e_sim = @MMatrix zeros(Float64,4,4) # tetrad matrix at the surface point in the sim frame
        inve_plot = @MMatrix zeros(Float64,4,4) # tetrad matrix at the surface point in the plot frame
        e_plot = @MMatrix zeros(Float64,4,4) # tetrad matrix at the surface point in the plot frame

        pos_sim = @MVector zeros(Float64,4) # position of the surface point in the sim frame
        pos_plot = @MVector zeros(Float64,4) # position of the surface point in the plot frame

        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad

        for index in eachindex(points_plot)

            # work out frequency and θ, ϕ angles of momentum vector for r_point_to_obs vector in terms of local stationary tetrad at the surface point
                r = points_surface_to_obs[index]
                r_v1 = points_vertices_surface_to_obs[index][1]
                r_v2 = points_vertices_surface_to_obs[index][2]
                r_v3 = points_vertices_surface_to_obs[index][3]
                r_v4 = points_vertices_surface_to_obs[index][4]
                r_norm = norm(r)
                pvec_plot[1] = p0
                pvec_plot[2] = p0 * r[1] / r_norm
                pvec_plot[3] = p0 * r[2] / r_norm
                pvec_plot[4] = p0 * r[3] / r_norm

            # work out frequency and θ, ϕ angles in sim tetrad frame by converting momentum vector

                pos_sim .= [1.0, points_sim[index][1], points_sim[index][2], points_sim[index][3]] # position of the surface point in sim frame

                pos_plot .= [1.0, points_plot[index][1], points_plot[index][2], points_plot[index][3]] # position of the surface point in plot frame

                #println("pos_sim: $pos_sim, points_sim: $(points_sim[index])")

                DiplodocusPlots.DiplodocusTransport.InverseTetradComponents!(pos_plot,inve_plot,metric,coordinates,DiplodocusPlots.DiplodocusTransport.StationaryObserverTetrad())
                DiplodocusPlots.DiplodocusTransport.InverseTetradComponents!(pos_sim,inve_sim,metric,sim_coordinates,tetrad)
                DiplodocusPlots.DiplodocusTransport.TetradComponents!(pos_sim,e_sim,metric,sim_coordinates,tetrad)
                DiplodocusPlots.DiplodocusTransport.TetradComponents!(pos_plot,e_plot,metric,coordinates,DiplodocusPlots.DiplodocusTransport.StationaryObserverTetrad())

                #display(inve_plot)
                #display(inve_sim)
                #display(e_sim)

                #pvec_sim =  e_sim * inve_plot * pvec_plot # convert momentum vector from plot frame to sim frame

                #println("p0: $((inve_plot * pvec_plot)[1]), p0comps: $(norm((inve_plot * pvec_plot)[2:4]))")

                #println("p0vec_plot: $pvec_plot")
                #println("pos_plot: $pos_plot")
                CoordinateBasisTransform!(e_plot' * pvec_plot,pvec_tmp,pos_plot,coordinates,sim_coordinates) # convert momentum vector from plot frame to sim frame
                #println("p0vec_plot: $pvec_tmp")

                g = @MMatrix zeros(Float64,4,4)
                DiplodocusPlots.DiplodocusTransport.MetricComponents!(pos_sim,g,metric,sim_coordinates)

                #println(pvec_tmp' * g * pvec_tmp) # check that the momentum vector is null in the sim frame

                pvec_sim =  inve_sim' * pvec_tmp # convert momentum vector from plot frame to sim frame

                #println(-pvec_sim[1] * pvec_sim[1] + dot(pvec_sim[2:4],pvec_sim[2:4])) # check that the momentum vector is null in the sim frame
                
                #println("p0sim: $(pvec_sim[1]), p0comps: $(norm(pvec_sim[2:4]))")
            # work out value of distribution function at the surface point for the sim frame momentum vector

                p_sim = pvec_sim[1] # momentum p in sim frame
                u_sim = pvec_sim[4]/pvec_sim[1] # momentum u=cos(θ) = pcos(θ)/p in sim frame
                h_sim = mod(atan(pvec_sim[3],pvec_sim[2]),2pi) # momentum ϕ in sim frame tan(ϕ) = psin(θ)sin(ϕ)/psin(θ)cos(ϕ) in sim frame

                # TODO: make these next steps general 
                p_loc = DiplodocusPlots.DC.location(PhaseSpace.Momentum.px_low_list[2],PhaseSpace.Momentum.px_up_list[2],PhaseSpace.Momentum.px_num_list[2],p_sim,PhaseSpace.Momentum.px_grid_list[2])
                u_loc = DiplodocusPlots.DC.location(-1.0,1.0,PhaseSpace.Momentum.py_num_list[2],u_sim,PhaseSpace.Momentum.py_grid_list[2])
                h_loc = DiplodocusPlots.DC.location(0.0,2*pi,PhaseSpace.Momentum.pz_num_list[2],h_sim,PhaseSpace.Momentum.pz_grid_list[2])
                # location of the surface point in the sim frame grid
                x = indices_sim[index][1]
                y = indices_sim[index][2]
                z = indices_sim[index][3]
                state_index = DiplodocusPlots.DiplodocusTransport.GlobalIndicesToStateIndex(PhaseSpace,x,y,z,p_loc,u_loc,h_loc,2) # TODO: add species index to this function call

                #println("x: $x, y: $y, z: $z, p_loc: $p_loc, u_loc: $u_loc, h_loc: $h_loc, p_sim: $p_sim, u_sim: $u_sim, h_sim: $h_sim")

                f = sol.f[end][state_index]

            # work out intensity emitted from surface point by f*ν^3 where ν is in observer frame

                # dN/p^2dpdΩ = I_ν / p^3
                # f = dN/dpdΩ
                Iν = f * p0^3 / p_sim^2 # intensity emitted from surface point in observer frame

                #println("pos_plot: $(points_plot[index]), pos_sim: $(points_sim[index])")
                
                #println("pvec_plot: $pvec_plot, pvec_sim: $pvec_sim")
                #println("f: $f, Iν: $Iν")

            # work out position on the screen place where r_point_to_obs intercepts it

                # form vector between screen plane origin and intercept point from surface point to observer
                avec_screen = screen_distance * (pos_obs/norm(pos_obs) - r .* norm(pos_obs) ./(dot(r,pos_obs)))
                # calculate the x and y coordinates of the intercept point on the screen plane
                x_screen = dot(avec_screen,screen_xvec)
                y_screen = dot(avec_screen,screen_yvec)
                
            # add intensity point and screen point to list to be plotted

                push!(points_x, x_screen)
                push!(points_y, y_screen)
                push!(points_I, Iν)

            # work out positions on the screen place where r_vertices_to_obs intercepts it

                # form vector between screen plane origin and intercept point from surface point to observer
                avec_screen_v1 = screen_distance * (pos_obs/norm(pos_obs) - r_v1 .* norm(pos_obs) ./(dot(r_v1,pos_obs)))
                avec_screen_v2 = screen_distance * (pos_obs/norm(pos_obs) - r_v2 .* norm(pos_obs) ./(dot(r_v2,pos_obs)))
                avec_screen_v3 = screen_distance * (pos_obs/norm(pos_obs) - r_v3 .* norm(pos_obs) ./(dot(r_v3,pos_obs)))
                avec_screen_v4 = screen_distance * (pos_obs/norm(pos_obs) - r_v4 .* norm(pos_obs) ./(dot(r_v4,pos_obs)))

                # calculate the x and y coordinates of the intercept point on the screen plane
                x_screen_v1 = dot(avec_screen_v1,screen_xvec)
                y_screen_v1 = dot(avec_screen_v1,screen_yvec)
                x_screen_v2 = dot(avec_screen_v2,screen_xvec)
                y_screen_v2 = dot(avec_screen_v2,screen_yvec)
                x_screen_v3 = dot(avec_screen_v3,screen_xvec)
                y_screen_v3 = dot(avec_screen_v3,screen_yvec)
                x_screen_v4 = dot(avec_screen_v4,screen_xvec)
                y_screen_v4 = dot(avec_screen_v4,screen_yvec)

                push!(points_x_vertices, (x_screen_v1,x_screen_v2,x_screen_v3,x_screen_v4))
                push!(points_y_vertices, (y_screen_v1,y_screen_v2,y_screen_v3,y_screen_v4))

                poly = Polygon(Point2f[(x_screen_v1,y_screen_v1),(x_screen_v2,y_screen_v2),(x_screen_v3,y_screen_v3),(x_screen_v4,y_screen_v4)])
                push!(points_poly, poly)

        end

end # constant coordinate index

    ############ Plot Image Surface and Intensity ###############

    Iν_max = maximum(points_I)
    Iν_min = max(minimum(points_I),1e-10*Iν_max)

    hm = poly!(ax,points_poly,color=log10.(points_I),
        colormap=p0_color,
        colorrange=(log10(Iν_min),log10(Iν_max)),
        strokewidth=0.0,
        fxaa = true,
        transparency = true,
        alpha=α
    )

    Colorbar(fig[freq_idx,2],hm;label="log10(Intensity)",labelcolor=:white,ticklabelcolor=:white,minortickcolor=:white)

    # SED plot 
    νFν = sum(GeometryOps.area.(points_poly) .* points_I) * p0 #TODO: fix

    scatter!(ax_SED,[log10(p0)],[log10(νFν)],markersize=10,color=p0_SED_color,marker=:circle)

end # freq plots

fig

DiplodocusPlots.save("M87_jet_images_new.png",fig)

    #f = Figure()
    #ax = Axis(f[1,1];title="Jet Intensity Image",aspect=DataAspect(),backgroundcolor = :black)
    #ax2 = Axis(f[1,2];title="Jet Intensity Image (log scale)",aspect=DataAspect(),backgroundcolor = :black)

    #f

    #length(points_x)
    
    #DiplodocusPlots.save("M87_jet_render.svg",f)
    #DiplodocusPlots.save("M87_jet_render.png",f)

#end


function CoordinateTransform!(pos_from::AbstractVector{T},pos_to::AbstractVector{T},::DiplodocusPlots.DiplodocusTransport.Cartesian,::DiplodocusPlots.DiplodocusTransport.Paraboloidal) where T

    #t = pos_from[1]
    #ϕ = pos_from[2]
    #u = pos_from[3]
    #v = pos_from[4]

    t = pos_from[1]
    x = pos_from[2]
    y = pos_from[3]
    z = pos_from[4]

    ρ = sqrt(z^2+x^2+y^2)
    ϕ = atan(y, x)
    u = sqrt(sqrt(ρ^2+z^2)+z)
    v = sqrt(sqrt(ρ^2+z^2)-z)

    pos_to[1] = t 
    pos_to[2] = ϕ
    pos_to[3] = u
    pos_to[4] = v

    return nothing

end

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::Paraboloidal,::Cartesian) where T

    t = pos_from[1]
    ϕ = pos_from[2]
    u = pos_from[3]
    v = pos_from[4]

    vt = vec_from[1]
    vϕ = vec_from[2]
    vu = vec_from[3]
    vv = vec_from[4]

    vt::T = vt
    vx::T = v*cos(ϕ) * vu + u*cos(ϕ) * vv - u*v*sin(ϕ) * vϕ
    vy::T = v*sin(ϕ) * vu + u*sin(ϕ) * vv + u*v*cos(ϕ) * vϕ
    vz::T = u * vu - v * vv

    vec_to[1] = vt 
    vec_to[2] = vx
    vec_to[3] = vy
    vec_to[4] = vz

    return nothing

end

function CoordinateBasisTransform!(vec_from::AbstractVector{T},vec_to::AbstractVector{T},pos_from::AbstractVector{T},::DiplodocusPlots.DiplodocusTransport.Cartesian,::DiplodocusPlots.DiplodocusTransport.Paraboloidal) where T

    t = pos_from[1]
    x = pos_from[2]
    y = pos_from[3]
    z = pos_from[4]

    r = sqrt(x^2+y^2+z^2)
    ρ = sqrt(x^2+y^2)

    #println("x: $x, y: $y, z: $z, r: $r, ρ: $ρ")

    vt = vec_from[1]
    vx = vec_from[2]
    vy = vec_from[3]
    vz = vec_from[4]

    vt = vt
    vϕ = -y/ρ^2 * vx + x/ρ^2 * vy
    vu = x*sqrt(r-z)/(2r*ρ) * vx + y*sqrt(r-z)/(2r*ρ) * vy + sqrt(r+z)/(2r) * vz
    vv = x*sqrt(r+z)/(2r*ρ) * vx + y*sqrt(r+z)/(2r*ρ) * vy - sqrt(r-z)/(2r) * vz

    #println("vt: $vt, vϕ: $vϕ, vu: $vu, vv: $vv")

    vec_to[1] = vt 
    vec_to[2] = vϕ
    vec_to[3] = vu
    vec_to[4] = vv

    return nothing

end


using CairoMakie
using CairoMakie.GeometryBasics
using GeometryOps


f = Figure()
Axis(f[1, 1]; backgroundcolor = :gray15)

# vector of polygons
ps = [Polygon(rand(Point2f, 4) .+ Point2f(i, j))
    for i in 1:5 for j in 1:10]

poly!(ps, color = rand(RGBf, length(ps)))

f