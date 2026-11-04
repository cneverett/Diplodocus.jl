#function JetIntensityImage(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,observerDistance::Float64,observerInclination::Float64,constantCoordinateIndex::Int64)


using DiplodocusPlots
using LinearAlgebra
using StaticArrays
using CairoMakie
using GeometryBasics
using GeometryOps
using ColorSchemes

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

fileName = "M87_t300_10x10_med_a1.jld2"
fileName = "M87_t300_10x10_med_a5.jld2"
fileName = "M87_t300_10x10_med_a9.jld2"

fileName = "M87_t300_10x10_med_am5.jld2"

fileName = "M87_t256_8x8_med_a1.jld2"
fileName = "M87_t256_8x8_med_a5.jld2"
fileName = "M87_t256_8x8_med_a9.jld2"

fileLocation = joinpath(pwd(),"..","..","Papers","Thesis","Thesis_julia","M87Sims");

(PhaseSpace,sol) = DiplodocusPlots.SolutionFileLoad(fileLocation,fileName);

    name_list = PhaseSpace.name_list
    photon_idx = findfirst(x->x=="Pho",name_list)

    Grids = PhaseSpace.Grids
    xr_sim = Grids.xr
    yr_sim = Grids.yr
    zr_sim = Grids.zr
    mx_sim = Grids.mx
    my_sim = Grids.my
    mz_sim = Grids.mz
    dp = Grids.dpx_list[photon_idx] # for photon
    du = Grids.dpy_list[photon_idx] # for photon 
    dh = Grids.dpz_list[photon_idx] # for photon
    mp = Grids.mpx_list[photon_idx] # for photon

CONST_c = DiplodocusPlots.DiplodocusTransport.CONST_c
CONST_h = DiplodocusPlots.DiplodocusTransport.CONST_ħ * (2pi)
CONST_mEle = DiplodocusPlots.DiplodocusTransport.CONST_mEle
CHAR_n = PhaseSpace.Characteristic.CHAR_number_density
CHAR_l = PhaseSpace.Characteristic.CHAR_length

constantCoordinateIndex_list = [3] # index of the constant coordinate (z=0 plane)
observerDistance = 16.8e6 * DiplodocusPlots.DiplodocusTransport.CONST_pcs / CHAR_l # distance from observer to origin in code units (16.8 Mpc for M87)
observerInclination = (17) * π / 180 # inclination angle of the observer in radians (0 = 0 looking down the z-axis, π/2 = looking along the x-y plane)
screen_rotate = pi/4+(0.0)
screen_distance = 1.0 # distance from observer to image plane
ϕsamples = 61 
usamples = 5 # increased number of points in between simulation grid points should be 2^n-1 values
two_jet = true # if true plots both approaching and receding jets

nsub = usamples
a = 1/sqrt(2)
ratio = sqrt(2)

DiplodocusPlots.DiplodocusTransport.CONST_pcs
6.5e9*DiplodocusPlots.DiplodocusTransport.CONST_rs_Sun


################# define frequencies for images ######################

# 10^11 Hz = 8.0e-10 m_e, roughly EHT (roughly 86 GHz)
p0_EHT = 8.0e-10 
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

white_blue_black = ColorScheme([ColorSchemes.Colors.RGB(1.0,1.0,1.0), ColorSchemes.Colors.RGB(17/255, 63/255, 250/255), ColorSchemes.Colors.RGB(0.0,0.0,0.0)]);

white_green_black = ColorScheme([ColorSchemes.Colors.RGB(1.0,1.0,1.0), ColorSchemes.Colors.RGB(0.0, 255/255, 31/255), ColorSchemes.Colors.RGB(0.0,0.0,0.0)]);

p0_list = 8 * 10 .^ range(-12,5,length=18)
p0_list_plot = [8e-10,8e-6,8e-3,8e3,8e5] # sample momenta for the intensity calculation
p0_color_list = [:gist_heat,ColorSchemes.ColorScheme(get(white_blue_black,range(1.0,0.0,length=256))),ColorSchemes.ColorScheme(get(white_purple_black,range(1.0,0.0,length=256))),ColorSchemes.ColorScheme(get(white_orange_black,range(1.0,0.0,length=256))),ColorSchemes.ColorScheme(get(white_green_black,range(1.0,0.0,length=256)))] # colors for the intensity plots
p0_list_names = [L"100\,\text{GHz, (Radio)}",L"10^{15}\,\text{Hz, (HST, JWST)}",L"10^{18}\,\text{Hz}=4.14\,\text{keV, (Chandra, NuSTAR)}",L"10^{24}\,\text{Hz}=4.14\,\text{GeV, (Fermi)}",L"10^{26}\,\text{Hz}=414\,\text{GeV, (HESS, MAGIC, VERITAS)}"]
p0_SED_color_list = [:red,:blue,:purple,:orange,:green]
p0_ColorLimits_list = [(5,10),(5,10),(5,10),(0,5),(-10,0)]

nuFnu = zeros(Float64,length(p0_list))


####### Figure setup ####################### 

fig = Figure(size=(720,800),backgroundcolor=:black);
ax_EHT = Axis(fig[1,1];title="$(p0_list_names[1])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_HST = Axis(fig[2,1];title="$(p0_list_names[2])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_Chandra = Axis(fig[3,1];title="$(p0_list_names[3])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_Fermi = Axis(fig[4,1];title="$(p0_list_names[4])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_HESS = Axis(fig[5,1];title="$(p0_list_names[5])",aspect=DataAspect(),backgroundcolor = :black,titlecolor=:white)
ax_list = [ax_EHT, ax_HST, ax_Chandra, ax_Fermi, ax_HESS]
#ax_SED = Axis(fig[6,1];title="SED",aspect=DataAspect())
hidedecorations!.(ax_list)
fill!(nuFnu,0.0)

####### Plots on figure ####################

for freq_idx in eachindex(p0_list)

    p0 = p0_list[freq_idx]

    points_Iv = Vector{Float64}(undef, 0) # Intensity points on the image plane to be plotted
    points_poly = Vector{Polygon}(undef, 0) # polygons for the surface to be plotted

    num_constantcoords = length(constantCoordinateIndex_list)
    α = 1-(1.0-0.95)^(1/num_constantcoords) # target an alpha value of 0.95 for fully overlapped elements

for coord_idx in eachindex(constantCoordinateIndex_list)

    constantCoordinateIndex = constantCoordinateIndex_list[coord_idx]

    ###########  Define Coordinates and Transforms ############

        sim_coordinates = PhaseSpace.Spacetime.coordinates
        coordinates = DiplodocusPlots.DiplodocusTransport.Cartesian()

    ########### Define Image Surface and Resolution ###########

        # ϕ grid for image surface
        xr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.UniformGrid(2pi,0.0,ϕsamples))
        mx_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(xr_points_sim)

        # u grid for image surface
        nsub = usamples
        yr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.StretchGrid(0.0,1.0,nsub,8*(nsub+1)/2))
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
                push!(indices_sim, [1, ceil(Int64,j/(length(my_points_sim)/length(my_sim))), constantCoordinateIndex])
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
            nsub = usamples
            zr_points_sim = DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.StretchGrid(0.0,1.0,nsub,8*(nsub+1)/2))
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
                println(length(mz_points_sim))
                push!(indices_sim, [1, constantCoordinateIndex,ceil(Int64,k/(length(mz_points_sim)/length(mz_sim)))]) 
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

        g = @MMatrix zeros(Float64,4,4)

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

                println("x: $x, y: $y, z: $z, p_loc: $p_loc, u_loc: $u_loc, h_loc: $h_loc, p_sim: $p_sim, u_sim: $u_sim, h_sim: $h_sim")

                f = sol.f[end][state_index]

            # work out intensity emitted from surface point by f*ν^3 where ν is in observer frame

                # dN/p^2dpdΩ = I_ν / ν^3
                # f = dN/dpdudh * dpdudh
                # therfore dN/p^2dpdΩ = f / p^2dpdudh, where p is the sim p and dp is the width of the bin
                Iν = (f/(mp[p_loc]^2*dp[p_loc]*du[u_loc]*dh[h_loc])) * p0^3 # intensity emitted from surface point in observer frame
                # convert to SI 
                Iν *= CHAR_n * (CONST_mEle * CONST_c)^3 / (CONST_mEle * CONST_c)^3 * CONST_c^3 / CONST_h^3 * CONST_h^4 / CONST_c^2 # convert to SI units of W/m^2/Hz/sr

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
                push!(points_Iv, Iν)

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

    println("p0: $p0")
    if p0 in p0_list_plot
        println("here2")
        plot_idx = findfirst(x->x==p0,p0_list_plot)
        p0_color = p0_color_list[plot_idx]
        p0_name = p0_list_names[plot_idx]
        p0_SED_color = p0_SED_color_list[plot_idx]
        #p0_ColorLimits = p0_ColorLimits_list[freq_idx]
        ax = ax_list[plot_idx]

        if plot_idx == 1 
            Iν_max = 10^-5
            Iν_min = 10^-11
        elseif plot_idx == 2
            Iν_max = 10^-9
            Iν_min = 10^-15
        elseif plot_idx == 3
            Iν_max = 10^-13
            Iν_min = 10^-19
        elseif plot_idx == 4
            Iν_max = 10^-21
            Iν_min = 10^-27
        elseif plot_idx == 5
            Iν_max = 10^-29
            Iν_min = 10^-35
        end

        hm = poly!(ax,points_poly,color=log10.(points_Iv),
            colormap=p0_color,
            colorrange=(log10(Iν_min),log10(Iν_max)),
            strokewidth=0.0,
            fxaa = true,
            transparency = false,
            alpha=α,
        )

        Colorbar(fig[plot_idx,2];label=L"\log_{10}\left(I_\nu \left[\text{W}/\text{m}^2\,\text{sr}\,\text{Hz}\right]\right)",colormap=hm.colormap,colorrange=hm.colorrange,labelcolor=:white,ticklabelcolor=:white,minortickcolor=:white,tickcolor=:white,topspinecolor=:white,bottomspinecolor=:white,leftspinecolor=:white,rightspinecolor=:white,lowclip = DiplodocusPlots.Makie.RGB{Float64}(0.0,0.0,0.0))

        line_length = screen_distance * tand(0.1/3600000) # 0.1 miliarcsecond length on the screen plane

        lines!(ax,[-line_length,0],[-line_length,-line_length],color=:white,linestyle=:solid,linewidth=2)
        text!(ax,(-line_length,-line_length),text=L"0.1~\text{mas}",color=:white,align=(:center,:top))

    end

    # SED plot 
    νFν = sum(GeometryOps.area.(points_poly) .* points_Iv) * p0 * CONST_c / CONST_h / num_constantcoords / (observerDistance #=* CHAR_l=#)^2 # SI units of W/m^2

    #scatter!(ax_SED,[log10(p0)],[log10(νFν)],markersize=10,color=p0_SED_color,marker=:circle)

    nuFnu[freq_idx] = νFν

end # freq plots

fig

DiplodocusPlots.save("M87/M87_jet_images_a1.svg",fig)
DiplodocusPlots.save("M87/M87_jet_images_a5.svg",fig)
DiplodocusPlots.save("M87/M87_jet_images_a9.svg",fig)

nuFnu1 = copy(nuFnu)
nuFnu5 = copy(nuFnu)
nuFnu9 = copy(nuFnu)

plot(log10.((CONST_mEle * CONST_c^2 / CONST_h) * p0_list),log10.(nuFnu1),color=:black,marker=:circle,markersize=6,label="SED")
plot(log10.((CONST_mEle * CONST_c^2 / CONST_h) * p0_list),log10.(nuFnu5),color=:black,marker=:circle,markersize=6,label="SED")
plot(log10.((CONST_mEle * CONST_c^2 / CONST_h) * p0_list),log10.(nuFnu9),color=:black,marker=:circle,markersize=6,label="SED")


############## Multi-wavelength SED plot ####################

nuFnu1 =  [1.4580812643621608e-16,9.68985637681044e-16,2.9051503711644477e-15, 6.0612714719756845e-15, 8.654718033999387e-15, 5.684310222008691e-15, 1.76562681941999e-15, 7.926269127225466e-16, 1.0900355345635896e-15, 1.432683557413196e-15, 1.3215485364104159e-15, 8.173109532627206e-16, 3.648004156722957e-16, 1.3632894281061916e-16, 1.5388920271712214e-17, 2.7514908877812787e-19, 7.095299576611425e-22, 3.1038160834539794e-27] 


nuFnu5 = [6.699108418090639e-16,6.426246920208894e-15,1.7660379247902053e-14, 2.6685251723098268e-14, 2.8262972196387527e-14, 1.7272219946329712e-14, 5.7008295021228784e-15, 3.04936949788703e-15, 4.0362398479438864e-15, 4.8692776587083576e-15, 4.2896340766207766e-15, 2.661405808524775e-15, 1.2227527605069505e-15, 4.715929810341174e-16, 5.241786962396879e-17, 9.999654808970028e-19, 3.408985088740251e-21, 3.173852200208438e-26]

nuFnu9 =  [8.900605348317896e-16,1.1283341394442055e-14,3.791678205813041e-14, 6.570522608800822e-14, 8.067844465775689e-14, 6.844702338976266e-14, 3.0436704512808556e-14, 1.0873067584465521e-14, 1.029492363529832e-14, 1.3156677055074721e-14, 1.3163335448111218e-14, 9.128800441449115e-15, 4.837043720635128e-15, 2.146753058760761e-15, 4.519794130081836e-16, 1.9572383028963127e-17, 2.021288979411133e-19, 4.9985244704340276e-23]

nu = (CONST_mEle * CONST_c^2 / CONST_h) * p0_list
eV = nu * CONST_h / DiplodocusPlots.DiplodocusTransport.CONST_q

ergpercm2s_to_Wperm2 = 1e-3
eV_to_nu = DiplodocusPlots.DiplodocusTransport.CONST_q / CONST_h

### M87 data from EHT 2017
VLBI_nu = [1.70E09,8.40E09,1.50E10,2.40E10,2.20E10,2.20E10,4.30E10,2.17E10,4.34E10,8.68E10,1.29E11,2.40E10,4.30E10,8.60E10,2.20E11,2.20E11,2.30E11]
VLBI_nuFnu = [1.28E-14,1.09E-13,1.70E-13,2.93E-13,2.90E-13,2.95E-13,4.73E-13,3.54E-13,5.69E-13,9.46E-13,1.07E-12,2.00E-13,3.62E-13,6.30E-13,2.87E-12,2.90E-12,1.51E-12]
VLBI_nuFnu_sigma = [0.14E-14
0.11E-13
0.17E-13
0.29E-13
0.29E-13
0.29E-13
0.47E-13
0.37E-13
0.56E-13
0.95E-13
0.32E-12
0.20E-13
0.36E-13
0.60E-13
0.29E-12
0.29E-12
0.37E-12]

HST_nu = [5.20E14,1.10E15]
HST_nuFnu = [2.55E-12,1.80E-12]
HST_nuFnu_sigma = [0.21E-12
0.26E-12]

Swift_nu = [5.47E14,6.82E14,8.65E14,1.15E15,1.33E15,1.55E15]
Swift_nuFnu = [4.33E-12,2.64E-12,1.29E-12,9.43E-13,2.55E-12,1.38E-12]
Swift_nuFnu_sigma = [1.56E-12
1.43E-12
0.70E-12
4.25E-13
0.39E-12
0.36E-12]
Swift_nu_upper = [1.08E18]
Swift_nuFnu_upper = [2.25E-12]
# Swift frequency bounds are only on the last data point
Swift_nu_low = [2.0] .* 1e3 * eV_to_nu 
Swift_nu_up = [10.0] .* 1e3 * eV_to_nu 

ChNu_nu = [1.09E17,1.40E17,1.79E17,2.29E17,2.93E17,3.75E17,4.80E17,6.13E17,7.85E17,1.00E18,1.28E18,1.64E18,2.10E18,2.69E18,3.44E18,4.40E18,5.62E18,7.19E18,9.20E18,1.18E19]
ChNu_nuFnu = [1.04E-12,1.02E-12,1.00E-12,9.89E-13,9.70E-13,9.53E-13,9.36E-13,9.22E-13,9.03E-13,8.85E-13,8.70E-13,8.58E-13,8.44E-13,8.31E-13,8.15E-13,8.05E-13,7.92E-13,7.77E-13,7.65E-13,7.55E-13]
ChNu_nuFnu_sigma = [0.09E-12
0.08E-12
0.06E-12
0.52E-13
0.41E-13
0.31E-13
0.25E-13
0.24E-13
0.28E-13
0.35E-13
0.43E-13
0.53E-13
0.60E-13
0.69E-13
0.78E-13
0.86E-13
0.94E-13
1.02E-13
1.10E-13
1.18E-13]
ChNu_nu_low = [0.4,0.51,0.65,0.84,1.07,1.37,1.75,2.24,2.87,3.67,4.69,6.00,7.67,9.82,12.6,16.1,20.5,26.3,33.6,43.0] .* 1e3 * eV_to_nu 
ChNu_nu_up = [0.51,0.65,0.84,1.07,1.37,1.75,2.24,2.87,3.67,4.69,6.00,7.67,9.82,12.6,16.1,20.5,26.3,33.6,43.0,55.0] .* 1e3 * eV_to_nu 

Fermi_nu = [7.65E22,7.65E23,7.65E24]
Fermi_nuFnu = [3.27E-12,1.87E-12,4.96E-12]
Fermi_nuFnu_sigma = [1.01E-12,0.75E-12,2.96E-12]
Fermi_nu_upper = [7.65E25]
Fermi_nuFnu_upper = [5.26E-11]
Fermi_nu_low = [0.1,1.0,10.0] * 1e9 * eV_to_nu
Fermi_nu_up = [1.0,10.0,100.0] .* 1e9 * eV_to_nu
Fermi_nu_low_upper = [100.0] * 1e9 * eV_to_nu
Fermi_nu_up_upper = [1000.0] .* 1e9 * eV_to_nu

HMV_nu = [3.82E25,6.06E25,9.60E25,1.52E26,2.41E26,3.82E26,1.14E26,1.61E26]
HMV_nuFnu = [1.40E-13,6.67E-13,3.83E-13,7.72E-13,4.70E-13,8.21E-13,7.32E-13,6.47E-13]
HMV_nuFnu_sigma = [0.82E-13,4.68E-13,3.54E-13,2.94E-13,2.41E-13,3.30E-13,3.90E-14,3.70E-13]

HMV_nu_low = [126.0,199.0,315.0,500.0,792.0,1256.0,398.0,562.0] .* 1e9 * eV_to_nu
HMV_nu_up = [199.0,315.0,500.0,792.0,1256.0,1991.0,562.0,793.0] .* 1e9 * eV_to_nu

HMV_nu_upper = [7.31E25,1.16E26,1.83E26,2.91E26,4.61E26,5.73E25,8.10E25,2.28E26,3.22E26]
HMV_nuFnu_upper = [3.22E-12,1.44E-12,9.74E-13,8.88E-13,7.70E-13,9.97E-12,1.70E-12,9.62E-13,1.13E-12]

HMV_nu_low_upper = [250.0,380.0,603.0,995.0,1512.0,200.0,281.0,793.0,1120.0] .* 1e9 * eV_to_nu
HMV_nu_up_upper = [380.0,603.0,995.0,1512.0,2397.0,281.0,398.0,1120.0,1583.0] .* 1e9 * eV_to_nu




function M87SEDPlot(;theme=DiplodocusLight())

    CairoMakie.activate!(inline=true)

    with_theme(theme) do 

    fig = Figure(size=(500,300))

    ax = Axis(fig[1,1],xlabel=L"\log_{10}\left(\nu\left[\text{Hz}\right]\right)",ylabel=L"\log_{10}\left(\nu F_\nu \left[\text{W}/\text{m}^2\right]\right)",xticks=8:4:28)

    markersize = 5


    # EHT 2017 data
    # VLBI 
    color = theme.colormap[][(1 - 1) / (6 - 1)]
    scatter!(ax,log10.(VLBI_nu),log10.(VLBI_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{Radio VLBI}")
    errorbars!(ax,log10.(VLBI_nu),log10.(VLBI_nuFnu * ergpercm2s_to_Wperm2),log10.((VLBI_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((VLBI_nuFnu .- VLBI_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((VLBI_nuFnu .+ VLBI_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((VLBI_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)

    # HST 
    color = theme.colormap[][(2 - 1) / (6 - 1)]
    scatter!(ax,log10.(HST_nu),log10.(HST_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{HST}")
    errorbars!(ax,log10.(HST_nu),log10.(HST_nuFnu * ergpercm2s_to_Wperm2),log10.((HST_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((HST_nuFnu .- HST_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((HST_nuFnu .+ HST_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((HST_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)

    # Swift
    color = theme.colormap[][(3 - 1) / (6 - 1)]
    scatter!(ax,log10.(Swift_nu),log10.(Swift_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{Swift}")
    errorbars!(ax,log10.(Swift_nu),log10.(Swift_nuFnu * ergpercm2s_to_Wperm2),log10.((Swift_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((Swift_nuFnu .- Swift_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((Swift_nuFnu .+ Swift_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((Swift_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)
    errorbars!(ax,log10.(Swift_nu_upper),log10.(Swift_nuFnu_upper * ergpercm2s_to_Wperm2),log10.(Swift_nu_upper) .- log10.(Swift_nu_low),log10.(Swift_nu_up) .- log10.(Swift_nu_upper),color=color,whiskerwidth = 0.0,direction=:x)
    scatter!(ax,log10.(Swift_nu_upper),log10.(Swift_nuFnu_upper * ergpercm2s_to_Wperm2),color=color,marker=:dtriangle,markersize=6)

    # Chandra+NuSTAR
    color = theme.colormap[][(4 - 1) / (6 - 1)]
    scatter!(ax,log10.(ChNu_nu),log10.(ChNu_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{Chandra+NuSTAR}")
    errorbars!(ax,log10.(ChNu_nu),log10.(ChNu_nuFnu * ergpercm2s_to_Wperm2),log10.((ChNu_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((ChNu_nuFnu .- ChNu_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((ChNu_nuFnu .+ ChNu_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((ChNu_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)
    errorbars!(ax,log10.(ChNu_nu),log10.(ChNu_nuFnu * ergpercm2s_to_Wperm2),log10.(ChNu_nu) .- log10.(ChNu_nu_low),log10.(ChNu_nu_up) .- log10.(ChNu_nu),color=color,whiskerwidth = 0.0,direction=:x)

    # Fermi
    color = theme.colormap[][(5 - 1) / (6 - 1)]
    scatter!(ax,log10.(Fermi_nu),log10.(Fermi_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{Fermi}")
    errorbars!(ax,log10.(Fermi_nu),log10.(Fermi_nuFnu * ergpercm2s_to_Wperm2),log10.((Fermi_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((Fermi_nuFnu .- Fermi_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((Fermi_nuFnu .+ Fermi_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((Fermi_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)
    errorbars!(ax,log10.(Fermi_nu),log10.(Fermi_nuFnu * ergpercm2s_to_Wperm2),log10.(Fermi_nu) .- log10.(Fermi_nu_low),log10.(Fermi_nu_up) .- log10.(Fermi_nu),color=color,whiskerwidth = 0.0,direction=:x)
    scatter!(ax,log10.(Fermi_nu_upper),log10.(Fermi_nuFnu_upper * ergpercm2s_to_Wperm2),color=color,marker=:dtriangle,markersize=6)
    errorbars!(ax,log10.(Fermi_nu_upper),log10.(Fermi_nuFnu_upper * ergpercm2s_to_Wperm2),log10.(Fermi_nu_upper) .- log10.(Fermi_nu_low_upper),log10.(Fermi_nu_up_upper) .- log10.(Fermi_nu_upper),color=color,whiskerwidth = 0.0,direction=:x)

    # HESS+Vertias+MAGIC
    color = theme.colormap[][(6 - 1) / (6 - 1)]
    scatter!(ax,log10.(HMV_nu),log10.(HMV_nuFnu * ergpercm2s_to_Wperm2),color=color,marker=:circle,markersize=markersize,label=L"\text{HESS+VERITAS+MAGIC}")
    errorbars!(ax,log10.(HMV_nu),log10.(HMV_nuFnu * ergpercm2s_to_Wperm2),log10.((HMV_nuFnu) * ergpercm2s_to_Wperm2) .- log10.((HMV_nuFnu .- HMV_nuFnu_sigma) * ergpercm2s_to_Wperm2),log10.((HMV_nuFnu .+ HMV_nuFnu_sigma) * ergpercm2s_to_Wperm2) .- log10.((HMV_nuFnu) * ergpercm2s_to_Wperm2),color=color,whiskerwidth = 0.0)
    errorbars!(ax,log10.(HMV_nu),log10.(HMV_nuFnu * ergpercm2s_to_Wperm2),log10.(HMV_nu) .- log10.(HMV_nu_low),log10.(HMV_nu_up) .- log10.(HMV_nu),color=color,whiskerwidth = 0.0,direction=:x)

    scatter!(ax,log10.(HMV_nu_upper),log10.(HMV_nuFnu_upper * ergpercm2s_to_Wperm2),color=color,marker=:dtriangle,markersize=6)
    errorbars!(ax,log10.(HMV_nu_upper),log10.(HMV_nuFnu_upper * ergpercm2s_to_Wperm2),log10.(HMV_nu_upper) .- log10.(HMV_nu_low_upper),log10.(HMV_nu_up_upper) .- log10.(HMV_nu_upper),color=color,whiskerwidth = 0.0,direction=:x)

    # my sims 
    scatterlines!(ax,log10.(nu),log10.(nuFnu1),marker=:circle,markersize=0.0,linewidth=2.0,color=theme.textcolor[],label=L"a=0.1",linestyle=:dash)
    scatterlines!(ax,log10.(nu),log10.(nuFnu5),marker=:circle,markersize=0.0,linewidth=2.0,color=theme.textcolor[],label=L"a=0.5",linestyle=:dot)
    scatterlines!(ax,log10.(nu),log10.(nuFnu9),marker=:circle,markersize=0.0,linewidth=2.0,color=theme.textcolor[],label=L"a=0.9",linestyle=:solid)

    axislegend(ax,position=:cb,framevisible=true,textsize=10,rowgap=-10,nbanks=3)

    ax.limits = (8,28,-20,-12)

    ax_t = Axis(fig[1,1],xlabel =L"\log_{10}\left(p\left[m_ec\right]\right)",xaxisposition=:top,xticks=-12:4:8)
    hidespines!(ax_t)
    hidexdecorations!(ax_t,ticklabels=false,ticks=false,label=false,minorgrid=false,minorticks=false)
    hideydecorations!(ax_t,ticklabels=false,ticks=false,minorgrid=false,minorticks=false)
    ax_t.limits = (log10(1e8 / (CONST_mEle*CONST_c^2/CONST_h)),log10(1e28 / (CONST_mEle*CONST_c^2/CONST_h)),-20,-12)
    scatterlines!(ax_t,log10.(nu ./ (CONST_mEle*CONST_c^2/CONST_h)),log10.(nuFnu9),color=:transparent,markersize=10.0)

    return fig

    end # theme
end

fig = M87SEDPlot(;theme=DiplodocusLight())
DiplodocusPlots.save("M87/M87_jet_SED.pdf",fig)

p0_list



