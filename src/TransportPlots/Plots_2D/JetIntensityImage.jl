#function JetIntensityImage(type::Static,PhaseSpace::PhaseSpaceStruct,sol::OutputStruct,time_idx::Int64,observerDistance::Float64,observerInclination::Float64,constantCoordinateIndex::Int64)

constantCoordinateIndex = 4 # index of the constant coordinate (z=0 plane)
observerDistance = 10000.0
observerInclination = 78 * π / 180 # inclination angle of the observer in radians (0 = looking down the z-axis, π/2 = looking along the x-y plane)
screen_rotate = pi/4
momenta_sample = [1e-5,1e0,1e5] # sample momenta for the intensity calculation
ϕsamples = 101 
u_sub_samples = 4
using DiplodocusPlots
using LinearAlgebra
using StaticArrays
using CairoMakie

fileName = "M87_test_10_Spin.jld2"
fileName = "M87_t10_10x10_med_a5.jld2"
fileLocation = joinpath(pwd(),"..","DiplodocusTransport","TestingGround","Data");
(PhaseSpace,sol) = DiplodocusPlots.SolutionFileLoad(fileLocation,fileName);

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

        x_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.UniformGrid(2pi,0.0,ϕsamples)))
        a = 1/sqrt(2)
        ratio = sqrt(2)
        nsub = 12#u_sub_samples
        cells = 10
        y_points_sim = DiplodocusPlots.DiplodocusTransport.meanVector(DiplodocusPlots.DiplodocusTransport.SpacetimeGrid(DiplodocusPlots.DiplodocusTransport.StretchGrid(0.0,a*ratio^((1-2^nsub)/2^nsub),ratio^(1/2^nsub),10*(2^nsub))))
        z_point_sim = mz_sim[constantCoordinateIndex]

        #### Observer Position and Orientation ####

        pos_obs = zeros(Float64,3)
        pos_obs[1] = observerDistance*sin(observerInclination)
        pos_obs[2] = 0.0
        pos_obs[3] = observerDistance*cos(observerInclination)

        screen_distance = 10.0 # distance from observer to image plane

        # form x-y basis of image plane perpendicular to the observer's line of sight
        screen_normvec = pos_obs/norm(pos_obs) # normal vector of the image plane
        screen_xvec = [0,0,1] - dot([0,0,1],screen_normvec)*screen_normvec # x basis vector of the image plane
        screen_xvec = screen_xvec/norm(screen_xvec) # normalize x basis vector
        screen_yvec = cross(screen_normvec,screen_xvec) # y basis vector of the image plane      

        @inline function normalToSurface!(pos_plot_normal::Vector{Float64},pos_sim::Vector{Float64})
            ϕ = pos_sim[1]
            u = pos_sim[2]
            v = pos_sim[3]
            norm = sqrt(u^2+v^2)
            pos_plot_normal[1] = u*cos(ϕ)/norm 
            pos_plot_normal[2] = u*sin(ϕ)/norm
            pos_plot_normal[3] = -v/norm
            return nothing
        end

        # calculate plotting points excluding ones that cannot be seen by the observer

        
        points_sim = Vector{Vector{Float64}}(undef, 0)
        points_plot = Vector{Vector{Float64}}(undef, 0)
        points_surface_to_obs = Vector{Vector{Float64}}(undef, 0)
        indices_sim = Vector{Vector{Int64}}(undef, 0)

        pos_plot = zeros(Float64,4)
        pos_plot_normal = zeros(Float64,3)
        for j in 1:length(y_points_sim), i in 1:length(x_points_sim),

            x_sim = x_points_sim[i]
            y_sim = y_points_sim[j]
            z_sim = z_point_sim
            
            DiplodocusPlots.DiplodocusTransport.CoordinateTransform!([1,x_sim,y_sim,z_sim],pos_plot,sim_coordinates,coordinates)

            x_plot = pos_plot[2]
            y_plot = pos_plot[3]
            z_plot = pos_plot[4]

            normalToSurface!(pos_plot_normal,[x_sim,y_sim,z_sim])

            pos_surface_to_obs = pos_obs - [x_plot,y_plot,z_plot]
            visable = dot(pos_surface_to_obs,pos_plot_normal) > 0.0 # check if the surface point is visible to the observer

            #println("j:$j, $(ceil(Int64,j/(cells))), $cells")
            
            if visable 
                push!(points_sim, [x_sim,y_sim,z_sim])
                push!(indices_sim, [1, ceil(Int64,j/(length(y_points_sim)/10)), constantCoordinateIndex])
                push!(points_plot, [x_plot,y_plot,z_plot])
                push!(points_surface_to_obs, pos_surface_to_obs)
            end
        end
  
    ########### Sample Image Surface and Compute Intensity ####

        points_x = Vector{Float64}(undef, 0) # x points on the image plane to be plotted
        points_y = Vector{Float64}(undef, 0) # y points on the image plane to be plotted
        points_I = Vector{Float64}(undef, 0) # Intensity points on the image plane to be plotted

        pvec_plot = zeros(Float64,4) # momentum vector of photon emitted to observer in plot frame
        pvec_sim = zeros(Float64,4) # momentum vector of photon emitted to observer in the sim frame
        pvec_tmp = zeros(Float64,4) # momentum vector of photon emitted to observer in the sim frame

        inve_sim = @MMatrix zeros(Float64,4,4) # inverse tetrad matrix at the surface point in the sim frame
        e_sim = @MMatrix zeros(Float64,4,4) # tetrad matrix at the surface point in the sim frame
        inve_plot = @MMatrix zeros(Float64,4,4) # tetrad matrix at the surface point in the plot frame

        pos_sim = @MVector zeros(Float64,4) # position of the surface point in the sim frame
        pos_plot = @MVector zeros(Float64,4) # position of the surface point in the plot frame

        metric = PhaseSpace.Spacetime.metric
        tetrad = PhaseSpace.Spacetime.tetrad

        for index in eachindex(points_plot)

            p0 = 1e-10 # example photon energy in observer frame

            # work out frequency and θ, ϕ angles of momentum vector for r_point_to_obs vector in terms of local stationary tetrad at the surface point
                r = points_surface_to_obs[index]
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

                #display(inve_plot)
                #display(inve_sim)
                #display(e_sim)

                #pvec_sim =  e_sim * inve_plot * pvec_plot # convert momentum vector from plot frame to sim frame

                #println("p0: $((inve_plot * pvec_plot)[1]), p0comps: $(norm((inve_plot * pvec_plot)[2:4]))")

                #println("p0vec_plot: $pvec_plot")
                #println("pos_plot: $pos_plot")
                CoordinateBasisTransform!(pvec_plot,pvec_tmp,pos_plot,coordinates,sim_coordinates) # convert momentum vector from plot frame to sim frame
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
                u_loc = DiplodocusPlots.DC.location(-1.0,1.0,9,u_sim,PhaseSpace.Momentum.py_grid_list[2])
                h_loc = DiplodocusPlots.DC.location(0.0,2*pi,4,h_sim,PhaseSpace.Momentum.pz_grid_list[2])
                # location of the surface point in the sim frame grid
                x = indices_sim[index][1]
                y = indices_sim[index][2]
                z = indices_sim[index][3]
                state_index = DiplodocusPlots.DiplodocusTransport.GlobalIndicesToStateIndex(PhaseSpace,x,y,z,p_loc,u_loc,h_loc,2) # TODO: add species index to this function call

                #println("x: $x, y: $y, z: $z, p_loc: $p_loc, u_loc: $u_loc, h_loc: $h_loc, p_sim: $p_sim, u_sim: $u_sim, h_sim: $h_sim")

                f = sol.f[end][state_index]

            # work out intensity emitted from surface point by f*ν^3 where ν is in observer frame

                Iν = f * p0^3 # intensity emitted from surface point in observer frame

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

        end

    ############ Plot Image Surface and Intensity ###############
        Iν_max = maximum(points_I)
        Iν_min = max(minimum(points_I),1e-40*Iν_max)

        #=fig = DiplodocusPlots.Figure()
        ax = DiplodocusPlots.Axis(fig[1,1];title="Jet Intensity Image",xlabel="x [Code Units]",ylabel="y [Code Units]",aspect=DiplodocusPlots.DataAspect())
        for point in eachindex(points_x)
            DiplodocusPlots.scatter!(ax,points_x[point],points_y[point],
                color=log10(points_I[point]),
                colormap=:plasma,
                colorrange=(log10(Iν_min),log10(Iν_max)),
                markersize=5.0)
        end=#
        using CairoMakie
        f,ax,sc = CairoMakie.scatter(points_x,points_y,
                color=log10.(points_I),
                colormap=:binary,
                colorrange=(log10(Iν_min),log10(Iν_max)),
                markersize=5.0)

        CairoMakie.hidedecorations!(ax)
        CairoMakie.aspect!(ax,DataAspect())
        
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