"""
    BuildFluxMatrices(PhaseSpace;debug_mode=false,Precision=Float32)

Function that builds the flux matrices associated with coordinate forces and regular forces. First the sparse arrays representing the fluxes are filled, then the flux matrices are returned as an immutable `FluxMatricesStruct`.

# Optional arguments
- `debug_mode::Bool=false`: If true, all flux matrices are built and returned. If false, only the essential flux matrices are built and returned to save memory.
"""
function BuildFluxMatrices(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce};debug_mode::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    n_space = x_num*y_num*z_num

    n = n_momentum*n_space

    Ap_Flux_I::Vector{Int32} = Int32[]
    Ap_Flux_J::Vector{Int32} = Int32[]
    Ap_Flux_V::Vector{Precision} = Precision[]
    Am_Flux_I::Vector{Int32} = Int32[]
    Am_Flux_J::Vector{Int32} = Int32[]
    Am_Flux_V::Vector{Precision} = Precision[]

    B_Flux_I::Vector{Int32} = Int32[]
    B_Flux_J::Vector{Int32} = Int32[]
    B_Flux_V::Vector{Precision} = Precision[]

    C_Flux_I::Vector{Int32} = Int32[]
    C_Flux_J::Vector{Int32} = Int32[]
    C_Flux_V::Vector{Precision} = Precision[]

    D_Flux_I::Vector{Int32} = Int32[]
    D_Flux_J::Vector{Int32} = Int32[]
    D_Flux_V::Vector{Precision} = Precision[]

    I_Flux_I::Vector{Int32} = Int32[]
    I_Flux_J::Vector{Int32} = Int32[]
    I_Flux_V::Vector{Precision} = Precision[]

    J_Flux_I::Vector{Int32} = Int32[]
    J_Flux_J::Vector{Int32} = Int32[]
    J_Flux_V::Vector{Precision} = Precision[]

    K_Flux_I::Vector{Int32} = Int32[]
    K_Flux_J::Vector{Int32} = Int32[]
    K_Flux_V::Vector{Precision} = Precision[]

    Vol::Vector{Precision} = zeros(Precision,n_space)

    X_Flux_I::Vector{Int32} = Int32[]
    X_Flux_J::Vector{Int32} = Int32[]
    X_Flux_V::Vector{Precision} = Precision[]

    P_Flux_I::Vector{Int32} = Int32[]
    P_Flux_J::Vector{Int32} = Int32[]
    P_Flux_V::Vector{Precision} = Precision[]

    if debug_mode

        println("Building time flux matrices...")
        FillTimeFlux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)
        Ap_Flux = sparse(Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        Am_Flux= sparse(Am_Flux_I,Am_Flux_J,Am_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32} 

        println("Building space flux matrices...")
        FillSpaceFlux!(PhaseSpace,B_Flux_I,B_Flux_J,B_Flux_V,C_Flux_I,C_Flux_J,C_Flux_V,D_Flux_I,D_Flux_J,D_Flux_V)
        B_Flux = sparse(B_Flux_I,B_Flux_J,B_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        C_Flux = sparse(C_Flux_I,C_Flux_J,C_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        D_Flux = sparse(D_Flux_I,D_Flux_J,D_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}

        println("Building momentum flux matrices...")
        FillMomentumFlux!(PhaseSpace,Forces,I_Flux_I,I_Flux_J,I_Flux_V,J_Flux_I,J_Flux_J,J_Flux_V,K_Flux_I,K_Flux_J,K_Flux_V)
        I_Flux = sparse(I_Flux_I,I_Flux_J,I_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        J_Flux = sparse(J_Flux_I,J_Flux_J,J_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        K_Flux = sparse(K_Flux_I,K_Flux_J,K_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}

        # when values cancel they can still leave a saved zero in arrays
        dropzeros!(B_Flux)
        dropzeros!(C_Flux)
        dropzeros!(D_Flux)
        dropzeros!(I_Flux)
        dropzeros!(J_Flux)
        dropzeros!(K_Flux)
        
        #=Fill_A_Flux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)
        Ap_Flux = sparse(Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}
        Am_Flux= sparse(Am_Flux_I,Am_Flux_J,Am_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64} 

        # Fill flux matrices
        println("Building B flux matrix...")
        Fill_B_Flux!(PhaseSpace,B_Flux_I,B_Flux_J,B_Flux_V)
        B_Flux = sparse(B_Flux_I,B_Flux_J,B_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building C flux matrix...")
        Fill_C_Flux!(PhaseSpace,C_Flux_I,C_Flux_J,C_Flux_V)
        C_Flux = sparse(C_Flux_I,C_Flux_J,C_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building D flux matrix...")
        Fill_D_Flux!(PhaseSpace,D_Flux_I,D_Flux_J,D_Flux_V)
        D_Flux = sparse(D_Flux_I,D_Flux_J,D_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building I flux matrix...")
        Fill_I_Flux!(PhaseSpace,Forces,I_Flux_I,I_Flux_J,I_Flux_V)
        I_Flux = sparse(I_Flux_I,I_Flux_J,I_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building J flux matrix...")
        Fill_J_Flux!(PhaseSpace,Forces,J_Flux_I,J_Flux_J,J_Flux_V)
        J_Flux = sparse(J_Flux_I,J_Flux_J,J_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building K flux matrix...")
        Fill_K_Flux!(PhaseSpace,Forces,K_Flux_I,K_Flux_J,K_Flux_V)
        K_Flux = sparse(K_Flux_I,K_Flux_J,K_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}=#
        
        println("Calculating volume elements...")
        FillVolume!(Vol,PhaseSpace)

        println("Building X and P flux matrix...")
        P_Flux = I_Flux + J_Flux + K_Flux
        X_Flux = B_Flux + C_Flux + D_Flux 
        
    else

        println("Building time flux matrices...")
        FillTimeFlux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)
        Ap_Flux = sparse(Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}
        Am_Flux= sparse(Am_Flux_I,Am_Flux_J,Am_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32} 

        println("Building space flux matrices...")
        FillSpaceFlux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V)
        X_Flux = sparse(X_Flux_I,X_Flux_J,X_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}

        println("Building momentum flux matrices...")
        FillMomentumFlux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V)
        P_Flux = sparse(P_Flux_I,P_Flux_J,P_Flux_V,n,n)::SparseMatrixCSC{Precision,Int32}

        # when values cancel they can still leave a saved zero in arrays
        # Remove any remaining small values
        X_Flux.nzval[abs.(X_Flux.nzval) .< eps(eltype(X_Flux.nzval))] .= zero(eltype(X_Flux.nzval))
        P_Flux.nzval[abs.(P_Flux.nzval) .< eps(eltype(P_Flux.nzval))] .= zero(eltype(P_Flux.nzval))
        dropzeros!(X_Flux)
        dropzeros!(P_Flux)
        
        #=println("Building time flux matrices...")
        Fill_A_Flux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)
        Ap_Flux = sparse(Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}
        Am_Flux = sparse(Am_Flux_I,Am_Flux_J,Am_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building space flux matrices...")
        Fill_B_Flux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V)
        Fill_C_Flux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V)
        Fill_D_Flux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V)
        X_Flux = sparse(X_Flux_I,X_Flux_J,X_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}

        println("Building momentum flux matrices...")
        Fill_I_Flux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V)
        Fill_J_Flux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V)
        Fill_K_Flux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V)
        P_Flux = sparse(P_Flux_I,P_Flux_J,P_Flux_V,n,n)::SparseMatrixCSC{Precision,Int64}=#

        println("Calculating volume elements...")
        FillVolume!(Vol,PhaseSpace)

        # space fluxes
        B_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32}
        C_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32}
        D_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32}
        # momentum fluxes
        I_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32}
        J_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32}
        K_Flux = spzeros(Precision,Int32,0,0)::SparseMatrixCSC{Precision,Int32} 
    end

    if isdiag(Ap_Flux)
        Ap_Flux = convert(Vector{Precision},diag(Ap_Flux))::Vector{Precision}
    else
        error("Aplus flux is not diagonal, which should be the case for stationary spacetimes")
    end
    if isdiag(Am_Flux)
        Am_Flux = convert(Vector{Precision},diag(Am_Flux))::Vector{Precision}
    else
        error("Aminus flux is not diagonal, which should be the case for stationary spacetimes")
    end

    size = round(Base.summarysize(X_Flux),sigdigits=3)

    if size > 1e9
        println("X_Flux is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("X_Flux is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("X_Flux is approx. $(size/1e3) KB in memory")
    else
        println("X_Flux is approx. $size bytes in memory")
    end


    size = round(Base.summarysize(P_Flux),sigdigits=3)

    if size > 1e9
        println("P_Flux is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("P_Flux is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("P_Flux is approx. $(size/1e3) KB in memory")
    else
        println("P_Flux is approx. $size bytes in memory")
    end


    return FluxMatricesStruct{Precision}(Ap_Flux,Am_Flux,X_Flux,P_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

end

"""
    Fill_A_Flux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `Ap_Flux` and `Am_Flux` matrices.
"""
function Fill_A_Flux!(PhaseSpace::PhaseSpaceStruct,Ap_Flux_I::Vector{Int32},Ap_Flux_J::Vector{Int32},Ap_Flux_V::Vector{T},Am_Flux_I::Vector{Int32},Am_Flux_J::Vector{Int32},Am_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
                b = a

                # integration sign introduced here
                A_plus = AFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                A_minus = -AFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                # normalisation
                norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # fill
                push!(Ap_Flux_I,Int32(a))
                push!(Ap_Flux_J,Int32(b))
                push!(Ap_Flux_V,convert(T,A_plus / norm))

                push!(Am_Flux_I,Int32(a))
                push!(Am_Flux_J,Int32(b))
                push!(Am_Flux_V,convert(T,A_minus / norm))

            end
        end
    end
end

"""
    FillTimeFlux!(PhaseSpace,Ap_Flux_I,Ap_Flux_J,Ap_Flux_V,Am_Flux_I,Am_Flux_J,Am_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `Ap_Flux` and `Am_Flux` matrices.
"""
function FillTimeFlux!(PhaseSpace::PhaseSpaceStruct,Ap_Flux_I::Vector{Int32},Ap_Flux_J::Vector{Int32},Ap_Flux_V::Vector{T},Am_Flux_I::Vector{Int32},Am_Flux_J::Vector{Int32},Am_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    metric = Spacetime.metric
    coordinates = Spacetime.coordinates
    tetrad = Spacetime.tetrad
    momentum_coords = Momentum.coordinates

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    CFSpaceVectorPlus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceVectorMinus = MVector{4,Float64}(zeros(Float64,4))
    CFMomentumVector = MVector{4,Float64}(zeros(Float64,4))

    @inline CoordinateFluxAIntegrand!(xyzt,A) = CoordinateFluxSpaceAIntegrand!(xyzt,A,metric,coordinates,tetrad)
    
    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        t0 = Grids.tr[1]
        t1 = Grids.tr[2]
        x0 = Grids.xr[x]
        x1 = Grids.xr[x+1]
        y0 = Grids.yr[y]
        y1 = Grids.yr[y+1]
        z0 = Grids.zr[z]
        z1 = Grids.zr[z+1]

        xyzt0::SVector{4,Float64} = [x0,y0,z0,t0]
        xyzt1::SVector{4,Float64} = [x1,y1,z1,t1]
        n::SVector{3,Int64} = [8,8,8] # number of points for integration, can be changed to increase accuracy

        # Integrate space part of force
        FluxSimpson3D!(CoordinateFluxAIntegrand!,CFSpaceVectorPlus,CFSpaceVectorMinus,xyzt0,xyzt1,n)

        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                CoordinateFluxMomentumIntegral!(CFMomentumVector,PhaseSpace,name,px,py,pz)

                Aplus = 0.0
                Aminus = 0.0

                # integration sign introduced here
                for a in 1:4
                    Aplus += CFSpaceVectorPlus[a]*CFMomentumVector[a]
                    Aminus -= CFSpaceVectorMinus[a]*CFMomentumVector[a]
                end

                a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
                b = a

                Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # fill
                push!(Ap_Flux_I,Int32(a))
                push!(Ap_Flux_J,Int32(b))
                push!(Ap_Flux_V,convert(T,Aplus / Norm))

                push!(Am_Flux_I,Int32(a))
                push!(Am_Flux_J,Int32(b))
                push!(Am_Flux_V,convert(T,Aminus / Norm))
            end
        end
    end
end

"""
    FillVolume!(Vol,PhaseSpace)

Generates volume elements for each space sub-domain.
"""
function FillVolume!(Vol::Vector{T},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        space = (x-1)*y_num*z_num+(y-1)*z_num+z

        Vol[space] = convert(T,VolFunction(PhaseSpace,1,x,y,z))

    end

end

function SchemeCoefficients(scheme::String,flux_plus::T,flux_minus::T) where T<:Union{Float32,Float64}

    if scheme == "upwind"
        if sign(flux_plus) == 1
            h_plus_right = 0
            h_plus_left = 1
        else
            h_plus_right = 1
            h_plus_left = 0
        end
        if sign(flux_minus) == 1
            h_minus_right = 1
            h_minus_left = 0
        else
            h_minus_right = 0
            h_minus_left = 1
        end
    elseif scheme == "central"
        h_plus_right = 0.5
        h_plus_left = 0.5
        h_minus_right = 0.5
        h_minus_left = 0.5
    else
        error("Unknown scheme")
    end

    return h_plus_right, h_plus_left, h_minus_right, h_minus_left

end
function RightBound(p::Int64,num::Int64,BC::AbstractBoundaryCondition) 

    if p+1 > num # right boundary
        if BC isa Periodic
            return 1
        else
            return num
        end
    else
        return p+1
    end

end
function LeftBound(p::Int64,num::Int64,BC::AbstractBoundaryCondition)

    if p-1 < 1 # left boundary
        if BC isa Periodic
            return num
        else
            return 1
        end
    else
        return p-1
    end

end

function AssignFlux!(PhaseSpace::PhaseSpaceStruct,flux::String,Flux_I::Vector{Int32},Flux_J::Vector{Int32},Flux_V::Vector{T},Flux_plus::Float64,Flux_minus::Float64,x::Int64,y::Int64,z::Int64,px::Int64,py::Int64,pz::Int64,name::Int64,scheme::String,BCp::AbstractBoundaryCondition,BCm::AbstractBoundaryCondition) where T<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    x_num = PhaseSpace.Spacetime.x_num
    y_num = PhaseSpace.Spacetime.y_num
    z_num = PhaseSpace.Spacetime.z_num
    px_num = PhaseSpace.Momentum.px_num_list[name]
    py_num = PhaseSpace.Momentum.py_num_list[name]
    pz_num = PhaseSpace.Momentum.pz_num_list[name]

    a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
    b = a 
    Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
    if flux == "B"
        xp = RightBound(x,x_num,BCp)
        xm = LeftBound(x,x_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,xp,y,z,px,py,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,xm,y,z,px,py,pz,name)
        Norm_p = Norm
        Norm_m = Norm
    elseif flux == "C"
        yp = RightBound(y,y_num,BCp)
        ym = LeftBound(y,y_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,x,yp,z,px,py,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,x,ym,z,px,py,pz,name)
        Norm_p = Norm
        Norm_m = Norm
    elseif flux == "D"
        zp = RightBound(z,z_num,BCp)
        zm = LeftBound(z,z_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,zp,px,py,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,zm,px,py,pz,name)
        Norm_p = Norm
        Norm_m = Norm
    elseif flux == "I"
        pxp = RightBound(px,px_num,BCp)
        pxm = LeftBound(px,px_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,pxp,py,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,pxm,py,pz,name)
        Norm_p = MomentumSpaceNorm(Grids,name,pxp,py,pz)
        Norm_m = MomentumSpaceNorm(Grids,name,pxm,py,pz)
    elseif flux == "J"
        pyp = RightBound(py,py_num,BCp)
        pym = LeftBound(py,py_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,pyp,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,pym,pz,name)
        Norm_p = MomentumSpaceNorm(Grids,name,px,pyp,pz)
        Norm_m = MomentumSpaceNorm(Grids,name,px,pym,pz)
    elseif flux == "K"
        pzp = RightBound(pz,pz_num,BCp)
        pzm = LeftBound(pz,pz_num,BCm)
        bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pzp,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pzm,name)
        Norm_p = MomentumSpaceNorm(Grids,name,px,py,pzp)
        Norm_m = MomentumSpaceNorm(Grids,name,px,py,pzm)
    end

    (h_plus_right, h_plus_left, h_minus_right, h_minus_left) = SchemeCoefficients(scheme,Flux_plus,Flux_minus)

    #=
    ________________________________
    a |  F_m   | F_m+F_p | F_p    |
    __|________|_________|________|_
          b-1       b        b+1  
    =#

    #= note on boundaries:
     b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
     b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)
    =#

    # normalised fluxes
    if b != bp
        push!(Flux_I,a)
        push!(Flux_J,bp)
        push!(Flux_V,convert(T,(Flux_plus * h_plus_right) / Norm_p))
        #I_Flux[a,bp] += convert(T,(I_plus * h_plus_right) / Mom_Normp)
        push!(Flux_I,a)
        push!(Flux_J,b)
        push!(Flux_V,convert(T,(Flux_plus * h_plus_left) / Norm))
        #I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
    elseif BCp isa Open # b=bp
        push!(Flux_I,a)
        push!(Flux_J,b)
        push!(Flux_V,convert(T,(Flux_plus * h_plus_left) / Norm))
        # I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
    end
    if b != bm
        push!(Flux_I,a)
        push!(Flux_J,b)
        push!(Flux_V,convert(T,(Flux_minus * h_minus_right) / Norm))
        #I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm)
        push!(Flux_I,a)
        push!(Flux_J,bm)
        push!(Flux_V,convert(T,(Flux_minus * h_minus_left) / Norm_m)) 
        #I_Flux[a,bm] += convert(T,(I_minus * h_minus_left) / Mom_Normm) 
    elseif BCm isa Open # b=bm
        push!(Flux_I,a)
        push!(Flux_J,b)
        push!(Flux_V,convert(T,(Flux_minus * h_minus_right) / Norm))
        #I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
    end

end

"""
    FillMomentumFlux!(PhaseSpace,Forces,I_Flux_I,I_Flux_J,I_Flux_V,J_Flux_I,J_Flux_J,J_Flux_V,K_Flux_I,K_Flux_J,K_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `I_Flux`, `J_Flux`, and `K_Flux` matrices.
"""
#=function FillMomentumFlux!(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce},I_Flux_I::Vector{Int32},I_Flux_J::Vector{Int32},I_Flux_V::Vector{T},J_Flux_I::Vector{Int32},J_Flux_J::Vector{Int32},J_Flux_V::Vector{T},K_Flux_I::Vector{Int32},K_Flux_J::Vector{Int32},K_Flux_V::Vector{T}) where T<:Union{Float32,Float64}
    
    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    metric = Spacetime.metric
    coordinates = Spacetime.coordinates
    tetrad = Spacetime.tetrad
    momentum_coords = Momentum.coordinates
    scheme = Momentum.scheme

    BC_pxp = momentum_coords.xp_BC
    BC_pxm = momentum_coords.xm_BC 
    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert BC_pxp isa Closed || BC_pxm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

    BC_pyp = momentum_coords.yp_BC
    BC_pym = momentum_coords.ym_BC
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert (BC_pyp isa Closed) && (BC_pym isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

    BC_pzp = momentum_coords.zp_BC
    BC_pzm = momentum_coords.zm_BC
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    @assert (BC_pzp isa Periodic) && (BC_pzm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    # setup for coordinate forces
    pos = MVector{4,Float64}(zeros(Float64,4))
    e = MMatrix{4,4,Float64,16}(zeros(Float64,4,4))
    inve = MMatrix{4,4,Float64,16}(zeros(Float64,4,4))
    Γ = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    ∂inve = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    CFSpaceArray = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    CFMomentumArray = MArray{Tuple{4,4,4,3,2},Float64,5,384}(zeros(Float64,4,4,4,3,2))

    @inline inv_tetrad_func!(inve_local,pos_local) = InverseTetradComponents!(pos_local,inve_local,metric,coordinates,tetrad)
    cfg::ForwardDiff.JacobianConfig = ForwardDiff.JacobianConfig(inv_tetrad_func!,inve,pos)
    @inline CoordinateForceSpaceIntegrand_func!(pos_local,CFSpaceArray_local) = CoordinateForceSpaceIntegrand!(pos_local,CFSpaceArray_local,metric,coordinates,tetrad,e,inve,Γ,∂inve,inv_tetrad_func!,cfg)

    non_dimensional_factor::Float64 = CONST_c * Characteristic.CHAR_time / Characteristic.CHAR_length

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        println("Calculating momentum fluxes for space cell ($x,$y,$z)")
    
        if CoordinateForce() in Forces

            t0 = Grids.tr[1]
            t1 = Grids.tr[2]
            x0 = Grids.xr[x]
            x1 = Grids.xr[x+1]
            y0 = Grids.yr[y]
            y1 = Grids.yr[y+1]
            z0 = Grids.zr[z]
            z1 = Grids.zr[z+1]

            a::SVector{4,Float64} = [t0,x0,y0,z0]
            b::SVector{4,Float64} = [t1,x1,y1,z1]
            n::SVector{4,Int64} = [2,16,16,16] # number of points for integration, can be changed to increase accuracy

            # Integrate space part of force
            Simpson4D!(CoordinateForceSpaceIntegrand_func!,CFSpaceArray,a,b,n)

            CFSpaceArray .*= non_dimensional_factor

        end
            
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]
    
            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                println("Calculating momentum fluxes for momentum cell ($px,$py,$pz)")

                I_plus = 0.0
                I_minus = 0.0
                J_plus = 0.0
                J_minus = 0.0
                K_plus = 0.0
                K_minus = 0.0

                for f in 1:length(Forces)

                    if Forces[f] isa CoordinateForce

                        CoordinateForceMomentumIntegral!(CFMomentumArray,PhaseSpace,name,px,py,pz)
                        # Integrate momentum part of force and calculate product

                        for a in 1:4, b in 1:4, c in 1:4

                            if isnan(CFMomentumArray[a,b,c,2,1])
                                println("$a,$b,$c,$px,$py,$pz")
                            end

                            # integration sign introduced here
                            I_plus += CFMomentumArray[a,b,c,1,1] * CFSpaceArray[a,b,c]
                            I_minus -= CFMomentumArray[a,b,c,1,2] * CFSpaceArray[a,b,c]
                            J_plus += CFMomentumArray[a,b,c,2,1] * CFSpaceArray[a,b,c]
                            J_minus -= CFMomentumArray[a,b,c,2,2] * CFSpaceArray[a,b,c]
                            K_plus += CFMomentumArray[a,b,c,3,1] * CFSpaceArray[a,b,c]
                            K_minus -= CFMomentumArray[a,b,c,3,2] * CFSpaceArray[a,b,c]

                        end
                        
                    else 

                        # integration sign introduced here
                        I_plus += IFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                        I_minus -= IFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz) 
                        J_plus += JFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                        J_minus -= JFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
                        K_plus += KFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                        K_minus -= KFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                    end

                end # force loop

                AssignFlux!(PhaseSpace,"I",I_Flux_I,I_Flux_J,I_Flux_V,I_plus,I_minus,x,y,z,px,py,pz,name,scheme,BC_pxp,BC_pxm)
                AssignFlux!(PhaseSpace,"J",J_Flux_I,J_Flux_J,J_Flux_V,J_plus,J_minus,x,y,z,px,py,pz,name,scheme,BC_pyp,BC_pym)
                AssignFlux!(PhaseSpace,"K",K_Flux_I,K_Flux_J,K_Flux_V,K_plus,K_minus,x,y,z,px,py,pz,name,scheme,BC_pzp,BC_pzm)

            end # momentum  loop
        end # name loop
    end # space loop
end=#


function FillMomentumFlux!(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce},I_Flux_I::Vector{Int32},I_Flux_J::Vector{Int32},I_Flux_V::Vector{T},J_Flux_I::Vector{Int32},J_Flux_J::Vector{Int32},J_Flux_V::Vector{T},K_Flux_I::Vector{Int32},K_Flux_J::Vector{Int32},K_Flux_V::Vector{T}) where T<:Union{Float32,Float64}
    
    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    metric = Spacetime.metric
    coordinates = Spacetime.coordinates
    tetrad = Spacetime.tetrad
    momentum_coords = Momentum.coordinates
    scheme = Momentum.scheme

    BC_pxp = momentum_coords.xp_BC
    BC_pxm = momentum_coords.xm_BC 
    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert BC_pxp isa Closed || BC_pxm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

    BC_pyp = momentum_coords.yp_BC
    BC_pym = momentum_coords.ym_BC
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert (BC_pyp isa Closed) && (BC_pym isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

    BC_pzp = momentum_coords.zp_BC
    BC_pzm = momentum_coords.zm_BC
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    @assert (BC_pzp isa Periodic) && (BC_pzm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    # setup for coordinate forces
    pos = MVector{4,Float64}(zeros(Float64,4))
    e = MMatrix{4,4,Float64,16}(zeros(Float64,4,4))
    inve = MMatrix{4,4,Float64,16}(zeros(Float64,4,4))
    Γ = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    ∂inve = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    CFSpaceArray = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    CFMomentumArray = MArray{Tuple{4,4,4,3,2},Float64,5,384}(zeros(Float64,4,4,4,3,2))
    @inline inv_tetrad_func!(inve_local,pos_local) = InverseTetradComponents!(pos_local,inve_local,metric,coordinates,tetrad)
    cfg::ForwardDiff.JacobianConfig = ForwardDiff.JacobianConfig(inv_tetrad_func!,inve,pos)
    @inline CoordinateForceSpaceIntegrand_func!(pos_local,CFSpaceArray_local) = CoordinateForceSpaceIntegrand!(pos_local,CFSpaceArray_local,metric,coordinates,tetrad,e,inve,Γ,∂inve,inv_tetrad_func!,cfg)

    # setup for space scalar forces
    MScaForce = MMatrix{3,2,Float64}(zeros(Float64,3,2))

    # setup for space vector forces
    SVecForce = MVector{4,Float64}(zeros(Float64,4))
    MVecForce = MArray{Tuple{3,4,2},Float64,3,24}(zeros(Float64,3,4,2))

    # setup for space matrix forces
    SMatForce = MMatrix{4,4,Float64}(zeros(Float64,4,4))
    MMatForce = MArray{Tuple{3,4,4,2},Float64,4,96}(zeros(Float64,3,4,4,2))

    # setup for space 3 tensor forces
    S3TenForce = MArray{Tuple{4,4,4},Float64,3,64}(zeros(Float64,4,4,4))
    M3TenForce = MArray{Tuple{3,4,4,4,2},Float64,5,384}(zeros(Float64,3,4,4,4,2))

    # non-dimensional normalisation
    non_dimensional_factor::Float64 = CONST_c * Characteristic.CHAR_time / Characteristic.CHAR_length

    for name in eachindex(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        I_plus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)
        I_minus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)
        J_plus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)
        J_minus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)
        K_plus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)
        K_minus_array = zeros(Float64,x_num,y_num,z_num,px_num,py_num,pz_num)

        for f in eachindex(Forces)

            println("Assigning force flux for force $(Forces[f]) and species $(name_list[name])")

            if Forces[f] isa CoordinateForce

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    if hasfield(typeof(Forces[f]), :Domain)
                        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
                        if !isnothing(Forces[f].Domain) && !(off_space in Forces[f].Domain)
                            continue
                        end
                    end

                    t0 = Grids.tr[1]
                    t1 = Grids.tr[2]
                    x0 = Grids.xr[x]
                    x1 = Grids.xr[x+1]
                    y0 = Grids.yr[y]
                    y1 = Grids.yr[y+1]
                    z0 = Grids.zr[z]
                    z1 = Grids.zr[z+1]

                    a::SVector{4,Float64} = [t0,x0,y0,z0]
                    b::SVector{4,Float64} = [t1,x1,y1,z1]
                    n::SVector{4,Int64} = [2,16,16,16] # number of points for integration, can be changed to increase accuracy

                    # Integrate space part of coordinate force
                    Simpson4D!(CoordinateForceSpaceIntegrand_func!,CFSpaceArray,a,b,n)

                    # put a tolerance on integration to avoid small numerical errors causing a non-zero flux when the flux should be zero
                    pos[1] = (t1+t0)/2
                    pos[2] = (x1+x0)/2
                    pos[3] = (y1+y0)/2
                    pos[4] = (z1+z0)/2
                    χ = VolumeElement(pos,metric,coordinates)
                    tolerance::T = 1e-8
                    for idx in eachindex(CFSpaceArray)
                        if abs(CFSpaceArray[idx]) < tolerance * χ
                            CFSpaceArray[idx] = zero(T)
                        end
                    end

                    CFSpaceArray .*= non_dimensional_factor

                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                        # Integrate momentum part of coordinate force
                        CoordinateForceMomentumIntegral!(CFMomentumArray,PhaseSpace,name,px,py,pz)

                        for a in 1:4, b in 1:4, c in 1:4
                            # integration sign introduced here
                            I_plus_array[x,y,z,px,py,pz] += CFMomentumArray[a,b,c,1,1] * CFSpaceArray[a,b,c]
                            I_minus_array[x,y,z,px,py,pz] -= CFMomentumArray[a,b,c,1,2] * CFSpaceArray[a,b,c]
                            J_plus_array[x,y,z,px,py,pz] += CFMomentumArray[a,b,c,2,1] * CFSpaceArray[a,b,c]
                            J_minus_array[x,y,z,px,py,pz] -= CFMomentumArray[a,b,c,2,2] * CFSpaceArray[a,b,c]
                            K_plus_array[x,y,z,px,py,pz] += CFMomentumArray[a,b,c,3,1] * CFSpaceArray[a,b,c]
                            K_minus_array[x,y,z,px,py,pz] -= CFMomentumArray[a,b,c,3,2] * CFSpaceArray[a,b,c]
                            
                        end
                    end
                end

            elseif Forces[f] isa SpaceScalarForce

                for x in 1:x_num, y in 1:y_num, z in 1:z_num
                    SScaForce = SpaceScalarForceFunction(Forces[f],PhaseSpace,x,y,z)
                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num
                        MomentumScalarForceFunction!(Forces[f],MScaForce,PhaseSpace,name,px,py,pz)
                        # M[flux_dir,pm]
                        # integration sign introduced here
                        I_plus_array[x,y,z,px,py,pz] += SScaForce * MScaForce[1,1]
                        I_minus_array[x,y,z,px,py,pz] -= SScaForce * MScaForce[1,2]
                        J_plus_array[x,y,z,px,py,pz] += SScaForce * MScaForce[2,1]
                        J_minus_array[x,y,z,px,py,pz] -= SScaForce * MScaForce[2,2]
                        K_plus_array[x,y,z,px,py,pz] += SScaForce * MScaForce[3,1]
                        K_minus_array[x,y,z,px,py,pz] -= SScaForce * MScaForce[3,2]
                    end
                end

            elseif Forces[f] isa SpaceVectorForce

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    if hasfield(typeof(Forces[f]), :Domain)
                        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
                        if !isnothing(Forces[f].Domain) && !(off_space in Forces[f].Domain)
                            continue
                        end
                    end
                    
                    SpaceVectorForceFunction!(Forces[f],SVecForce,PhaseSpace,x,y,z)
                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num
                        MomentumVectorForceFunction!(Forces[f],MVecForce,PhaseSpace,name,px,py,pz)
                        # M[flux_dir,a,pm]
                        # integration sign introduced here
                        @inbounds for a in 1:4
                            I_plus_array[x,y,z,px,py,pz] += SVecForce[a] * MVecForce[1,a,1]
                            I_minus_array[x,y,z,px,py,pz] -= SVecForce[a] * MVecForce[1,a,2]
                            J_plus_array[x,y,z,px,py,pz] += SVecForce[a] * MVecForce[2,a,1]
                            J_minus_array[x,y,z,px,py,pz] -= SVecForce[a] * MVecForce[2,a,2]
                            K_plus_array[x,y,z,px,py,pz] += SVecForce[a] * MVecForce[3,a,1]
                            K_minus_array[x,y,z,px,py,pz] -= SVecForce[a] * MVecForce[3,a,2]
                        end
                    end
                end

            elseif Forces[f] isa SpaceMatrixForce

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    if hasfield(typeof(Forces[f]), :Domain)
                        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
                        if !isnothing(Forces[f].Domain) && !(off_space in Forces[f].Domain)
                            continue
                        end
                    end
                    
                    SpaceMatrixForceFunction!(Forces[f],SMatForce,PhaseSpace,x,y,z)
                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num
                        MomentumMatrixForceFunction!(Forces[f],MMatForce,PhaseSpace,name,px,py,pz)
                        # M[flux_dir,a,pm]
                        # integration sign introduced here
                        @inbounds for a in 1:4, b in 1:4
                            I_plus_array[x,y,z,px,py,pz] += SMatForce[a,b] * MMatForce[1,a,b,1]
                            I_minus_array[x,y,z,px,py,pz] -= SMatForce[a,b] * MMatForce[1,a,b,2]
                            J_plus_array[x,y,z,px,py,pz] += SMatForce[a,b] * MMatForce[2,a,b,1]
                            J_minus_array[x,y,z,px,py,pz] -= SMatForce[a,b] * MMatForce[2,a,b,2]
                            K_plus_array[x,y,z,px,py,pz] += SMatForce[a,b] * MMatForce[3,a,b,1]
                            K_minus_array[x,y,z,px,py,pz] -= SMatForce[a,b] * MMatForce[3,a,b,2]
                        end
                    end
                end

            elseif Forces[f] isa Space3TensorForce

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    if hasfield(typeof(Forces[f]), :Domain)
                        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
                        if !isnothing(Forces[f].Domain) && !(off_space in Forces[f].Domain)
                            continue
                        end
                    end
                    
                    Space3TensorForceFunction!(Forces[f],S3TenForce,PhaseSpace,x,y,z)
                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num
                        Momentum3TensorForceFunction!(Forces[f],M3TenForce,PhaseSpace,name,px,py,pz)
                        # M[flux_dir,a,pm]
                        # integration sign introduced here
                        @inbounds for a in 1:4, b in 1:4, c in 1:4
                            I_plus_array[x,y,z,px,py,pz] += S3TenForce[a,b,c] * M3TenForce[1,a,b,c,1]
                            I_minus_array[x,y,z,px,py,pz] -= S3TenForce[a,b,c] * M3TenForce[1,a,b,c,2]
                            J_plus_array[x,y,z,px,py,pz] += S3TenForce[a,b,c] * M3TenForce[2,a,b,c,1]
                            J_minus_array[x,y,z,px,py,pz] -= S3TenForce[a,b,c] * M3TenForce[2,a,b,c,2]
                            K_plus_array[x,y,z,px,py,pz] += S3TenForce[a,b,c] * M3TenForce[3,a,b,c,1]
                            K_minus_array[x,y,z,px,py,pz] -= S3TenForce[a,b,c] * M3TenForce[3,a,b,c,2]
                        end
                    end
                end

            elseif Forces[f] isa AnalyticForce 

                for x in 1:x_num, y in 1:y_num, z in 1:z_num

                    if hasfield(typeof(Forces[f]), :Domain)
                        off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1
                        if !isnothing(Forces[f].Domain) && !(off_space in Forces[f].Domain)
                            continue
                        end
                    end
                    
                    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num
                    # integration sign introduced here
                    I_plus_array[x,y,z,px,py,pz] += IFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    I_minus_array[x,y,z,px,py,pz] -= IFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz) 
                    J_plus_array[x,y,z,px,py,pz] += JFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    J_minus_array[x,y,z,px,py,pz] -= JFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
                    K_plus_array[x,y,z,px,py,pz] += KFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    K_minus_array[x,y,z,px,py,pz] -= KFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
                    end
                end

            else
                error("Unknown force type of $(Forces[f])")
            end

        end # force loop

        # assign flux values
        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num
            I_plus = I_plus_array[x,y,z,px,py,pz]
            I_minus = I_minus_array[x,y,z,px,py,pz]
            J_plus = J_plus_array[x,y,z,px,py,pz]
            J_minus = J_minus_array[x,y,z,px,py,pz]
            K_plus = K_plus_array[x,y,z,px,py,pz]
            K_minus = K_minus_array[x,y,z,px,py,pz]
            AssignFlux!(PhaseSpace,"I",I_Flux_I,I_Flux_J,I_Flux_V,I_plus,I_minus,x,y,z,px,py,pz,name,scheme,BC_pxp,BC_pxm)
            AssignFlux!(PhaseSpace,"J",J_Flux_I,J_Flux_J,J_Flux_V,J_plus,J_minus,x,y,z,px,py,pz,name,scheme,BC_pyp,BC_pym)
            AssignFlux!(PhaseSpace,"K",K_Flux_I,K_Flux_J,K_Flux_V,K_plus,K_minus,x,y,z,px,py,pz,name,scheme,BC_pzp,BC_pzm)
        end

    end # name loop

    

end
FillMomentumFlux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V) = FillMomentumFlux!(PhaseSpace,Forces,P_Flux_I,P_Flux_J,P_Flux_V,P_Flux_I,P_Flux_J,P_Flux_V,P_Flux_I,P_Flux_J,P_Flux_V)


"""
    FillSpaceFlux!(PhaseSpace,B_Flux_I,B_Flux_J,B_Flux_V,C_Flux_I,C_Flux_J,C_Flux_V,D_Flux_I,D_Flux_J,D_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `B_Flux`, `C_Flux`, and `D_Flux` matrices.
"""
function FillSpaceFlux!(PhaseSpace::PhaseSpaceStruct,B_Flux_I::Vector{Int32},B_Flux_J::Vector{Int32},B_Flux_V::Vector{T},C_Flux_I::Vector{Int32},C_Flux_J::Vector{Int32},C_Flux_V::Vector{T},D_Flux_I::Vector{Int32},D_Flux_J::Vector{Int32},D_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    metric = Spacetime.metric
    coordinates = Spacetime.coordinates
    tetrad = Spacetime.tetrad
    momentum_coords = Momentum.coordinates
    scheme = Spacetime.scheme

    BC_xp = coordinates.xp_BC
    BC_xm = coordinates.xm_BC 
    BC_yp = coordinates.yp_BC
    BC_ym = coordinates.ym_BC
    BC_zp = coordinates.zp_BC
    BC_zm = coordinates.zm_BC

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    CFSpaceBPlus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceBMinus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceCPlus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceCMinus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceDPlus = MVector{4,Float64}(zeros(Float64,4))
    CFSpaceDMinus = MVector{4,Float64}(zeros(Float64,4))
    CFMomentumVector = MVector{4,Float64}(zeros(Float64,4))

    @inline CoordinateFluxBIntegrand!(yztx,B) = CoordinateFluxSpaceBIntegrand!(yztx,B,metric,coordinates,tetrad)
    @inline CoordinateFluxCIntegrand!(ztxy,C) = CoordinateFluxSpaceCIntegrand!(ztxy,C,metric,coordinates,tetrad)
    @inline CoordinateFluxDIntegrand!(txyz,D) = CoordinateFluxSpaceDIntegrand!(txyz,D,metric,coordinates,tetrad)

    non_dimensional_factor::Float64 = CONST_c * Characteristic.CHAR_time / Characteristic.CHAR_length

    for x in 1:x_num, y in 1:y_num, z in 1:z_num

        t0 = Grids.tr[1]
        t1 = Grids.tr[2]
        x0 = Grids.xr[x]
        x1 = Grids.xr[x+1]
        y0 = Grids.yr[y]
        y1 = Grids.yr[y+1]
        z0 = Grids.zr[z]
        z1 = Grids.zr[z+1]

        aB::SVector{4,Float64} = [y0,z0,t0,x0]
        bB::SVector{4,Float64} = [y1,z1,t1,x1]
        nB::SVector{3,Int64} = [16,16,2] 

        aC::SVector{4,Float64} = [z0,t0,x0,y0]
        bC::SVector{4,Float64} = [z1,t1,x1,y1]
        nC::SVector{3,Int64} = [16,2,16] 

        aD::SVector{4,Float64} = [t0,x0,y0,z0]
        bD::SVector{4,Float64} = [t1,x1,y1,z1]
        nD::SVector{3,Int64} = [2,16,16] 

        # Integrate space part of flux
        FluxSimpson3D!(CoordinateFluxBIntegrand!,CFSpaceBPlus,CFSpaceBMinus,aB,bB,nB)
        FluxSimpson3D!(CoordinateFluxCIntegrand!,CFSpaceCPlus,CFSpaceCMinus,aC,bC,nC)
        FluxSimpson3D!(CoordinateFluxDIntegrand!,CFSpaceDPlus,CFSpaceDMinus,aD,bD,nD)

        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                B_plus = 0.0
                B_minus = 0.0
                C_plus = 0.0
                C_minus = 0.0
                D_plus = 0.0
                D_minus = 0.0

                CoordinateFluxMomentumIntegral!(CFMomentumVector,PhaseSpace,name,px,py,pz)

                for a in 1:4

                    # integration sign introduced here
                    B_plus += CFSpaceBPlus[a] * CFMomentumVector[a]
                    B_minus -= CFSpaceBMinus[a] * CFMomentumVector[a]
                    C_plus += CFSpaceCPlus[a] * CFMomentumVector[a]
                    C_minus -= CFSpaceCMinus[a] * CFMomentumVector[a]
                    D_plus += CFSpaceDPlus[a] * CFMomentumVector[a]
                    D_minus -= CFSpaceDMinus[a] * CFMomentumVector[a]

                end

                B_plus *= non_dimensional_factor
                B_minus *= non_dimensional_factor
                C_plus *= non_dimensional_factor
                C_minus *= non_dimensional_factor
                D_plus *= non_dimensional_factor
                D_minus *= non_dimensional_factor

                # fill
                AssignFlux!(PhaseSpace,"B",B_Flux_I,B_Flux_J,B_Flux_V,B_plus,B_minus,x,y,z,px,py,pz,name,scheme,BC_xp,BC_xm)
                AssignFlux!(PhaseSpace,"C",C_Flux_I,C_Flux_J,C_Flux_V,C_plus,C_minus,x,y,z,px,py,pz,name,scheme,BC_yp,BC_ym)
                AssignFlux!(PhaseSpace,"D",D_Flux_I,D_Flux_J,D_Flux_V,D_plus,D_minus,x,y,z,px,py,pz,name,scheme,BC_zp,BC_zm)

            end
        end
    end
end
FillSpaceFlux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V) = FillSpaceFlux!(PhaseSpace,X_Flux_I,X_Flux_J,X_Flux_V,X_Flux_I,X_Flux_J,X_Flux_V,X_Flux_I,X_Flux_J,X_Flux_V)


"""
    Fill_I_Flux!(PhaseSpace,Forces,I_Flux_I,I_Flux_J,I_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `I_Flux` matrix.
"""
function Fill_I_Flux!(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce},I_Flux_I::Vector{Int32},I_Flux_J::Vector{Int32},I_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme

    BCp = momentum_coords.xp_BC
    BCm = momentum_coords.xm_BC 

    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert BCp isa Closed || BCm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            pxp = px+1
            pxm = px-1
            if pxp > px_num # right boundary
                if BCp isa Closed || BCp isa Open
                    pxp = px_num
                elseif BCp isa Periodic
                    pxp = 1
                end
            end
            if pxm < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    pxm = 1
                elseif BCm isa Periodic
                    pxm = px_num
                end
            end

            bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,pxp,py,pz,name)
            bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,pxm,py,pz,name)

            for f in 1:length(Forces)
                # integration sign introduced here
                I_plus = IFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                I_minus = -IFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
        
                # scheme
                if scheme == "upwind"
                    if sign(I_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(I_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 0.5
                    h_plus_left = 0.5
                    h_minus_right = 0.5
                    h_minus_left = 0.5
                else
                    error("Unknown scheme")
                end

                
                #=
                ________________________________
                a |  I_m   | I_m+I_p | I_p    |
                __|________|_________|________|_
                      b-1       b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                Mom_Normp = MomentumSpaceNorm(Grids,name,pxp,py,pz)
                Mom_Normm = MomentumSpaceNorm(Grids,name,pxm,py,pz)

                # normalised fluxes
                if b != bp
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(bp))
                    push!(I_Flux_V,convert(T,(I_plus * h_plus_right) / Mom_Normp))
                    #I_Flux[a,bp] += convert(T,(I_plus * h_plus_right) / Mom_Normp)
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(b))
                    push!(I_Flux_V,convert(T,(I_plus * h_plus_left) / Mom_Norm))
                    #I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
                elseif BCp isa Open # b=bp
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(b))
                    push!(I_Flux_V,convert(T,(I_plus * h_plus_left) / Mom_Norm))
                   # I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
                end
                if b != bm
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(b))
                    push!(I_Flux_V,convert(T,(I_minus * h_minus_right) / Mom_Norm))
                    #I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm)
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(bm))
                    push!(I_Flux_V,convert(T,(I_minus * h_minus_left) / Mom_Normm)) 
                    #I_Flux[a,bm] += convert(T,(I_minus * h_minus_left) / Mom_Normm) 
                elseif BCm isa Open # b=bm
                    push!(I_Flux_I,Int32(a))
                    push!(I_Flux_J,Int32(b))
                    push!(I_Flux_V,convert(T,(I_minus * h_minus_right) / Mom_Norm))
                    #I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
                end

            end # Forces loop
        end # end coordinates loop
    end
end

"""
    Fill_J_Flux!(PhaseSpace,Forces,J_Flux_I,J_Flux_J,J_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `J_Flux` matrix.
"""
function Fill_J_Flux!(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce},J_Flux_I::Vector{Int32},J_Flux_J::Vector{Int32},J_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.yp_BC
    BCm = momentum_coords.ym_BC 
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert (BCp isa Closed) && (BCm isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            pyp = py+1
            pym = py-1
            if pyp > py_num # right boundary
                if BCp isa Closed || BCp isa Open
                    pyp = py_num
                elseif BCp isa Periodic
                    pyp = 1
                end
            end
            if pym < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    pym = 1
                elseif BCm isa Periodic
                    pym = py_num
                end
            end

            bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,pyp,pz,name)
            bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,pym,pz,name)

            for f in 1:length(Forces)
                # integration sign introduced here
                J_plus = JFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                J_minus = -JFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                # scheme
                if scheme == "upwind"
                    if sign(J_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(J_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 0.5
                    h_plus_left = 0.5
                    h_minus_right = 0.5
                    h_minus_left = 0.5
                else
                    error("Unknown scheme")
                end
    
                #=
                ________________________________
                a |  J_m   | J_m+J_p | J_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                Mom_Normp = MomentumSpaceNorm(Grids,name,px,pyp,pz)
                Mom_Normm = MomentumSpaceNorm(Grids,name,px,pym,pz)

                # normalised fluxes
                if b != bp
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(bp))
                    push!(J_Flux_V,convert(T,(J_plus * h_plus_right) / Mom_Normp))
                    #J_Flux[a,bp] += convert(T,(J_plus * h_plus_right) / Mom_Normp)
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(b))
                    push!(J_Flux_V,convert(T,(J_plus * h_plus_left) / Mom_Norm))
                    #J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
                elseif BCp isa Open # b=bp
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(b))
                    push!(J_Flux_V,convert(T,(J_plus * h_plus_left) / Mom_Norm))
                    #J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
                end
                if b != bm
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(b))
                    push!(J_Flux_V,convert(T,(J_minus * h_minus_right) / Mom_Norm))
                    #J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(bm))
                    push!(J_Flux_V,convert(T,(J_minus * h_minus_left) / Mom_Normm))
                    #J_Flux[a,bm] += convert(T,(J_minus * h_minus_left) / Mom_Normm)
                elseif BCm isa Open # b=bm
                    push!(J_Flux_I,Int32(a))
                    push!(J_Flux_J,Int32(b))
                    push!(J_Flux_V,convert(T,(J_minus * h_minus_right) / Mom_Norm))
                    #J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
                end

            end # Force loop
        end
    end
end

"""
    Fill_K_Flux!(PhaseSpace,Forces,K_Flux_I,K_Flux_J,K_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `K_Flux` matrix.
"""
function Fill_K_Flux!(PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce},K_Flux_I::Vector{Int32},K_Flux_J::Vector{Int32},K_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.zp_BC
    BCm = momentum_coords.zm_BC 
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    @assert (BCp isa Periodic) && (BCm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            pzp = pz+1
            pzm = pz-1
            if pzp > pz_num # right boundary
                if BCp isa Closed || BCp isa Open
                    pzp = pz_num
                elseif BCp isa Periodic
                    pzp = 1
                end
            end
            if pzm < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    pzm = 1
                elseif BCm isa Periodic
                    pzm = pz_num
                end
            end

            bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pzp,name)
            bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pzm,name)

            for f in 1:length(Forces)
                # integration sign introduced here
                K_plus = KFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                K_minus = -KFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                # scheme
                if scheme == "upwind"
                    if sign(K_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(K_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 1/2
                    h_plus_left = 1/2
                    h_minus_right = 1/2
                    h_minus_left = 1/2
                else
                    error("Unknown scheme")
                end

                #=
                ________________________________
                a |  K_m   | K_m+K_p | K_p    |
                __|________|_________|________|_
                    b-1        b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                Mom_Normp = MomentumSpaceNorm(Grids,name,px,py,pzp)
                Mom_Normm = MomentumSpaceNorm(Grids,name,px,py,pzm)

                if b != bp
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(bp))
                    push!(K_Flux_V,convert(T,(K_plus * h_plus_right) / Mom_Normp))
                    #K_Flux[a,bp] += convert(T,(K_plus * h_plus_right) / Mom_Normp)
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(b))
                    push!(K_Flux_V,convert(T,(K_plus * h_plus_left) / Mom_Norm))
                    #K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
                elseif BCp isa Open # b=bp
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(b))
                    push!(K_Flux_V,convert(T,(K_plus * h_plus_left) / Mom_Norm))
                    #K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
                end
                if b != bm
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(b))
                    push!(K_Flux_V,convert(T,(K_minus * h_minus_right) / Mom_Norm))
                    #K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(bm))
                    push!(K_Flux_V,convert(T,(K_minus * h_minus_left) / Mom_Normm))
                    #K_Flux[a,bm] += convert(T,(K_minus * h_minus_left) / Mom_Normm)
                elseif BCm isa Open # b=bm
                    push!(K_Flux_I,Int32(a))
                    push!(K_Flux_J,Int32(b))
                    push!(K_Flux_V,convert(T,(K_minus * h_minus_right) / Mom_Norm))
                    #K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
                end
                
            end # force loop
        end
    end
end


"""
    Fill_B_Flux!(PhaseSpace,B_Flux_I,B_Flux_J,B_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `B_Flux` matrix.
"""
function Fill_B_Flux!(PhaseSpace::PhaseSpaceStruct,B_Flux_I::Vector{Int32},B_Flux_J::Vector{Int32},B_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme

    BCp = space_coords.xp_BC
    BCm = space_coords.xm_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            xp = x+1
            xm = x-1
            if xp > x_num # right boundary
                if BCp isa Closed || BCp isa Open || BCp isa Escape
                    xp = x_num
                elseif BCp isa Periodic
                    xp = 1
                end
            end
            if xm < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    xm = 1
                elseif BCm isa Periodic
                    xm = x_num
                end
            end

            bp = GlobalIndicesToStateIndex(PhaseSpace,xp,y,z,px,py,pz,name)
            bm = GlobalIndicesToStateIndex(PhaseSpace,xm,y,z,px,py,pz,name)

            # integration sign introduced here
            B_plus = BFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            B_minus = -BFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(B_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(B_minus) == 1
                    h_minus_right = 1
                    h_minus_left = 0
                else
                    h_minus_right = 0
                    h_minus_left = 1
                end
            elseif scheme == "central"
                h_plus_right = 0.5
                h_plus_left = 0.5
                h_minus_right = 0.5
                h_minus_left = 0.5
            else
                error("Unknown scheme")
            end

            
            #=
            ________________________________
            a |  B_m   | B_m+B_p | B_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(bp))
                push!(B_Flux_V,convert(T,(B_plus * h_plus_right) / Mon_Norm))
                #B_Flux[a,bp] += convert(T,(B_plus * h_plus_right) / Mon_Norm) 
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(b))
                push!(B_Flux_V,convert(T,(B_plus * h_plus_left) / Mon_Norm))
                #B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
            elseif BCp isa Open || BCp isa Escape # b=bp
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(b))
                push!(B_Flux_V,convert(T,(B_plus * h_plus_left) / Mon_Norm))
                #B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
            end
            if b != bm
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(b))
                push!(B_Flux_V,convert(T,(B_minus * h_minus_right) / Mon_Norm))
                #B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(bm))
                push!(B_Flux_V,convert(T,(B_minus * h_minus_left) / Mon_Norm))
                #B_Flux[a,bm] += convert(T,(B_minus * h_minus_left) / Mon_Norm) 
            elseif BCm isa Open # b=bm
                push!(B_Flux_I,Int32(a))
                push!(B_Flux_J,Int32(b))
                push!(B_Flux_V,convert(T,(B_minus * h_minus_right) / Mon_Norm))
                #B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    Fill_C_Flux!(PhaseSpace,C_Flux_I,C_Flux_J,C_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `C_Flux` matrix.
"""
function Fill_C_Flux!(PhaseSpace::PhaseSpaceStruct,C_Flux_I::Vector{Int32},C_Flux_J::Vector{Int32},C_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = space_coords.yp_BC
    BCm = space_coords.ym_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            yp = y+1
            ym = y-1
            if yp > y_num # right boundary
                if BCp isa Closed || BCp isa Open
                    yp = y_num
                elseif BCp isa Periodic
                    yp = 1
                end
            end
            if ym < 1 # left boundary
                if BCm isa Closed || BCm isa Open
                    ym = 1
                elseif BCm isa Periodic
                    ym = y_num
                end
            end

            bp = GlobalIndicesToStateIndex(PhaseSpace,x,yp,z,px,py,pz,name)
            bm = GlobalIndicesToStateIndex(PhaseSpace,x,ym,z,px,py,pz,name)

            # integration sign introduced here
            C_plus = CFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            C_minus = -CFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(C_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(C_minus) == 1
                    h_minus_right = 1
                    h_minus_left = 0
                else
                    h_minus_right = 0
                    h_minus_left = 1
                end
            elseif scheme == "central"
                h_plus_right = 0.5
                h_plus_left = 0.5
                h_minus_right = 0.5
                h_minus_left = 0.5
            else
                error("Unknown scheme")
            end

            
            #=
            ________________________________
            a |  C_m   | C_m+C_p | C_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(bp))
                push!(C_Flux_V,convert(T,(C_plus * h_plus_right) / Mom_Norm))
                #C_Flux[a,bp] += convert(T,(C_plus * h_plus_right) / Mom_Norm) 
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(b))
                push!(C_Flux_V,convert(T,(C_plus * h_plus_left) / Mom_Norm))
                #C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
            elseif BCp isa Open # b=bp
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(b))
                push!(C_Flux_V,convert(T,(C_plus * h_plus_left) / Mom_Norm))
                #C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
            end
            if b != bm
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(b))
                push!(C_Flux_V,convert(T,(C_minus * h_minus_right) / Mom_Norm))
                #C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(bm))
                push!(C_Flux_V,convert(T,(C_minus * h_minus_left) / Mom_Norm))
                #C_Flux[a,bm] += convert(T,(C_minus * h_minus_left) / Mom_Norm) 
            elseif BCm isa Open # b=bm
                push!(C_Flux_I,Int32(a))
                push!(C_Flux_J,Int32(b))
                push!(C_Flux_V,convert(T,(C_minus * h_minus_right) / Mom_Norm))
                #C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    Fill_D_Flux!(PhaseSpace,D_Flux_I,D_Flux_J,D_Flux_V)

Generates the vectors `...I` of row indices, `...J` of column indices, and `...V` of values for the `D_Flux` matrix.
"""
function Fill_D_Flux!(PhaseSpace::PhaseSpaceStruct,D_Flux_I::Vector{Int32},D_Flux_J::Vector{Int32},D_Flux_V::Vector{T}) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids
    Characteristic = PhaseSpace.Characteristic

    space_coords = Space.space_coordinates
    momentum_coords = Momentum.momentum_coordinates
    scheme = Space.scheme

    BCp = space_coords.zp_BC
    BCm = space_coords.zm_BC 

    name_list = PhaseSpace.name_list
    x_num = Space.x_num
    y_num = Space.y_num
    z_num = Space.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list
    
    for name in 1:length(name_list)

        px_num = px_num_list[name]
        py_num = py_num_list[name]
        pz_num = pz_num_list[name]

        for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

            a = GlobalIndicesToStateIndex(PhaseSpace,x,y,z,px,py,pz,name)
            b = a 
            zp = z+1
            zm = z-1
            if zp > z_num # right boundary
                if BCp isa Closed || BCp isa Open || BCp isa Reflective
                    zp = z_num
                elseif BCp isa Periodic
                    zp = 1
                end
            end
            if zm < 1 # left boundary
                if BCm isa Closed || BCm isa Open || BCm isa Reflective
                    zm = 1
                elseif BCm isa Periodic
                    zm = z_num
                end
            end

            if BCp isa Reflective && (z+1 > z_num)
                @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                # mirror u momentum at boundary
                py_ref = py_num - py + 1
                bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,zp,px,py_ref,pz,name)
            else
                bp = GlobalIndicesToStateIndex(PhaseSpace,x,y,zp,px,py,pz,name)
            end
            
            if BCm isa Reflective && (z-1 < 1)
                @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                # mirror u momentum at boundary
                py_ref = py_num - py + 1
                bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,zm,px,py_ref,pz,name)
            else
                bm = GlobalIndicesToStateIndex(PhaseSpace,x,y,zm,px,py,pz,name)
            end

            # integration sign introduced here
            D_plus = DFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
            D_minus = -DFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
    
            # scheme
            if scheme == "upwind"
                if sign(D_plus) == 1
                    h_plus_right = 0
                    h_plus_left = 1
                else
                    h_plus_right = 1
                    h_plus_left = 0
                end
                if sign(D_minus) == 1
                    h_minus_right = 1
                    h_minus_left = 0
                else
                    h_minus_right = 0
                    h_minus_left = 1
                end
            elseif scheme == "central"
                h_plus_right = 0.5
                h_plus_left = 0.5
                h_minus_right = 0.5
                h_minus_left = 0.5
            else
                error("Unknown scheme")
            end

            
            #=
            ________________________________
            a |  D_m   | D_m+D_p | D_p    |
            __|________|_________|________|_
                  b-1       b        b+1  
            =#

            ## note on boundaries:
            # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
            # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

            Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

            # normalised fluxes
            if b != bp
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(bp))
                push!(D_Flux_V,convert(T,(D_plus * h_plus_right) / Mon_Norm))
                #D_Flux[a,bp] += convert(T,(D_plus * h_plus_right) / Mon_Norm)
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(b))
                push!(D_Flux_V,convert(T,(D_plus * h_plus_left) / Mon_Norm)) 
                #D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
            elseif BCp isa Open # b=bp
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(b))
                push!(D_Flux_V,convert(T,(D_plus * h_plus_left) / Mon_Norm))
                #D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
            end
            if b != bm
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(b))
                push!(D_Flux_V,convert(T,(D_minus * h_minus_right) / Mon_Norm))
                #D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(bm))
                push!(D_Flux_V,convert(T,(D_minus * h_minus_left) / Mon_Norm))
                #D_Flux[a,bm] += convert(T,(D_minus * h_minus_left) / Mon_Norm) 
            elseif BCm isa Open # b=bm
                push!(D_Flux_I,Int32(a))
                push!(D_Flux_J,Int32(b))
                push!(D_Flux_V,convert(T,(D_minus * h_minus_right) / Mon_Norm))
                #D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
            end

        end # end coordinates loop
    end
end

"""
    MomentumSpaceNorm(Grids,species_index,px_idx,py_idx,pz_idx)

Returns the normalization factor for a species distribution function in given momentum space sub-domain.
"""
function MomentumSpaceNorm(Grids::GridsStruct,species_index,px_idx,py_idx,pz_idx)

    dpx = Grids.dpx_list[species_index][px_idx]
    dpy = Grids.dpy_list[species_index][py_idx]
    dpz = Grids.dpz_list[species_index][pz_idx]

    return dpx*dpy*dpz

end



#===== OLD CODE =====#
#=  

    """
        Allocate_Flux(PhaseSpace,Precision,debug_mode)

    Allocates arrays for fluxes and volume elements.
    """
    function Allocate_Flux(PhaseSpace::PhaseSpaceStruct,Precision::DataType,debug_mode::Bool)

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum

        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        n_space = x_num*y_num*z_num

        n = n_momentum*n_space

        size = n*n*sizeof(Precision)

        if size > 1e9
            println("Flux will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("Flux will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("Flux will be approx. $(size/1e3) KB in memory if dense")
        else
            println("Flux will be approx. $size bytes in memory if dense")
        end


        # boundary terms included in arrays
        # fluxes allocated with all zeros
        # time fluxes
        Ap_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n) 
        Am_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
        # volume element 
        Vol::Vector{Precision} = zeros(Precision,n_space)
        # sum of space and momentum fluxes 
        X_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
        P_Flux::SparseMatrixCSC{Precision,Int64} = spzeros(Precision,n,n)
        if debug_mode
            # space fluxes
            B_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
            C_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
            D_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
            # momentum fluxes
            I_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
            J_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
            K_Flux = spzeros(Precision,n,n)::SparseMatrixCSC{Precision,Int64}
        else
            # space fluxes
            B_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
            C_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
            D_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
            # momentum fluxes
            I_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
            J_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64}
            K_Flux = spzeros(Precision,0,0)::SparseMatrixCSC{Precision,Int64} 
        end

        return (Ap_Flux,Am_Flux,X_Flux,P_Flux,Vol,B_Flux,C_Flux,D_Flux,I_Flux,J_Flux,K_Flux)

    end

    """
        Fill_A_Flux!(Ap_Flux,Am_Flux,PhaseSpace)

    Generates `Ap_Flux` and `Am_Flux` terms.
    """
    function Fill_A_Flux!(Ap_Flux::SparseMatrixCSC{T,Int64},Am_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)


            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num

                for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                    a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                    b = a

                    # integration sign introduced here
                    A_plus = AFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    A_minus = -AFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                    # normalisation
                    norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                    # fill
                    Ap_Flux[a,b] += convert(T,A_plus / norm)
                    Am_Flux[a,b] += convert(T,A_minus / norm)

                end
            end
        end
    end

    """
    Fill_I_Flux!(I_Flux,PhaseSpace,Forces)

    Generates `I_Flux` terms.
    """
    function Fill_I_Flux!(I_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce}) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Momentum.scheme

        BCp = momentum_coords.xp_BC
        BCm = momentum_coords.xm_BC 

        #= Boundary Conditions:
            Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
        =#
        @assert BCp isa Closed || BCm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                pxp = px+1
                pxm = px-1
                if pxp > px_num # right boundary
                    if BCp isa Closed || BCp isa Open
                        pxp = px_num
                    elseif BCp isa Periodic
                        pxp = 1
                    end
                end
                if pxm < 1 # left boundary
                    if BCm isa Closed || BCm isa Open
                        pxm = 1
                    elseif BCm isa Periodic
                        pxm = px_num
                    end
                end

                bp = GlobalIndices_To_StateIndex(x,y,z,pxp,py,pz,name,PhaseSpace)
                bm = GlobalIndices_To_StateIndex(x,y,z,pxm,py,pz,name,PhaseSpace)

                for f in 1:length(Forces)
                    # integration sign introduced here
                    I_plus = IFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    I_minus = -IFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
            
                    # scheme
                    if scheme == "upwind"
                        if sign(I_plus) == 1
                            h_plus_right = 0
                            h_plus_left = 1
                        else
                            h_plus_right = 1
                            h_plus_left = 0
                        end
                        if sign(I_minus) == 1
                            h_minus_right = 1
                            h_minus_left = 0
                        else
                            h_minus_right = 0
                            h_minus_left = 1
                        end
                    elseif scheme == "central"
                        h_plus_right = 0.5
                        h_plus_left = 0.5
                        h_minus_right = 0.5
                        h_minus_left = 0.5
                    else
                        error("Unknown scheme")
                    end

                    
                    #=
                    ________________________________
                    a |  I_m   | I_m+I_p | I_p    |
                    __|________|_________|________|_
                        b-1       b        b+1  
                    =#

                    ## note on boundaries:
                    # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                    # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                    Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                    Mom_Normp = MomentumSpaceNorm(Grids,name,pxp,py,pz)
                    Mom_Normm = MomentumSpaceNorm(Grids,name,pxm,py,pz)

                    # normalised fluxes
                    if b != bp
                        I_Flux[a,bp] += convert(T,(I_plus * h_plus_right) / Mom_Normp)
                        I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
                    elseif BCp isa Open # b=bp
                        I_Flux[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
                    end
                    if b != bm
                        I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
                        I_Flux[a,bm] += convert(T,(I_minus * h_minus_left) / Mom_Normm) 
                    elseif BCm isa Open # b=bm
                        I_Flux[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
                    end

                end # Forces loop
            end # end coordinates loop
        end
    end

    """
        Fill_J_Flux!(J_Flux,PhaseSpace,Forces)

    Generates `J_Flux` terms.
    """
    function Fill_J_Flux!(J_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce}) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Momentum.scheme
        offset = PhaseSpace.Grids.momentum_species_offset

        BCp = momentum_coords.yp_BC
        BCm = momentum_coords.ym_BC 
        #= Boundary Conditions:
            Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
        =#
        @assert (BCp isa Closed) && (BCm isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                pyp = py+1
                pym = py-1
                if pyp > py_num # right boundary
                    if BCp isa Closed || BCp isa Open
                        pyp = py_num
                    elseif BCp isa Periodic
                        pyp = 1
                    end
                end
                if pym < 1 # left boundary
                    if BCm isa Closed || BCm isa Open
                        pym = 1
                    elseif BCm isa Periodic
                        pym = py_num
                    end
                end

                bp = GlobalIndices_To_StateIndex(x,y,z,px,pyp,pz,name,PhaseSpace)
                bm = GlobalIndices_To_StateIndex(x,y,z,px,pym,pz,name,PhaseSpace)

                for f in 1:length(Forces)
                    # integration sign introduced here
                    J_plus = JFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    J_minus = -JFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                    # scheme
                    if scheme == "upwind"
                        if sign(J_plus) == 1
                            h_plus_right = 0
                            h_plus_left = 1
                        else
                            h_plus_right = 1
                            h_plus_left = 0
                        end
                        if sign(J_minus) == 1
                            h_minus_right = 1
                            h_minus_left = 0
                        else
                            h_minus_right = 0
                            h_minus_left = 1
                        end
                    elseif scheme == "central"
                        h_plus_right = 0.5
                        h_plus_left = 0.5
                        h_minus_right = 0.5
                        h_minus_left = 0.5
                    else
                        error("Unknown scheme")
                    end
        
                    #=
                    ________________________________
                    a |  J_m   | J_m+J_p | J_p    |
                    __|________|_________|________|_
                        b-1        b        b+1  
                    =#

                    ## note on boundaries:
                    # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                    # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                    Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                    Mom_Normp = MomentumSpaceNorm(Grids,name,px,pyp,pz)
                    Mom_Normm = MomentumSpaceNorm(Grids,name,px,pym,pz)

                    # normalised fluxes
                    if b != bp
                        J_Flux[a,bp] += convert(T,(J_plus * h_plus_right) / Mom_Normp)
                        J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
                    elseif BCp isa Open # b=bp
                        J_Flux[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm)
                    end
                    if b != bm
                        J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
                        J_Flux[a,bm] += convert(T,(J_minus * h_minus_left) / Mom_Normm)
                    elseif BCm isa Open # b=bm
                        J_Flux[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm)
                    end

                end # Force loop
            end
        end
    end

    """
        Fill_K_Flux!(K_Flux,PhaseSpace)

    Generates `K_Flux` terms.
    """
    function Fill_K_Flux!(K_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct,Forces::Vector{AbstractForce}) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Momentum.scheme
        offset = PhaseSpace.Grids.momentum_species_offset

        BCp = momentum_coords.zp_BC
        BCm = momentum_coords.zm_BC 
        #= Boundary Conditions:
            Flux on K boundaries should always be periodic
        =#
        @assert (BCp isa Periodic) && (BCm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                pzp = pz+1
                pzm = pz-1
                if pzp > pz_num # right boundary
                    if BCp isa Closed || BCp isa Open
                        pzp = pz_num
                    elseif BCp isa Periodic
                        pzp = 1
                    end
                end
                if pzm < 1 # left boundary
                    if BCm isa Closed || BCm isa Open
                        pzm = 1
                    elseif BCm isa Periodic
                        pzm = pz_num
                    end
                end

                bp = GlobalIndices_To_StateIndex(x,y,z,px,py,pzp,name,PhaseSpace)
                bm = GlobalIndices_To_StateIndex(x,y,z,px,py,pzm,name,PhaseSpace)

                for f in 1:length(Forces)
                    # integration sign introduced here
                    K_plus = KFluxFunction(Forces[f],PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                    K_minus = -KFluxFunction(Forces[f],PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

                    # scheme
                    if scheme == "upwind"
                        if sign(K_plus) == 1
                            h_plus_right = 0
                            h_plus_left = 1
                        else
                            h_plus_right = 1
                            h_plus_left = 0
                        end
                        if sign(K_minus) == 1
                            h_minus_right = 1
                            h_minus_left = 0
                        else
                            h_minus_right = 0
                            h_minus_left = 1
                        end
                    elseif scheme == "central"
                        h_plus_right = 1/2
                        h_plus_left = 1/2
                        h_minus_right = 1/2
                        h_minus_left = 1/2
                    else
                        error("Unknown scheme")
                    end

                    #=
                    ________________________________
                    a |  K_m   | K_m+K_p | K_p    |
                    __|________|_________|________|_
                        b-1        b        b+1  
                    =#

                    ## note on boundaries:
                    # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                    # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                    Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)
                    Mom_Normp = MomentumSpaceNorm(Grids,name,px,py,pzp)
                    Mom_Normm = MomentumSpaceNorm(Grids,name,px,py,pzm)

                    if b != bp
                        K_Flux[a,bp] += convert(T,(K_plus * h_plus_right) / Mom_Normp)
                        K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
                    elseif BCp isa Open # b=bp
                        K_Flux[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm)
                    end
                    if b != bm
                        K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
                        K_Flux[a,bm] += convert(T,(K_minus * h_minus_left) / Mom_Normm)
                    elseif BCm isa Open # b=bm
                        K_Flux[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm)
                    end
                    
                end # force loop
            end
        end
    end


    """
        Fill_B_Flux!(B_Flux,PhaseSpace)

    Generates `B_Flux` terms.
    """
    function Fill_B_Flux!(B_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Space.scheme

        BCp = space_coords.xp_BC
        BCm = space_coords.xm_BC 

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                xp = x+1
                xm = x-1
                if xp > x_num # right boundary
                    if BCp isa Closed || BCp isa Open
                        xp = x_num
                    elseif BCp isa Periodic
                        xp = 1
                    end
                end
                if xm < 1 # left boundary
                    if BCm isa Closed || BCm isa Open
                        xm = 1
                    elseif BCm isa Periodic
                        xm = x_num
                    end
                end

                bp = GlobalIndices_To_StateIndex(xp,y,z,px,py,pz,name,PhaseSpace)
                bm = GlobalIndices_To_StateIndex(xm,y,z,px,py,pz,name,PhaseSpace)

                # integration sign introduced here
                B_plus = BFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                B_minus = -BFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
        
                # scheme
                if scheme == "upwind"
                    if sign(B_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(B_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 0.5
                    h_plus_left = 0.5
                    h_minus_right = 0.5
                    h_minus_left = 0.5
                else
                    error("Unknown scheme")
                end

                
                #=
                ________________________________
                a |  B_m   | B_m+B_p | B_p    |
                __|________|_________|________|_
                    b-1       b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # normalised fluxes
                if b != bp
                    B_Flux[a,bp] += convert(T,(B_plus * h_plus_right) / Mon_Norm) 
                    B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
                elseif BCp isa Open # b=bp
                    B_Flux[a,b] += convert(T,(B_plus * h_plus_left) / Mon_Norm) 
                end
                if b != bm
                    B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
                    B_Flux[a,bm] += convert(T,(B_minus * h_minus_left) / Mon_Norm) 
                elseif BCm isa Open # b=bm
                    B_Flux[a,b] += convert(T,(B_minus * h_minus_right) / Mon_Norm) 
                end

            end # end coordinates loop
        end
    end

    """
        Fill_C_Flux!(C_Flux,PhaseSpace)

    Generates `C_Flux` terms.
    """
    function Fill_C_Flux!(C_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Space.scheme
        offset = PhaseSpace.Grids.momentum_species_offset

        BCp = space_coords.yp_BC
        BCm = space_coords.ym_BC 

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                yp = y+1
                ym = y-1
                if yp > y_num # right boundary
                    if BCp isa Closed || BCp isa Open
                        yp = y_num
                    elseif BCp isa Periodic
                        yp = 1
                    end
                end
                if ym < 1 # left boundary
                    if BCm isa Closed || BCm isa Open
                        ym = 1
                    elseif BCm isa Periodic
                        ym = y_num
                    end
                end

                bp = GlobalIndices_To_StateIndex(x,yp,z,px,py,pz,name,PhaseSpace)
                bm = GlobalIndices_To_StateIndex(x,ym,z,px,py,pz,name,PhaseSpace)

                # integration sign introduced here
                C_plus = CFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                C_minus = -CFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
        
                # scheme
                if scheme == "upwind"
                    if sign(C_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(C_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 0.5
                    h_plus_left = 0.5
                    h_minus_right = 0.5
                    h_minus_left = 0.5
                else
                    error("Unknown scheme")
                end

                
                #=
                ________________________________
                a |  C_m   | C_m+C_p | C_p    |
                __|________|_________|________|_
                    b-1       b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mom_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # normalised fluxes
                if b != bp
                    C_Flux[a,bp] += convert(T,(C_plus * h_plus_right) / Mom_Norm) 
                    C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
                elseif BCp isa Open # b=bp
                    C_Flux[a,b] += convert(T,(C_plus * h_plus_left) / Mom_Norm) 
                end
                if b != bm
                    C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
                    C_Flux[a,bm] += convert(T,(C_minus * h_minus_left) / Mom_Norm) 
                elseif BCm isa Open # b=bm
                    C_Flux[a,b] += convert(T,(C_minus * h_minus_right) / Mom_Norm) 
                end

            end # end coordinates loop
        end
    end

    """
        Fill_D_Flux!(D_Flux,PhaseSpace)

    Generates `D_Flux` terms.
    """
    function Fill_D_Flux!(D_Flux::SparseMatrixCSC{T,Int64},PhaseSpace::PhaseSpaceStruct) where T<:Union{Float32,Float64}

        Spacetime = PhaseSpace.Spacetime
        Momentum = PhaseSpace.Momentum
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic

        space_coords = Space.space_coordinates
        momentum_coords = Momentum.momentum_coordinates
        scheme = Space.scheme

        BCp = space_coords.zp_BC
        BCm = space_coords.zm_BC 

        name_list = PhaseSpace.name_list
        x_num = Space.x_num
        y_num = Space.y_num
        z_num = Space.z_num
        px_num_list = Momentum.px_num_list
        py_num_list = Momentum.py_num_list
        pz_num_list = Momentum.pz_num_list
        
        for name in 1:length(name_list)

            px_num = px_num_list[name]
            py_num = py_num_list[name]
            pz_num = pz_num_list[name]

            for x in 1:x_num, y in 1:y_num, z in 1:z_num, px in 1:px_num, py in 1:py_num, pz in 1:pz_num

                a = GlobalIndices_To_StateIndex(x,y,z,px,py,pz,name,PhaseSpace)
                b = a 
                zp = z+1
                zm = z-1
                if zp > z_num # right boundary
                    if BCp isa Closed || BCp isa Open || BCp isa Reflective
                        zp = z_num
                    elseif BCp isa Periodic
                        zp = 1
                    end
                end
                if zm < 1 # left boundary
                    if BCm isa Closed || BCm isa Open || BCm isa Reflective
                        zm = 1
                    elseif BCm isa Periodic
                        zm = z_num
                    end
                end

                if BCp isa Reflective && (z+1 > z_num)
                    @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                    # mirror u momentum at boundary
                    py_ref = py_num - py + 1
                    bp = GlobalIndices_To_StateIndex(x,y,zp,px,py_ref,pz,name,PhaseSpace)
                else
                    bp = GlobalIndices_To_StateIndex(x,y,zp,px,py,pz,name,PhaseSpace)
                end
                
                if BCm isa Reflective && (z-1 < 1)
                    @assert space_coords isa Cartesian "Reflective BCs only implemented for Cartesian coordinates"
                    # mirror u momentum at boundary
                    py_ref = py_num - py + 1
                    bm = GlobalIndices_To_StateIndex(x,y,zm,px,py_ref,pz,name,PhaseSpace)
                else
                    bm = GlobalIndices_To_StateIndex(x,y,zm,px,py,pz,name,PhaseSpace)
                end

                # integration sign introduced here
                D_plus = DFluxFunction(PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
                D_minus = -DFluxFunction(PhaseSpace,name,"minus",1,x,y,z,px,py,pz)
        
                # scheme
                if scheme == "upwind"
                    if sign(D_plus) == 1
                        h_plus_right = 0
                        h_plus_left = 1
                    else
                        h_plus_right = 1
                        h_plus_left = 0
                    end
                    if sign(D_minus) == 1
                        h_minus_right = 1
                        h_minus_left = 0
                    else
                        h_minus_right = 0
                        h_minus_left = 1
                    end
                elseif scheme == "central"
                    h_plus_right = 0.5
                    h_plus_left = 0.5
                    h_minus_right = 0.5
                    h_minus_left = 0.5
                else
                    error("Unknown scheme")
                end

                
                #=
                ________________________________
                a |  D_m   | D_m+D_p | D_p    |
                __|________|_________|________|_
                    b-1       b        b+1  
                =#

                ## note on boundaries:
                # b = bp means at right boundary, therefore no plus flux from either direction if closed boundary (no particles leave/enter) or only plus flux from the left if open boundary (particles may only leave)
                # b = bm means at left boundary, therefore no minus flux from either direction if closed boundary (no particles leave/enter) or only minus flux from the right if open boundary (particles may only leave)

                Mon_Norm = MomentumSpaceNorm(Grids,name,px,py,pz)

                # normalised fluxes
                if b != bp
                    D_Flux[a,bp] += convert(T,(D_plus * h_plus_right) / Mon_Norm) 
                    D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
                elseif BCp isa Open # b=bp
                    D_Flux[a,b] += convert(T,(D_plus * h_plus_left) / Mon_Norm) 
                end
                if b != bm
                    D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
                    D_Flux[a,bm] += convert(T,(D_minus * h_minus_left) / Mon_Norm) 
                elseif BCm isa Open # b=bm
                    D_Flux[a,b] += convert(T,(D_minus * h_minus_right) / Mon_Norm) 
                end

            end # end coordinates loop
        end
    end
    =#