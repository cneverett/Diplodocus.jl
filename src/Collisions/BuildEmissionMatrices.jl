"""
    BuildEmissionMatrices(PhaseSpace,DataDirectory;loading_check)

Function that builds the emission matrices associated with emissive interactions and their corresponding reaction forces if applicable. If there are such interactions, first space is allocated for the arrays, then data is loaded into these arrays from the desired `DataDirectory` location and finally the big matrices are returned as an immutable `EmissionMatricesStruct`.

For spatial indices within `Binary_Domain` (if set to `nothing` all cells are in the Binary_Domain) the emission matrices will be defined as dense matrices for compatibility with GPU `copyto!`. Otherwise they will be defined as sparse matrices.
"""
function BuildEmissionMatrices(PhaseSpace::PhaseSpaceStruct,Emission_list::Vector{EmissiveInteraction},DataDirectory::String,BinM::BinaryMatricesStruct;Emi_corrected::Bool=true)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum
    Spacetime = PhaseSpace.Spacetime

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = Spacetime.x_num*Spacetime.y_num*Spacetime.z_num

    n = n_momentum*n_space

    M_Emi = Vector{Union{Matrix{Precision},SparseMatrixCSC{Precision,Int32}}}(undef,n_space)

    size = (n)^2*sizeof(Precision)

    if !isempty(Emission_list)
        LoadMatrices_Emi!(M_Emi,Emission_list,DataDirectory,PhaseSpace,BinM;Emi_corrected=Emi_corrected)
    end

    size = Base.summarysize(M_Emi)

    if size > 1e9
        println("M_Emi is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("M_Emi is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("M_Emi is approx. $(size/1e3) KB in memory")
    else
        println("M_Emi is approx. $size bytes in memory")
    end

    return EmissionMatricesStruct{Precision}(M_Emi,Emission_list)

end


"""
    Allocate_M_Emi(Momentum::MomentumStruct)

Allocates a big matrix `M_Emi` which stores the interaction rates for emission interactions between all particles in the simulation.
"""
function Allocate_M_Emi(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType,Emi_sparse::Bool)

    Momentum = PhaseSpace.Momentum
    Space = PhaseSpace.Space

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)
    n_space = Space.x_num*Space.y_num*Space.z_num

    n = n_momentum*n_space

    size = (n)^2*sizeof(Precision)

    if loading_check

        if size > 1e9
            println("M_Emi will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Emi will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Emi will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Emi will be approx. $size bytes in memory if dense")
        end

        #println("Do you want to proceed? (y/n)")
        #proceed = readline()
        #proceed = "y" ## FIX LATER
        #if proceed != "y"
        #    error("Error: User aborted emission interaction matrix loading")
        #end

    else

        println("Building emission interaction matrix")
        if size > 1e9
            println("M_Emi will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Emi will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Emi will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Emi will be approx. $size bytes in memory if dense")
        end

    end

    if Emi_sparse
        M_Emi = spzeros(Precision,n,n)
    else 
        M_Emi = zeros(Precision,n,n)
    end

    return M_Emi

end

function Fill_M_Emi!(PhaseSpace::PhaseSpaceStruct,names::Tuple{Int64,Int64,Int64},x::Int64,y::Int64,z::Int64;GainMatrix2=nothing,GainMatrix3=nothing,LossMatrix1=nothing,mode::AbstractMode=Ani(),M_Emi::Union{Nothing,Matrix{F}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids

    (name1,name2,name3) = names

    # interaction is name1 -> name2 + name3
    # absorption of name1 and emission of name2 not implemented with name1 == name2

    #GainMatrix_to_M_Emi!(M_Emi,PhaseSpace,GainMatrix2,name2,name1,x,y,z,mode;M_Emi=M_Emi,M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
    GainMatrix_to_M_Emi!(PhaseSpace,GainMatrix3,name3,name1,x,y,z,mode;M_Emi=M_Emi,M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
    #LossMatrix_to_M_Emi!(M_Emi,PhaseSpace,LossMatrix1,name1,x,y,z,mode;M_Emi=M_Emi,M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)

    #GainMatrix2 = nothing
    GainMatrix3 = nothing
    #LossMatrix1 = nothing

    GC.gc()

end


function GainMatrix_to_M_Emi!(PhaseSpace::PhaseSpaceStruct,GainMatrix::AbstractArray{Float64,6},name2::Int64,name1::Int64,x::Int64,y::Int64,z::Int64,mode::AbstractMode;M_Emi::Union{Nothing,Matrix{F}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    vol = VolFunction(PhaseSpace,1,x,y,z)

    dpy1 = PhaseSpace.Grids.dpy_list[name1]
    dpy2 = PhaseSpace.Grids.dpy_list[name2]
    dpz1 = PhaseSpace.Grids.dpz_list[name1]
    dpz2 = PhaseSpace.Grids.dpz_list[name2]

    is_sparse = isnothing(M_Emi)
    
    for px1 in axes(GainMatrix,4), px2 in axes(GainMatrix,1)

        #if px1 == 1 
        #    continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        #end

        if mode isa Iso
            val = 0.0
            w = 1.0/(sum(dpy1)*sum(dpz1)*sum(dpy2)*sum(dpz2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in axes(GainMatrix,5), pz1 in axes(GainMatrix,6), py2 in axes(GainMatrix,2), pz2 in axes(GainMatrix,3)
                val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpy2[py2] * dpz2[pz2] * dpy1[py1] * dpz1[pz1]
            end
        end

        for py1 in axes(GainMatrix,5), py2 in axes(GainMatrix,2)
            
            if mode isa Axi
                val = 0.0
                w = 1.0/(sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3)
                    val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpz2[pz2] * dpz1[pz1]
                end
            end

            for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3)

                if mode isa Ani 
                    val = GainMatrix[px2,py2,pz2,px1,py1,pz1]
                    w = 1.0
                end

                a = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px2,py2,pz2,name2) # M_Emi now stored in vector of M_Emi for each location
                b = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px1,py1,pz1,name1)

                if val == 0.0
                    continue
                end

                if is_sparse
                    push!(M_Emi_I,Int32(a))
                    push!(M_Emi_J,Int32(b))
                    push!(M_Emi_V,convert(F,val*w*vol))
                else
                    M_Emi[a,b] += convert(F,val*w*vol) 
                end

            end

        end

    end

end

function LossMatrix_to_M_Emi!(PhaseSpace::PhaseSpaceStruct,LossMatrix::Array{Float64,2},name1::Int64,x::Int64,y::Int64,z::Int64,mode::AbstractMode;M_Emi::Union{Nothing,Matrix{F}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    vol = VolFunction(PhaseSpace,1,x,y,z)

    dpy1 = PhaseSpace.Grids.dpy_list[name1]
    dpz1 = PhaseSpace.Grids.dpz_list[name1]

    is_sparse = isnothing(M_Emi)

    for px1 in axes(LossMatrix,1)

        if px1 == 1 
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end

        if mode isa Iso

            val = 0.0 
            w = 1.0 

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in axes(LossMatrix,2), pz1 in axes(LossMatrix,3)
                val += LossMatrix[px1,py1,pz1] * dpy1[py1] * dpz1[pz1]
            end

            w /= (sum(dpy1) * sum(dpz1))

        end

        for py1 in axes(LossMatrix,2)

            if mode isa Axi

                val = 0.0 
                w = 1.0 

                # average over incoming phi angles (pz)
                for pz1 in axes(LossMatrix,3)
                    val += LossMatrix[px1,py1,pz1] * dpz1[pz1]
                end

                w /= sum(dpz1)

            end

            for pz1 in axes(LossMatrix,3)

                if mode isa Ani 

                    val = 0.0 
                    w = 1.0 

                    val += LossMatrix[px1,py1,pz1]

                end

                a = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px1,py1,pz1,name1) # M_Emi now stored in vector of M_Emi for each location
                b = a

                if val == 0.0 
                    continue
                end 

                if is_sparse
                    push!(M_Emi_I,Int32(a))
                    push!(M_Emi_J,Int32(b))
                    push!(M_Emi_V,-convert(F,val*w*vol))
                else
                    M_Emi[a,b] -= convert(F,val*w*vol)
                end 

            end

        end
        
    end

end


"""
    Fill_I_Emi!(PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx;M_Emi,M_Emi_I,M_Emi_J,M_Emi_V)

Generates `I_Flux` terms in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_I_Emi!(PhaseSpace::PhaseSpaceStruct,Force::AbstractForce,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64;M_Emi::Union{Nothing,Matrix{T}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{T}}=nothing) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    momentum_coords = Momentum.coordinates
    scheme = Momentum.scheme

    BCp = momentum_coords.xp_BC
    BCm = momentum_coords.xm_BC 

    #= Boundary Conditions:
        Flux on I boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert BCp isa Closed || BCm isa Closed "I flux boundaries incorrectly defined, i.e. not closed"

    is_sparse = isnothing(M_Emi)

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx

    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,name) # M_Emi now stored in vector of M_Emi for each location
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

        bp = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,pxp,py,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,pxm,py,pz,name)

        # integration sign introduced here, opposite sign compared to IFluxFunction as I_Flux is on left hand side of transport equation but M_Emi is on the right hand side.
        I_plus = -IFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        I_minus = IFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(I_plus) == -1 # sign change as above
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(I_minus) == -1 # sign change as above
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
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bp))
                push!(M_Emi_V,convert(T,(I_plus * h_plus_right) / Mom_Normp))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(I_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,bp] += convert(T,(I_plus * h_plus_right) / Mom_Normp)
                M_Emi[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
            end
        elseif BCp isa Open # b=bp
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(I_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(I_plus * h_plus_left) / Mom_Norm) 
            end
        end
        if b != bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bm))
                push!(M_Emi_V,convert(T,(I_minus * h_minus_left) / Mom_Normm))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(I_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
                M_Emi[a,bm] += convert(T,(I_minus * h_minus_left) / Mom_Normm) 
            end
        elseif BCm isa Open # b=bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(I_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(I_minus * h_minus_right) / Mom_Norm) 
            end
        end

    end # end coordinates loop

end

"""
    Fill_J_Emi!(PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx;M_Emi,M_Emi_I,M_Emi_J,M_Emi_V)

Generates `J_Flux` term in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_J_Emi!(PhaseSpace::PhaseSpaceStruct,Force::AbstractForce,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64;M_Emi::Union{Nothing,Matrix{T}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{T}}=nothing) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    momentum_coords = Momentum.coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.yp_BC
    BCm = momentum_coords.ym_BC 
    #= Boundary Conditions:
        Flux on J boundaries should always be closed i.e. no particles leave/enter from the domain bound
    =#
    @assert (BCp isa Closed) && (BCm isa Closed) "J flux boundaries incorrectly defined, i.e. not closed"

    is_sparse = isnothing(M_Emi)

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx
    
    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,name) # M_Emi now stored in vector of M_Emi for each location
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

        bp = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,pyp,pz,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,pym,pz,name)

        # integration sign introduced here, opposite sign compared to JFluxFunction as J_Flux is on left hand side of transport equation but M_Emi is on the right hand side.
        J_plus = -JFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        J_minus = JFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(J_plus) == -1 # sign change as above
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(J_minus) == -1 # sign change as above
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
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bp))
                push!(M_Emi_V,convert(T,(J_plus * h_plus_right) / Mom_Normp))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(J_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,bp] += convert(T,(J_plus * h_plus_right) / Mom_Normp)
                M_Emi[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm) 
            end
        elseif BCp isa Open # b=bp
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(J_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(J_plus * h_plus_left) / Mom_Norm) 
            end
        end
        if b != bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bm))
                push!(M_Emi_V,convert(T,(J_minus * h_minus_left) / Mom_Normm))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(J_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm) 
                M_Emi[a,bm] += convert(T,(J_minus * h_minus_left) / Mom_Normm) 
            end
        elseif BCm isa Open # b=bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(J_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(J_minus * h_minus_right) / Mom_Norm) 
            end
        end

    end
end

"""
    Fill_K_Emi!(PhaseSpace,Force,x_idx,y_idx,z_idx,species_idx;M_Emi,M_Emi_I,M_Emi_J,M_Emi_V)

Generates `K_Flux` terms in the Emission matrix `M_Emi` if the emission interaction has an associated `Force`.
"""
function Fill_K_Emi!(PhaseSpace::PhaseSpaceStruct,Force::AbstractForce,x_idx::Int64,y_idx::Int64,z_idx::Int64,species_idx::Int64;M_Emi::Union{Nothing,Matrix{T}}=nothing,M_Emi_I::Union{Nothing,Vector{Int32}}=nothing,M_Emi_J::Union{Nothing,Vector{Int32}}=nothing,M_Emi_V::Union{Nothing,Vector{T}}=nothing) where T<:Union{Float32,Float64}

    Spacetime = PhaseSpace.Spacetime
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    momentum_coords = Momentum.coordinates
    scheme = Momentum.scheme
    offset = PhaseSpace.Grids.momentum_species_offset

    BCp = momentum_coords.zp_BC
    BCm = momentum_coords.zm_BC 
    #= Boundary Conditions:
        Flux on K boundaries should always be periodic
    =#
    @assert (BCp isa Periodic) && (BCm isa Periodic) "K flux boundaries incorrectly defined, i.e. not periodic"

    is_sparse = isnothing(M_Emi)

    name_list = PhaseSpace.name_list
    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    x = x_idx
    y = y_idx
    z = z_idx
    name = species_idx

    px_num = px_num_list[name]
    py_num = py_num_list[name]
    pz_num = pz_num_list[name]

    for px in 1:px_num, py in 1:py_num, pz in 1:pz_num

        a = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pz,name) # M_Emi now stored in vector of M_Emi for each location
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

        bp = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pzp,name)
        bm = GlobalIndicesToStateIndex(PhaseSpace,1,1,1,px,py,pzm,name)

        # integration sign introduced here, opposite sign compared to KFluxFunction as K_Flux is on left hand side of transport equation but M_Emi is on the right hand side.
        K_plus = -KFluxFunction(Force,PhaseSpace,name,"plus",1,x,y,z,px,py,pz)
        K_minus = KFluxFunction(Force,PhaseSpace,name,"minus",1,x,y,z,px,py,pz)

        # scheme
        if scheme == "upwind"
            if sign(K_plus) == -1 # sign change as above
                h_plus_right = 0
                h_plus_left = 1
            else
                h_plus_right = 1
                h_plus_left = 0
            end
            if sign(K_minus) == -1 # sign change as above
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
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bp))
                push!(M_Emi_V,convert(T,(K_plus * h_plus_right) / Mom_Normp))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(K_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,bp] += convert(T,(K_plus * h_plus_right) / Mom_Normp) 
                M_Emi[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm) 
            end
        elseif BCp isa Open # b=bp
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(K_plus * h_plus_left) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(K_plus * h_plus_left) / Mom_Norm) 
            end
        end
        if b != bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(bm))
                push!(M_Emi_V,convert(T,(K_minus * h_minus_left) / Mom_Normm))
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(K_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm) 
                M_Emi[a,bm] += convert(T,(K_minus * h_minus_left) / Mom_Normm) 
            end
        elseif BCm isa Open # b=bm
            if is_sparse
                push!(M_Emi_I,Int32(a))
                push!(M_Emi_J,Int32(b))
                push!(M_Emi_V,convert(T,(K_minus * h_minus_right) / Mom_Norm))
            else
                M_Emi[a,b] += convert(T,(K_minus * h_minus_right) / Mom_Norm) 
            end
        end
                
    end
end
