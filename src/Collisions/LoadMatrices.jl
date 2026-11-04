"""
    LoadMatrices_Binary(...,Binary_list,DataDirectory,PhaseSpace,mode,Ani(),corrected=true)

Loads the binary interaction matrices from the specified `DataDirectory` for each interaction in `Binary_list` and fills the big matrix `M_Bin` either directly if `M_Bin` is dense or the vectors of rows, columns and values `M_Bin_I`, `M_Bin_J``, `M_Bin_V` if sparse, with the interaction rates for binary interactions between all particles in the simulation. The `corrected` argument specifies whether to load the number and energy corrected matrices. The `mode` argument specifies whether to assume these arrays are anisotropic, axisymmetric or isotropic.
"""
function LoadMatrices_Binary(Binary_list::Vector{BinaryInteraction},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,mode::AbstractMode=Ani(),corrected::Bool=true;symmetric::Bool=false,M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int32}}=nothing,M_Bin_J::Union{Nothing,Vector{Int32}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing,gain_loss::Union{Nothing,String}=nothing) where F<:AbstractFloat

    Bin_Norm = PhaseSpace.Characteristic.Bin_Norm
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = PhaseSpace.Grids.n_momentum

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        Indistinguishable_12 = name1 == name2

        println("Filling M_Bin for binary interaction: $(name1) + $(name2) -> $(name3) + $(name4)")

        # Memory optimisation by allowing Ele and Pos populations to be modelled as identical thus only requiring one to be defined
        if !isnothing(findfirst(==("Pos"),[name1,name2,name3,name4])) && isnothing(findfirst(==("Pos"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho") # "Pos" compton scattering
                continue # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            elseif (name1,name2,name3,name4) == ("Ele","Pos","Pho","Pho") # "Ele" "Pos" annihilation
                GainScale = 1.0/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                LossScale = 1.0/4.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                LossScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            LossScale = 1.0
            name1_loc = findfirst(==(name1),name_list)
            name2_loc = findfirst(==(name2),name_list)
            name3_loc = findfirst(==(name3),name_list)
            name4_loc = findfirst(==(name4),name_list)
        end

        # ele pos swap for compton i.e. use ele compton matrices for pos compton
        if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho")
            name1 = "Ele"
            name3 = "Ele"
        end

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        px4_grid::String = px_grid_list[name4_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        py4_grid::String = py_grid_list[name4_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
        pz4_grid::String = pz_grid_list[name4_loc]

        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        px4_num::Int64 = px_num_list[name4_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        py4_num::Int64 = py_num_list[name4_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]
        pz4_num::Int64 = pz_num_list[name4_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]
        px4_low::Float64 = px_low_list[name4_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]
        px4_up::Float64 = px_up_list[name4_loc]

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]
        m4::Float64 = PhaseSpace.Grids.mass_list[name4_loc]

        Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)

        filename::String = BinaryFileName(Parameters)

        println(filename)

        Output = BinaryFileLoad_Matrix(DataDirectory,filename,corrected=corrected)
        
        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainScale *= Bin_Norm 
        LossScale *= Bin_Norm

        Parameters = Output[1]
        GainMatrix3 = Output[2]
        GainMatrix4 = Output[3]
        LossMatrix = Output[4]

        if !isnothing(gain_loss) # used for timescale plots where it is useful to only load the gain or loss matrices to generate loss and gain rates
            if gain_loss == "gain"
                LossScale = 0.0  
            elseif gain_loss == "loss"
                GainScale = 0.0
            end
        end

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        #DoesConserve(Output) # print conversion statistic
        M_Bin_I, M_Bin_J, M_Bin_V = Fill_M_Bin!(name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix,GainScale,LossScale,Indistinguishable_12;mode=mode,symmetric=symmetric,M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

        println("size in memory of M_Bin_I: $(Base.summarysize(M_Bin_I)/1e9) GB")
        println("size in memory of M_Bin_J: $(Base.summarysize(M_Bin_J)/1e9) GB")
        println("size in memory of M_Bin_V: $(Base.summarysize(M_Bin_V)/1e9) GB")

    end # for

end

#=
function LoadMatrices_BinaryGraphLaplacian(Binary_list::Vector{BinaryInteraction},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,mode::AbstractMode=Ani(),corrected::Bool=true;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:AbstractFloat

    Bin_Norm = PhaseSpace.Characteristic.Bin_Norm
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        println("Filling M_Bin for binary interaction: $(name1) + $(name2) -> $(name3) + $(name4)")

        # Memory optimisation by allowing Ele and Pos populations to be modelled as identical thus only requiring one to be defined
        if !isnothing(findfirst(==("Pos"),[name1,name2,name3,name4])) && isnothing(findfirst(==("Pos"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho") # "Pos" compton scattering
                continue # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            elseif (name1,name2,name3,name4) == ("Ele","Pos","Pho","Pho") # "Ele" "Pos" annihilation
                GainScale = 1.0/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                LossScale = 1.0/4.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                LossScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            LossScale = 1.0
            name1_loc = findfirst(==(name1),name_list)
            name2_loc = findfirst(==(name2),name_list)
            name3_loc = findfirst(==(name3),name_list)
            name4_loc = findfirst(==(name4),name_list)
        end

        # ele pos swap for compton i.e. use ele compton matrices for pos compton
        if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho")
            name1 = "Ele"
            name3 = "Ele"
        end

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        px4_grid::String = px_grid_list[name4_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        py4_grid::String = py_grid_list[name4_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
        pz4_grid::String = pz_grid_list[name4_loc]

        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        px4_num::Int64 = px_num_list[name4_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        py4_num::Int64 = py_num_list[name4_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]
        pz4_num::Int64 = pz_num_list[name4_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]
        px4_low::Float64 = px_low_list[name4_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]
        px4_up::Float64 = px_up_list[name4_loc]

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]
        m4::Float64 = PhaseSpace.Grids.mass_list[name4_loc]

        Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)

        filename::String = BinaryFileName(Parameters)

        println(filename)

        Output = BinaryFileLoad_Matrix(DataDirectory,filename,corrected=corrected)
        
        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainScale *= Bin_Norm 
        LossScale *= Bin_Norm

        Parameters = Output[1]
        GainMatrix3 = Output[2] .* GainScale
        GainMatrix4 = Output[3] .* GainScale
        LossMatrix1 = Output[4] .* LossScale
        LossMatrix2 = Output[5] .* LossScale

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        DoesConserve(Output) # print conversion statistic
        Fill_M_BinGraphLaplacian!(name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,n_momentum;mode=mode,M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)


    end # for

end

function LoadMatrices_BinaryPatankar(Binary_list::Vector{BinaryInteraction},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,mode::AbstractMode=Ani(),corrected::Bool=true;Gijk::Union{Nothing,Matrix{F}}=nothing,Gijk_I::Union{Nothing,Vector{Int64}}=nothing,Gijk_J::Union{Nothing,Vector{Int64}}=nothing,Gijk_V::Union{Nothing,Vector{F}}=nothing,Lij::Union{Nothing,Matrix{F}}=nothing,Lij_I::Union{Nothing,Vector{Int64}}=nothing,Lij_J::Union{Nothing,Vector{Int64}}=nothing,Lij_V::Union{Nothing,Vector{F}}=nothing) where F<:AbstractFloat

    Bin_Norm = PhaseSpace.Characteristic.Bin_Norm
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        println("Filling M_Bin for binary interaction: $(name1) + $(name2) -> $(name3) + $(name4)")

        # Memory optimisation by allowing Ele and Pos populations to be modelled as identical thus only requiring one to be defined
        if !isnothing(findfirst(==("Pos"),[name1,name2,name3,name4])) && isnothing(findfirst(==("Pos"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho") # "Pos" compton scattering
                continue # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            elseif (name1,name2,name3,name4) == ("Ele","Pos","Pho","Pho") # "Ele" "Pos" annihilation
                GainScale = 1.0/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                LossScale = 1.0/4.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                LossScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            LossScale = 1.0
            name1_loc = findfirst(==(name1),name_list)
            name2_loc = findfirst(==(name2),name_list)
            name3_loc = findfirst(==(name3),name_list)
            name4_loc = findfirst(==(name4),name_list)
        end

        # ele pos swap for compton i.e. use ele compton matrices for pos compton
        if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho")
            name1 = "Ele"
            name3 = "Ele"
        end

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        px4_grid::String = px_grid_list[name4_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        py4_grid::String = py_grid_list[name4_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
        pz4_grid::String = pz_grid_list[name4_loc]

        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        px4_num::Int64 = px_num_list[name4_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        py4_num::Int64 = py_num_list[name4_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]
        pz4_num::Int64 = pz_num_list[name4_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]
        px4_low::Float64 = px_low_list[name4_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]
        px4_up::Float64 = px_up_list[name4_loc]

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]
        m4::Float64 = PhaseSpace.Grids.mass_list[name4_loc]

        Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)

        filename::String = BinaryFileName(Parameters)

        println(filename)

        Output = BinaryFileLoad_Matrix(DataDirectory,filename,corrected=corrected)
        
        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainScale *= Bin_Norm 
        LossScale *= Bin_Norm

        Parameters = Output[1]
        GainMatrix3 = Output[2] .* GainScale
        GainMatrix4 = Output[3] .* GainScale
        LossMatrix1 = Output[4] .* LossScale
        LossMatrix2 = Output[5] .* LossScale

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        DoesConserve(Output) # print conversion statistic
        Fill_M_BinPatankar!(name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,n_momentum;mode=mode,Gijk=Gijk,Gijk_I=Gijk_I,Gijk_J=Gijk_J,Gijk_V=Gijk_V,Lij=Lij,Lij_I=Lij_I,Lij_J=Lij_J,Lij_V=Lij_V)


    end # for

end

function LoadMatrices_BinaryPatankarSymmetric(Binary_list::Vector{BinaryInteraction},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,mode::AbstractMode=Ani(),corrected::Bool=true;Gijk::Union{Nothing,Matrix{F}}=nothing,Gijk_I::Union{Nothing,Vector{Int64}}=nothing,Gijk_J::Union{Nothing,Vector{Int64}}=nothing,Gijk_V::Union{Nothing,Vector{F}}=nothing,Lijk::Union{Nothing,Matrix{F}}=nothing,Lijk_I::Union{Nothing,Vector{Int64}}=nothing,Lijk_J::Union{Nothing,Vector{Int64}}=nothing,Lijk_V::Union{Nothing,Vector{F}}=nothing,Liij::Union{Nothing,Matrix{F}}=nothing,Liij_I::Union{Nothing,Vector{Int64}}=nothing,Liij_J::Union{Nothing,Vector{Int64}}=nothing,Liij_V::Union{Nothing,Vector{F}}=nothing,Liik::Union{Nothing,Matrix{F}}=nothing,Liik_I::Union{Nothing,Vector{Int64}}=nothing,Liik_J::Union{Nothing,Vector{Int64}}=nothing,Liik_V::Union{Nothing,Vector{F}}=nothing) where F<:AbstractFloat

    Bin_Norm = PhaseSpace.Characteristic.Bin_Norm
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    for i in eachindex(Binary_list)

        interaction = Binary_list[i]

        name1 = interaction.name1
        name2 = interaction.name2
        name3 = interaction.name3
        name4 = interaction.name4

        println("Filling M_Bin for binary interaction: $(name1) + $(name2) -> $(name3) + $(name4)")

        # Memory optimisation by allowing Ele and Pos populations to be modelled as identical thus only requiring one to be defined
        if !isnothing(findfirst(==("Pos"),[name1,name2,name3,name4])) && isnothing(findfirst(==("Pos"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho") # "Pos" compton scattering
                continue # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            elseif (name1,name2,name3,name4) == ("Ele","Pos","Pho","Pho") # "Ele" "Pos" annihilation
                GainScale = 1.0/4.0 # ele and pos populations are half the total ele population so scale gain matrix by 1/4
                LossScale = 1.0/4.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==("Ele"),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==(name4),name_list)
            elseif (name1,name2,name3,name4) == ("Pho","Pho","Ele","Pos") # "Ele" "Pos" pair production
                GainScale = 1.0
                LossScale = 1.0
                name1_loc = findfirst(==(name1),name_list)
                name2_loc = findfirst(==(name2),name_list)
                name3_loc = findfirst(==(name3),name_list)
                name4_loc = findfirst(==("Ele"),name_list)
            end
        else
            GainScale = 1.0
            LossScale = 1.0
            name1_loc = findfirst(==(name1),name_list)
            name2_loc = findfirst(==(name2),name_list)
            name3_loc = findfirst(==(name3),name_list)
            name4_loc = findfirst(==(name4),name_list)
        end

        # ele pos swap for compton i.e. use ele compton matrices for pos compton
        if (name1,name2,name3,name4) == ("Pos","Pho","Pos","Pho")
            name1 = "Ele"
            name3 = "Ele"
        end

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        px4_grid::String = px_grid_list[name4_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        py4_grid::String = py_grid_list[name4_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
        pz4_grid::String = pz_grid_list[name4_loc]

        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        px4_num::Int64 = px_num_list[name4_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        py4_num::Int64 = py_num_list[name4_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]
        pz4_num::Int64 = pz_num_list[name4_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]
        px4_low::Float64 = px_low_list[name4_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]
        px4_up::Float64 = px_up_list[name4_loc]

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]
        m4::Float64 = PhaseSpace.Grids.mass_list[name4_loc]

        Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,px4_low,px4_up,px4_grid,px4_num,py4_grid,py4_num,pz4_grid,pz4_num)

        filename::String = BinaryFileName(Parameters)

        println(filename)

        Output = BinaryFileLoad_Matrix(DataDirectory,filename,corrected=corrected)
        
        # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
        GainScale *= Bin_Norm 
        LossScale *= Bin_Norm

        Parameters = Output[1]
        GainMatrix3 = Output[2] .* GainScale
        GainMatrix4 = Output[3] .* GainScale
        LossMatrix1 = Output[4] .* LossScale
        LossMatrix2 = Output[5] .* LossScale

        name_locs = (name1_loc,name2_loc,name3_loc,name4_loc)

        DoesConserve(Output) # print conversion statistic
        Fill_M_BinPatankarSymmetric!(name_locs,PhaseSpace,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2,n_momentum;mode=mode,Gijk=Gijk,Gijk_I=Gijk_I,Gijk_J=Gijk_J,Gijk_V=Gijk_V,Lijk=Lijk,Lijk_I=Lijk_I,Lijk_J=Lijk_J,Lijk_V=Lijk_V,Liij=Liij,Liij_I=Liij_I,Liij_J=Liij_J,Liij_V=Liij_V,Liik=Liik,Liik_I=Liik_I,Liik_J=Liik_J,Liik_V=Liik_V)


    end # for

end
=#


function LoadMatrices_Emi!(M_Emi::Vector{Union{Matrix{F},SparseMatrixCSC{F,Int32}}},Emission_list::Vector{EmissiveInteraction},DataDirectory::String,PhaseSpace::PhaseSpaceStruct,BinM::BinaryMatricesStruct;Emi_corrected::Bool=true) where F<:AbstractFloat

    Emi_Norm = PhaseSpace.Characteristic.Emi_Norm

    # to define if MEmi should be sparse or dense for each spatial cell
    Binary_Domain = BinM.Domain
    Binary_list = BinM.Binary_list
    
    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Spacetime = PhaseSpace.Spacetime

    px_up_list = Momentum.px_up_list
    px_low_list = Momentum.px_low_list
    px_grid_list = Momentum.px_grid_list
    px_num_list = Momentum.px_num_list

    py_up_list = Momentum.py_up_list
    py_low_list = Momentum.py_low_list
    py_grid_list = Momentum.py_grid_list
    py_num_list = Momentum.py_num_list

    pz_up_list = Momentum.pz_up_list
    pz_low_list = Momentum.pz_low_list
    pz_grid_list = Momentum.pz_grid_list
    pz_num_list = Momentum.pz_num_list

    x_num = Spacetime.x_num
    y_num = Spacetime.y_num
    z_num = Spacetime.z_num

    n_momentum = sum(px_num_list.*py_num_list.*pz_num_list)

    M_Emi_D = zeros(F,n_momentum,n_momentum)
    M_Emi_I::Vector{Int32} = Int32[]
    M_Emi_J::Vector{Int32} = Int32[]
    M_Emi_V::Vector{F} = F[]

    for i in eachindex(Emission_list)

        interaction = Emission_list[i]

        name1::String = interaction.name1
        name2::String = interaction.name2
        name3::String = interaction.name3
        type::String = interaction.EmiName
        Ext_sampled::Vector{Float64} = interaction.Ext_sampled
        mode::AbstractMode = interaction.mode
        Force::Bool = interaction.Force
        Domain::Union{Vector{Int64},Nothing} = interaction.Domain

        println("Loading emission matrices for $type emission: $(name1) -> $(name2) + $(name3)")

        if !isnothing(findfirst(==("Pos"),[name1,name2,name3])) && isnothing(findfirst(==("Pos"),name_list)) # if "Pos" is in interactions but not in "name_list" ∴ "Ele" population is taken to be "Ele"+"Pos" with identical populations of each particle
            if (name1,name2,name3) == ("Pos","Pos","Pho") # "Pos" synchrotron 
                continue # skip loading pos compton matrices as this is correctly accounted by ele population including pos population
            end
        end

        # ele pos swap for synchrotron
        if name1 == "Pos" && name2 == "Pos" && name3 == "Pho"
            name1 = "Ele"
            name2 = "Ele"
            name3 = "Pho"
        end

        name1_loc = findfirst(==(name1),name_list)
        name2_loc = findfirst(==(name2),name_list)
        name3_loc = findfirst(==(name3),name_list)

        m1::Float64 = PhaseSpace.Grids.mass_list[name1_loc]
        m2::Float64 = PhaseSpace.Grids.mass_list[name2_loc]
        m3::Float64 = PhaseSpace.Grids.mass_list[name3_loc]

        z1::Float64 = PhaseSpace.Grids.charge_list[name1_loc]
        z2::Float64 = PhaseSpace.Grids.charge_list[name2_loc]
        z3::Float64 = PhaseSpace.Grids.charge_list[name3_loc]

        px1_grid::String = px_grid_list[name1_loc]
        px2_grid::String = px_grid_list[name2_loc]
        px3_grid::String = px_grid_list[name3_loc]
        py1_grid::String = py_grid_list[name1_loc]
        py2_grid::String = py_grid_list[name2_loc]
        py3_grid::String = py_grid_list[name3_loc]
        pz1_grid::String = pz_grid_list[name1_loc]
        pz2_grid::String = pz_grid_list[name2_loc]
        pz3_grid::String = pz_grid_list[name3_loc]
    
        px1_num::Int64 = px_num_list[name1_loc]
        px2_num::Int64 = px_num_list[name2_loc]
        px3_num::Int64 = px_num_list[name3_loc]
        py1_num::Int64 = py_num_list[name1_loc]
        py2_num::Int64 = py_num_list[name2_loc]
        py3_num::Int64 = py_num_list[name3_loc]
        pz1_num::Int64 = pz_num_list[name1_loc]
        pz2_num::Int64 = pz_num_list[name2_loc]
        pz3_num::Int64 = pz_num_list[name3_loc]

        px1_low::Float64 = px_low_list[name1_loc]
        px2_low::Float64 = px_low_list[name2_loc]
        px3_low::Float64 = px_low_list[name3_loc]

        px1_up::Float64 = px_up_list[name1_loc]
        px2_up::Float64 = px_up_list[name2_loc]
        px3_up::Float64 = px_up_list[name3_loc]

        # build Vector of Emission Gain matrices for each sampled value of Ext
        GainMatrix3_All = Vector{Array{Float64,6}}(undef,length(Ext_sampled))

        for Ext_idx in eachindex(Ext_sampled)

            Parameters = (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,px1_low,px1_up,px1_grid,px1_num,py1_grid,py1_num,pz1_grid,pz1_num,px2_low,px2_up,px2_grid,px2_num,py2_grid,py2_num,pz2_grid,pz2_num,px3_low,px3_up,px3_grid,px3_num,py3_grid,py3_num,pz3_grid,pz3_num,Ext_sampled[Ext_idx])

            filename = EmissionFileName(Parameters)
            println(filename)

            #Parameters = DC.fload_Matrix_Sync(DataDirectory,filename)[1] # 1 is Parameters
            #matrix = DC.fload_Matrix_Sync(DataDirectory,filename)[2] # 1 is Parameters
            GainMatrix3_All[Ext_idx] = EmissionFileLoad_Matrix(DataDirectory,filename)[2] # remove later

            # apply correct scaling from DiplodocusCollisions to non-dimensionalisation of DiplodocusTransport
            GainMatrix3_All[Ext_idx] .*= Emi_Norm 

            # apply energy correction
            if Emi_corrected
                EmissionCorrection!(PhaseSpace,GainMatrix3_All[Ext_idx],Parameters)
            end

        end   

        println("Filling M_Emi for $type emission: $(name1) -> $(name2) + $(name3)")

        name_locs = (name1_loc,name2_loc,name3_loc)

        for x in 1:x_num, y in 1:y_num, z in 1:z_num

            off_space = (x-1)*y_num*z_num+(y-1)*z_num+z-1

            in_Binary = !isempty(Binary_list) && (isnothing(Binary_Domain) || in(off_space,Binary_Domain)) 

            if isnothing(Domain) || in(off_space,Domain)

                # reset temporary arrays
                if in_Binary 
                    fill!(M_Emi_D,zero(F))
                else
                    empty!(M_Emi_I)
                    empty!(M_Emi_J)
                    empty!(M_Emi_V)
                end

                B_field = PhaseSpace.Grids.B_field[x,y,z]
                Ext_idx = findmin(abs.(Ext_sampled .- B_field))[2]

                println("x: $x, y: $y, z: $z, Ext_idx: $Ext_idx, Ext_val: $(Ext_sampled[Ext_idx])")

                if type=="Sync" && Force 
                    
                    force = SyncRadReact(mode=mode,B=Ext_sampled[Ext_idx])

                    if in_Binary
                    Fill_I_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi=M_Emi_D)
                    Fill_J_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi=M_Emi_D)
                    Fill_K_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi=M_Emi_D)
                    else
                    Fill_I_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
                    Fill_J_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
                    Fill_K_Emi!(PhaseSpace,force,x,y,z,name1_loc;M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
                    end

                end

                GainMatrix3 = GainMatrix3_All[Ext_idx]

                # Fill_M_Emi! is called for each spatial grid point as the emission correction is dependent on space through the electromagnetic fields
                if in_Binary
                    Fill_M_Emi!(PhaseSpace,name_locs,x,y,z;GainMatrix3=GainMatrix3,mode=mode,M_Emi=M_Emi_D)
                else
                    Fill_M_Emi!(PhaseSpace,name_locs,x,y,z;GainMatrix3=GainMatrix3,mode=mode,M_Emi_I=M_Emi_I,M_Emi_J=M_Emi_J,M_Emi_V=M_Emi_V)
                end

                # add values to Vector of Matrices for each spatial grid point
                if isassigned(M_Emi,off_space+1)
                    if in_Binary
                        M_Emi[off_space+1] += M_Emi_D
                    else
                        M_Emi[off_space+1] += sparse(M_Emi_I,M_Emi_J,M_Emi_V,n_momentum,n_momentum)
                    end
                else
                    if in_Binary
                        M_Emi[off_space+1] = zeros(F,n_momentum,n_momentum)
                        copyto!(M_Emi[off_space+1], M_Emi_D)
                    else
                        M_Emi[off_space+1] = sparse(M_Emi_I,M_Emi_J,M_Emi_V,n_momentum,n_momentum)
                    end
                end

            else
                continue
            end

        end # space loop

    end # for

end

