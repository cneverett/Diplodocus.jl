function UserBinaryParameters()

    name1::String = getfield(Main,Symbol("name1"))
    name2::String = getfield(Main,Symbol("name2"))
    name3::String = getfield(Main,Symbol("name3"))
    name4::String = getfield(Main,Symbol("name4"))

    p1_low::Float64 = getfield(Main,Symbol("p_low_"*name1))
    p1_up::Float64 = getfield(Main,Symbol("p_up_"*name1))
    p1_grid::String = getfield(Main,Symbol("p_grid_"*name1))
    p1_num::Int64 = getfield(Main,Symbol("p_num_"*name1))

    u1_grid::String = getfield(Main,Symbol("u_grid_"*name1))
    u1_num::Int64 = getfield(Main,Symbol("u_num_"*name1))
    
    h1_grid::String = getfield(Main,Symbol("h_grid_"*name1))
    h1_num::Int64 = getfield(Main,Symbol("h_num_"*name1))

    p2_low::Float64 = getfield(Main,Symbol("p_low_"*name2))
    p2_up::Float64 = getfield(Main,Symbol("p_up_"*name2))
    p2_grid::String = getfield(Main,Symbol("p_grid_"*name2))
    p2_num::Int64 = getfield(Main,Symbol("p_num_"*name2))

    u2_grid::String = getfield(Main,Symbol("u_grid_"*name2))
    u2_num::Int64 = getfield(Main,Symbol("u_num_"*name2))

    h2_grid::String = getfield(Main,Symbol("h_grid_"*name2))
    h2_num::Int64 = getfield(Main,Symbol("h_num_"*name2))

    p3_low::Float64 = getfield(Main,Symbol("p_low_"*name3))
    p3_up::Float64 = getfield(Main,Symbol("p_up_"*name3))
    p3_grid::String = getfield(Main,Symbol("p_grid_"*name3))
    p3_num::Int64 = getfield(Main,Symbol("p_num_"*name3))

    u3_grid::String = getfield(Main,Symbol("u_grid_"*name3))
    u3_num::Int64 = getfield(Main,Symbol("u_num_"*name3))

    h3_grid::String = getfield(Main,Symbol("h_grid_"*name3))
    h3_num::Int64 = getfield(Main,Symbol("h_num_"*name3))

    p4_low::Float64 = getfield(Main,Symbol("p_low_"*name4))
    p4_up::Float64 = getfield(Main,Symbol("p_up_"*name4))
    p4_grid::String = getfield(Main,Symbol("p_grid_"*name4))
    p4_num::Int64 = getfield(Main,Symbol("p_num_"*name4))

    u4_grid::String = getfield(Main,Symbol("u_grid_"*name4))
    u4_num::Int64 = getfield(Main,Symbol("u_num_"*name4))

    h4_grid::String = getfield(Main,Symbol("h_grid_"*name4))
    h4_num::Int64 = getfield(Main,Symbol("h_num_"*name4))

    m1::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name1))
    m2::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name2))
    m3::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name3))
    m4::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name4))

    scale::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = getfield(Main,Symbol("scale"))

    p1loc_low::Int64 = try 
        getfield(Main,Symbol("ploc_low_"*name1))
        catch
        1
        end
    p1loc_up::Int64 = try 
        getfield(Main,Symbol("ploc_up_"*name1))
        catch
        p1_num
        end
    p2loc_low::Int64 = try 
        getfield(Main,Symbol("ploc_low_"*name2))
        catch
        1
        end
    p2loc_up::Int64 = try 
        getfield(Main,Symbol("ploc_up_"*name2))
        catch
        p2_num
        end

    Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num)

    bins = (p1loc_low,p1loc_up,p2loc_low,p2loc_up)

    numLoss::Int64 = getfield(Main,Symbol("numLoss"))   
    numGain::Int64 = getfield(Main,Symbol("numGain"))  
    numThreads::Int64 = getfield(Main,Symbol("numThreads"))

    fileLocation::String = try 
        getfield(Main,Symbol("fileLocation"))
    catch
        joinpath(pwd(),"Data")
    end
    fileName::String = BinaryFileName(Parameters)

    Setup::Tuple = (Parameters,bins,scale, numLoss, numGain, numThreads, fileLocation, fileName)

    return (Setup, fileName)
end

"""
    FileName(Parameters)

Generates a file name for a binary interaction based on user provided parameters
"""
function BinaryFileName(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64})
    
    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    @inline function stringval(val::Float64)
        if val < 0.0
            return "m"*string(round(Int64,abs(val)))
        else
            return "p"*string(round(Int64,val))
        end
    end
    
    fileName = name1*name2*name3*name4
    fileName *= "_"*stringval(p1_low)*stringval(p1_up)*p1_grid*string(p1_num)
    fileName *= u1_grid*string(u1_num)
    fileName *= h1_grid*string(h1_num)

    fileName *= "_"*stringval(p2_low)*stringval(p2_up)*p2_grid*string(p2_num)
    fileName *= u2_grid*string(u2_num)
    fileName *= h2_grid*string(h2_num)

    fileName *= "_"*stringval(p3_low)*stringval(p3_up)*p3_grid*string(p3_num)
    fileName *= u3_grid*string(u3_num)
    fileName *= h3_grid*string(h3_num)

    fileName *= "_"*stringval(p4_low)*stringval(p4_up)*p4_grid*string(p4_num)
    fileName *= u4_grid*string(u4_num)
    fileName *= h4_grid*string(h4_num)

    return fileName
end

"""
    EmissionFileName(Parameters)

Generates a file name for an emissive interaction based on user provided parameters
"""
function EmissionFileName(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64, Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64})
    
    (name1,name2,name3,type,mu1,mu2,mu3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext) = Parameters

    @inline function stringval(val::Float64)
        if val < 0.0
            return "m"*string(round(Int64,abs(val)))
        else
            return "p"*string(round(Int64,val))
        end
    end
    
    fileName = type*name1*name2*name3
    fileName *= "_"*stringval(p1_low)*stringval(p1_up)*p1_grid*string(p1_num)
    fileName *= u1_grid*string(u1_num)
    fileName *= h1_grid*string(h1_num)

    fileName *= "_"*stringval(p2_low)*stringval(p2_up)*p2_grid*string(p2_num)
    fileName *= u2_grid*string(u2_num)
    fileName *= h2_grid*string(h2_num)

    fileName *= "_"*stringval(p3_low)*stringval(p3_up)*p3_grid*string(p3_num)
    fileName *= u3_grid*string(u3_num)
    fileName *= h3_grid*string(h3_num)

    fileName *= "_"*string(Ext)
    
    fileName *= ".jld2";

    return fileName
end
#

"""
    UserEmissionParameters()

Returns the integration `Setup` and `fileName` for an emissive integration based on user provided parameters defined in `Main`.
"""
function UserEmissionParameters(Ext::Float64)

    type::String = getfield(Main,Symbol("type"))
    name1::String = getfield(Main,Symbol("name1"))
    name2::String = getfield(Main,Symbol("name2"))
    name3::String = getfield(Main,Symbol("name3"))

    p1_low::Float64 = getfield(Main,Symbol("p_low_"*name1))
    p1_up::Float64 = getfield(Main,Symbol("p_up_"*name1))
    p1_grid::String = getfield(Main,Symbol("p_grid_"*name1))
    p1_num::Int64 = getfield(Main,Symbol("p_num_"*name1))

    u1_grid::String = getfield(Main,Symbol("u_grid_"*name1))
    u1_num::Int64 = getfield(Main,Symbol("u_num_"*name1))
    
    h1_grid::String = getfield(Main,Symbol("h_grid_"*name1))
    h1_num::Int64 = getfield(Main,Symbol("h_num_"*name1))

    p2_low::Float64 = getfield(Main,Symbol("p_low_"*name2))
    p2_up::Float64 = getfield(Main,Symbol("p_up_"*name2))
    p2_grid::String = getfield(Main,Symbol("p_grid_"*name2))
    p2_num::Int64 = getfield(Main,Symbol("p_num_"*name2))

    u2_grid::String = getfield(Main,Symbol("u_grid_"*name2))
    u2_num::Int64 = getfield(Main,Symbol("u_num_"*name2))

    h2_grid::String = getfield(Main,Symbol("h_grid_"*name2))
    h2_num::Int64 = getfield(Main,Symbol("h_num_"*name2))

    p3_low::Float64 = getfield(Main,Symbol("p_low_"*name3))
    p3_up::Float64 = getfield(Main,Symbol("p_up_"*name3))
    p3_grid::String = getfield(Main,Symbol("p_grid_"*name3))
    p3_num::Int64 = getfield(Main,Symbol("p_num_"*name3))

    u3_grid::String = getfield(Main,Symbol("u_grid_"*name3))
    u3_num::Int64 = getfield(Main,Symbol("u_num_"*name3))

    h3_grid::String = getfield(Main,Symbol("h_grid_"*name3))
    h3_num::Int64 = getfield(Main,Symbol("h_num_"*name3))

    m1::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name1))
    m2::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name2))
    m3::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name3))

    z1::Float64 = getfield(DiplodocusCollisions,Symbol("z"*name1))
    z2::Float64 = getfield(DiplodocusCollisions,Symbol("z"*name2))
    z3::Float64 = getfield(DiplodocusCollisions,Symbol("z"*name3))

    scale::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = getfield(Main,Symbol("scale"))

    println("Ext = $Ext")

    p1loc_low::Int64 = try 
        getfield(Main,Symbol("ploc_low_"*name1))
        catch
        1
        end
    p1loc_up::Int64 = try 
        getfield(Main,Symbol("ploc_up_"*name1))
        catch
        p1_num
        end
    p3loc_low::Int64 = try 
        getfield(Main,Symbol("ploc_low_"*name3))
        catch
        1
        end
    p3loc_up::Int64 = try 
        getfield(Main,Symbol("ploc_up_"*name3))
        catch
        p3_num
        end

    Parameters = (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext)

    bins = (p1loc_low,p1loc_up,p3loc_low,p3loc_up)

    numLoss::Int64 = getfield(Main,Symbol("numLoss"))   
    numGain::Int64 = getfield(Main,Symbol("numGain"))  
    numThreads::Int64 = getfield(Main,Symbol("numThreads"))

    fileLocation::String = try 
        getfield(Main,Symbol("fileLocation"))
    catch
        joinpath(pwd(),"Data")
    end
    fileName::String = EmissionFileName(Parameters)

    Setup::Tuple = (Parameters,bins,scale, numLoss, numGain, numThreads, fileLocation, fileName)

    return (Setup, fileName)
end