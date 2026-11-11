abstract type BinaryInteraction end
struct SphSphSphSph <: BinaryInteraction end

struct BinaryParameters

    name1::String
    name2::String
    name3::String
    name4::String

    mu1::Float64
    mu2::Float64
    mu3::Float64
    mu4::Float64

    p1_low::Float64
    p1_up::Float64
    p1_grid::String
    p1_num::Int64

    u1_grid::String
    u1_num::Int64

    h1_grid::String
    h1_num::Int64

    p2_low::Float64
    p2_up::Float64
    p2_grid::String
    p2_num::Int64

    u2_grid::String
    u2_num::Int64

    h2_grid::String
    h2_num::Int64

    p3_low::Float64
    p3_up::Float64
    p3_grid::String
    p3_num::Int64

    u3_grid::String
    u3_num::Int64

    h3_grid::String
    h3_num::Int64

    p4_low::Float64
    p4_up::Float64
    p4_grid::String
    p4_num::Int64

    u4_grid::String
    u4_num::Int64

    h4_grid::String
    h4_num::Int64

end

mutable struct BinaryUserInput <: Function

    Parameters::BinaryParameters
    numTiter::Int64
    numSiter::Int64
    fileLocation::String
    fileName::String
    numThreads::Int64
    MinMax::Bool

    function BinaryUserInput(Parameters::BinaryParameters,fileLoc::String,numT::Int64,numS::Int64;numThreads=1,MinMax=false)
        self = new()

        self.Parameters = Parameters
        self.numTiter = numT
        self.numSiter = numS
        self.fileName = FileName(Parameters)
        self.fileLocation = fileLoc
        self.numThreads = numThreads
        self.MinMax = MinMax

        return self

    end

end

mutable struct MonteCarloArrays <: Function

    GainTotal3::Array{Float64,9}
    GainTotal4::Array{Float64,9}

    GainTally3::Array{UInt32,9}
    GainTally4::Array{UInt32,9}

    LossTotal::Array{Float64,6}
    LossTally::Array{UInt32,6}

    GainMatrix3::Array{Float64,9}
    GainMatrix4::Array{Float64,9}
    LossMatrix1::Array{Float64,6}
    LossMatrix2::Array{Float64,6}

    function MonteCarloArrays(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64})

        self = new()

        #=filePath = UserInput.fileLocation*"\\"*UserInput.fileName
        fileExist = isfile(filePath)
        MinMax = UserInput.MinMax
        P = UserInput.Parameters=#
  
        (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters


        self.GainTotal3 = zeros(Float64,(p3_num+1),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num); 
        self.GainTotal4 = zeros(Float64,(p4_num+1),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num); 
        self.LossTotal = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        # GainTally have first dimension elements [k1,k2,k3,...,k(n+1),N]
        self.GainTally3 = zeros(UInt32,(p3_num+2),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        self.GainTally4 = zeros(UInt32,(p4_num+2),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        self.LossTally = zeros(UInt32,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);

        self.GainMatrix3 = zeros(Float64,(p3_num+1),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        self.GainMatrix4 = zeros(Float64,(p4_num+1),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        self.LossMatrix1 = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        self.LossMatrix2 = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);

        return self

    end

end

mutable struct OldMonteCarloArrays <: Function

    GainTally3::Array{UInt32,9}
    GainTally4::Array{UInt32,9}

    GainMatrix3::Array{Float64,9}
    GainMatrix4::Array{Float64,9}

    LossTally::Array{UInt32,6}
    LossMatrix1::Array{Float64,6}
    LossMatrix2::Array{Float64,6}

    function OldMonteCarloArrays(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},filePath::String)

        self = new()

        #=filePath = UserInput.fileLocation*"\\"*UserInput.fileName
        fileExist = isfile(filePath)
        MinMax = UserInput.MinMax
        P = UserInput.Parameters=#
  
        (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters
        #filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            self.GainTally3 = f["GainTally3"];
            self.GainMatrix3 = f["GainMatrix3"];
            self.GainTally4 = f["GainTally4"];
            self.GainMatrix4 = f["GainMatrix4"];
            self.LossTally = f["LossTally"];
            self.LossMatrix1 = f["LossMatrix1"];
            self.LossMatrix2 = f["LossMatrix2"];
            close(f)
        else
            #=p1_num = P.p1_num
            u1_num = P.u1_num
            h1_num = P.h1_num
            p2_num = P.p2_num
            u2_num = P.u2_num
            h2_num = P.h2_num
            p3_num = P.p3_num
            u3_num = P.u3_num
            h3_num = P.h3_num
            p4_num = P.p4_num
            u4_num = P.u4_num
            h4_num = P.h4_num=#

            # GainTally have first dimension elements [k1,k2,k3,...,k(n+1)], N is not needed to be saved
            self.GainTally3 = zeros(UInt32,(p3_num+1),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.GainTally4 = zeros(UInt32,(p4_num+1),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.LossTally = zeros(UInt32,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.GainMatrix3 = zeros(Float64,(p3_num+1),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.GainMatrix4 = zeros(Float64,(p4_num+1),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.LossMatrix1 = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
            self.LossMatrix2 = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
        end

        return self

    end

end