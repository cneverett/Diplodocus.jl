using Test
using Pkg
using Diplodocus

@testset "Diplodocus and Sub-packages" begin 

    @testset "Diplodocus" begin
        # Write your tests here.
        include("Aqua.jl")
    end

    @testset "DiplodocusCollisions" begin 
        Pkg.activate(joinpath(@__DIR__, "..", "lib", "DiplodocusCollisions"))
        Pkg.test()
    end

    @testset "DiplodocusTransport" begin 
        Pkg.activate(joinpath(@__DIR__, "..", "lib", "DiplodocusTransport"))
        Pkg.test()
    end

    @testset "DiplodocusPlots" begin 
        Pkg.activate(joinpath(@__DIR__, "..", "lib", "DiplodocusPlots"))
        Pkg.test()
    end

end
