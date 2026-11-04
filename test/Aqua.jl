using Aqua 

@testset "Aqua.jl" begin
    Aqua.test_all(
        DiplodocusPlots;
        deps_compat = (check_extras = false,)
    )
end