using Aqua 

@testset "Aqua.jl" begin
    Aqua.test_all(
        Diplodocus;
        deps_compat = (check_extras = false,)
    )
end