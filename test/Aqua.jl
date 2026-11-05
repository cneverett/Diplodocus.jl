using Aqua 

@testset "Aqua.jl" begin
    Aqua.test_all(
        DiplodocusCollisions;
        deps_compat = (check_extras = false,)
    )
end