using Aqua 

@testset "Aqua.jl" begin
    Aqua.test_all(
        DiplodocusTransport;
        deps_compat = (check_extras = false,)
    )
end