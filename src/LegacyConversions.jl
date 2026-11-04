function Base.convert(::Type{AbstractTetrad},x::JLD2.ReconstructedMutable{Symbol("ParabolicForceFreeFieldTetrad{Float64}"),(:B0, :Ω, :Bfunction),Tuple{Float64, Any, Any}})
    Ω_old = x.Ω
    Bfun_old = x.Bfunction
    return ParabolicForceFreeFieldTetrad{Float64,typeof(Ω_old),typeof(Bfun_old)}(Float64(x.B0),Ω_old,Bfun_old)
end