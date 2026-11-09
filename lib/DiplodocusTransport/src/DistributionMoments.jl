"""
    FourFlow(state,PhaseSpace,species)

Returns the contravariant (upper indices) four-flow vector Nᵃ `AbstractVector{T}` for a given particle `species` given a `state` vector. The four-flow is calculated as the first moment of the distribution function over momentum space. 

Takes additional arguments `x_idx`, `y_idx`, and `z_idx` to specify the spatial location in the phase space to evaluate the four-flow.
"""
function FourFlow(state::AbstractVector{T},PhaseSpace::PhaseSpaceStruct,species::S;x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where T<:AbstractFloat where S<:Union{String,Int64}

    if typeof(species) == Int64
        species_index = species
    else
        name_list = PhaseSpace.name_list
        species_index = findfirst(==(species),name_list)
    end

    f1D = copy(LocationSpeciesToStateVector(state,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

    m = PhaseSpace.Grids.mass_list[species_index]
    p_num = PhaseSpace.Momentum.px_num_list[species_index]
    u_num = PhaseSpace.Momentum.py_num_list[species_index]
    h_num = PhaseSpace.Momentum.pz_num_list[species_index]
    pr = PhaseSpace.Grids.pxr_list[species_index]
    ur = PhaseSpace.Grids.pyr_list[species_index]
    hr = PhaseSpace.Grids.pzr_list[species_index]

    f3D = T.(reshape(f1D,(p_num,u_num,h_num)))

    Nᵃ = zeros(T,4)
    
    for p in 1:p_num, u in 1:u_num, h in 1:h_num

        Nᵃ[1] += f3D[p,u,h]
        Nᵃ[2] += f3D[p,u,h]*((sqrt(m^2 + pr[p]^2) - sqrt(m^2 + pr[p+1]^2)) * (ur[u]*sqrt(1 - ur[u]^2) - ur[u+1]*sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1])) * (sinpi(hr[h]/pi) - sinpi(hr[h+1]/pi))) / (2(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
        Nᵃ[3] += f3D[p,u,h]*((sqrt(m^2 + pr[p]^2) - sqrt(m^2 + pr[p+1]^2)) * (ur[u]*sqrt(1 - ur[u]^2) - ur[u+1]*sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1])) * (cospi(hr[h]/pi) - cospi(hr[h+1]/pi))) / (2(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
        Nᵃ[4] += f3D[p,u,h]*(sqrt(pr[p+1]^2+m^2) - sqrt(pr[p]^2+m^2)) * (ur[u+1]+ur[u]) / (pr[p+1]-pr[p]) / 2

    end

    return Nᵃ

end

function HydroFourVelocity(Nᵃ::AbstractVector{T}) where T

    Uᵃ = zeros(T,4)

    metric = zeros(T,4,4)
    metric[1,1] = T(-1)
    metric[2,2] = T(1)
    metric[3,3] = T(1)
    metric[4,4] = T(1)

    sqrtNᵃ2 = sqrt(abs(Nᵃ'*metric*Nᵃ))
    Uᵃ .= Nᵃ/sqrtNᵃ2

    return Uᵃ

end

function HydroThreeVelocity(Uᵃ::AbstractVector{T}) where T

    Vⁱ = zeros(T,3)

    Vⁱ .= Uᵃ[2:4] / Uᵃ[1]
    return Vⁱ

end

function ProjectionTensor(nₐ::AbstractVector{T}) where T

    Δab = zeros(T,4,4)

    metric = zeros(T,4,4)
    metric[1,1] = T(-1)
    metric[2,2] = T(1)
    metric[3,3] = T(1)
    metric[4,4] = T(1)

    for i in 1:4, j in 1:4
        Δab[i,j] = metric[i,j] + nₐ[i]*nₐ[j]
    end

    return Δab

end


"""
    StressEnergyTensor(state,PhaseSpace,species)

Returns the contravariant (upper indices) stress-energy tensor Tab for a given particle `species` given a `state` vector. The stress-energy is calculated as the second moment of the distribution function over momentum space.
"""
function StressEnergyTensor(state::Vector{T},PhaseSpace::PhaseSpaceStruct,species::S;x_idx::Int64=1,y_idx::Int64=1,z_idx::Int64=1) where T<:AbstractFloat where S<:Union{String,Int64}

    if typeof(species) == Int64
        species_index = species
    else
        name_list = PhaseSpace.name_list
        species_index = findfirst(==(species),PhaseSpace.name_list)
    end

    f1D = copy(LocationSpeciesToStateVector(state,PhaseSpace,species_index=species_index,x_idx=x_idx,y_idx=y_idx,z_idx=z_idx))

    m = PhaseSpace.Grids.mass_list[species_index]
    p_num = PhaseSpace.Momentum.px_num_list[species_index]
    u_num = PhaseSpace.Momentum.py_num_list[species_index]
    h_num = PhaseSpace.Momentum.pz_num_list[species_index]
    pr = PhaseSpace.Grids.pxr_list[species_index]
    ur = PhaseSpace.Grids.pyr_list[species_index]
    hr = PhaseSpace.Grids.pzr_list[species_index]

    f3D = T.(reshape(f1D,(p_num,u_num,h_num)))

    Tᵃᵇ = zeros(T,4,4)

    #=for i in 1:u_num
        du[i] = ur[i+1]-ur[i]
        du2[i] = (ur[i+1]^2-ur[i]^2)/2
        du3[i] = (ur[i+1]^3-ur[i]^3)/3
        duplusu3[i] = ur[i+1]-ur[i] - (ur[i+1]^3-ur[i]^3)/3
    end
    for i in 1:p_num 
        dp[i] = (pr[i+1]-pr[i])
        dp2[i] = (pr[i+1]^2-pr[i]^2)/2
        if m == 0e0
            dpfunc[i] = pr64[i+1]^2 - pr64[i]^2
        else
            if pr[i+1]/m < 1e-3
                dpfunc[i] = 2/3 * pr[i+1]^3/m - 1/5 * pr[i+1]^5/m^3
                dpfunc[i] -= 2/3 * pr[i]^3/m - 1/5 * pr[i]^5/m^3
            else
                dpfunc[i] = pr64[i+1]*sqrt(pr64[i+1]^2+m^2) - m^2*atanh(pr64[i+1]/sqrt(pr64[i+1]^2+m^2))
                dpfunc[i] -= pr64[i]*sqrt(pr64[i]^2+m^2) - m^2*atanh(pr64[i]/sqrt(pr64[i]^2+m^2))
            end
        end

    end
    =#

    for p in 1:p_num, u in 1:u_num, h in 1:h_num

        if m != 0e0

            Tᵃᵇ[1,1] += f3D[p,u,h]*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2*(asinh(pr[p]/m) - asinh(pr[p+1]/m)))/(2(pr[p] - pr[p+1]))

            Tᵃᵇ[1,2] += f3D[p,u,h] * (pr[p] + pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(sin(hr[h]) - sin(hr[h+1])) / (4(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
            Tᵃᵇ[1,3] += f3D[p,u,h] * (-pr[p] - pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(cos(hr[h]) - cos(hr[h+1])) / (4(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
            Tᵃᵇ[1,4] += f3D[p,u,h] * (1/4 * (pr[p] + pr[p+1])*(ur[u] + ur[u+1]))
            Tᵃᵇ[2,2] += f3D[p,u,h] * (3 - ur[u]^2 - ur[u]ur[u+1] - ur[u+1]^2) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) * (hr[h] - hr[h+1] + cos(hr[h])sin(hr[h]) - cos(hr[h+1])sin(hr[h+1])) / (12(pr[p] - pr[p+1])*(hr[h] - hr[h+1]))

            Tᵃᵇ[2,3] += f3D[p,u,h] * (-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) * (cos(2hr[h]) - cos(2hr[h+1])) / (24(pr[p] - pr[p+1])*(hr[h] - hr[h+1]))
            Tᵃᵇ[2,4] += f3D[p,u,h] * (-sqrt(1 - ur[u]^2) + ur[u]^2 * sqrt(1 - ur[u]^2) + sqrt(1 - ur[u+1]^2) - ur[u+1]^2 * sqrt(1 - ur[u+1]^2)) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) * (sin(hr[h]) - sin(hr[h+1])) / (6(pr[p] - pr[p+1]) * (ur[u] - ur[u+1]) * (hr[h] - hr[h+1]))

            Tᵃᵇ[3,3] += f3D[p,u,h] * (3 - ur[u]^2 - ur[u]ur[u+1] - ur[u+1]^2) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) * (hr[h] - hr[h+1] - cos(hr[h])sin(hr[h]) + cos(hr[h+1])sin(hr[h+1])) / (12(pr[p] - pr[p+1])*(hr[h] - hr[h+1]))

            Tᵃᵇ[3,4] += f3D[p,u,h] * (sqrt(1 - ur[u]^2) - ur[u]^2 * sqrt(1 - ur[u]^2) - sqrt(1 - ur[u+1]^2) + ur[u+1]^2 * sqrt(1 - ur[u+1]^2)) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) * (cos(hr[h]) - cos(hr[h+1])) / (6(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))

            Tᵃᵇ[4,4] += f3D[p,u,h] * (ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2) * (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m))) / (6(pr[p] - pr[p+1]))

        else
            
            Tᵃᵇ[1,1] += f3D[p,u,h]*(pr[p] + pr[p+1])/2

            Tᵃᵇ[1,2] += f3D[p,u,h]*((pr[p] + pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(sin(hr[h]) - sin(hr[h+1])))/(4(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
            Tᵃᵇ[1,3] += f3D[p,u,h]*(-((pr[p] + pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(cos(hr[h]) - cos(hr[h+1])))/(4(ur[u] - ur[u+1])*(hr[h] - hr[h+1])))
            Tᵃᵇ[1,4] += f3D[p,u,h]*(1/4 * (pr[p] + pr[p+1])*(ur[u] + ur[u+1]))
            Tᵃᵇ[2,2] += f3D[p,u,h]*(1/(12(hr[h] - hr[h+1])))*(pr[p] + pr[p+1])*(-3 + ur[u]^2 + ur[u]*ur[u+1] + ur[u+1]^2)*(-hr[h] + hr[h+1] - cos(hr[h])sin(hr[h]) + cos(hr[h+1])sin(hr[h+1]))

            Tᵃᵇ[2,3] += f3D[p,u,h]*((pr[p] + pr[p+1])*(-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)*(cos(2hr[h]) - cos(2hr[h+1])))/(24(hr[h] - hr[h+1]))
            Tᵃᵇ[2,4] += f3D[p,u,h]*(-(((pr[p] + pr[p+1]) * ((1 - ur[u]^2)^(3/2) - (1 - ur[u+1]^2)^(3/2)) * (sin(hr[h]) - sin(hr[h+1])))/(6(ur[u] - ur[u+1]) * (hr[h] - hr[h+1]))))

            Tᵃᵇ[3,3] += f3D[p,u,h]*(-(1/(12(hr[h] - hr[h+1])))*(pr[p] + pr[p+1]) * (-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2) * (hr[h] - hr[h+1] - cos(hr[h])sin(hr[h]) + cos(hr[h+1])sin(hr[h+1])))

            Tᵃᵇ[3,4] += f3D[p,u,h]*(((pr[p] + pr[p+1]) * ((1 - ur[u]^2)^(3/2) - (1 - ur[u+1]^2)^(3/2)) * (cos(hr[h]) - cos(hr[h+1])))/(6(ur[u] - ur[u+1]) * (hr[h] - hr[h+1])))

            Tᵃᵇ[4,4] += f3D[p,u,h]*(1/6) * (pr[p] + pr[p+1]) * (ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)

        end

    end

    Tᵃᵇ[2,1] = Tᵃᵇ[1,2]
    Tᵃᵇ[3,1] = Tᵃᵇ[1,3]
    Tᵃᵇ[4,1] = Tᵃᵇ[1,4]

    Tᵃᵇ[3,2] = Tᵃᵇ[2,3]
    Tᵃᵇ[4,2] = Tᵃᵇ[2,4]

    Tᵃᵇ[4,3] = Tᵃᵇ[3,4]

    # unscale by dp*du
    #for i in axes(f2D,1), j in axes(f2D,2)
    #    f2D[i,j] /= dp[i] * du[j]
    #end

    #Tab[1,1] = dE' * f2D * du
    #Tab[2,2] = 1/4 * dpfunc' * f2D * duplusu3
    #Tab[3,3] = Tab[2,2]
    #Tab[4,4] = 1/2 * dpfunc' * f2D * du3 
    #Tab[1,4] = dp2' * f2D * du2
    #Tab[4,1] = Tab[1,4]

    # unscale by dp*du
    #Tab[1,1] = (dE./dp)' * f2D * (du./du)
    #Tab[2,2] = 1/4 * (dpfunc./dp)' * f2D * (duplusu3./du)
    #Tab[3,3] = Tab[2,2]
    #Tab[4,4] = 1/2 * (dpfunc./dp)' * f2D * (du3./du) 
    #Tab[1,4] = (dp2./dp)' * f2D * (du2./du)
    #Tab[4,1] = Tab[1,4]
    
    return Tᵃᵇ
end

function ScalarNumberDensity(Nᵃ::AbstractVector{T},Uₐ::AbstractVector{T}) where T

    n::T = 0.0

    if Nᵃ[1] == 0
        n = T(0.0)
    else
        n = - dot(Nᵃ,Uₐ)
    end

    return n

end

function ScalarMassDensity(n::AbstractVector{T},PhaseSpace::PhaseSpaceStruct,species::S) where S<:Union{String,Int64} where T

    if typeof(species) == Int64
        species_index = species
    else
        name_list = PhaseSpace.name_list
        species_index = findfirst(==(species),PhaseSpace.name_list)
    end
    m = PhaseSpace.Grids.mass_list[species_index]

    ρ::T = n*m

    return ρ

end

function ScalarEnergyDensity(Tᵃᵇ::AbstractMatrix{T},Uₐ::AbstractVector{T},n::T;perparticle=false) where T

    e::T = 0.0
    en::T = Uₐ' * Tᵃᵇ * Uₐ

    if perparticle
        if n == 0.0
            e = 0.0
        else
            e = en/n
        end

        return e
    else 
        return en
    end

end

function ScalarPressure(Tᵃᵇ::AbstractMatrix{T},Δab::AbstractMatrix{T}) where T

    p::T = 0.0

    metric = zeros(T,4,4)
    metric[1,1] = T(-1)
    metric[2,2] = T(1)
    metric[3,3] = T(1)
    metric[4,4] = T(1)

    #p = (1/3) * sum(Tab .* (metric * Δab * metric))
    p = (1/3) * sum(Tᵃᵇ .* Δab)

    return p

end

function ScalarTemperature(p::T,n::T) where T

    Temp::T = 0.0

    kb = CONST_kb
    c = CONST_c
    mEle = CONST_mEle
    
    Temp = p/n / kb * mEle * c^2 # Kelvin

    return Temp

end


