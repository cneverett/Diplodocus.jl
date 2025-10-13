# Adding New Collisions

## Adding New Binary Collisions

Users may add their own binary interaction cross sections to the `/src/commom/DifferentialCrossSectionFunctions.jl` file of `DiplodocusCollisions.jl`. Functions should be named in the following format: `"sigma_name1name2name3name4"` and `"dsigmadt_name1name2name3name4"` where the names are three letter abbreviations of the particles involved (see [Particles, Grids and Units](@ref)). The named pairs `name1name2` and `name3name4` should be in alphabetical order. Both the total cross section and differential cross sections must be provided for a single interaction.    

All cross sections are to be defined in terms of the Mandelstram variables $s=(p_1^\mu+p_2^\mu)^2$, $t=(p_1^\mu-p_3^\mu)^2$ and $u=(p_2^\mu-p_3^\mu)^2$. To maintain accuracy of cross sections and avoid DivZero issues when momenta is small compared to the mass of the particles (at `Float64` precision), each Mandelstram variable in the cross sections should be split into two components:
- `s=sSmol+sBig` where $sBig = (m_1+m_2)^2$
- `t=tSmol+tBig` where $tBig = (m_1-m_3)^2$
- `u=uSmol+uBig` where $uBig = (m_2-m_3)^2$
The "Big" part typically cancels with terms in the cross sections, leading to better accuracy. Therefore the function should defined as follows:
```julia
function sigma_name1name2name3name4(sSmol::Float64,sBig::Float64,tSmol::Float64,tBig::Float64,uSmol::Float64,uBig::Float64)
    ... 
end

function dsigmadt_name1name2name3name4(sSmol::Float64,sBig::Float64,tSmol::Float64,tBig::Float64,uSmol::Float64,uBig::Float64)
    ... 
end
```
where all of `sSmol`, `sBig`, `tSmol`, `tBig`, `uSmol`, and `uBig` must be included in the function definition irrespective of if they actually appear in the cross section. It is also important to define a normalisation for the cross sections to avoid emission and absorption terms becoming small compared to the `Float64` minimum.

::: warning

Total cross sections may typically be defined/derived in textbooks and other sources to include division by $1/2$ if output states are identical. This factor should NOT be included here as this factor is included separably in the code.

:::

## Adding New Emissive Interactions

Coming soon