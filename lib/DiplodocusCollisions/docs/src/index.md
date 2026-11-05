# BoltzmannCollisionIntegral.jl

`BoltzmannCollisionIntegral.jl` is a package for the evaluation of the Boltzmann collision integral for binary interactions of arbitrary momentum and anisotropic particle distributions using a method of momentum discretisation and Monte-Carlo integration. 

The code now also includes a module for the evaluation of synchrotron emissions, with details found in [Synchrotron Radiation](@ref).

The evaluation of the collision integral is of great use when studying kinetic system. The evolution of particle distributions ``f(x^\mu,\boldsymbol{p})`` within such a system is dictated by the collisional relativistic Boltzmann equation ([Everett+Cotter 2024](https://doi.org/10.1093/rasti/rzae036)):
```math
p^\nu\partial_\nu f(x^\mu,\boldsymbol{p})+\partial_{p^\nu}\left(F^\nu f(x^\mu,\boldsymbol{p})\right)=C(x^\mu,\boldsymbol{p}),
```
where ``p^\mu=(p^0,\boldsymbol{p})`` is the 4-Momentum and ``C(x^\mu,\boldsymbol{p})`` is the collision integral at a particular point in spacetime ``x^\mu``. From here we shall neglect the terms involving advection in space and advection on momentum-space (external forcing), as such the dependence on ``x^\mu`` is also dropped. What is left can be written as 
```math
    \partial_t f(\boldsymbol{p}) = \frac{1}{p^0}C(\boldsymbol{p})
```
The collision integral for a particular particle ``3`` in the interaction ``12\rightleftharpoons34`` is given by 
```math
    C(\boldsymbol{p}_3)=\int\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\delta_{12}}W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)- \frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\delta_{12}}W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu)\right],
```
with the transition rate ``W`` being given by:
```math
W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)=\frac{sp_{\text{in}}^{*2}}{\pi}\frac{\mathrm{d}\sigma_{12|34}}{\mathrm{d}t}(s,t)\delta^{(4)}(p_1^\mu+p_2^\mu-p_3^\mu-p_4^\mu),
```
and ``\frac{\mathrm{d}\sigma_{12|34}}{\mathrm{d}t}(s,t)`` is the Lorentz invariant differential cross section of the forward reaction ``12\to34``.

The Boltzmann equation can then be discretised in momentum space by integrating over a discrete volume element ``\Delta^3\boldsymbol{p}_3`` to yield: 
```math 
\int_{\Delta^3\boldsymbol{p}_3}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}p_3^\mu\partial_\mu f(x^\mu,\boldsymbol{p}_3)=\int_{\Delta^3\boldsymbol{p}_3}\mathrm{d}^3\boldsymbol{p}_3~S_{spe}(x^\mu,\boldsymbol{p}_3)-T_{spe}(x^\mu,\boldsymbol{p}_3),
```
where the collision integral has been split into two parts, an emission spectrum
```math
S_{\text{spe}}(\boldsymbol{p}_3)=\frac{1}{1+\delta_{12}}\frac{1}{p^0_3}\int\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu),
```math
dictating the rate of gain of particles of a specific type at a specific momentum from a binary interaction, and an absorption spectrum 
```math
T_{\text{spe}}(\boldsymbol{p}_3)=\frac{1}{1+\delta_{12}}\frac{1}{p^0_3}\int\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0} f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu),
```
dictating the rate of loss of particles of a specific type at a specific momentum from the same interaction.

## Assumptions
In order to enable the evaluation of the emission and absorption terms, we make the following set of assumptions (dropping dependence on space-time for notational simplicity):  
- The system is axisymmetric in momentum space coordinates i.e. only dependant on ``p=|\boldsymbol{p}|`` and the cosine of the angle to the axis of symmetry (``\hat{z}``) ``u`` where ``\hat{z}\cdot\boldsymbol{p}=p\cos\theta=pu`` 
- The distribution functions are redefined in an axisymmetric form i.e. ``f(\boldsymbol{p})\equiv\frac{f(p,u)}{2\pi p^2}``. This is done such that a single particle is described by the distribution ``f(p,u)=\delta(p-p_i)\delta(u-u_i)``.
- The axisymmetric distribution functions are averaged over phase space intervals ``\Delta p\Delta u``, to generate discrete values i.e. ``f(t,p_i, u_j)=\frac{1}{\Delta p_i\Delta u_j}\int_{\Delta p_i\Delta u_j} \mathrm{d}p \mathrm{d} u~f(p, u)``, such that they act as if they are constant over that interval of phase space.  
- Particle masses, momenta and energies are normalised by a factor of the rest mass of the electron ``m_\text{Ele}`` and the speed of light ``c=1``.

## Discrete Form of the Boltzmann Equation
Following the assumptions laid out in the previous section. The axisymmetric, momentum-discrete Boltzmann equation is given by (remember the neglect of space and momentum-space advection terms)
```math
\partial_t f_{3,ij}=f_{1,kl}f_{2,mn}S_{12|34,ijklmn}-f_{3,ij}f_{4,kl}T_{34|12,ijkl}
```
The discrete emission spectrum term is a 6D array, given by
```math
f_{1,kl}f_{2,mn}S_{12|34,ijklmn}=f(p_{1,k}, u_{1,l})f(p_{2,m}, u_{2,n})\frac{\Delta p_{1,k}\Delta  u_{1,l}\Delta p_{2,m}\Delta  u_{2,n}}{\Delta p_{3,i}}\frac{1}{N}\sum^{N}_{a=1}\left[S_{val}\right](\{\boldsymbol{p}_1,\boldsymbol{p}_2,\boldsymbol{p}_3\}_a), 
```
with 
```math
S_{val}=\frac{1}{1+\delta_{12}}\sum_\pm\frac{2p_\pm^2\mathcal{F}_{12}^2}{p_1^0p_2^0}\frac{\mathrm{d}\sigma_{12|34}}{\mathrm{d}t}(s,t_\pm)\frac{1}{p_1^0p_\pm-p_\pm^0p_1\cos\Theta_{\pm1}+p_2^0p_\pm-p_\pm^0p_2\cos\Theta_{\pm2}},
```
where ``p_{\pm}`` are the two roots of ``(s+t-m_3^2-m_2^2-2p^ u_1p_{4 u})``. 

The discrete absorption term is a 4D Array, given by
```math
f_{3,ij}f_{4,kl}T_{34|12,ijkl}=f(p_{3,i}, u_{3,j})f(p_{4,k}, u_{4,l})\Delta\boldsymbol{p}_{4,k}\Delta u_{4,l}\frac{1}{N}\sum^{N}_{a=1}\left[T_{val}\right](\{\boldsymbol{p}_3,\boldsymbol{p}_4\}_a),
```
with
```math
T_{val}=\frac{1}{1+\delta_{12}}\frac{\mathcal{F}_{34}(s)\sigma_{34|12}(s)}{p_3^0p_4^0}.
```
and ``\sigma_{34|12}(s)`` the lorentz invariant total cross section of the reverse reaction ``34\to12``.

## Evaluation
Given a user input of interaction to evaluate and what discretisation to use, the discrete emission spectra ``S_{12|34,ijklmn}``, ``S_{12|43,ijklmn}`` and absorption spectra ``T_{12|34,ijkl}``, ``T_{21|34,ijkl}`` are solved simultaneously by Monte-Carlo integration. This provides all matricies required to self consistently model the evolution of particles due to the forward interaction ``12\to34``. 

!!! note
    The order of the interaction in the absorption spectrum (``34|12 \to 12|34``) is changed for easier evaluation as sampled points can be used for both emission and absorption terms. Rather than solving for the forward and reverse reaction, only the forward reaction is solved for. Hence two sets of emission and absorption spectra are need to be generated in general to describe a single reversible reaction.  

The evaluation is performed by the functions [SpectraEvaluateSerial](@ref) or [SpectraEvaluateMultiThread](@ref) depending on if the code is to run in serially on a single core or multi-threaded on multiple cores of a single CPU. The code will then save the resulting multidimensional arrays.

!!! warning 
    The evaluation does not guarantee convergence of the Monte-Carlo integration, instead it will sample only a user defined number of points. A rough measure of convergence is calculated (based on comparison to the output of the previous run) and stored as `SConverge` and `TConverge` in the output file.

Evaluation is performed using `Float64` precision, with all [Internal Functions](@ref) designed to maintain accuracy to this precision.

An in-depth setup guide is provided in the [Getting Started](@ref) section.  