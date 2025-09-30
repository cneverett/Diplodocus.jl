# DIPLODOCUS Framework

::: tip **DIPLODOCUS**

**D**istribution **I**n **PL**ateaux meth**ODO**logy for the comp**U**tation of transport equation**S**.

:::

Diplodocus is a novel framework for evolving the a set of particle distribution functions ``f(\boldsymbol{x},\boldsymbol{p})`` through the seven dimensions of phase space: one time, three space and three momentum. Evolution includes advection of particles, continuous forcing and discrete interactions between particles.   

## Transport Equations
Transport equations refers to a set of equations that dictate the evolution of particles through phase space. In Diplodocus, that takes the form of an integrated Boltzmann equation: 
```math
\int_{\partial Q} f(\boldsymbol{x},\boldsymbol{p}) \boldsymbol{\omega} = \int_{Q} \boldsymbol{C}(\boldsymbol{x},\boldsymbol{p}).
```
Terms on the left-hand-side dictate the continuous transport of particle through phase space including advection through space and forcing through momentum. The right-hand-side, described termination and beginning of particle worldlines in a volume of phase space due to discrete interaction between particles. 

## Distribution-In-Plateaux
To solve the evolution described by the transport equations computationally, a method called "Distribution-In-Plateaux" is used to discretise particle distribution functions. In brief, a continuous distribution over a surface in phase space is divided into a number of plateaux over sub-areas of the surface. The "height" of this plateaux is then taken to be the average value of the continuous distribution over that sub-area. 

*INSERT IMAGE OF CONTINIOUS SURFACE TO GRID*

This discretisation allows the effects of discrete interactions, collisions, to be pre-computed while simultaneously providing a flexible and conservative numerical scheme for the transport of distribution function through phase space. 

Pre-computation of collision terms is implemented in the `DiplodocusCollisions.jl` package. 

Transport of the distribution functions is handled by the `DiplodocusTransport.jl` package.

Plotting of results is split into a separate package `DiplodocusPlots.jl`.

