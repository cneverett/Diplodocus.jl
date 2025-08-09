# DIPLODOCUS: Distribution In Plateaux Methodology for the Computation of transport EquationS

[![Build Status](https://github.com/cneverett/DiplodocusDocs.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/cneverett/DiplodocusDocs.jl/actions/workflows/CI.yml?query=branch%3Amaster)

> [!CAUTION]
> This project is currently under development

*Diplodocus.jl* is a Julia package for the transport of particle distribution functions through phase space with the inclusion of continous forces and discrete interactions between particle. It contains three sub-packages: [*DiplodocusCollisions.jl*](https://github.com/cneverett/DiplodocusCollisions.jl) generates collision arrays via Monte-Carlo integration that describe discrete particle interactions; these pre-computable arrays are then fed into [*DiplodocusTransport.jl*](https://github.com/cneverett/DiplodocusTransport.jl), where advection terms are added and the evolution of particle distribution function can be simuilated; results of these simulations can then be plotted using function from [*DiplodocusPlots.jl*](https://github.com/cneverett/DiplodocusPlots.jl).

Documentation, installation instructions and tutorials can be found [here](https://cneverett.github.io/Diplodocus.jl/). 

## Basic Installation 
  Coming soon!
