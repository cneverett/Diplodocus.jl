# DIPLODOCUS: Distribution In Plateaux Methodology for the Computation of transport EquationS

[![][docs-latest-img]][docs-latest-url]
[![Build Status][gha-img]][gha-url]

[docs-latest-img]: https://img.shields.io/badge/Docs-Stable-lightgrey.svg
[docs-latest-url]: https://cneverett.github.io/Diplodocus.jl/dev/

[gha-img]: https://github.com/cneverett/Diplodocus.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gha-url]: https://github.com/cneverett/Diplodocus.jl/actions/workflows/CI.yml?query=branch%3Amain

*Diplodocus.jl* is a Julia package for the transport of particle distribution functions through phase space with the inclusion of continous forces and discrete interactions between particle. It contains three sub-packages: [*DiplodocusCollisions.jl*](https://github.com/cneverett/DiplodocusCollisions.jl) generates collision arrays via Monte-Carlo integration that describe discrete particle interactions; these pre-computable arrays are then fed into [*DiplodocusTransport.jl*](https://github.com/cneverett/DiplodocusTransport.jl), where advection terms are added and the evolution of particle distribution function can be simuilated; results of these simulations can then be plotted using function from [*DiplodocusPlots.jl*](https://github.com/cneverett/DiplodocusPlots.jl).

Documentation, installation instructions and tutorials can be found [here](https://cneverett.github.io/Diplodocus.jl/). 
