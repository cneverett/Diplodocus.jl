<img align="left" width="64" height="64" src="./docs/src/assets/favicon.svg" alt="Diplodocus icon">

# Diplodocus.jl
**D**istribution **I**n **PL**ateaux meth**ODO**logy for the **C**omp**U**tation of transport equation**S**

[![][docs-latest-img]][docs-latest-url]
[![][docs-dev-img]][docs-dev-url]
[![Build Status][gha-img]][gha-url]
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![ascl](https://img.shields.io/badge/ascl-2605.023-blue.svg?colorB=262255)](https://ascl.net/2605.023)

[docs-latest-img]: https://img.shields.io/badge/Docs-Stable-lightgrey.svg
[docs-latest-url]: https://cneverett.github.io/Diplodocus.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/Docs-Dev-lightgrey.svg
[docs-dev-url]: https://cneverett.github.io/Diplodocus.jl/dev/

[gha-img]: https://github.com/cneverett/Diplodocus.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gha-url]: https://github.com/cneverett/Diplodocus.jl/actions/workflows/CI.yml?query=branch%3Amain

*Diplodocus.jl* is a Julia package for the transport of particle distribution functions through phase space with the inclusion of continuous forces and discrete interactions between particle. It contains three sub-packages: [*DiplodocusCollisions.jl*](https://github.com/cneverett/DiplodocusCollisions.jl) generates collision arrays via Monte-Carlo integration that describe discrete particle interactions; these pre-computable arrays are then fed into [*DiplodocusTransport.jl*](https://github.com/cneverett/DiplodocusTransport.jl), where advection terms are added and the evolution of particle distribution function can be simulated; results of these simulations can then be plotted using function from [*DiplodocusPlots.jl*](https://github.com/cneverett/DiplodocusPlots.jl).

Documentation, installation instructions and tutorials can be found [here](https://cneverett.github.io/Diplodocus.jl/). 

This framework is being published as the DIPLODOCUS series of papers, the first two covering the theoretical background, numerical implementation, and demonstration its capability via test cases relevant to micro-physical problems, can be found at [Paper I](https://doi.org/10.48550/arXiv.2508.13296) and [Paper II](https://doi.org/10.48550/arXiv.2510.12505).
