# Tutorial 1: Evolution of a Population of Hard Spheres

---

In this tutorial we will consider a population of hard spheres which are homogenous in space. These spheres only interact via perfectly elastic collisions (a binary interaction), the cross section for which is given in [Implemented Collisions](@ref). As a result of which an initially non-thermal and anisotropic population is expected to tend towards thermal and isotropic over time, as can be seen in the resulting animation below:

```@raw html
<video autoplay loop muted playsinline controls src="./assets/HardSphereMomentumComboAnimation.mp4" style="max-height: 60vh;"/>
```

This tutorial is split into three parts. First we will generate the collision matrix needed to examine the effects of the elastic collisions between hard spheres. Second we will setup and run the phase space evolution of theses spheres. And third we will plot the results.

Raw julia scripts for each of these parts can be found in the Diplodocus.jl GitHub repository under the location `src/Tutorials/Tutorials1/` and file names `Tutorial1_collisions.jl`, `Tutorial1_transport.jl`, and `Tutorial1_plots.jl`. Each can be run directly inside a Julia REPL using, e.g.

```julia

juila> include("Tutorial1_collisions.jl")

```