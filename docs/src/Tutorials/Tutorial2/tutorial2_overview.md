# Tutorial 2: Non-Linear Electron Cooling

---

In this tutorial we will consider an *inverse-Compton catastrophe*! This refers to the rapid rise in the luminosity of inverse-Compton scattered photons cased by a rapidly cooling electron population [Longair_2011,Had_Code_Comp_2026](@citep).

If electrons emit synchrotron photons via radiation reaction with a magnetic field, this catastrophe can occur if the energy density of the magnetic field is lower than the energy density of the synchrotron photons. Therefore a runaway process can occur: low-energy synchrotron photons are up-scattered to medium energies by the electron population population. If the energy density of this medium-energy population is again greater than the low-energy photons, they can once again be up-scattered to high-energies and so on.... This process will continue until the highest-energy Compton scatterings begin to take place in the Klein-Nishina regime where scattering begins to be suppressed. 

As radiation reaction, synchrotron, and Compton scattering are all included in Diplodocus (see [Implemented Collisions](@ref) and [Implemented Coordinates and Forces](@ref)), we can examine this exact scenario. To do this we will use the setup of [Had_Code_Comp_2026](@citep) Sec.6.1. We will consider a single, spherical region with radius ``R=10^16\,\mathrm{m}`` with a tangled magnetic field of strength ``B=10\,\rm{G}``. In this region we will inject a power-law distribution of electrons with index ``\alpha=2`` spanning from a minimum Lorentz factor of ``\gamma_{e,\text{min}}=10^{1.9}`` to a maximum ``\gamma_{e,\text{min}}=10^{2.1}``. This injected distribution will have a compactness ``l_e^\text{inj}=1`` and luminosity ``L_e^\text{inj}=4.6\times10^{45}\,\mathrm{erg}\mathrm{s}^{-1}``.

This tutorial is split into four parts. First we will generate the *binary* collision matrix needed for the Compton scattering. Second we will generate the *emission* matrix needed for the synchrotron photons. Third we will set up the problem and run the solver, and finally we will plot the results.

Raw julia scripts for each of these parts can be found in the Diplodocus.jl GitHub repository under the location `src/Tutorials/Tutorials2/` and file names `Tutorial2_binary_collisions.jl`, `Tutorial2_emission_collisions.jl`, `Tutorial2_transport.jl`, and `Tutorial2_plots.jl`. Each can be run directly inside a Julia REPL using, e.g.

```julia

juila> include("Tutorial2_binary_collisions.jl")

```