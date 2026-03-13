```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
    name: "Diplodocus.jl"
    text: ""
    tagline: Particle Transport and Emissions for Astrophysics
    image:
        src: /test1.gif
        alt: Diplodocus
        width: 480
        height: 480
    actions:
        -   theme: brand
            text: Get Started
            link: /Get Started/overview
        -   theme: alt
            text: Tutorials
            link: /Tutorials/overview
        -   theme: alt
            text: View on Github
            link: https://github.com/cneverett/Diplodocus

---

```

# What is Diplodocus.jl?

Diplodocus.jl is a Julia package for modelling of high energy astrophysical sources such as jets from Active Galactic Nuclei, Gamma-Ray Bursts, and X-Ray Binaries. 

It acts as a bridge between macro-scale simulations of spatial structure and dynamics (e.g. GRMHD) and micro-scale modelling of particle interactions an emissive processes (e.g. one-/multi-zone approaches); capable of including both regimes! Diplodocus.jl solves the phase-space transport of particle distribution functions, including: arbitrary forces and asymmetric interactions on a curved spacetime background with options to add complex electromagnetic field structures. 

By evolving particles in space, momentum, and time, Diplodocus.jl captures of all the these multi-scale effects simultaneously and self-consistently, allowing the generation of multi-wavelength emission spectra and imaging through which the true physics of a wide range of astrophysical objects can be studied.      

## Citing Diplodocus.jl

If you use Diplodocus.jl, any of its components or are just feeling kind, please cite one or all of the main Diplodocus series papers (Paper III coming soon):

>  C. N. Everett and G. Cotter, DIPLODOCUS I: Framework for the evaluation of relativistic transport equations with continuous forcing and discrete particle interactions (2026), The Open Journal of Astrophysics, 9, https://doi.org/10.33232/001c.155822

>  C. N. Everett, M. Klinger-Plaisier and G. Cotter, DIPLODOCUS II: Implementation of Transport Equations and Test Cases Relevant to Micro-Scale Physics of Jetted Astrophysical Sources (2026), The Open Journal of Astrophysics, 9, https://doi.org/10.33232/001c.157823

::: details Show BibTeX

```
@ARTICLE{DiplodocusI,
       author = {{Everett}, Christopher N. and {Cotter}, Garret},
        title = "{DIPLODOCUS I: Framework for the evaluation of relativistic transport equations with continuous forcing and discrete particle interactions}",
      journal = {The Open Journal of Astrophysics},
     keywords = {High Energy Astrophysical Phenomena},
         year = 2026,
        month = jan,
       volume = {9},
        pages = {55822},
          doi = {10.33232/001c.155822},
archivePrefix = {arXiv},
       eprint = {2508.13296},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026OJAp....955822E},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{DiplodocusII,
       author = {{Everett}, Christopher N. and {Klinger-Plaisier}, Marc and {Cotter}, Garret},
        title = "{DIPLODOCUS II: Implementation of transport equations and test cases relevant to micro-scale physics of jetted astrophysical sources}",
      journal = {The Open Journal of Astrophysics},
     keywords = {High Energy Astrophysical Phenomena},
         year = 2026,
        month = feb,
       volume = {9},
        pages = {57823},
          doi = {10.33232/001c.157823},
archivePrefix = {arXiv},
       eprint = {2510.12505},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026OJAp....957823E},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

:::
