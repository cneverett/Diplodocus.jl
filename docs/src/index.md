```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
    name: "Diplodocus.jl"
    text: "Particle transport for Astrophysics"
    tagline: A Julia framework for evaluating particle transport through phase space for emission modelling of astrophysical sources
    image:
        src: assets/diplodocus_sphere.png
        alt: Diplodocus
    actions:
        -   theme: alt
            text: View on Github
            link: https://github.com/cneverett/Diplodocus

features:
    -   title: Overview
        details: Overview of the Diplodocus framework, installation and examples
        link: /Overview/overview
    -   title: Collisions
        details: Evaluation of collision integrals for discrete interactions between particles using the DiplodocusCollisions.jl package
        link: /DiplodocusCollisions/overview
    -   title: Transport
        details: Evaluation of particle transport through phasespace, including forcing and discrete interactions using the DiplodocusTransport.jl package
        link: /DiplodocusTransport/overview
    -   title: Plotting
        details: Plotting the evolution of particle distributions due to transport using the DiplodocusPlots.jl package
        link: /DiplodocusPlots/overview

---
```