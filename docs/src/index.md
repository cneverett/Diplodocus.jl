```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
    name: "Diplodocus.jl"
    text: "Particle transport for Astrophysics"
    tagline: A Julia framework for evaluating transport of particle distributions through phase space for mesoscopic modelling of astrophysical sources
    image:
        src: /test1.gif
        alt: Diplodocus
        width: 480
        height: 480
    actions:
        -   theme: alt
            text: View on Github
            link: https://github.com/cneverett/Diplodocus

features:
    -   title: Get Started
        details: Overview of the Diplodocus framework, installation guide and the codes conventions
        link: /Get Started/overview
    -   title: Tutorials
        details: Worked examples of how to use Diplodocus
        link: /Tutorials/tutorials
    -   title: Collisions
        details: Evaluation of collision integrals for discrete interactions between particles using the DiplodocusCollisions.jl package
        link: /DiplodocusCollisions/overview
    -   title: Transport
        details: Evaluation of particle transport through phase space, including forcing and discrete interactions using the DiplodocusTransport.jl package
        link: /DiplodocusTransport/overview
    -   title: Plotting
        details: Plotting the evolution of particle distributions due to transport using the DiplodocusPlots.jl package
        link: /DiplodocusPlots/overview

---
```
