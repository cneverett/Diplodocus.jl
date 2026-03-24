# Installation

To use *Diplodocus* you will first need to install *Julia* on your machine (see guide for installing Julia [here](https://julialang.org/install/))

Then inside a Julia REPL you can install *Diplodocus.jl* using the julia package manager 
```julia-repl
juila> ] # ']' should be pressed
pkg> add Diplodocus
```

That's it! For guides on setting up simulations see the [Tutorials](@ref).

::: warning

The main branch of Diplodocus.jl is currently locked to v0.1 while v0.2 is being tested. To use the dev branch with v0.2 features use instead:  

    ```julia-repl
    juila> ] # ']' should be pressed
    pkg> add Diplodocus#dev
    ```

But this version is not guaranteed to be stable.

:::