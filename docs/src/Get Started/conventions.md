# Particles, Grids and Units

## Particles
Below is a table of the currently implemented particles (i.e. their particle properties are defined within the code)

| Particle | Abr. String | Mass (units of ``m_\text{Ele}``)     |  Charge (units of ``q`` the fundamental charge) |
| -------- | ----------- | --------------------- | ---- |
| Sphere   | `"Sph"`     |  1836.1528 (mass of Proton)     |  0.0 |
| Electron | `"Ele"`     |                    1.0          | -1.0 |
| Positron | `"Pos"`     |                    1.0          | 1.0  |
| Proton   | `"Pro"`     |                    1836.1528    | 1.0  |

## Grids
Three types of grids for coordinates in phase space have been implemented

| Grid Spacing | Abr. String | Notes                                        | 
| -------- | ----------- | -------------------------------------------- |
| Uniform   | `"u"`     |  Uniform grid spacing between the upper and lower bounds    | 
| Log10 | `"l"`     |  Exponentially increasing grid spacing, uniform in Log10 space from upper to lower bounds (bounds should be given as Log10(bound value))                                            | 
| Binary | `"b"`     |  Fractionally decreasing grid spacing for angular grids. Bin widths half as they move away from u=0 in both directions to u=+-1.                                            | 

## Units
Within Diplodocus, the speed of light ``c`` is taken to be unity and all masses, momenta and energies are normalised by the rest mass of the electron ``m_\text{Ele}``. Further all interaction rates are normalised by ``\sigma_\text{T}c``, where ``\sigma_\text{T}`` is the Thompson scattering cross section. Therefore time ``t`` is given by ``t [s] = t [\text{code units}] / \sigma_\text{T}c \approx 5\times10^{19}[s]*t [\text{code units}]``, with ``t [code units]`` corresponding to the characteristic time for an interaction with cross section ``\sigma_\text{T}`` and particle density ``1 [m^{-3}]``. 

