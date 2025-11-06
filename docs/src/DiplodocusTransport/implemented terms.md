# Implemented Spatial Coordinates and Forces

## Implemented Spatial Coordinates
These structs define what spatial coordinates are used when solving transport equations.

The currently implemented spatial coordinates include:

| Coordinates    | Coordinate Struct | Comment |
| ----------------- | ------------ | ------------------------------- |
| Cartesian Minkowski   | `Cartesian()`     |  Standard ``x,y,z`` flat Minkowski spacetime | 
| Cylindrical Minkowski      | `Cylindrical()`     |     Standard ``\rho,\vartheta,z`` flat Minkowski spacetime. Allows three arguments `Cylindrical(a,b,c)` for local "z-x-z" rotation of the momentum-space ``z`` axis from the spatial cylindrical ``z`` axis  |
| Spherical Minkowski | `Spherical()`   | Standard ``r,\theta,\psi`` flat Minkowski spacetime   | 


## Implemented Forces 
These structs define which forces to include when solving the transport equations and depend on the spatial coordinates of the system.

The currently implemented forces include:

| External Force    | Force Struct | Implemented Spatial Coordinates | Implemented Momentum Space Anisotropies |
| ----------------- | ------------ | ------------------------------- | --------------------------------------- |
| Ricci Rotation coefficients   | `CoordinateForce()`     |   `Cartesian()`,`Cylindrical()`,`Spherical()`   |  N/A | 
| Lorentz Force (Uniform and Orthogonal ``E`` and ``B`` Fields)      | `ExB()`     |     `Cartesian()`   |  N/A | 
| Synchrotron Radiation Reaction | `SyncRadReact()`   |             `Cartesian()`,`Cylindrical()`,`Spherical()`   | `Ani()`,`Axi()`,`Iso()` | 
