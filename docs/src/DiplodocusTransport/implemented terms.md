# Implemented Coordinates and Forces

## Implemented Spatial Coordinates
These structs define what spatial coordinates are used when solving transport equations.

The currently implemented spatial coordinates include:

| Coordinates    | Coordinate Struct | Comment |
| ----------------- | ------------ | ------------------------------- |
| Cartesian Minkowski   | `Cartesian()`     |  Standard ``x,y,z`` flat Minkowski spacetime | 
| Cylindrical Minkowski      | `Cylindrical()`     |     Standard ``\rho,\vartheta,z`` flat Minkowski spacetime  |
| Spherical Minkowski | `Spherical()`   | Standard ``r,\theta,\psi`` flat Minkowski spacetime   | 

### Global to Local Rotation

By default the local ortho-normal basis of momentum is aligned to the local coordinate directions. For `Cartesian` and `Cylindrical` spatial coordinates, the momentum-space polar axis is aligned to the ``z`` spatial coordinate direction, with the ``\phi=0`` direction aligned to the spatial ``x`` (``rho``) direction. For `Spherical` spatial coordinates, the polar axis is instead aligned with the ``x`` (``r``) direction with ``\phi=0`` aligned to the ``y`` (``\theta``) direction.

These default directions may be edited by defining `global_to_local_rotation` a `NTuple{3,Float64}` containing three angles (radians normalised by ``\pi``) ``(\alpha,\beta,\gamma)``, which give an ``z-x-z`` rotation of the local bases vectors. This can be used to specify the local magnetic field direction with respect to the spatial coordinate bases vectors.

## Implemented Forces 
These structs define which forces to include when solving the transport equations and depend on the spatial coordinates of the system.

The currently implemented forces include:

| External Force    | Force Struct | Implemented Spatial Coordinates | Implemented Momentum Space Anisotropies |
| ----------------- | ------------ | ------------------------------- | --------------------------------------- |
| Ricci Rotation coefficients   | `CoordinateForce()`     |   `Cartesian()`, `Cylindrical()`, `Spherical()`   |  N/A | 
| Lorentz Force (Uniform and Orthogonal ``E`` and ``B`` Fields)      | `ExB()`     |     `Cartesian()`   |  N/A | 
| Synchrotron Radiation Reaction | `SyncRadReact()`   |             `Cartesian()`, `Cylindrical()`, `Spherical()`   | `Ani()`,`Axi()`,`Iso()` | 
| Grad B. Force due to 1/z B Field | `InvGradB()`   |            `Cylindrical()`   | N/A | 
