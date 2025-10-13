```@meta
CurrentModule = DiplodocusCollisions
```
# Implemented Collisions

## Binary Collisions

These binary interactions have currently been implemented:

| Binary Interaction    | Abbreviation | Internal Functions (cross sections) |
| --------------------- | ------------ | ----------------------------------- | 
| Elastic collision of spheres                                    | `SphSphSphSph`     |   [`dsigmadt_SphSphSphSph`](@ref) [`sigma_SphSphSphSph`](@ref)             |  
| Photon pair production from electron-positron annihilation      | `ElePosPhoPho`     |    [`dsigmadt_ElePosPhoPho`](@ref) [`sigma_ElePosPhoPho`](@ref)            | 
| Electron-Positron pair production from photon pair annihilation | `PhoPhoElePos`     |             [`dsigmadt_PhoPhoElePos`](@ref) [`sigma_PhoPhoElePos`](@ref)   |
| Electron(or Positron)-Photon scattering (Compton scattering)    | `ElePhoElePho`     |          [`dsigmadt_ElePhoElePho`](@ref) [`sigma_ElePhoElePho`](@ref)      |

## Emission Interactions

These emissive interactions have currently been implemented:

| Emissive Interaction  | Abbreviation | Internal Functions (emissivity kernel) |
| --------------------- | ------------ | ----------------------------------- | 
| Synchrotron(cyclostron) emission of photons by a charged particle (`Name1`)  | `SyncName1Name1Pho`     |   [`SyncKernel`](@ref)           |  

# Internal Collision Functions

```@docs
sigma_SphSphSphSph
dsigmadt_SphSphSphSph

sigma_ElePosPhoPho
dsigmadt_ElePosPhoPho

sigma_PhoPhoElePos
dsigmadt_PhoPhoElePos

sigma_ElePhoElePho
dsigmadt_ElePhoElePho
```

```@docs

SyncKernel

```