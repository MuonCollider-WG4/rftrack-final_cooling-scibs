# Muon Collider Final Cooling — Space Charge and Intrabeam Scattering Studies

**Authors:** Paula Desire, Andrea Latina

## Overview

This repository contains simulations of the **final cooling section of the Muon Collider**, with a focus on the effects of **Space Charge (SC)** and **Intrabeam Scattering (IBS)** on the muon beam.

The current implementation considers the **first three cells** of the final cooling lattice.

The lattice design is based on the Geant4Beamline model developed by **R. Taylor et al.**, available in the [Muon Collider Final Cooling repository](https://github.com/MuonCollider-WG4/muon_final_cooling/tree/main).

## Objective

The objective of this study is to simulate and compare muon tracking through the first three cells of the final cooling lattice under three different conditions:

1. Without Space Charge or Intrabeam Scattering.
2. Including Space Charge only.
3. Including both Space Charge and Intrabeam Scattering.

This allows the impact of these collective effects on the beam dynamics during final cooling to be investigated.

## Simulation Tools

The simulations are performed using **RF-Track**.

RF-Track can be used through either its:

* Python interface
* Octave interface

Both implementations are provided where applicable.

## Simulation Cases

### `threecells`

Tracks the muon bunch through the three cooling cells without collective effects.

* Space Charge: **No**
* Intrabeam Scattering: **No**

### `threecells_SC`

Tracks the muon bunch through the three cooling cells including Space Charge.

* Space Charge: **Yes**
* Intrabeam Scattering: **No**

### `threecells_SC_IBS`

Tracks the muon bunch through the three cooling cells including both Space Charge and Intrabeam Scattering.

* Space Charge: **Yes**
* Intrabeam Scattering: **Yes**

## Output

Simulation results are stored in the following directories:

```text
Phase_space_Bunch6d/
Phase_space_Bunch6d_SC/
Phase_space_Bunch6d_SC_IBS/
```

These correspond respectively to:

| Simulation     | Output directory              |
| -------------- | ----------------------------- |
| No SC / No IBS | `Phase_space_Bunch6d/`        |
| SC             | `Phase_space_Bunch6d_SC/`     |
| SC + IBS       | `Phase_space_Bunch6d_SC_IBS/` |

## Analysis and Plotting

The `plotter_cuts` scripts can be used to analyze the resulting phase-space distributions and generate plots with **sigma cuts** applied to the beam distribution.

Implementations are available for both Python and Octave:

```text
plotter_cuts.py    # Python
plotter_cuts.m     # Octave
```

## RF-Track Documentation

For additional information about RF-Track, its physics models, available elements, and interfaces, refer to the [RF-Track Reference Manual](https://gitlab.cern.ch/rf-track/rf-track-reference-manual/-/raw/master/RF_Track_reference_manual.pdf?ref_type=heads).

## References

* R. Taylor et al., *Muon Collider Final Cooling lattice*, [MuonCollider-WG4/muon_final_cooling](https://github.com/MuonCollider-WG4/muon_final_cooling/tree/main).
* A. Latina, *RF-Track Reference Manual*, CERN. [RF-Track Reference Manual](https://gitlab.cern.ch/rf-track/rf-track-reference-manual/-/raw/master/RF_Track_reference_manual.pdf?ref_type=heads)
