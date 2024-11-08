# A Structure for Describing Macroscopic Structure

## Context

The main job of spinW is to calculate the energy levels of a magnetic system over a range of q values. This
task is relatively self-contained and well defined by just the material's cystallographic and magnetic properties.

However, there are a number of things that are more properly understood as a property of a particular experiment. 
The main features are to include here is whether or not the crystal is twinned or is a powder.

There is a question about whether it should be the location for the option to apply neutron scattering selection rules.


## Decision

We will create a base `Sample` class with derived `TwinnedSample` and `PowderSample` classes which will have methods for preprocessing the input **Q**-vector list and post-processing the output spectra from `spinwave()` for the case of twinned crystals and powder (polycrystalline) samples. The base class will be used by default for untwinned crystals.

## Status

Accepted

## Consequences

This will add another new object that users are unfamiliar with. (see backwards compatability proposal) 