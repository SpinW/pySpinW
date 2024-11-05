# A Structure for Describing Macroscopic Structure

## Context

The main job of spinW is to calculate the energy levels of a magnetic system over a range of q values. This
task is relatively self-contained and well defined by just the material's cystallographic and magnetic properties.

However, there are a number of things that are more properly understood as a property of a particular experiment. 
The main features are to include here is whether or not the crystal is twinned or is a powder.

There is a question about whether it should be the location for the option to apply neutron scattering selection rules.


## Decision

Create a class representing the macroscopic structure that is responsible for working out the appropriate inputs for the
main calculation, and constructing the measurable output.

## Status

Proposed


## Consequences

This will add another new object that users are unfamiliar with. (see backwards compatability proposal) 