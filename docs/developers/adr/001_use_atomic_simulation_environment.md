# PySpinW will use the Atomic Simulation Environment as a dependency

## Context

The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) is a Python library for setting up, running, visualizing and analyzing atomistic simulations.
Primarily, it supports `calculators` based on density function theory or force fields.
In most cases these codes are not written in Python but `ase` provides an interface to drive the calculations from Python.
In particular it has an `Atoms` [class](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) which embodies a collection of atoms which can be the lattice which SpinW calculations need.
Each `Atom` in an `Atoms` object has a `position` and magnetic moment `magmom`, as well as `charge` and `mass` which SpinW does not need.
Most importantly, the `Atoms` can be visualized using the `ase.visualize.view` [function](https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html)
which provides an internal viewer (`ase.gui`) as well as interface to a variety of external viewers.
Finally, there is a spacegroup package with a `crystal` [constructor](https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html) which creates an `Atoms` object with a particular space group symmetry.
Thus much of the functionality provided by the current `genlattice` and `addatom` methods could be outsourced to `ase` reducing the amount of code needed to be written in PySpinW.


## Decision

We will use the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) as a dependency in PySpinW to handle constructing and manipulating a lattice.


## Status

Proposed


## Consequences

This adds a significant dependency to the project. The possible downsides are:

* Needing to adapt spinw code to the `Atoms` class which may not be flexible enough for our needs.
* If `ase` changes its code this might break anything we might build on top of it.
* Risk that the `ase` project becomes abandoned (considered low as it has a large user base and an STFC staff member is a contributor).

The advantages is:

* Saving a large amount of coding in order to write a `Lattice` class or the lattice handling part of the `spinw` class.
