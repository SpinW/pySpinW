# Classes for different kinds of magnetic structures

## Context

SpinW has different methods for dealing with different kinds of simulation, some simulations only work for systems decribed by a single (possibly complex) propagation vector, some (e.g. "biquadratic") only work
when the propagation vector is (0,0,0). It seems that a useful distinction is between systems where every unique magnetic state is explicitly represented, and ones where it is implicitly defined. 
The notion of commensurability maps on to this. Commensurate propagation vectors permit a finite (though potentially large) explicit representation.

The plan is to make separate classes for these descriptions, provide some degree of conversion between the two, and associate appropriate calculations with each one.
1) The explicit structure - every unique site is described explicity - explcit supercell, only (0,0,0), supports biquadratic
2) The implicit structure - described by a unit cell, and a single propagation vector - implicit supercell, supports incommensurate and true helical structures. The impilcit structure requires the "trick" used in spinw to calculate intensities at +/- the propagation vector.

In both cases the classes should implement a method to return the set of `zed` and `eta` vectors needed for the spin wave calculations as detailed in [this doc](../design/001_linear_spinwave_theory.md#local-spin-directions).

In addition, it would be useful if the user can specify a set of propagation vector(s) and basis(es) in the explicit representation rather than every spin orientation in the unit cell, but this is an implementation detail to be discussed later.

## Decision

We will create two classes for handling the magnetic structure, a `CommensurateStructure` for the "explicit" representation and a `RotatingFrameStructure` for the "implicit" representation.

## Status

Accepted

## Advantages

* Clear distinction between two kinds of description
* Separation of code
