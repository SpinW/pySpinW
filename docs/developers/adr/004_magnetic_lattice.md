# Classes for different kinds of magnetic structures

## Context

SpinW has different methods for dealing with different kinds of simulation, some simulations only work for systems decribed by a single (possibly complex) propagation vector, some (e.g. "biquadratic") only work
when the propagation vector is (0,0,0). It seems that a useful distinction is between systems where every unique magnetic state is explicitly represented, and ones where it is implicitly defined. 
The notion of commensurability maps on to this. Commensurate propagation vectors permit a finite (though potentially large) explicit representation.

The plan is to make separate classes for these descriptions, provide some degree of conversion between the two, and associate appropriate calculations with each one.
1) The explicit structure - every unique site is described explicity - explcit supercell, only (0,0,0), supports biquadratic
2) The implicit structure - decribed by a unit cell, and a propagation vector - implicit supercell, supports incommensuate and true helical structures. The implcit structure requires the "trick" used in spinw to calculate intensities at +/- the propagation vector.

Possible names could be ExplicitMagneticStructure/ImplicitMagneticStructure, MagneticLattice/MagneticPsudolattice, FullMagneticStructure/PartialMagneticStructure, etc.

## Status

Proposed

## Advantages

* Clear distinction between two kinds of description
* Separation of code
