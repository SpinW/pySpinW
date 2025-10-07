Decomposable Systems
====================

A decomposable system is one where the couplings are such that the set of sites can be split into multiple sets where there are no couplings between them.
This might be because we have two separate spinwave systems, or, it might be that we just have a single, non-interacting ion. The two cases are the same from this perspective,
and it is possible that we can decompose into more than just two systems.

This means that the physics of these systems will be independent, and as such, we can potentially separate out the systems and perfrom the calculations on them independently.
This might be worth doing, fixing some potential edge cases.


Questions
=========

* Is it worth doing this decomposition?
* Once decomposed how do we work out the energy and scattering intensity
  * Energy is trivial, we just add them
  * Scattering intensity is less clear??? Though it should follow simply from the scattering intensity of an non-decomposed system.

Multi-magnon Systems
====================

If the calculation that results in the scattering intensity is the same as that for calculating multiple magnons, then there is another motivation for doing this.
