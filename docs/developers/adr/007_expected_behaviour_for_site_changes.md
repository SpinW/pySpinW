Expected behaviour for coupling and symmetry when sites change
==============================================================

There are many different ways one can approach the interaction of user site and coupling changes with symmetry,
and changes thereof.

It's a bit much to articulate all the possibilities for how this can work, instead, here is some desirata, and
a description of a system that fits

Desirata
--------

* Symmetry considerations should apply to sites and couplings
* The user should be able to work in a symmetry lower than the one of their system
* Sites that cause conflicts due to a symmetry constraint should be flagged
* Couplings that are impossible due to a symmetry constraint should be flagged
* A set of couplings defined in a high symmetry system should also be definable in a lower symmetry system with as much similarity as possible


The coupling group object
-------------------------

The `CouplingGroup` object defines the coupling between sites in an abstract form.

This object, or it's usage, should be such that it obeys the following transitivity relationship, defined schematically as 

   couplings(symmetry(sites), P1) <=> couplings(P1(sites), symmetry)

That is, the couplings of sites defined by a symmetry group (the sites, and those implied by symmetry) should be the
same as between a list of sites defined in P1 symmetry.

The implementation of this is rather trivial, but it suggests a design choice. That is that the `CouplingGroup` object
acts a list of all equivalent sites according to the symmetry group - as opposed to, for example,
creating couplings then applying a symmetry. 

Implications for couplings under symmetry
-----------------------------------------

There is then a question of what we do with couplings that are invalid due to symmetry. For example, a DM coupling
from A->B along with the same coupling but B->A would have zero contribution to the Hamiltonian (by the anti-symmetry 
of the cross product characteristic of DM couplings).

Do we check this, or just ignore it?

