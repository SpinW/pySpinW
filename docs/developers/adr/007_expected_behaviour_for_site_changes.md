Expected behaviour for coupling and symmetry when sites change
==============================================================

When the sites are modified the associated couplings.

Sites are split into two kinds
1) real/explicit/independent sites, these are set by the user, and correspond to magnetic freedoms in the solver
2) implied/virtual sites, these are implied by the symmetry and are completely determined real sites and the symmetry




There are different ways that you can approach symmetry in the system
1) Use the correct symmetry for your system
2) Use lower symmetry
3)


Adding a site
-------------

This is trivial

Moving a site
-------------

Moving a site can affect the symmetry, so, we need to check if this is now a duplicate of an existing site (i.e. if it
has been moved to a higher symmetry point with a corresponding )

Changing the magnetism on a site
--------------------------------



Reifying a site
---------------

Virtualising a site
-------------------

Removing a site
---------------

  * What happens when an explicit site is removed, but it is replaced with a virtual site