Handling Higher Symmetry
========================

We will take the approach that one does not need to represent every component of the system to run a calculation. 
This allows a user to represent only the minimal number of atoms for performing a calculation, and also has 
the advantage of simplifying the system.

However, this leads to questions about how we treat symmetry. For example, 
a P1 lattice with only one site will have much higher symmetry than P1 in the sense that
there will be more groups that can leave it invariant than P1.

Once we make the decision that the symmetry group we specify is a "lower bound", 
there a problem arrises - one that exists in crystallographic group symmetry 
already - but that becomes more of a general issue: the setting problem. That is to say,
depending on how you translate the crystallographic origin within the unit cell, 
the definitions of groups can be different. Outside of the abstract categorisation 
of symmetries, the choice of origin is more significant as we are no longer choosing it
for the sole purpose of finding unique groups.  

In other words, the categorisation done in group theory relies on being able to choose the origin
to be a high symmetry point, but when we don't represent the full detail of our system, 
this might not correspond to anything we have explicitly represented.

Consequence
===========

Ultimately, there is an obvious and natural design choice that accompanies this. 
That is, we assume the origin used for symmetry is the cartesian origin we're using for
building the system. 

This is something that should be remembered when choosing groups, and may have consequences
when interacting with other software and file formats.



