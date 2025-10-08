Decomposable Systems
====================

A decomposable system is one where the couplings are such that the set of sites can be split into multiple sets where there are no couplings between them.
This might be because we have two separate spinwave systems, or, it might be that we just have a single, non-interacting ion. The two cases are the same from this perspective,
and it is possible that we can decompose into more than just two systems. Clearly, single, uncoupled spins won't exhibit spinwaves, but more complex decomposable systems will.

By definion, a system that is decomposable has multiple, uncoupled subsystems, therefore, it would be possible to excite one without exciting the other. Or, perhaps more relevant,
we can excite one with a given wavevector and energy ($q_1$ and $E_1$) and another with a different wavevector and energy ($q_2$ and $E_2$) with a single neutron with $q = q_1 + q_2$ 
and $E = E_1 + E_2$. This is in-fact a two magnon excitation. If our systems are decomposable, proper treatmeant appears to require doing a multi-magnon calculation, which is problematic 
as whilst achievable in the 2 magnon case with some computational expense, will have an expense that grows as $n^{3m}$, where $n$ is the number of $q$ samples in one dimension, and $m$ is 
the number of components in the system.

Multi-magnon Systems
====================

If the calculation that results in the scattering intensity is the same as that for calculating multiple magnons, then there is another motivation for doing this.

Proposal
========

We should be checking for connectivity of all non-zero moments, and warn the user or error if this is the case.
We should not be trying to automatically work out scattering for these complex cases.


