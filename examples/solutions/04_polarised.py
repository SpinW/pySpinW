"""
In this worksheet we will explore doing a simulation of an XYZ polarised experiment
"""

from pyspinw import *

# Set up a J1-J2 square antiferromagnet structrure with nearest and next-nearest
# antiferromagnetic Heisenberg exchange interaction, with J2 less than 0.5J1
# and a propagation vector of (0.5, 0.5, 0) and moments in the a-b plane

struc = generate_structure(UnitCell(3,3,6), positions=[[0,0,0]], spins=[[1,0,0]],
          propagation_vectors=[[0.5,0.5,0]])
exchanges = generate_exchanges(struc, bond=1, j=1) + generate_exchanges(struc, bond=2, j=0.2)
ham = Hamiltonian(struc, exchanges)

# Now plot the (unpolarised) dispersion along the (100) and (110) directions

ham.spaghetti_plot(Path([[1,0,0],[0,0,0],[1,1,0]]))

# Now plot the same dispersion with Pxx, Pyy, Pzz polarisation,
# assuming that the sample was mounted with the a-b plane horizontal
# These are the components we can measure with longitudinal polarisation analysis

ham.spaghetti_plot(Path([[1,0,0],[0,0,0],[1,1,0]]), components='Pxx')

ham.spaghetti_plot(Path([[1,0,0],[0,0,0],[1,1,0]]), components='Pyy')

ham.spaghetti_plot(Path([[1,0,0],[0,0,0],[1,1,0]]), components='Pzz')

# QUESTION: What do you notice about the spectra compared to the unpolarised?
# ANSWER:   In an XYZ polarisation analysis experiment you can only measure
#           the spin component perpendicular to the polarisation *and* Q.
#           Since we are using the Blume-Maleev coordinate system where x||Q
#           we cannot see any precesion components parallel to x.
#           Furthermore, since we are only looking at the non-spin-flip
#           components here, we cannot see any magnetic signal in Pxx, and
#           Pyy is only sensitive to precession along z, whilst
#           Pzz is only sensitive to precession along y.
#           
#           We see only signal in Pzz along (110) implying precession in-plane.
#           Along (100) we see no signal at all, implying precession along (100)
