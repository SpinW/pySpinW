"""
In this worksheet we will explore using the rotating frame calculations
on a helical system with a Dzyaloshinskii-Moriya interaction
"""

from pyspinw import *
import time

# Create a helical magnet with the helix propagating in the a-direction 
# in a tetragonal structure with a pitch of 30 degrees.
# Make the structure have a short a-axis, and long b- and  c-axes.

struc = generate_helical_structure(????

# Now, we want to investigate how the Heisenberg and Dzyaloshinskii-Moriya
# interactions compete to generate a helical structure, and for this we
# need to define a supercell rather than use the rotating-frame system.
# This is because the energy minimisation routine to relax the spin
# structure requires all the spins to be independent (whereas for the 
# rotating frame the spins would be fixed by their pitch / propagation vector)

# This line converts the structure to a supercell
supercell = Structure(sites=struc.sites, unit_cell=struc.unit_cell,
                supercell=struc.supercell.approximant())

view(supercell)

#%%

# Now add two nearest neighbour interaction along the a-direction:
#   1. An antiferromagnetic Heisenberg interaction
#   2. A DM interaction with the DM-vector along the normal of the plane of the helix.
# (note that you have to specify all d_x, d_y, d_z components else they default to 1.0)

exchanges = generate_exchanges(supercell, ???

# We have to use an expanded supercell in order to make the spins independent for the
# ground state calculation below
ham = Hamiltonian(supercell, exchanges).expanded()

optstruc = ham.ground_state(initial_randomisation='randomised')

view(optstruc)
optstruc.print_summary()

# QUESTION: What optimised structure / ground state to do you get when the
#           DM interaction is small or zero?
# ANSWER:   When D=0 you should get an antiferromagnetic chain
# QUESTION: What happens as you increase D w.r.t. the Heisenberg J?
# ANSWER:   As D approaches around 50% of J, you should get a helix.
# QUESTION: What happens if you change the direction of the DM-vector?
# ANSWER:   The plane of the helix changes to be perpendicular to D

#%%

# Now let's plot the dispersion along the (100) direction.
# Note that because we have used an expanded supercell the reciprocal lattice now 
# extends to 1/k where k is the propagation vector chosen at the start of the 
# worksheet.

# First plot the *un*-optimised supercell

ham.spaghetti_plot(???

# Now let's compare this dispersion to what we would get with a rotating frame
# calculation. We should also time the calculations to see which is faster.
# (You can use time.time() to record the current system time)
# Note that although we can use the original structure defined at the start
# we need to redefine the exchange interactions w.r.t. this structure
# And likewise the Hamiltonian needs to be defined w.r.t. these new exchanges.

ex0 = generate_exchanges(struc, ???

ham0 = Hamiltonian(struc, ex0)
ham0.spaghetti_plot(???

# QUESTION: How does the dispersion compare in the two different methods?
# ANSWER:   They should look almost the same
# QUESTION: Which calculation is faster? (What happens if you decrease the pitch angle 
#           or reduce the propagation vector? [e.g. increase the supercell size])
# ANSWER:   The rotating frame calculation is shorter as it needs to solve a smaller
#           system with fewer magnetic atoms; but there are overhead costs in starting
#           the calculation such that for small-ish supercell (~10 atoms) it's not
#           significantly faster.

# Finally let's plot the *optimised* (ground-state) supercell structure's dispersion

optstruc.spaghetti_plot(???

# QUESTION: What do you notice about this dispersion compared to the previous two?
# ANSWER:   It could be very different depending on the value of J and D. This is 
#           because the rotating frame calculation locks in a helix but in real life
#           you might not have such a structure, depending on the values of J/D etc.
# QUESTION: What happens to the dispersion of the three models if you change the
#           values of J and/or D?
# ANSWER:   You should notice that as J tends to zero the dispersion of the three
#           models should converge.
