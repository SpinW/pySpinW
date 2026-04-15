"""
In this worksheet we will explore using the rotating frame calculations
on a helical system with a Dzyaloshinskii-Moriya interaction
"""

from pyspinw import *

# Create a helical magnet with the helix propagating in the a-direction 
# in an orthorhombic structure with a pitch of 30 degrees
# Make the orthrhombic structure have a short a-axis, intermediate b-axis
# and long c-axis.

struc = generate_helical_structure(UnitCell(3, 4, 5),
            positions=[[0,0,0]], spins=[[1,0,0]],
            propagation_vector=[10/180., 0, 0], perpendicular=[0,0,1])

# Note that "view" currently doesn't plot the rotating frame structure correctly 
# So we must first convert it into a supercell before plotting
supercell = Structure(sites=struc.sites, unit_cell=struc.unit_cell,
                supercell=struc.supercell.approximant())

#view(supercell)

#%%

# Now add a nearest neighbour Dzyaloshinskii-Moriya interaction along 
# the a-direction with the DM vector along the normal of the plane of the helix.
# Then add a next-nearest-neighbour antiferromagnetic Heisenberg 

exchanges = generate_exchanges(supercell, bond=1, j=1)# + \
#       generate_exchanges(supercell, bond=1, exchange_type=DMExchange, d_z=1.)

ham = Hamiltonian(supercell, exchanges)

#optstruc = ham.ground_state(initial_randomisation='randomised', planar_axis=[0,0,1])
optstruc = ham.ground_state(fixed=[supercell.sites[0]], step_size=0.0,
                            initial_randomisation='randomised')

view(optstruc)
