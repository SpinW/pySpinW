
from pyspinw import *
import time

# Create a helical magnet with the helix propagating in the a-direction
# in a tetragonal structure with a pitch of 30 degrees.
# Make the structure have a short a-axis, and long b- and  c-axes.

struc = generate_helical_structure(UnitCell(3, 4, 5),
            positions=[[0,0,0]], spins=[[1,0,0]],
            propagation_vector=[30/180., 0, 0], perpendicular=[0,0,1])

# Now, we want to investigate how the Heisenberg and Dzyaloshinskii-Moriya
# interactions compete to generate a helical structure, and for this we
# need to define a supercell rather than use the rotating-frame system.
# This is because the energy minimisation routine to relax the spin
# structure requires all the spins to be independent (whereas for the
# rotating frame the spins would be fixed by their pitch / propagation vector)

# This line converts the structure to a supercell
supercell = Structure(sites=struc.sites, unit_cell=struc.unit_cell,
                supercell=struc.supercell.approximant())

exchanges = generate_exchanges(supercell, bond=1, j=1) + \
       generate_exchanges(supercell, bond=1, exchange_type=DMExchange, d_x=0.0, d_y=0.0, d_z=1.0)

# We have to use an expanded supercell in order to make the spins independent for the
# ground state calculation below
ham = Hamiltonian(supercell, exchanges).expanded()

optstruc = ham.ground_state(initial_randomisation='randomised')

optstruc.print_summary()

path = Path([[0,0,0],[10,0,0]])
# ham.spaghetti_plot(path)
ham.new_spaghetti_plot(path)
