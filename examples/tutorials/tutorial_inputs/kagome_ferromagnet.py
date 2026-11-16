""" Kagome Ferromagnet example

Reproduces Tutorial 5: https://spinw.org/tutorials/05tutorial
"""
## title: Ferromagnetic Kagome lattice
## reproduces: 5, 6

## subtitle: Define the lattice

# We define the kagome lattice using the $P\overline{3}$ space group (Nr. 147).
# Using the 3 fold symmetry, we need to define only one magnetic atom.
# In this space group all first neighbour couplings will be equivalent (related by symmetry).
# Symmetry equivalent positions of the magnetic atoms are automatically generated when we specify the spacegroup.
# The magnetic atoms are Cu$^+$ with $S=1$ spin.

## capture-stdout
from pyspinw import *
unit_cell = UnitCell(6,6,10, gamma=120)
site = Site(0.5, 0, 0, 0, 0, 1, name="X")
struc = Structure([site], unit_cell=unit_cell, spacegroup=get_spacegroup('P -3'))
struc.print_summary()
## skip
view(struc)
## end-capture-stdout

## image
## snapshot(struc, "kagome_ferromagnetic_atoms.png", view_point=(0,0,15))

# Note that unlike in the equivalent Matlab tutorial here we place the ferromagnetic spins along the $c$-axis.
# This is because PySpinW enforces the symmetry operations on the spins as well as atoms, and ferromagnetic
# spins in any other direction will not satisfy the 3-fold symmetry axis along $c$.

## subtitle: Nearest-neighbour dispersion

# We new define a ferromagnetic nearest-neighbour only and create a Hamiltonian from this
# then use it to compute the spin wave dispersion
## capture-stdout
exchanges = generate_exchanges(struc,
                               bond=1,
                               exchange_type=HeisenbergExchange,
                               j=-1)

hamiltonian = Hamiltonian(struc, exchanges)

hamiltonian.print_summary()
## end-capture-stdout
## skip
view(hamiltonian)

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])

## skip
hamiltonian.energy_plot(path)

## image
## snapshot(hamiltonian, "kagome_ferromagnetic_structure.png", view_point=(0,-10,10))

## image
## fig = hamiltonian.energy_plot(path, show=False)
## fig.savefig("kagome_ferromagnet_nearest_disp.png")

## subtitle: Nearest-neighbour powder

# Now compute and plot the powder spectrum. ``n_samples`` is the number of $\mathbf{Q}$-points
# per $|Q|$-shell to use and ``method`` is how to sample a spherical shell to obtain the powder
# average. The default is ``method='fibonacci'``, and another possiblity is ``method='geodesic'``

sample = Powder(hamiltonian)
path1D = Path1D(0.1, 2.5, n_q=100)
## skip:2
sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
                     method='random', vmin=0, vmax=0.1)
## image
## fig = sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
##                            method='random', vmin=0, vmax=0.1, show=False)
## fig.savefig("kagome_ferromagnet_nearest_powder.png")

# We can see that the flat mode at 6 meV (a van Hove singularity) shows up the in powder spectrum as an extremely sharp band.
# Next we plot again but convolve the calculation with a Gaussian with a FWHM of 0.25 meV in energy.

## skip:2
sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
                     energy_stddev=0.25/2.35, method='random', vmin=0, vmax=1)
## image
## fig = sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
##                            energy_stddev=0.25/2.35, method='random', vmin=0, vmax=1, show=False)
## fig.savefig("kagome_ferromagnet_nearest_powder_dE.png")

## subtitle: Further neighbours

# Now, we create the kagome lattice with up to 3rd neighbor interactions.
# The symmetry related atoms are denoted by MCu1($i$)$_j$, where $i$ is the index of independent atomic positions,
# $j$ is the index of the generated atomic positions of the $i$-th independent position.

## capture-stdout
from pyspinw import *
unit_cell = UnitCell(6,6,10, gamma=120)
site = Site(0.5, 0, 0, 0, 0, 1, name="Cu1")
struc = Structure([site], unit_cell=unit_cell, spacegroup=get_spacegroup('P -3'))
struc.print_summary()
## end-capture-stdout

# We add couplings up to 3rd neighbor interactions.
# If the generation of the bond tables would depend on distance, J3a, J3b and J3c would be equivalent.
# However using the $P\bar{3}$ space group the three type of bonds are inequivalent, as physically expected in real systems
# (J3a goes through an intermediate magnetic atom, while the other two bonds are not).
# To see these three inequivalent sets which we use the ``use_sym`` option to ``generate_exchanges``.
# We then set values for these different sets of bonds and plot them:

## capture-stdout
j3sym = generate_exchanges(struc, bond=3, j=-1, use_sym=True)
for n, exchanges in enumerate(j3sym):
    print(f'Set {n}:')
    print('\n'.join([str(ex) for ex in exchanges]))
j3a = [HeisenbergExchange(j.site_1, j.site_2, 0.0, j.cell_offset, color=[1, 0.6, 0]) for j in j3sym[0]] 
j3b = [HeisenbergExchange(j.site_1, j.site_2, 0.17, j.cell_offset, color=[0, 0, 1]) for j in j3sym[1]] 
j3c = [HeisenbergExchange(j.site_1, j.site_2, 0.0, j.cell_offset, color=[0.5, 0, 0.5]) for j in j3sym[2]] 
h3 = Hamiltonian(struc, j3a + j3b + j3c)
## skip
view(h3, copies=(2,2,1))
## end-capture-stdout
## image
## snapshot(struc, "kagome_j3.png", view_point=(0,-10,10))

# We can see that the blue bonds do not pass through any intermediate Cu atoms but the other two (orange and purple do).
# Now we define the nearest and next nearest neighbours and then the full Hamiltonian
j1 = generate_exchanges(struc, bond=1, j=-1, color=[0, 1, 0])
j2 = generate_exchanges(struc, bond=2, j=0.1, color=[1, 0, 0])
hamiltonian = Hamiltonian(struc, j1 + j2 + j3a + j3b + j3c)

## subtitle: 3rd neighbour spin wave dispersion

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
## skip
hamiltonian.energy_plot(path)
## image
## fig = hamiltonian.energy_plot(path, show=False)
## fig.savefig("kagome_3nn_disp.png")

# And finally the powder spectrum

sample = Powder(hamiltonian)
path1D = Path1D(0.1, 2.5, n_q=100)
## skip:2
sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
                     method='random')
## image
## fig = sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
##                            method='random', vmin=0, vmax=0.1, show=False)
## fig.savefig("kagome_3nn_powder.png")

