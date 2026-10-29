""" Kagome Ferromagnet example

Reproduces Tutorial 5: https://spinw.org/tutorials/05tutorial
"""
## title: Ferromagnetic Kagome lattice
## reproduces: 5

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
