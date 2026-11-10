""" Kagome Ferromagnet with 3 nearest neighbours

Reproduces Tutorial 6: https://spinw.org/tutorials/06tutorial
"""
## title: Kagome ferromagnet with 3 nearest neighbours
## reproduces: 6

# We create the kagome lattice with up to 3rd neighbor interactions.
# The symmetry related atoms are denoted by MCu1($i$)$_j$, where $i$ is the index of independent atomic positions,
# $j$ is the index of the generated atomic positions of the $i$-th independent position.

## subtitle: Define the lattice

## capture-stdout
from pyspinw import *
unit_cell = UnitCell(6,6,10, gamma=120)
site = Site(0.5, 0, 0, 0, 0, 1, name="Cu1")
struc = Structure([site], unit_cell=unit_cell, spacegroup=get_spacegroup('P -3'))
struc.print_summary()
## end-capture-stdout

## subtitle: Define Hamiltonian

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

## We can see that the blue bonds do not pass through any intermediate Cu atoms but the other two (orange and purple do).
## Now we define the nearest and next nearest neighbours and then the full Hamiltonian
j1 = generate_exchanges(struc, bond=1, j=-1, color=[0, 1, 0])
j2 = generate_exchanges(struc, bond=2, j=0.1, color=[1, 0, 0])
hamiltonian = Hamiltonian(struc, j1 + j2 + j3a + j3b + j3c)

## subtitle: Spin wave dispersion

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
## skip
hamiltonian.energy_plot(path)
## image
## fig = hamiltonian.energy_plot(path, show=False)
## fig.savefig("kagome_3nn_disp.png")

## subtitle: Powder averaged spectrum

sample = Powder(hamiltonian)
path1D = Path1D(0.1, 2.5, n_q=100)
## skip:2
sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
                     method='random')
## image
## fig = sample.show_spectrum(path1D, n_samples=1000, max_energy=7, n_energy_bins=200,
##                            method='random', vmin=0, vmax=0.1, show=False)
## fig.savefig("kagome_3nn_powder.png")
