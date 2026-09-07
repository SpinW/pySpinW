""" Ferromagnetic chain example """

## title: Ferromagnetic Chain
## reproduces: 1

## subtitle: Introduction
# A ferromagnetic chain is the simplest system that we can simulate with pySpinW

## subtitle: A basic calculation using
# Import the main pySpinW module, doing this will give access to the majority of pySpinW methods and classes
# that are needed for most tasks.
from pyspinw import *

unit_cell = UnitCell(1,1,1)

only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

s = Structure([only_site], unit_cell=unit_cell)

exchanges = generate_exchanges(sites=[only_site],
                               unit_cell=unit_cell,
                               max_distance=1.1,
                               exchange_type=HeisenbergExchange,
                               j=-1,
                               direction_filter=filter([1,0,0]))

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [1,0,0]])

hamiltonian.spaghetti_plot(path, dE=0.4)
