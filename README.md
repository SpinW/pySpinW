
![](./docs/branding/logo_with_background_small.png)

pySpinW is a tool for calculating magnon energies and simulating scattering experiments   

Details of its use can be found at www.spinw.org/pyspinw

## Example

Here is an example of calculating dispersion curves for simple ferromagnetic chain in python:

```python

from pyspinw import *

# Define the unit cell
unit_cell = UnitCell(1,1,1)

# Specify a magetic atom at (0,0,0), with moment (0,0,1)
x = LatticeSite(0, 0, 0, 0, 0, 1, name="X")

# Create a magnetic structure (this could include a spacegroup, or 
# supercell structure, but we don't do so here)
structure = Structure([x], unit_cell=unit_cell)

# A single exchange between atoms in neighbouring unit cells (in the x direction)
exchanges = [HeisenbergExchange(x, x, cell_offset=(1,0,0), j=-1)]

# Create a Hamiltonian object that holds things together 
hamiltonian = Hamiltonian(structure, exchanges)

# Define a path through q-space for out dispersion curve
path = Path([[0,0,0], [1,0,0]])

# Show the spaghetti plot, with constant energy smearing of 0.4
hamiltonian.spaghetti_plot(path, dE=0.4)

```



## New features in the python version

There are a number of features in pySpinW that did not exist in the MATLAB version. 

 * A viewer that lets you see magnetic structures and their data simultaneously
 * Parameterisation of hamiltonians
 * Saving of objects (note: we are still in alpha, so expect there to be changes in file format in the near future)
 * Approximation of incommensurate structures
 * Tools for working with spacegroups

## Functional parity with SpinW

pySpinW is still in alpha, this means that a small number of features of SpinW are not present or complete
in pySpinW and that some features are not considered stable.

Not fully ported:
 * Fitting of incommensurate k-vectors (work in progress)
 * Powder fitting tools (basic powder fitting available, but not all tools are ported yet)

