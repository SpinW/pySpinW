#  Antiferromagnetic chain with rotating-frame calculation

Source: [python (.py)](antiferromagnetic_chain.py) | [jupyter (.ipynb)](antiferromagnetic_chain.ipynb) 


_This tutorial mirrors MATLAB spinW Tutorial 2_

 In this tutorial we will then take three different approaches to calculating the magnon spectrum
  - Creating unit cell containing with two spins pointing in different directions.
  - Using a unit cell with a single spin, performing a calculation that *explicitly* creates two unit cells
    containing opposing spins.
  - Using the "rotating frame" optimisation, which *implicitly* defines two cells. This has advantages
    in terms of speed and clarity for systems where it can be applied.

 As in the previous tutorial we start by importing pyspinw

```python
from pyspinw import *
```
##  Two spins per unit cell
 We create a magnetic structure with a 2x1x1 unit cell and two atoms - at 0.5 and 1.5 angstroms in $x$,
 facing in different directions in $z$.

```python
two_site_unit_cell = UnitCell(2, 1, 1)
two_site_site_1 = LatticeSite(0.25, 0.5, 0.5, sz=1, name="a")
two_site_site_2 = LatticeSite(0.75, 0.5, 0.5, sz=-1, name="b")
two_site_structure = Structure([two_site_site_1, two_site_site_2], two_site_unit_cell)
```
 Like in the previous example, there are different ways we can set up the exchanges, if we do it explicitly we
 have

```python
two_site_exchanges = [
    HeisenbergExchange(two_site_site_1, two_site_site_2, j=1),
    HeisenbergExchange(two_site_site_2, two_site_site_1, cell_offset=(1,0,0), j=1)
]
```
 Or using `generate_exchanges`, which should give the same result

```python
two_site_exchanges = generate_exchanges([two_site_site_1, two_site_site_2],
                               unit_cell=two_site_unit_cell,
                               max_distance=1.1,
                               direction_filter=filter([1,0,0], symmetric=True))
```
 If we check to see if it is the same we find that it is not identical to the hand-crafted version
 but it is equivalent, we have exchanges between a and b, one of which is offset by one cell in the
 $x$ (lattice a) direction

```python
for exchange in two_site_exchanges:
    print(exchange)
```

```text
HeisenbergExchange('J1:a-b', a <-> b, offset=(-1, 0, 0), j=1)
HeisenbergExchange('J1:a-b', a <-> b, offset=(0, 0, 0), j=1)
```

 Now we build a hamiltonian

```python
two_site_hamiltonian = Hamiltonian(two_site_structure, two_site_exchanges)
```
 View it

```python
view(two_site_hamiltonian)
```

![]( two_site_hamiltonian.png
)

 And plot the dispersion, here we specify the path in angstroms to make it comparable with different unit cells.

```python
path = Path([(0, 0, 0), (1, 0, 0)], convert_to_lattice_units_with=two_site_unit_cell)
```

```python
two_site_hamiltonian.spaghetti_plot(path)
```

![]( two_site_dispersion.png
)

##  One spin per unit cell (normal frame)
 We'll now do calculations on the same system, but using a magnetic cell that is different from the unit cell.
 We call magnetic cells supercells, and there are different options for this. First we'll look at the case
 which simply creates a 2x1x1 system explicitly during the calculation. There are a few options for this,
 we will use a `TransformationSupercell`. This applies transformations to the original spin based on the relative
 positions of the cells. We will tell it to use a rotation around the $b$ axis ($y$), with a 2 cell
 propagation vector in the $a$ axis ($x$), i.e. (1/2,0,0).

```python
one_spin_unit_cell = UnitCell(1,1,1)
one_spin_site = LatticeSite(0.5, 0.5, 0.5, sz=1, name="S")
supercell = rotation_supercell(directions=[(0.5, 0, 0)], axes=[(0, 1, 0)])
structure = Structure([one_spin_site], one_spin_unit_cell, supercell=supercell)
```
 We only need to specify one exchange in this case, as it is implicit in the supercell description that
  exchanges are the same in each repetition. Again, $j=1$ for an antiferromagnet.

```python
one_site_exchange = HeisenbergExchange(one_spin_site, one_spin_site, cell_offset=(1,0,0), j=1)
hamiltonian = Hamiltonian(structure, [one_site_exchange])
hamiltonian.print_summary()
```
 View it

```python
view(hamiltonian)
```

![]( one_site_hamiltonian.png
)


```python
path = Path([(0,0,0), (1,0,0)], convert_to_lattice_units_with=one_spin_unit_cell)
```

```python
hamiltonian.spaghetti_plot(path)
```

![]( one_site_dispersion_1.png
)

 We can view this too

```python
view(structure)
```

![]( first_one_spin_structure
)

 The following generates a 2x1x1 supercell in which the spins alternate in y (period 2 rotation around z)

```python
structure = generate_helical_structure(one_spin_unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])
```
 Generate Heisenberg exchanges based on distance, this will only be in x because of the shape of the unit cell

```python
exchanges = generate_exchanges(sites=structure,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)
```
 Build the Hamiltonian

```python
hamiltonian = Hamiltonian(structure, exchanges)
```
 And create a spaghetti plot. We will be comparing structures with different unit cells, so,

```python
path = Path([(0,0,0), (1,0,0)])
```
 You can get a summary of the Hamiltonian using `print_summary` (also works on `Structure`)

```python
hamiltonian.print_summary()
```

```text
Unit Cell: (1, 1, 1), alpha=90, beta=90, gamma=90
Spacegroup: P1
Supercell: Rotation
  size = (1, 1, 1)
  PropagationVector(0.5, 0, 0)
  Perpendicular=[0, 0, 1]
Sites:
  Site(MCu1, 0, 0, 0, spin=[0. 1. 0.])
Exchanges:
  HeisenbergExchange('J1:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 0, 1), j=1)
  HeisenbergExchange('J1:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 1, 0), j=1)
  HeisenbergExchange('J1:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 0, 0), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, -1, 1), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 0, 1), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 0, 1), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 1, 1), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 1, 0), j=1)
  HeisenbergExchange('J2:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 1, 0), j=1)
  HeisenbergExchange('J3:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, -1, 1), j=1)
  HeisenbergExchange('J3:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, -1, 1), j=1)
  HeisenbergExchange('J3:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 1, 1), j=1)
  HeisenbergExchange('J3:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 1, 1), j=1)
  HeisenbergExchange('J4:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 0, 2), j=1)
  HeisenbergExchange('J4:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 2, 0), j=1)
  HeisenbergExchange('J4:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 0, 0), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, -2, 1), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, -1, 2), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 0, 1), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 0, 2), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 0, 2), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 0, 1), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 1, 2), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 2, 1), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 1, 0), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 1, 0), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 2, 0), j=1)
  HeisenbergExchange('J5:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 2, 0), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, -2, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, -2, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, -1, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, -1, 2), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, -1, 2), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, -1, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 1, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 1, 2), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 1, 2), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 1, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 2, 1), j=1)
  HeisenbergExchange('J6:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 2, 1), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, -2, 2), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 0, 2), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 0, 2), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 2, 2), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 2, 0), j=1)
  HeisenbergExchange('J7:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 2, 0), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, -2, 1), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, -2, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, -2, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, -2, 1), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, -1, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, -1, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 0, 3), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 1, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 1, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-2, 2, 1), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(-1, 2, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(1, 2, 2), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(2, 2, 1), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(0, 3, 0), j=1)
  HeisenbergExchange('J8:MCu1-MCu1', MCu1 <-> MCu1, offset=(3, 0, 0), j=1)
```

 View it using view

```python
view(hamiltonian)
```

![]( structure.png
)

 Now we can plot a spaghetti diagram, first we generate a path though the lattice
 from $q=0$ to $1$ reciprocal lattice unit in $x$.

```python
path = Path([[0,0,0], [1,0,0]])
```
 Show a plot

```python
hamiltonian.spaghetti_plot(path, scale='log')
```

![]( spaghetti_plot.png
)

