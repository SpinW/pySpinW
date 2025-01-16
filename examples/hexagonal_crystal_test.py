from ase.spacegroup import crystal
from ase.visualize import view
import numpy as np

data = crystal("CNS", basis=np.array([
    [0,1,0],
    # [0,1,1], # Top of cell same as bottom
    [0.5, 0.5, 0.5],
    [0, 0, 0.5]
]),
        cellpar=[10,10,10,90,90,90],
        spacegroup=168,
        size=(3,3,3))



view(data)