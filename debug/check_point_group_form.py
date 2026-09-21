""" Check for the form of point group matrices

 We assume that it integer valued, is that reasonable
 """
import numpy as np
import spglib

for hall_number in range(1, 531):
    # Get the relevant data

    print(f"Hall {hall_number}:")

    op_data = spglib.get_symmetry_from_database(hall_number)

    rotations = op_data["rotations"]

    for i, rotation in enumerate(rotations):
        if np.allclose(np.array(rotation, dtype=int), rotation):
            print(f"  {i}: OK")
        else:
            print(f"  {i}:FAILED")
            print(rotation)
            assert False
