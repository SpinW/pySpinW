import numpy as np

from pyspinw.basis import angle_axis_rotation_matrix
from pyspinw.checks import check_sizes
from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.supercell import CommensurateSupercell, CommensuratePropagationVector, TransformationSupercell


@check_sizes(k=(3,), force_numpy=True)
def rotation_supercell(k: np.ndarray,
                       axis: np.ndarray | None =None,
                       initial_moment: np.ndarray | None = None,
                       orthogonal: bool=False) -> CommensurateSupercell:

    # Make sure k is a numpy array
    k = np.array(k)

    # Make sure axis is a numpy array
    if axis is None:
        axis = np.array([0,0,1])
    else:
        axis = np.array(axis)

    # Handle different options
    if orthogonal:
        k_vectors = [CommensuratePropagationVector(*v) for v in np.diag(k)]
        frequencies = [(2*np.pi)*component for component in k]
    else:
        k_vectors = [CommensuratePropagationVector(*k)]
        frequencies = [2*np.pi*np.sqrt(np.sum(k**2))]

    transforms = [angle_axis_rotation_matrix(frequency, axis) for frequency in frequencies]

    supercell = TransformationSupercell(propagation_vectors=k_vectors, transforms=transforms)

    if initial_moment is None:
        return supercell
    else:
        return supercell.summation_form(initial_moment)


def demo_rotation_supercell():
    import matplotlib.pyplot as plt

    m = np.array([1,0,0])

    supercell = rotation_supercell(k=[1/3, 1/2, 0], axis=[0, 0, 1], orthogonal=True)

    for i in range(4):
        for j in range(3):
            cell_offset = CellOffset(i=i, j=j, k=0)

            moment = supercell.evaluate(cell_offset, m)

            plt.arrow(i, j, moment[0], moment[1])

    plt.show()

if __name__ == "__main__":
    demo_rotation_supercell()
