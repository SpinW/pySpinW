import numpy as np

from pyspinw.gui.rendermodel import rotation_from_z
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite

class Parameterisation:
    def __init__(self, hamiltonian: Hamiltonian, constraints: list[np.ndarray | None]):
        """

        The parameterisation we use is spherical or angular coordinates, centred on the current moment *at every step*
        i.e. we make a new parameterisation each time, where the coordinates are (0,0) or (0).

        We can reuse the algorithm to rotate x,y,z to the z axis, then changes in
        this coordinate system will be a rotation matrix from z, i.e.

        M = [[ cos(a)    sin(a) sin(b)    cos(b) sin(a) ]
             [  0          cos(b)           -sin(b)     ]
             [-sin(a)    cos(a) sin(b)    cos(a) cos(b) ]]

        The total rotation being m = MRz

        The derivative of this with respect to a and b is

        dM/da = [[ 0  0  1 ]
                 [ 0  0  0 ]
                 [-1  0  0 ]]

        dM/db = [[ 0  0  0 ]
                 [ 0  0 -1 ]
                 [ 0  1  0 ]]

        The derivative of the hamiltonian with respect to the parameters will be

        dH/da = dH/m dM/da R z

        and similarly for b.

        For the version constrained to a plane, we want

        """

        self.hamiltonian = Hamiltonian

        sites = hamiltonian.structure.sites
        self.moments = np.array([site.base_moment for site in sites])
        self.magnitudes = np.sqrt(np.sum(self.moments**2, axis=1))
        self.constraints = constraints

        self.directions = self.moments / self.magnitudes

        self._site_uid_to_index = {site.unique_id: i for i, site in enumerate(sites)}
        self.n_sites = len(sites)

    def iterate(self):

        self.rotation_matrices = [rotation_from_z(direction) for direction in self.directions]

        dH_dm = np.zeros((self.n_sites,3))

        # Couplings have different gradients depending on whether or not they're between the same site
        #  Anisotropies always have the same gradient

