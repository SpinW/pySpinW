"""Calculates polarised cross-sections from Sab matrix"""
import numpy as np
from numpy._typing import ArrayLike
from enum import Enum
import re

# ruff: noqa: E702
class PolarisationType(Enum):
    """Types of Polarisation components"""

    PXX = 'Pxx'; PXY = 'Pxy'; PXZ = 'Pxz'  # Full polarisation
    PYX = 'Pyx'; PYY = 'Pyy'; PYZ = 'Pyz'
    PZX = 'Pzx'; PZY = 'Pzy'; PZZ = 'Pzz'
    PX  = 'Px';  PY  = 'Py';  PZ  = 'Pz'   # Half polarisation (incident only)
    MXX = 'Mxx'; MXY = 'Mxy'; MXZ = 'Mxz'  # Projection of Sab in the Blume-Maleev
    MYX = 'Myx'; MYY = 'Myy'; MYZ = 'Myz'  #   cartesian coordinate system
    MZX = 'Mzx'; MZY = 'Mzy'; MZZ = 'Mzz'

_STRFLOATDIC = {'':1., '+':1., '-':-1.}
_ABTOIND = {'x':0, 'y':1, 'z':2}

# The Levi-Civita symbol in 3 dimensions
LC = np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 0]],
               [[0, 0, -1], [0, 0, 0], [1, 0, 0]],
               [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]])

def _str_to_vec(in_str):
    ma = re.match(r'\[([0-9]*)\ *([0-9]*)\ *([0-9]*)\]', in_str)
    if ma:
        try:
            return np.array([float(v) for v in ma.groups()])
        except ValueError:
            raise RuntimeError(f'Could not parse string "{in_str}" as a vector')
    raise RuntimeError(f'Could not parse string "{in_str}" as a vector')

def parse_components(components: str, rlu_to_cart: ArrayLike):
    """Parses a polarised components string

    :param components: the component string which is a mathematical expression involving polarisation
                       components Pab, Pa or Mab, and optionally the horizontal scattering plane
                       defined after a comma. If not specified, we assume the a-b plane is horizontal.
    """
    if components == 'all':
        return [], np.array([0., 0., 1.]), True
    # First determine if horizontal plane is given
    components = components.split(',', maxsplit=1)
    if len(components) > 1:
        planestr = re.findall(r'(\[[0-9\ ,]*\])', components[1].replace(',', ' '))
        if len(planestr) == 1:
            normal = _str_to_vec(planestr[0])
        else:
            u, v = (_str_to_vec(s) for s in planestr)
            normal = np.cross(rlu_to_cart @ u, rlu_to_cart @ v)
            normal = normal / np.sqrt(np.sum(normal**2))
    else:
        normal = np.array([0., 0., 1.])
    components = re.sub(r'[\s;]', '', components[0].lower()).replace('p', 'P').replace('m', 'M')
    rv = []
    for comp in components.replace('+', ';+').replace('-', ';-').split(';'):
        pt = re.match(r'([+-]*[0-9\.]*)[\*]*([PM][xyz]*)', comp).groups()
        prefac = _STRFLOATDIC[pt[0]] if pt[0] in _STRFLOATDIC.keys() else float(pt[0])
        poltyp = PolarisationType(pt[1])
        if not rv:
            typeletter = pt[1][0]
            typelen = len(pt[1][0])
        elif pt[1][0] != typeletter:
            raise RuntimeError(f'Cannot mix polarisation types in component string: {components}')
        elif len(pt[1][0]) != typelen:
            raise RuntimeError(f'Cannot mix full and half polarisation in component string: {components}')
        rv.append((prefac, poltyp))
    return rv, normal, False

def _sum_components(mat: ArrayLike, components: list[tuple[float, PolarisationType]], nMode: int):
    intensity = np.zeros(nMode, dtype=complex)
    if len(components[0][1].value) == 2:
        for comp in components:
            intensity += comp[0] * mat[_ABTOIND[comp[1].value[1]], :].conj()
    else:
        for comp in components:
            intensity += comp[0] * mat[_ABTOIND[comp[1].value[1]], _ABTOIND[comp[1].value[2]], :].conj()
    return intensity

def calculate_polarised_intensity(Sab: ArrayLike,
                                  Sperp: ArrayLike,
                                  q_vectors: ArrayLike,
                                  components: str,
                                  rlu_to_cart: ArrayLike):
    """Calculates the polarised intensity for a given component string from the spin-spin correlation matrix

    :param Sab: the spin-spin correlation matrices computed previously
    :param q_vectors: a numpy array of the q-vectors corresponding to the set of Sab inputted
    :param components: the component string which is a mathematical expression involving polarisation
                       components Pab, Pa or Mab, and optionally the horizontal scattering plane
                       defined after a comma. If not specified, we assume the a-b plane is horizontal.

    Pab and Pa indicate the polarised neutron intensity in either full (Pab) or half (Pa) polarised modes.
    Here ab indicate the x, y, or z directions in the Blume-Maleev cartesian coordinate system where
    x is parallel to Q, y is perpendicular to Q in the horizontal plane and z is vertical.
    Mab will return the projection of the Sab matrix in the Blume-Maleev coordinate system.

    The horizontal plane used to define the Blume-Maleev system can be specified after a comma in the
    component string, either as a perpendicular vector (in the XYZ coordinate system), or two Q-vector
    (in the lattice units system) defining the plane. Note that these vectors defining the plane should
    be delimited by square brackets, and components can be delimited by commas or spaces

    Examples: 'Pxy+0.33*Pyz, [0,1,0]'  # y direction is normal (does not necessarily correspond to b-axis)
              'Mxy, [1 1 0]-[0 0 1]'   # [110]-[001] is defined as horizontal plane, returns projected Sab
              'Px'                     # a-b plane is horizontal (default), return half-polarised intensity
    """
    Sab, Sperp, q_vectors = (np.array(Sab), np.array(Sperp), np.array(q_vectors))
    components, zdir, return_all = parse_components(components, rlu_to_cart)
    # divide Sab into symmetric and anti-symmetric components
    SabA = (Sab - np.transpose(Sab,(1, 0, 2, 3))) / 2
    SabS = (Sab + np.transpose(Sab,(1, 0, 2, 3))) / 2
    nMode = Sab.shape[2]
    intensity = np.empty((nMode, q_vectors.shape[0]), dtype=complex)
    Mabs, intPs, Pabs = ([], [], [])
    for ii, q in enumerate(q_vectors):
        xdir = q @ rlu_to_cart
        xdir = xdir / np.sqrt(np.sum(xdir ** 2))
        ydir = -np.cross(xdir, zdir)
        ydir = ydir / np.sqrt(np.sum(ydir ** 2))
        Pi = np.array([xdir, ydir, zdir]).T
        invP = np.linalg.inv(Pi)
        if return_all or components[0][1].value.startswith('M'):
            # Just a projection of Sab in the Blume-Maleev xyz coordinate system
            Mab = np.einsum("ij,jkl,nk->inl", invP, Sab[:,:,:,ii], invP) # invP @ Sab @ invP'
            if not return_all:
                intensity[:,ii] = _sum_components(Mab, components, nMode)
                continue
            else:
                Mabs.append(Mab)
        # TODO: work out what code in sw_neutron is in relation to B-M equations...
        # TODO: For now just try to copy the code...
        Qdelta = np.tile(np.expand_dims(xdir, (0,1)), (3, 3, 1))
        mat1 = np.tile(np.expand_dims(np.sum(LC * Qdelta, axis=2).T, axis=2), 3) * Qdelta
        mat3 = np.sum(np.sum(np.transpose(np.tile(np.expand_dims(2*mat1 + LC, 3), nMode), (1,2,3,0))
                 * np.tile(np.expand_dims(SabA[:,:,:,ii], 3), 3),0),0).T
        if return_all or len(components[0][1].value) == 2:
            # Half-polarisation calculations
            intP = 1j*np.einsum('ij,ik', Pi, mat3) + np.tile(np.expand_dims(Sperp[ii,:], 1), 3).T
            if not return_all:
                intensity[:,ii] = _sum_components(intP, components, nMode)
                continue
            else:
                intPs.append(intP)
        # Full polarisation calculations
        mat4 = np.transpose(np.tile(np.expand_dims(Pi - np.outer(np.einsum('ij,i', Pi, xdir), xdir).T, 2), 3), (2,0,1))
        mat5 = -mat4 * np.tile(np.expand_dims(np.tile(np.expand_dims(xdir, 1), 3), 2), 3)
        Pab1 = np.transpose(np.sum(np.transpose(np.tile(np.expand_dims(mat4,3), nMode), (0,1,3,2))
                 * np.tile(np.expand_dims(SabS[:,:,:,ii],3), 3), 1), (0,2,1))
        Pab2 = np.transpose(np.tile(np.expand_dims(np.sum(np.sum(
                np.transpose(np.tile(np.expand_dims(mat5,3), nMode), (0,1,3,2))
                 * np.tile(np.expand_dims(SabS[:,:,:,ii],3), 3), 0), 0), 2), 3), (2,1,0)) \
                 * np.tile(np.expand_dims(xdir, (1,2)), (3,nMode))
        Pab3 = np.einsum('ij,k', Pi, Sperp[ii,:])
        Pab4 = np.tile(np.expand_dims(mat3, 1), (1,3,1))
        Pab = np.einsum('ij,jkm', invP, 2*(Pab1 + Pab2) - Pab3 - 1j*Pab4)
        if not return_all:
            intensity[:,ii] = _sum_components(Pab, components, nMode)
        else:
            Pabs.append(Pab)
    if return_all:
        return np.stack(Pabs, -1), np.stack(intPs, -1), np.stack(Mabs, -1)
    return intensity.T
