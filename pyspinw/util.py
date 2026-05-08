"""Generally helpful functions that don't obviously live anywhere else in particular"""
from typing import Callable

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes
from pyspinw.site import LatticeSite
from pyspinw.tolerances import tolerances


@check_sizes(v=(3,), force_numpy=True)
def triple_product_matrix(v: np.ndarray):
    """Find the matrix, M, such that for all vectors X, Y

    X^T M Y = V . (X x Y)
    """
    x, y, z = v

    return np.array([
        [ 0,  z, -y],
        [-z,  0,  x],
        [ y, -x,  0]])


def demo_triple_product_matrix():
    """Show an example of making a matrix that does the triple product"""
    v = [1, 2, 3]
    m = triple_product_matrix(v)
    print(m)

def problematic_sites(sites: list[LatticeSite],
                      implied_sites: list[LatticeSite],
                      site_to_implied_site: list[list[int]]):
    """ Find sites that are contradictory in terms of magnetism

    These sites have an image with same position but negated spin

    :param sites: Asymmetric cell sites
    :param implied_sites: Sites inferred by symmetry
    :param site_to_implied_site: mapping from site index to implied site index
    :returns: Indices of sites with this problem

    """
    bad_sites = []

    for i, site in enumerate(sites):
        for j in site_to_implied_site[i]:
            implied_site = implied_sites[j]
            if np.all(np.abs(site.ijk - implied_site.ijk) < tolerances.SAME_SITE_ABS_TOL):
                # Same position

                if np.all(np.abs(site.m + implied_site.m) < tolerances.SAME_SITE_ABS_TOL):
                    # Opposite spins

                    bad_sites.append(i)
                    break

    return bad_sites

@check_sizes(axis=(3,), force_numpy=True)
def rotation_matrix(angle, axis):
    """ Create a rotation matrix """
    mag = np.sqrt(np.sum(axis**2))

    if mag == 0:
        raise ValueError("Rotation matrix cannot be made for axis (0,0,0) ")

    axis = axis/mag  # Don't use /= because of dtype error possibility

    c = np.cos(angle)
    s = np.sin(angle)

    part1 = (axis.reshape(-1, 1) * axis.reshape(1, -1)) * (1-c)
    part2 = c * np.eye(3)
    part3 = s * triple_product_matrix(-axis)

    return part1 + part2 + part3

def connected_components(adjacency_matrix: np.ndarray) -> list[list[int]]:
    """ Get the connected components of a graph specied by an adjacency matrix

    :param adjacency_matrix: n-by-n numpy array of dtype bool representing the adjacency matrix
    :returns: list of subgraphs, themselves lists of indices for the adjacency matrix
    """
    assert len(adjacency_matrix.shape) == 2
    assert adjacency_matrix.shape[0] == adjacency_matrix.shape[1]

    n = adjacency_matrix.shape[0]
    visited = np.zeros((n, ), dtype=bool)
    components = []

    def search(current_node, component):
        visited[current_node] = True
        component.append(current_node)
        for i in range(n):
            if adjacency_matrix[current_node][i] and not visited[i]:
                search(i, component)

    for i in range(n):
        if not visited[i]:
            component = []
            search(i, component)
            components.append(component)

    return components

def arraylike_equality(array_1: ArrayLike, array_2: ArrayLike, abs_tol=1e-8):
    """ Check for approximate equality of arrays that may or may not have the same shape """
    array_1 = np.array(array_1)
    array_2 = np.array(array_2)

    if array_1.shape != array_2.shape:
        return False

    return np.all(np.abs(array_1 - array_2) < abs_tol)



def uniquetol(values: ArrayLike, tol: float=1e-5):
    """ Returns floating point unique values within a given tolerance """

    p1_norm = np.real(values) + np.imag(values) # This is not (supposed to be?) the absolute value
    order = np.argsort(p1_norm)
    is_bigger_than_tolerance = np.append(True, np.diff(p1_norm[order]) > tol)

    return values[order[is_bigger_than_tolerance]]


def remove_degenerate_and_ghost(energy: ArrayLike, intensity: ArrayLike, tol: float=1e-5, zeroint: int=0, is_series: bool=True):
    """ Removes degenerate and ghost (zero-intensity) modes from spectrum """

    energy, intensity = np.array(energy), np.array(intensity)

    output_energy, output_intensity = energy * np.nan, intensity * np.nan

    # Strategy is to round energies so that they correspond to "bins"
    if tol > 0:
        energy = np.round(energy / tol) * tol

    if zeroint > 0:
        energy[np.where(np.abs(np.real(intensity)) < zeroint)] = np.nan

    # Iterate over q
    for i in range(energy.shape[0]):
        eu = uniquetol(energy[i,:])
        output_energy[i,:len(eu)] = eu
        # Iterate over energy
        for j in range(len(eu)):
            output_intensity[i, j] = np.sum(intensity[i, np.where(np.abs(energy[i,:] - output_energy[i, j]) < tol)])

    nans = np.isnan(output_energy)
    nanmodes = np.sum(nans, axis=0)

    if is_series:
        # Ensure there are the same number of modes throughout
        for col in set([int(idx) for accidentals in np.where((nanmodes > 0) * (nanmodes < output_energy.shape[0]))[0]
                        for idx in np.where(nans[:,accidentals])[0]]):

            idnxt = 1 if col < output_energy.shape[0]-1 else -1
            idnan = np.where(~np.isnan(output_energy[col + idnxt,:]))[0]
            idx = np.array([np.nanargmin(np.abs(output_energy[col + idnxt,i] - output_energy[col])) for i in idnan])
            int_idx, bc = (output_intensity[col, np.where(~np.isnan(output_intensity[col,:]))], np.bincount(idx))
            if int_idx.shape[1] == len(bc) and all(bc > 0):
                output_energy[col,idnan] = np.take(output_energy[col, :], idx)
                output_intensity[col,idnan] = np.take(int_idx / bc, idx)

    # Removes columns which are all NaNs
    output_energy = np.delete(output_energy, np.where(nanmodes == output_energy.shape[0])[0], axis=1)
    output_intensity = np.delete(output_intensity, np.where(nanmodes == output_energy.shape[0])[0], axis=1)

    return output_energy, output_intensity


def energy_grid(energy: ArrayLike, intensity: ArrayLike, output_energies: ArrayLike, energy_fwhm: float | Callable[[float], float]):
    """Bins a set of energy/intensity into a spectrum with Gaussian broadening"""

    energy, intensity = np.real(np.array(energy)), np.array(intensity)

    fwhm = energy_fwhm(output_energies) if isinstance(energy_fwhm, Callable) else np.ones(output_energies.shape) * energy_fwhm
    sigma = fwhm / np.sqrt(8 * np.log(2))  # Convert from FWHM

    output_spectrum = np.zeros((energy.shape[0], len(output_energies)))

    # Get the normalisation prefactor for a Gaussian, based on the numerically obtained output delta E

    output_deltas = output_energies[1:] - output_energies[:-1]
    output_deltas = np.concatenate((output_deltas, output_deltas[-1:]))
    
    gaussian_normalisation_factor = output_deltas / (np.sqrt(2 * np.pi) * sigma)

    # Iterate over "input" energies
    for i in range(energy.shape[1]):
        output_spectrum += intensity[:,i,np.newaxis] * gaussian_normalisation_factor * \
            np.exp(-0.5 * ((output_energies[np.newaxis, :] - energy[:,i,np.newaxis]) / sigma) ** 2)

    return output_spectrum


def rotation_from_z(target_vector):
    """ Rotation matrix from (0,0,1) to the target vector direction """
    mag_sq = np.sum(target_vector**2)

    if mag_sq < 1e-9:
        return np.eye(3)

    v = target_vector / np.sqrt(mag_sq)

    x,y,z = v

    if z < 1e-9 - 1:
        return np.array([[-1,0,0],[0,1,0],[0,0,-1]])

    return np.array([
        [1-x**2 / (1+z), -x*y / (1+z), x],
        [-x*y / (1+z), 1 - y**2 / (1+z), y],
        [-x, -y, z]
    ])


if __name__ == "__main__":
    demo_triple_product_matrix()