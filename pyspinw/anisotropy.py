"""Different kinds of anisotropy."""

import numpy as np
from numpy._typing import ArrayLike

from pyspinw import UnitCell
from pyspinw.structure import Structure
from pyspinw.checks import check_sizes
from pyspinw.serialisation import (
    SPWSerialisationContext,
    SPWSerialisable,
    numpy_serialise,
    numpy_deserialise,
    SPWDeserialisationContext,
)
from pyspinw.site import LatticeSite
from pyspinw.symmetry.operations import SpaceOperation
from pyspinw.tolerances import tolerances


_anisotropy_id_counter = -1


def _generate_unique_anisotropy_id():
    """Generate a unique ID for each anisotropy currently loaded."""
    global _anisotropy_id_counter  # noqa: PLW0603
    _anisotropy_id_counter += 1
    return _anisotropy_id_counter


class Anisotropy(SPWSerialisable):
    """Represent the anisotropy term at a given site.

    Parameters
    ----------
    site
        Lattice site to which the anisotropy term applies.
    anisotropy_matrix
        3x3 matrix defining the anisotropy contribution for the site.
    """

    #: Type name used when serialising anisotropy terms to SPW data.
    serialisation_name = "anisotropy"

    #: Names of scalar fields that can be varied by parameter definitions.
    scalar_parameters = []

    @check_sizes(anisotropy_matrix=(3, 3), force_numpy=True)
    def __init__(self, site: LatticeSite, anisotropy_matrix: ArrayLike):
        self._site = site
        self._anisotropy_matrix = np.array(anisotropy_matrix)
        self._unique_id = _generate_unique_anisotropy_id()

    @property
    def unique_id(self):
        """Unique ID for this anisotropy."""
        return self._unique_id

    @property
    def site(self):
        """Get the site for this anisotropy."""
        return self._site

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """Matrix specifying the anisotropy - `A` term in the Hamiltonian."""
        if self._anisotropy_matrix is None:
            raise ValueError("Anisotropy matrix not initialised - this shouldn't happen")
        else:
            return self._anisotropy_matrix

    def _obeys_symmetry(self, unit_cell: UnitCell, operations: list[SpaceOperation]):
        """ Main part of the check for whether the anisotropy matrix obeys symmetry """
        # We only care if it gives the same quadratic expression for the energy, not if the actual matrices
        # are the same, so, look at the symmetrised version (also multiples by 2, but that doesn't matter)
        m = self.anisotropy_matrix + self.anisotropy_matrix.T

        for operation in operations:
            transform = operation.point_operation_in_cartesian(unit_cell)
            transformed = transform @ m @ transform.T

            if not np.allclose(transformed, m):
                return False

        return True

    def obeys_symmetry(self, structure: Structure):
        """Check that this anisotropy is consistent with the symmetry group"""
        spacegroup = structure.spacegroup
        unit_cell = structure.unit_cell

        # Checking is easier than finding the list of symmetry groups
        operations = spacegroup.operations_between_sites(self._site, self._site)
        operations = [operation for operation in operations if operation.symmorphic] # TODO - does this matter?

        return self._obeys_symmetry(unit_cell, operations)


    def symmetry_copy(self,
                      structure: "Structure",
                      site: LatticeSite):
        """ Copy this anisotropy using symmetry operations """
        # We want to copy the anisotropy under symmetry operations
        # There might be more than one symmetry operation that maps the pair of sites
        #  however, the effect on the anisotropy should be the same for all these operations,
        #  this means we can just pick an arbitrary one.
        # If this turns out not to be the case, then anisotropy itself does not need to obey the
        # symmetry constraints

        spacegroup = structure.spacegroup
        unit_cell = structure.unit_cell

        if self.obeys_symmetry(spacegroup):
            # find the operations that map the sites
            operations = spacegroup.operations_between_sites(self.site, site)

            if len(operations) == 0:
                raise ValueError("New points are not related to the original by symmetry")

            # Pick one element for the transformation
            op = next(iter(operations))
            transform = op.point_operation_in_cartesian(unit_cell)

            new_anisotropy_matrix = transform @ self.anisotropy_matrix @ transform.T

            return Anisotropy(site,
                            anisotropy_matrix=new_anisotropy_matrix,
                            #name=f"{self.name} [{op.text_form}]"
                              )

        else:
            raise ValueError("Anisotropy does not obey symmetry constraints, cannot use symmetry to copy")

    def symmetry_fill(self,
                      structure: "Structure",
                      include_original=False,
                      specialisation_rounding_exponent: int | None = tolerances.SPECIALISE_ROUNDING_EXPONENT):
        """ Make multiple copies of this anisotropy so that symmetry is satisfied """
        # Get the symmetry related sites

        if not self.obeys_symmetry(structure):
            raise ValueError(f"{self} does not obey symmetry constraints of {structure.spacegroup}, "
                             f"cannot use symmetry to copy")

        unit_cell: UnitCell = structure.unit_cell

        related = structure.symmetry_related(self.site)


        new_anisotropies = []
        for site, operations in related:

            # Generate new anisotropy
            for operation in operations:

                # Get the transformed matrix
                op = operation.point_operation_in_cartesian(unit_cell)
                new_matrix = op @ self._anisotropy_matrix @ op.T

                name = self.name + " " + ", ".join([f"({operation.text_form})" for operation in operations])
                new_anisotropies.append(Anisotropy(site,
                                              anisotropy_matrix=new_matrix,
                                              #name = name
                                                   ))

        # Hacky way of excluding the original, add to the list at the start, then remove first
        #  element after the duplicate removal.
        # We also want to include self instead of any copy
        new_anisotropies = [self] + new_anisotropies

        # Now we want to remove copies of the same anisotropies
        for i in range(len(new_anisotropies)): # Basically a while loop, but we've made sure it is finite
            if i >= len(new_anisotropies):
                break

            to_keep: Anisotropy = new_anisotropies[i]
            to_remove: list[int] = []
            for j, to_check in enumerate(new_anisotropies[i+1:]):

                # One anisotropy per site
                if to_keep.site == to_check.site:
                    to_remove.append(j + i + 1)


            # Remove the ones we need to remove, as its ordered we can do
            to_remove.reverse() # in place reverse
            for j in to_remove:
                del new_anisotropies[j]

        # Remove the original input anisotropy if not included, it will always still be at the front
        if not include_original:
            new_anisotropies = new_anisotropies[1:]

        return [specialise_anisotropy(anisotropy) for anisotropy in new_anisotropies]

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"site": self._site._serialise(context), "anisotropy_matrix": numpy_serialise(self._anisotropy_matrix)}

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        site = LatticeSite._deserialise(json["site"], context)
        anisotropy_matrix = numpy_deserialise(json["anisotropy_matrix"])
        return Anisotropy(site, anisotropy_matrix)

    @property
    def parameter_string(self):
        """A string representation of the parameters."""
        m = self.anisotropy_matrix.reshape(-1)
        return f"[[{m[0]}, {m[1]}, {m[2]}], [{m[3]}, {m[4]}, {m[5]}], [{m[6]}, {m[7]}, {m[8]}]]"

    def __repr__(self):
        return f"Anisotropy({self.site.name}, {self.parameter_string})"

    def updated(self, site: LatticeSite | None = None, anisotropy_matrix: ArrayLike | None = None):
        """Return a copy of this anisotropy term with variables replaced.

        Parameters
        ----------
        site
            Replacement lattice site. If omitted, the current site is reused.
        anisotropy_matrix
            Replacement anisotropy matrix. If omitted, the current matrix is reused.
        """
        return Anisotropy(
            site=self.site if site is None else site,
            anisotropy_matrix=self.anisotropy_matrix if anisotropy_matrix is None else np.array(anisotropy_matrix),
        )

    def generalise(self):
        """ Returns the most general version of this class """
        return Anisotropy(site=self.site, anisotropy_matrix=self.anisotropy_matrix)

    def specialise(self):
        """ Convert this to the most specific anisotropy class possible """
        return specialise_anisotropy(self)

    @staticmethod
    def _specialise(anisotropy: "Anisotropy") -> "Anisotropy | None":
        """ Try to get an anisotropy of this class

        As this is the most general case, it will always succeed
        """
        return anisotropy




class AxisMagnitudeAnisotropy(Anisotropy):
    """Represent a uniaxial anisotropy term.

    Parameters
    ----------
    site
        Lattice site to which the anisotropy term applies.
    a
        Anisotropy constant.
    direction
        Principal anisotropy direction. The vector is normalized before use.
    """

    #: Names of scalar fields that can be varied by parameter definitions.
    scalar_parameters = ["a"]

    @check_sizes(direction=(3,), force_numpy=True)
    def __init__(self, site: LatticeSite, a: float, direction: ArrayLike = np.array([0, 0, 1])):
        direction = np.array(direction)
        mag = np.sqrt(np.sum(direction**2))

        if np.isclose(mag, 0, atol=tolerances.VECTOR_TOL):
            self._direction = np.zeros((3,), dtype=float)
        else:
            self._direction = direction / mag

        anisotropy_matrix = a * self._direction.reshape(3, 1) * self._direction.reshape(1, 3)

        super().__init__(site, anisotropy_matrix)
        self._a = a

    @property
    def a(self):
        """Amount of anisotropy (anisotropy constant)."""
        return self._a

    @property
    def constant(self):
        """Amount of anisotropy (anisotropy constant), alias for `a`."""
        return self._a

    @property
    def direction(self):
        """The principal direction of the anisotropy."""
        return self._direction

    @property
    def parameter_string(self):
        """A string representation of the parameters."""
        return f"a={self.constant}, axis={self.direction}"

    def __repr__(self):
        return f"Anisotropy({self.site.name}, {self.parameter_string})"

    def _serialise(self, context: SPWSerialisationContext):
        return {
            "site": self._site._serialise(context),
            "direction": numpy_serialise(self._direction),
            "a": float(self._a),
        }

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return AxisMagnitudeAnisotropy(
            LatticeSite._deserialise(json["site"], context),
            json["a"],
            numpy_deserialise(json["direction"]),
        )

    def updated(self, site: LatticeSite | None = None, a: float | None = None, direction: ArrayLike | None = None):
        """Return a copy of this anisotropy term with variables replaced.

        Parameters
        ----------
        site
            Replacement lattice site. If omitted, the current site is reused.
        a
            Replacement anisotropy constant. If omitted, the current constant is reused.
        direction
            Replacement principal anisotropy direction. If omitted, the current direction is reused.
        """
        return AxisMagnitudeAnisotropy(
            site=self.site if site is None else site,
            a=self.a if a is None else a,
            direction=self.direction if direction is None else np.array(direction),
        )

    @staticmethod
    def _specialise(anisotropy: "Anisotropy", tol=1e-14) -> "Anisotropy | None":
        """ Try to get an anisotropy of this class """

        # If it is already this form, return it
        if isinstance(anisotropy, AxisMagnitudeAnisotropy):
            return anisotropy

        # The form will be
        #   a (V outer V)
        # as the diagonal entries of the matrix are squared,
        # we can check the diagonal components for signs, which need to be consistent

        m = anisotropy.anisotropy_matrix
        m = 0.5*(m + m.T) # Force it to be symmetric

        vec_sq = np.diagonal(m)
        signs = np.sign(vec_sq)

        all_pos = np.all(signs >= -tol)
        all_neg = np.all(signs <= tol)

        if not (all_pos or all_neg):
            # Neither: Cant write as axis-magnitude
            return None

        if all_pos and all_neg:
            # Must be zeros along diagonal, which is either a magnitude zero axis-angle
            return None

        if all_neg:
            m *= -1

        # Check that the matrix is consistent with being V outer V, e.g. (xy)^2 = (x^2) (y^2)

        if not np.isclose(m[0,0]*m[1,1], m[0,1]**2):
            return None

        if not np.isclose(m[0,0] * m[2,2], m[0,2] ** 2):
            return None

        if not np.isclose(m[2,2] * m[1,1], m[2,1] ** 2):
            return None

        # The diagonals should be squared version of the vector, further scaled by a
        # a is actually just x^2 + y^2 + z^2 because the vector has magnitude 1
        a_constant = np.sum(vec_sq)

        # Work out the signs of the anisotropy based on the xy, yz and zx parts,
        # we can only do this up to an overall sign
        #  All positive: all x,y,z positive, or all negative
        #  Two positive: impossible
        #  One positive: one of x,y,z is negative (the one shared by the negatives), or
        #               two are negative (the two in the positive one)
        # Zero positive: impossible
        #
        # If any entry is zero, things are more complicated, it could be positive or negative
        #  but, if one entry is zero, at least one other must also be zero, e.g.
        #   xy = 0 => x = 0 | y = 0 => zx = 0 | yz = 0
        #  we know which ones should be zero though, from the diagonal, and we've checked things match up
        #
        # As we only have a choice of x,y,z up to an overall sign, we can choose the sign of the
        #  first non-zero entry freely, then based on the signs of off-diagonals, we can assign
        #  a corresponding sign to the others

        # Find first non-zero diagonal
        reference = 0
        for i in range(3):
            if not np.isclose(m[i,i], 0):
                reference = i
                break

        vec = np.sqrt(vec_sq / a_constant) # should normalise it

        # Then, for all of them after, apply sign of the cross term
        for i in range(reference+1, 3):
            vec[i] *= np.sign(m[reference, i]) # Could be zero, but that would be fine

        if all_neg:
            a_constant *= -1

        return AxisMagnitudeAnisotropy(anisotropy.site, a=float(a_constant), direction=vec)

_specialisation_search = [AxisMagnitudeAnisotropy, Anisotropy]

def specialise_anisotropy(anisotropy: Anisotropy) -> Anisotropy:
    """ Find the narrowest anisotropy subclass to fit the exchange """

    for An in _specialisation_search:
        specialised = An._specialise(anisotropy)
        if specialised is not None:
            return specialised

    return anisotropy