"""Functions and classes for finding allowed, disallowed and restricted symmetries."""

import numpy as np

from pyspinw.util import rotation_matrix


def rref(A, tol=1e-12):
    """Reduced row echelon form - not in scipy or numpy.

    Parameters
    ----------
    A : numpy.ndarray
        2D matrix to reduce.
    tol : float, optional
        Numerical tolerance below which entries are treated as zero.
        Default is ``1e-12``.

    Returns
    -------
    numpy.ndarray
        Reduced row echelon form of ``A``.
    """
    A = A.astype(float).copy()
    rows, cols = A.shape
    r = 0

    for c in range(cols):
        if r >= rows:
            break

        pivot = np.argmax(np.abs(A[r:, c])) + r

        if abs(A[pivot, c]) < tol:
            continue

        A[[r, pivot]] = A[[pivot, r]]

        A[r] /= A[r, c]

        for i in range(rows):
            if i != r:
                A[i] -= A[i, c] * A[r]

        r += 1

    A[np.abs(A) < tol] = 0
    return A

_transpose_matrix = np.array([
    [ 1, 0, 0, 0, 0, 0, 0, 0, 0 ],
    [ 0, 0, 0, 1, 0, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
    [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
    [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1 ]])

_symantisym_matrix = np.array([
    [1, 0, 0, 0, 0, 0,  0,  0,  0],
    [0, 1, 0, 0, 0, 0,  0,  0,  1],
    [0, 0, 1, 0, 0, 0,  0, -1,  0],
    [0, 1, 0, 0, 0, 0,  0,  0, -1],
    [0, 0, 0, 1, 0, 0,  0,  0,  0],
    [0, 0, 0, 0, 1, 0,  1,  0,  0],
    [0, 0, 1, 0, 0, 0,  0,  1,  0],
    [0, 0, 0, 0, 1, 0, -1,  0,  0],
    [0, 0, 0, 0, 0, 1,  0,  0,  0]
], dtype=float)

_sym_matrix = _symantisym_matrix[:,:6]

_symantisym_matrix_inv = np.linalg.inv(_symantisym_matrix)

def free_and_zero(reduced: np.ndarray, tol: float):
    """Calculate which variables are free and which are zero according to the equation ``reduced = 0``.

    Also creates a list of constraints.

    Parameters
    ----------
    reduced : numpy.ndarray
        Matrix in reduced row echelon form.
    tol : float
        Numerical tolerance for treating matrix entries as zero.

    Returns
    -------
    free : numpy.ndarray
        Boolean array indicating which variables are completely free.
    zero : numpy.ndarray
        Boolean array indicating which variables are constrained to be zero.
    constraints : list of numpy.ndarray
        List of row vectors representing linear constraint equations between
        non-zero, non-free variables.

    Notes
    -----
    The reduced form is interpretable as a list of constraints on the matrix.
    If a row has entries starting at position ``i``, then the next row will have
    entries starting at position ``i`` or greater. Scanning the rows:

    * If a column is skipped (no row starting with it), it is completely free.
    * If a row has a single non-zero entry, that entry must be zero.
    * Otherwise, the row defines a linear relationship between entries that
      must be obeyed.
    """
    size = reduced.shape[1]

    free = np.ones((size,), dtype=bool)  # Assume free by default
    zero = np.zeros((size,), dtype=bool)  # ... and nonzero
    constraints = []

    for row in range(size):
        # Which entries in this row are zeros?
        non_zeros = np.abs(reduced[row, :]) > tol
        zeros = ~non_zeros

        if np.all(zeros):
            # The row is all zeros, do nothing
            continue

        # Get first non-zero entry
        col = np.argmax(non_zeros)

        # check remaining entries to the right
        if np.all(zeros[col + 1:]):
            # If it's a single entry, then we have the form x_i=0, and we set to zero
            zero[col] = True
            free[col] = False
        else:
            # Otherwise, we have a more complex constraint, and it counts as non-free
            constraints.append(reduced[row, :])
            free[non_zeros] = False

    # None of the constraints should contain zero or free entries
    for constraint in constraints:
        nonzero = np.abs(constraint) > tol
        assert np.all(~zero[nonzero])
        assert np.all(~free[nonzero]), (constraint, free)

    return free, zero, constraints

class MatrixConstraints:
    """Base class for constraint data.

    Parameters
    ----------
    n : int
        Number of independent matrix parameters being constrained.
    """

    def __init__(self, n: int):
        self.free = None
        self.zero = None
        self.constraints = None
        self.n = n

    def _matrix_form_strings(self) -> tuple[list[str], list[str]]:
        """Get the strings used to build matrices and constraint equations.

        Returns
        -------
        strings : list of str
            Formatted symbol strings for each of the ``n`` matrix parameters.
        eqns : list of str
            Formatted linear constraint equations relating parameters.
        """
        strings = [" 0" for _ in range(self.n)]
        characters = list(reversed("abcdefghi"))

        for i, is_free in enumerate(self.free):
            if is_free:
                strings[i] = " "+characters.pop()

        eqns = []

        for constraint in self.constraints:

            # Check if the constraint is just +- another one
            nonzero = np.abs(constraint) > self.tol
            nonzero_entries = np.array(np.where(nonzero)).reshape(-1)

            if len(nonzero_entries) == 2:
                first, second = nonzero_entries

                # First value should be 1
                if np.abs(constraint[second] - 1) < self.tol:
                    # second is +1
                    chr = characters.pop()
                    strings[int(first)] = " " + chr
                    strings[int(second)] = "-" + chr
                    continue

                elif np.abs(constraint[second] + 1) < self.tol:
                    # second is -1
                    chr = characters.pop()
                    strings[int(first)] = " "+chr
                    strings[int(second)] = " "+chr
                    continue

            # It's not a simple pair, assign letters to each one
            # and write an equation
            eqn = ""
            for index, position in enumerate(nonzero_entries):
                entry_value = constraint[int(position)]

                chr = characters.pop()
                strings[int(position)] = " " + chr

                if np.abs(entry_value - 1) < self.tol:
                    # +1
                    if index == 0:
                        eqn += chr
                    else:
                        eqn += " + "+chr
                    continue

                elif np.abs(entry_value + 1) < self.tol:
                    # +1
                    if index == 0:
                        eqn += "-" + chr
                    else:
                        eqn += " - " + chr
                    continue

                # Should not be the first entry
                sign = "+" if entry_value > 0 else "-"

                if abs((entry_value - round(entry_value))) < 1e-12:
                    entry_value = int(round(entry_value))

                value = str(abs(entry_value))

                eqn += " " + sign + " " + value + chr

            eqn += " = 0"
            eqns.append(eqn)

        return strings, eqns

    def text_summary(self, simplify_output=True, unicode=True):
        """Summary of matrix constraints as a string (empty in base class).

        Parameters
        ----------
        simplify_output : bool, optional
            Whether to replace all-zero or all-free matrices with a short phrase.
            Default is ``True``.
        unicode : bool, optional
            Whether to use Unicode box-drawing characters for matrices.
            Default is ``True``.

        Returns
        -------
        str
            Formatted summary string.
        """
        return ""

    def print_summary(self, simplify_output=True, unicode=True):
        """Print out a summary of the constraint properties.

        Parameters
        ----------
        simplify_output : bool, optional
            Whether to replace all-zero or all-free matrices with a short phrase.
            Default is ``True``.
        unicode : bool, optional
            Whether to use Unicode box-drawing characters for matrices.
            Default is ``True``.
        """
        print(self.text_summary(simplify_output, unicode))

    def __repr__(self):
        return self.text_summary()




class ExchangeMatrixConstraints(MatrixConstraints):
    """Object representing the constraints on the exchange matrix based on symmetry transformations.

    Parameters
    ----------
    transformations : list of numpy.ndarray
        List of 3x3 point-operation matrices leaving the bond sites unswapped.
    swapped_transformations : list of numpy.ndarray or None, optional
        List of 3x3 point-operation matrices that swap the two bond sites.
        Default is ``None``.
    tol : float, optional
        Numerical tolerance for row reduction and zero checks. Default is ``1e-12``.

    Raises
    ------
    TypeError
        If ``transformations`` or ``swapped_transformations`` is not a list.
    """

    def __init__(self,
                 transformations: list[np.ndarray],
                 swapped_transformations: list[np.ndarray] | None = None,
                 tol=1e-12):

        super().__init__(9)

        swapped_transformations = [] if swapped_transformations is None else swapped_transformations

        if not isinstance(transformations, list):
            raise TypeError("Expected transformations to be a list")

        if not isinstance(swapped_transformations, list):
            raise TypeError("Expected swapped_transformations to be a list")

        self.transformations = transformations
        self.swapped_transformations = swapped_transformations
        self.tol = tol

        unswapped = [(np.kron(transformation, transformation) - np.eye(9)) @ _symantisym_matrix
                     for transformation in transformations]

        swapped = [(np.kron(transformation, transformation) - _transpose_matrix) @ _symantisym_matrix
                   for transformation in swapped_transformations]

        matrix = np.vstack(unswapped + swapped)

        reduced = rref(matrix)

        self.free, self.zero, self.constraints = free_and_zero(reduced, tol)

    def text_summary(self, simplify_output=True, unicode=True):
        """Textual summary of symmetry allowed exchange properties.

        Parameters
        ----------
        simplify_output : bool, optional
            Whether to simplify all-zero or all-free components. Default is ``True``.
        unicode : bool, optional
            Whether to use Unicode matrix brackets. Default is ``True``.

        Returns
        -------
        str
            Formatted textual summary of allowed symmetric and antisymmetric exchange.
        """
        m, equations = self._matrix_form_strings()

        if unicode:
            symmetric_output = (
                f"Symmetry allowed symmetric matrices are of the form:\n\n"
                f"  \u23A1{m[0]} {m[1]} {m[2]} \u23A4\n"
                f"  \u23A2{m[1]} {m[3]} {m[4]} \u23A5\n"
                f"  \u23A3{m[2]} {m[4]} {m[5]} \u23A6"
            )
        else:
            symmetric_output = (
                f"Symmetry allowed symmetric matrices are of the form:\n\n"
                f"  /{m[0]} {m[1]} {m[2]} \\\n"
                f"  |{m[1]} {m[3]} {m[4]} |\n"
                f"  \\{m[2]} {m[4]} {m[5]} /"
            )

        antisymmetric_output = ("Symmetry allowed antisymmetric exchanges have DM vectors of the form:\n\n"
                             f"  [{m[6]} {m[7]} {m[8]} ]")

        # print("Zero:", self.zero)
        # print("Free:", self.free)
        # print("Constraints:", self.constraints)
        if simplify_output:
            if np.all(self.zero[:6]):
                symmetric_output = "No symmetric exchange allowed"

            if np.all(self.free[:6]):
                symmetric_output = "Any symmetric exchange allowed"

            if np.all(self.zero[6:]):
                antisymmetric_output = "No antisymmetric exchange allowed"

            if np.all(self.free[6:]):
                antisymmetric_output = "Any antisymmetric exchange allowed"

        if len(equations) > 0:
            constraints = ["\n".join(["Constraints are:\n"] + ["  " + eqn for eqn in equations])]
        else:
            constraints = []


        return "\n\n".join([symmetric_output, antisymmetric_output] + constraints)


class AnisotropyMatrixConstraints(MatrixConstraints):
    """Object representing the constraints on anisotropy matrices based on symmetry transformations.

    Parameters
    ----------
    transformations : list of numpy.ndarray
        List of 3x3 point-operation matrices leaving the site invariant.
    tol : float, optional
        Numerical tolerance for row reduction and zero checks. Default is ``1e-12``.

    Raises
    ------
    TypeError
        If ``transformations`` is not a list.
    """

    def __init__(self,
                 transformations: list[np.ndarray],
                 tol=1e-12):

        super().__init__(6)

        if not isinstance(transformations, list):
            raise TypeError(f"Expected transformations to be of type list, got {type(transformations)}")

        self.transformations = transformations
        self.tol = tol

        equations = [(np.kron(transformation, transformation) - np.eye(9)) @ _sym_matrix
                     for transformation in transformations]

        matrix = np.vstack(equations)

        reduced = rref(matrix)

        # print(reduced)

        self.free, self.zero, self.constraints = free_and_zero(reduced, tol)


    def text_summary(self, simplify_output=True, unicode=True):
        """Textual summary of symmetry allowed anisotropy properties.

        Parameters
        ----------
        simplify_output : bool, optional
            Whether to simplify all-zero or all-free anisotropy matrices.
            Default is ``True``.
        unicode : bool, optional
            Whether to use Unicode matrix brackets. Default is ``True``.

        Returns
        -------
        str
            Formatted textual summary of allowed anisotropy matrices.
        """
        m, equations = self._matrix_form_strings()

        if unicode:
            main_output = (
                f"Symmetry allowed anisotropy matrices are equivalent to those of the form:\n\n"
                f"  \u23A1{m[0]} {m[1]} {m[2]} \u23A4\n"
                f"  \u23A2{m[1]} {m[3]} {m[4]} \u23A5\n"
                f"  \u23A3{m[2]} {m[4]} {m[5]} \u23A6"
            )
        else:
            main_output = (
                f"Symmetry allowed anisotropy matrices are equivalent to those of the form:\n\n"
                f"  /{m[0]} {m[1]} {m[2]} \\\n"
                f"  |{m[1]} {m[3]} {m[4]} |\n"
                f"  \\{m[2]} {m[4]} {m[5]} /"
            )


        if simplify_output:
            if np.all(self.zero):
                main_output = "No anisotropy allowed"

            if np.all(self.free):
                main_output = "Any anisotropy allowed"

        if len(equations) > 0:
            constraints = ["\n".join(["Constraints are:\n"] + ["  " + eqn for eqn in equations])]
        else:
            constraints = []


        return "\n\n".join([main_output] + constraints)


if __name__ == "__main__":
    transformation = rotation_matrix(2*np.pi/3, (0,1,0))

    AnisotropyMatrixConstraints([transformation]).print_summary()




