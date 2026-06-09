import numpy as np

def rref(A, tol=1e-12):
    """ Reduced row echelon form - not in scipy or numpy"""

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

_symasym_matrix_inv = np.linalg.inv(_symantisym_matrix)

class PopString:
    def __init__(self, characters: str):
        self.characters = characters

    def pop(self) -> str:
        if len(self.characters) > 0:
            x = self.characters[0]
            self.characters = self.characters[1:]
            return x
        else:
            return ""

class ExchangeMatrixConstraints:
    """ Object representing the constraints on the exchange matrix based on symmetry transformations """
    def __init__(self,
                 transformations: list[np.ndarray],
                 swapped_transformations: list[np.ndarray] | None = None,
                 tol=1e-12):

        swapped_transformations = [] if swapped_transformations is None else swapped_transformations

        if not isinstance(transformations, list):
            raise TypeError("Expected transformations to be a list")

        if not isinstance(swapped_transformations, list):
            raise TypeError("Expected swapped_transformations to be a list")

        self.transformations = transformations
        self.swapped_transformations = swapped_transformations
        self.tol = tol

        unswapped = [(np.kron(transformation, transformation) - np.eye(9)) @ _symasym_matrix_inv.T
                            for transformation in transformations]

        swapped = [(np.kron(transformation, transformation) - _transpose_matrix) @ _symasym_matrix_inv.T
                            for transformation in transformations]

        matrix = np.vstack(unswapped + swapped)

        reduced = rref(matrix)

        # print(reduced)

        # The reduced form should be interpretable as a list of constraints on the matrix
        # It will not have non-zero entries in unless the row directly above has entries in the same column or to the left,
        # or, in other words, if a row has entries starting at position i, then the next row will have entries starting at
        #  position i or greater
        #
        # We can scan the rows,
        #   there will be at most 9 non-zero of them
        #   each will be an equation for one of the entries, and that will be the leftmost
        #   * if a column is "skipped" (i.e. no row starting with it) then it is completely free
        #   * if it has a single non-zero entry, then that entry must be zero
        #   * otherwise, the row will define a relationship between entries that must be obeyed

        self.free = np.ones((9,), dtype=bool)   # Assume free by default
        self.zero = np.zeros((9, ), dtype=bool) # ... and nonzero
        self.constraints = []

        for row in range(9):
            # Which entries in this row are zeros?
            non_zeros = np.abs(reduced[row, :]) > tol
            zeros = ~non_zeros

            if np.all(zeros):
                # The row is all zeros, do nothing
                continue

            # Get first non-zero entry
            col = np.argmax(non_zeros)

            # check remaining entries to the right
            if np.all(zeros[col+1:]):
                # If it's a single entry, then we have the form x_i=0, and we set to zero
                self.zero[col] = True
                self.free[col] = False
            else:
                # Otherwise, we have a more complex constraint, and it counts as non-free
                self.constraints.append(reduced[row, :])
                self.free[non_zeros] = False

        # None of the constraints should contain zero or free entries
        for constraint in self.constraints:
            nonzero = np.abs(constraint) > tol
            assert np.all(~self.zero[nonzero])
            assert np.all(~self.free[nonzero]), (constraint, self.free)

        # print("Free:", self.free)
        # print("Zero:", self.zero)
        # print("Constraints:")
        # for constraint in self.constraints:
        #     print(constraint)

    def _matrix_form_strings(self) -> tuple[list[str], list[str]]:
        """ Gets the strings used to build matrices and constraint equations"""
        strings = [" 0" for _ in range(9)]
        characters = PopString("abcdefghi")
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
                value = str(entry_value)

                eqn += sign + " " + sign + " " + value + chr

            eqn += " = 0"
            eqns.append(eqn)

        return strings, eqns

    def text_summary(self, simplify_output=True, unicode=True):
        """ Textual summary of symmetry allowed exchange properties"""
        m, equations = self._matrix_form_strings()

        if unicode:
            symmetric_output = (
                f"Symmetry allowed symmetric matrices are of the form:\n\n"
                f"  \u23A1{m[0]} {m[1]} {m[2]} \u23A4\n"
                f"  \u23A2{m[1]} {m[4]} {m[3]} \u23A5\n"
                f"  \u23A3{m[2]} {m[3]} {m[5]} \u23A6"
            )
        else:
            symmetric_output = (
                f"Symmetry allowed symmetric matrices are of the form:\n\n"
                f"  /{m[0]} {m[1]} {m[2]} \\\n"
                f"  |{m[1]} {m[4]} {m[3]} |\n"
                f"  \\{m[2]} {m[3]} {m[5]} //"
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

    def print_summary(self, simplify_output=True, unicode=True):
        """ Print out a summary of the exchange properties """
        print(self.text_summary(simplify_output, unicode))

    def __repr__(self):
        return self.text_summary()

if __name__ == "__main__":
    p0 = np.eye(3)
    p1 = np.array([
        [1,0,0],
        [0,0,1],
        [0,1,0]])
    pr1 = np.array([
        [0,0,0],
        [0,1,0],
        [0,0,1]])
    pr2 = np.array([
        [0,0,0],
        [0,0,1],
        [0,1,0]
    ])

    s2 = np.sqrt(1 / 2)

    pr3 = np.array([
        [0,0,0],
        [0, s2, s2],
        [0, -s2, s2]
    ])

    s43 = np.sqrt(4/3)
    pr4 = np.array([
        [0,   0,   0],
        [0, 0.5, s43],
        [0,-s43, 0.5]
    ])

    ExchangeMatrixConstraints([p0]).print_summary()
    ExchangeMatrixConstraints([p1]).print_summary()

    ExchangeMatrixConstraints([p0, p1]).print_summary()

    ExchangeMatrixConstraints([p0, pr1]).print_summary()

    ExchangeMatrixConstraints([pr1, pr2]).print_summary()

    ExchangeMatrixConstraints([pr3]).print_summary()
    ExchangeMatrixConstraints([pr4]).print_summary()

