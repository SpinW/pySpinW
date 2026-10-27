"""Operations are the members of space groups and magnetic space groups."""

from fractions import Fraction

import numpy as np
from numpy.typing import ArrayLike

from pyspinw.symmetry.symbolic import evaluate as evaluate_text

PointOperationType = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
TranslationType = tuple[Fraction, Fraction, Fraction]

# Used for parsing text
_parsing_zero = {"x": 0, "y": 0, "z": 0}
_parsing_variables = {
    "x": np.array([1,0,0], dtype=float),
    "y": np.array([0,1,0], dtype=float),
    "z": np.array([0,0,1], dtype=float),
}

class Operation:
    """Base class for SpaceOperation and MagneticOperation."""

    def transform_positions(self, positions: ArrayLike):
        """Transform an array of positions with this operation.

        Parameters
        ----------
        positions : array_like
            Array of fractional positions to transform.

        Returns
        -------
        numpy.ndarray
            Transformed fractional positions.

        Raises
        ------
        NotImplementedError
            Always raised in the base class.
        """
        raise NotImplementedError("Transformation of position not implemented in base class")

    def transform_positions_and_spins(self, positions: ArrayLike, spins: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Transform an array of positions and an array of spins with this operation.

        Parameters
        ----------
        positions : array_like
            Array of fractional positions to transform.
        spins : array_like
            Array of spin vectors to transform.

        Returns
        -------
        new_positions : numpy.ndarray
            Transformed fractional positions.
        new_spins : numpy.ndarray
            Transformed spin vectors.

        Raises
        ------
        NotImplementedError
            Always raised in the base class.
        """
        raise NotImplementedError("Transformation of positions and spins not implemented in base class")

class SpaceOperation(Operation):
    """Spacegroup operation.

    Parameters
    ----------
    point_operation : PointOperationType
        3x3 integer matrix representing the point group operation in lattice coordinates.
    translation : TranslationType
        3-tuple of ``Fraction`` values representing the fractional translation in ``[0, 1)``.
    name : str or None, optional
        Optional name or label for the operation. Default is ``None``.

    Raises
    ------
    ValueError
        If ``point_operation`` is not invertible or does not consist of 3x3 entries in
        ``{-1, 0, 1}``, or if ``translation`` does not have 3 values in ``[0, 1)``.
    """

    def __init__(self, point_operation: PointOperationType, translation: TranslationType, name: str | None = None):
        self.name = name

        self.point_operation = point_operation
        self.point_operation_matrix = np.array(point_operation, dtype=float)

        try:
            self.inverse_point_operation_matrix = np.linalg.inv(self.point_operation_matrix)
        except np.linalg.LinAlgError as e:
            raise ValueError(f"Point operation matrix should be invertible, got: {point_operation}") from e

        self.translation = translation
        self.translation_vector = np.array(self.translation, dtype=float)

        self._symmorphic = np.allclose(self.translation, 0.0)

        # Validate data in the point_operation field
        if len(point_operation) != 3:
            raise ValueError("Rotation must have 3 rows")

        for i in range(3):
            col = point_operation[i]
            if len(col) != 3:
                raise ValueError("Rotation columns must have length 3")

            for j in range(3):
                if col[j] not in (-1, 0, 1):
                    raise ValueError("Rotation entries must be -1, 0 or 1")


        # Validator for translation component
        if len(translation) != 3:
            raise ValueError("Translation must have 3 values")

        for i in range(3):
            if translation[i] < 0 or translation[i] >= 1:
                raise ValueError("Translation entries must take values in the half open interval [0, 1),"
                                 f" got translation of {translation}")

    @property
    def symmorphic(self):
        """Whether this operation is symmorphic (i.e. has non-integer translation)."""
        return self._symmorphic

    def __lt__(self, other: "BaseMagneticOperation") -> bool:
        return (self.point_operation, self.translation) < \
            (other.point_operation, other.translation)


    def __eq__(self, other: "BaseMagneticOperation"):
        return (self.point_operation, self.translation) == \
            (other.point_operation, other.translation)

    def __hash__(self):
        return hash((self.point_operation, self.translation))


    @staticmethod
    def _from_numpy(point_operation: np.ndarray, translation: np.ndarray) -> \
                    tuple[PointOperationType, TranslationType]:
        """Convert numpy float data into ints and fractions.

        Parameters
        ----------
        point_operation : numpy.ndarray
            3x3 float array representing the point operation.
        translation : numpy.ndarray
            Length-3 float array representing the translation vector.

        Returns
        -------
        point_operation : PointOperationType
            3x3 tuple of integers.
        translation : TranslationType
            3-tuple of ``Fraction`` objects.
        """
        point_operation = tuple(tuple(int(round(point_operation[i,j])) for i in range(3)) for j in range(3))
        translation = tuple(Fraction(float(translation[i])).limit_denominator() for i in range(3))

        return point_operation, translation


    @property
    def text_form(self) -> str:
        """Representation of this operation as a triple of equations in xyz."""
        strings = []
        for i in range(3):
            string = ""
            first = True
            for j, symbol in enumerate("xyz"):

                match int(self.point_operation[i][j]):
                    case 1:
                        string += symbol if first else " + " + symbol
                        first = False
                    case 0:
                        pass
                    case -1:
                        string += "-" + symbol if first else " - " + symbol
                        first = False
                    case _:
                        raise ValueError(f"Magnetic space group rotation matrices should only contain "
                                         f"-1, 0 or 1, got {self.rotation[i, j]}")

            if self.translation[i] < 0:
                string += f" - {-self.translation[i]}"
            elif self.translation[i] > 0:
                string += f" + {self.translation[i]}"

            strings.append(string)

        return ", ".join(strings)

    @staticmethod
    def from_text(operation_string: str):
        """Parse a string form of an expression, e.g. 'x,x-y,z+1/2'.

        Parameters
        ----------
        operation_string : str
            Comma-separated coordinate transformation string in ``x, y, z``.

        Returns
        -------
        SpaceOperation
            Parsed space operation.

        Raises
        ------
        ValueError
            If ``operation_string`` does not contain three comma-separated expressions.
        """
        parts = [part.strip() for part in operation_string.split(",")]
        if len(parts) != 3:
            raise ValueError(f"Expected three comma separated values, e.g. 'x,y,z', got {operation_string}")

        translation = np.array([evaluate_text(part, _parsing_zero) for part in parts])

        matrix_components = np.array([evaluate_text(part, _parsing_variables) for part in parts])
        matrix_components -= translation.reshape(-1, 1)
        matrix_components = matrix_components.T

        return SpaceOperation.from_numpy(matrix_components, translation, operation_string)


    def __repr__(self):
        return self.text_form


    def and_then(self, other: "SpaceOperation"):
        """Compose this operation with another operation.

        Parameters
        ----------
        other : SpaceOperation
            Operation to apply after this operation.

        Returns
        -------
        SpaceOperation
            Composite operation representing ``self`` followed by ``other``.

        Notes
        -----
        Notationally, it is common in group theory to use implicit multiplication
        for left composition, and dot for right composition elsewhere. Because
        ``__mul__`` would be ambiguous, ``and_then`` makes the order of
        application explicit.
        """
        # Is this going to be slower than converting to matrices and back? Maybe
        new_point_operation = [[0,0,0],[0,0,0],[0,0,0]]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    new_point_operation[i][j] += other.point_operation[i][k] * self.point_operation[k][j]

        # Tupleise
        new_point_operation = tuple(tuple(x) for x in new_point_operation)

        # Translation
        new_translation = tuple(
            (sum(a*b for a, b in zip(point_op_row, self.translation)) + other_trans) % 1
                for point_op_row, other_trans in zip(other.point_operation, other.translation))


        return SpaceOperation(
            point_operation=new_point_operation,
            translation=new_translation)



    @staticmethod
    def from_numpy(point_operation: np.ndarray, translation: np.ndarray, name: str | None = None) -> "SpaceOperation":
        """Create from a description in terms of numpy float array.

        Parameters
        ----------
        point_operation : numpy.ndarray
            3x3 array representing the point operation.
        translation : numpy.ndarray
            Length-3 array representing the translation vector.
        name : str or None, optional
            Optional name for the operation. Default is ``None``.

        Returns
        -------
        SpaceOperation
            Constructed space operation.
        """
        point_operation, translation = \
            SpaceOperation._from_numpy(point_operation, translation)

        return SpaceOperation(point_operation=point_operation,
                              translation=translation,
                              name=name)

    @staticmethod
    def from_transformation_matrix(matrix: np.ndarray, name: str | None = None):
        """Get operation from a 4x4 transformation matrix (currently unused, but useful for debugging).

        Parameters
        ----------
        matrix : numpy.ndarray
            4x4 augmented transformation matrix.
        name : str or None, optional
            Optional name for the operation. Default is ``None``.

        Returns
        -------
        SpaceOperation
            Constructed space operation.
        """
        point_operation = tuple(tuple(int(matrix[j, i]) for i in range(3)) for j in range(3))
        translation = tuple(Fraction(matrix[i, 3]).limit_denominator(100) % 1 for i in range(3))

        return SpaceOperation(
            point_operation=point_operation,
            translation=translation,
            name=name)

    def apply_without_mod(self, points):
        """Apply this operation to a list of points without moving back to the unit cell.

        Parameters
        ----------
        points : array_like
            Array of fractional coordinates of shape ``(N, 3)``.

        Returns
        -------
        numpy.ndarray
            Transformed fractional coordinates without modulo wrapping.
        """
        point_operation = np.array(self.point_operation, dtype=float)
        translation = np.array([float(f) for f in self.translation]).reshape(-1, 1)

        points = np.array(points)

        new_points = (point_operation @ points.T + translation).T

        return new_points

    def point_operation_in_cartesian(self, cell: "UnitCell"):
        """Get the point group in cartesian coordinates.

        Parameters
        ----------
        cell : UnitCell
            Unit cell defining the coordinate transformation.

        Returns
        -------
        numpy.ndarray
            3x3 point operation matrix in Cartesian coordinates.
        """
        return cell._xyz.T @ self.point_operation @ cell._xyz_inv.T


    def transform_positions(self, positions: ArrayLike) -> np.ndarray:
        """Apply this operation to a list of points.

        Parameters
        ----------
        positions : array_like
            Array of fractional coordinates to transform.

        Returns
        -------
        numpy.ndarray
            Transformed fractional coordinates wrapped into ``[0, 1)``.
        """
        positions = np.array(positions)

        point_operation = self.point_operation_matrix
        translation = self.translation_vector.reshape(-1, 1)

        new_positions = (point_operation @ positions.T + translation).T

        return new_positions % 1

    def transform_positions_and_spins(self, positions: ArrayLike, spins: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Apply this operation to a list of points and spins.

        IMPORTANT: This expects spins to be in lattice coordinates (but with preserved magnitude).

        Parameters
        ----------
        positions : array_like
            Array of fractional coordinates to transform.
        spins : array_like
            Array of spin vectors in lattice coordinates.

        Returns
        -------
        new_positions : numpy.ndarray
            Transformed fractional coordinates wrapped into ``[0, 1)``.
        new_spins : numpy.ndarray
            Transformed spin vectors.
        """
        positions = np.array(positions, dtype=float).reshape(-1, 3)
        spins = np.array(spins, dtype=float).reshape(-1, 3)

        point_operation = self.point_operation_matrix
        translation = self.translation_vector.reshape(-1, 1)

        new_positions = (point_operation @ positions.T + translation).T
        new_spins = (point_operation @ spins.T).T

        # print(positions)
        # print(spins)
        # print(new_positions)
        # print(new_spins)

        return new_positions % 1, new_spins



class MagneticOperation(Operation):
    """Magnetic operation.

    Parameters
    ----------
    point_operation : PointOperationType
        3x3 integer matrix representing the point group operation in lattice coordinates.
    translation : TranslationType
        3-tuple of ``Fraction`` values representing the fractional translation in ``[0, 1)``.
    time_reversal : int
        Time reversal factor, either ``1`` or ``-1``.
    name : str or None, optional
        Optional name or label for the operation. Default is ``None``.

    Raises
    ------
    ValueError
        If ``time_reversal`` is not ``1`` or ``-1``, ``point_operation`` does not consist
        of 3x3 entries in ``{-1, 0, 1}``, or ``translation`` does not have 3 values in ``[0, 1)``.
    """

    def __init__(self,
                point_operation: PointOperationType,
                translation: TranslationType,
                time_reversal: int,
                name: str | None = None):

        self.point_operation = point_operation
        self.translation = translation
        self.time_reversal = time_reversal
        self.name = name

        # Validate time reversal field
        if time_reversal not in (1, -1):
            raise ValueError("Time inversion must be either 1 or -1.")


        # Validate point operation field
        if len(point_operation) != 3:
            raise ValueError("Rotation must have 3 rows")

        for i in range(3):
            col = point_operation[i]
            if len(col) != 3:
                raise ValueError("Rotation columns must have length 3")

            for j in range(3):
                if col[j] not in (-1, 0, 1):
                    raise ValueError("Rotation entries must be -1, 0 or 1")

        # Validate translation field
        if len(translation) != 3:
            raise ValueError("Translation must have 3 values")

        for i in range(3):
            if translation[i] < 0 or translation[i] >= 1:
                raise ValueError("Translation entries must take values in the half open interval [0, 1)")


    def __lt__(self, other: "BaseMagneticOperation") -> bool:
        return (self.point_operation, self.translation, self.time_reversal) < \
            (other.point_operation, other.translation, other.time_reversal)


    def __eq__(self, other: "BaseMagneticOperation"):
        return (self.point_operation, self.translation, self.time_reversal) == \
            (other.point_operation, other.translation, other.time_reversal)


    @staticmethod
    def _from_numpy(point_operation: np.ndarray, translation: np.ndarray, time_reversal: np.ndarray) -> \
                    tuple[PointOperationType, TranslationType, int]:
        """Convert numpy float data into ints and fractions.

        Parameters
        ----------
        point_operation : numpy.ndarray
            3x3 float array representing the point operation.
        translation : numpy.ndarray
            Length-3 float array representing the translation vector.
        time_reversal : numpy.ndarray or int or float
            Time reversal factor (``1`` or ``-1``).

        Returns
        -------
        point_operation : PointOperationType
            3x3 tuple of integers.
        translation : TranslationType
            3-tuple of ``Fraction`` objects.
        time_reversal : int
            Time reversal integer (``1`` or ``-1``).
        """
        point_operation = tuple(tuple(int(point_operation[i,j]) for i in range(3)) for j in range(3))
        translation = tuple(Fraction(float(translation[i])).limit_denominator() for i in range(3))
        time_reversal = int(time_reversal)

        return point_operation, translation, time_reversal


    @property
    def text_form(self) -> str:
        """Representation of this generator as a triple of equations in xyz."""
        strings = []
        for i in range(3):
            string = ""
            first = True
            for j, symbol in enumerate("xyz"):

                match int(self.point_operation[i][j]):
                    case 1:
                        string += symbol if first else " + " + symbol
                        first = False
                    case 0:
                        pass
                    case -1:
                        string += "-" + symbol if first else " - " + symbol
                        first = False
                    case _:
                        raise ValueError(f"Magnetic space group rotation matrices should only contain "
                                         f"-1, 0 or 1, got {self.rotation[i, j]}")

            if self.translation[i] < 0:
                string += f" - {-self.translation[i]}"
            elif self.translation[i] > 0:
                string += f" + {self.translation[i]}"

            strings.append(string)


        strings.append(str(self.time_reversal))

        return ", ".join(strings)


    def and_then(self, other: "MagneticOperation"):
        """Compose this operation with another operation.

        Parameters
        ----------
        other : MagneticOperation
            Operation to apply after this operation.

        Returns
        -------
        MagneticOperation
            Composite operation representing ``self`` followed by ``other``.

        Notes
        -----
        Notationally, it is common in group theory to use implicit multiplication
        for left composition, and dot for right composition elsewhere. Because
        ``__mul__`` would be ambiguous, ``and_then`` makes the order of
        application explicit.
        """
        # Is this going to be slower than converting to matrices and back? Maybe
        new_point_operation = [[0,0,0],[0,0,0],[0,0,0]]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    new_point_operation[i][j] += other.point_operation[i][k] * self.point_operation[k][j]

        # Tupleise
        new_point_operation = tuple(tuple(x) for x in new_point_operation)

        # Translation
        new_translation = tuple(
            (sum(a*b for a, b in zip(point_op_row, self.translation)) + other_trans) % 1
                for point_op_row, other_trans in zip(other.point_operation, other.translation))

        new_time_reversal = self.time_reversal * other.time_reversal

        return MagneticOperation(
            point_operation=new_point_operation,
            translation=new_translation,
            time_reversal=new_time_reversal)



    @staticmethod
    def from_numpy(point_operation: np.ndarray, translation: np.ndarray,
                   time_reversal: np.ndarray | float | int, name: str | None = None) -> "MagneticOperation":
        """Create magnetic operation from data as numpy arrays.

        Parameters
        ----------
        point_operation : numpy.ndarray
            3x3 array representing the point operation.
        translation : numpy.ndarray
            Length-3 array representing the translation vector.
        time_reversal : numpy.ndarray or float or int
            Time reversal factor (``1`` or ``-1``).
        name : str or None, optional
            Optional name for the operation. Default is ``None``.

        Returns
        -------
        MagneticOperation
            Constructed magnetic operation.
        """
        point_operation, translation, time_reversal = \
            MagneticOperation._from_numpy(point_operation, translation, time_reversal)

        return MagneticOperation(point_operation=point_operation,
                                 translation=translation,
                                 time_reversal=time_reversal,
                                 name=name)

    def transform_positions(self, positions: ArrayLike):
        """Apply this operation to a list of positions.

        Parameters
        ----------
        positions : array_like
            Array of fractional coordinates to transform.

        Returns
        -------
        numpy.ndarray
            Transformed fractional coordinates wrapped into ``[0, 1)``.
        """
        positions = np.array(positions, dtype=float).reshape(-1, 3)

        point_operation = np.array(self.point_operation, dtype=float)
        translation = np.array([float(f) for f in self.translation]).reshape(-1, 1)

        new_positions = (point_operation @ positions.T + translation).T

        return new_positions % 1


    def transform_positions_and_spins(self, positions: ArrayLike, spins: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        """Apply this operation to a list of positions, and a potentially different sized list of spins.

        Parameters
        ----------
        positions : array_like
            Array of fractional coordinates to transform.
        spins : array_like
            Array of spin vectors to transform.

        Returns
        -------
        new_positions : numpy.ndarray
            Transformed fractional coordinates wrapped into ``[0, 1)``.
        new_spins : numpy.ndarray
            Transformed spin vectors including time reversal.
        """
        positions = np.array(positions, dtype=float).reshape(-1, 3)
        spins = np.array(spins, dtype=float).reshape(-1, 3)

        point_operation = np.array(self.point_operation, dtype=float)
        translation = np.array([float(f) for f in self.translation]).reshape(-1, 1)

        new_positions = (point_operation @ positions.T + translation).T
        new_spins = self.time_reversal * (point_operation @ spins.T).T

        return new_positions % 1, new_spins

    @staticmethod
    def from_text(operation_string: str):
        """Parse a string form of an expression, e.g. 'x,x-y,z+1/2'.

        Parameters
        ----------
        operation_string : str
            Comma-separated coordinate transformation string and time reversal factor.

        Returns
        -------
        MagneticOperation
            Parsed magnetic operation.

        Raises
        ------
        ValueError
            If ``operation_string`` does not contain four comma-separated values.
        """
        parts = [part.strip() for part in operation_string.split(",")]
        if len(parts) != 4:
            raise ValueError(f"Expected three comma separated values, e.g. 'x,y,z', got {operation_string}")

        translation = np.array([evaluate_text(part, _parsing_zero) for part in parts[:3]])

        matrix_components = np.array([evaluate_text(part, _parsing_variables) for part in parts[:3]])
        matrix_components -= translation.reshape(-1, 1)
        matrix_components = matrix_components.T

        time_reversal = int(parts[3])

        return MagneticOperation.from_numpy(matrix_components, translation, time_reversal, operation_string)

    def space_operation(self):
        """Get the corresponding pure space operation for this magnetic operation.

        Returns
        -------
        SpaceOperation
            Pure spatial operation without time reversal.
        """
        return SpaceOperation(self.point_operation, self.translation)

    def __repr__(self):
        return f"MagneticOperation({self.text_form})"

    def __hash__(self):
        return hash((self.point_operation, self.translation, self.time_reversal))


