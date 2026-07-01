""" Cell offsets (i.e. reference to a particular unit cell in a crystal)"""
import numpy as np

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext



class CellOffset(SPWSerialisable):
    """ Representation of the relative position of individual unit cells within the lattice """

    serialisation_name = "cell_offset"

    def __init__(self, i: int, j: int, k: int):
        self._i = int(i)
        self._j = int(j)
        self._k = int(k)

        self._vector = np.array([i, j, k], dtype=int)

    @property
    def i(self):
        """ First cell offset component"""
        return self._i

    @property
    def j(self):
        """ Second cell offset component"""
        return self._j

    @property
    def k(self):
        """ Third cell offset component"""
        return self._k

    @property
    def as_tuple(self):
        """ As a tuple of ints"""
        return self._i, self._j, self._k

    @property
    def vector(self):
        """ Vector of index (int array)"""
        return self._vector

    def __repr__(self):
        return str(self.as_tuple)

    def position_in_supercell(self, supercell_size: tuple[int, int, int]):
        """ Get the position within a supercell, rather than an infinite crystal (do mods)"""
        si, sj, sk = supercell_size

        return CellOffset(
            i=self.i % si,
            j=self.j % sj,
            k=self.k % sk)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "i": self._i,
            "j": self._j,
            "k": self._k }

    @staticmethod
    def _deserialise(data: dict, context: SPWSerialisationContext):
        return CellOffset(data["i"], data["j"], data["k"])

    @staticmethod
    def coerce(cell_offset_or_data):
        """ Convert a cell offset specification to a CellOffset object

        Input could be a tuple, None or CellOffset
        """
        if isinstance(cell_offset_or_data, CellOffset):
            return cell_offset_or_data


        elif cell_offset_or_data is None:
            return CellOffset(0,0,0)

        elif isinstance(cell_offset_or_data, (tuple, list, np.ndarray)):
            numpyed = np.array(cell_offset_or_data).reshape(-1)
            if numpyed.shape == (3,):
                intified = [int(round(x)) for x in numpyed]
                if np.allclose(cell_offset_or_data, intified):
                    return CellOffset(*intified)
                else:
                    raise ValueError(f"Could not convert {cell_offset_or_data} to cell offset, "
                                     f"it is not approximately an integer list")
            else:
                raise ValueError(f"Could not convert {cell_offset_or_data} to cell offset, "
                                     f"it cannot be interpreted as length 3")



        raise TypeError(f"Could not convert {cell_offset_or_data} to cell offset, should be "
                        f"CellOffset, ArrayLike[Number,Number,Number] or None,"
                        f"got {type(cell_offset_or_data)}")

    def __neg__(self):
        return CellOffset(-self.i, -self.j, -self.k)

CellOffsetCoercible = CellOffset | tuple[int, int, int] | None
