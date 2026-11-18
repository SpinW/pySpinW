from typing import overload

import numpy as np
from numpy._typing import ArrayLike


class Component:
    """ Base class for the representation of scattering intensities / polarisation """

    def __add__(self, other: "Component"):
        return Sum(self, other)

    def __neg__(self):
        return Negation(self)

    def __sub__(self, other: "Component | float | int"):
        return self.__add__(other.__neg__())

    def __mul__(self, other: "Component | float | int"):
        if isinstance(other, (complex, float, int)):
            return Product(self, Constant(other))
        return Product(self, other)

    def __rmul__(self, other: "Component | float | int"):
        if isinstance(other, (complex, float, int)):
            return Product(Constant(other), self)
        return Product(other, self)

    def __truediv__(self, other: "Component | float | int"):
        if isinstance(other, (float, int)):
            return Division(self, Constant(other))
        return Division(other, self)

    def __rtruediv__(self, other: "Component | float | int"):
        if isinstance(other, (float, int)):
            return Division(Constant(other), self)
        return Division(other, self)

    def __call__(self,
                 s_ab: np.ndarray,
                 q_vectors: np.ndarray,
                 rlu_to_cart: np.ndarray,
                 multimagnon: np.ndarray | None = None):
        raise NotImplementedError("Cannot call base component")

    def __repr__(self):
        return "Component Base Class"

    @property
    def allows_negative(self) -> bool:
        raise NotImplementedError("allows_negative is not implemented")

class Division(Component):
    """ Component division"""
    def __init__(self, numerator: Component, denominator: Component):
        self.numerator = numerator
        self.denominator = denominator

    def __call__(self, *data):
        return self.numerator(*data) / self.denominator(*data)

    def __repr__(self):
        return f"({self.numerator}/{self.denominator})"

    @property
    def allows_negative(self) -> bool:
        return self.numerator.allows_negative or self.denominator.allows_negative

class Product(Component):
    """ Component product """
    def __init__(self, left: Component, right: Component):
        self.left = left
        self.right = right

    def __call__(self, *data):
        return self.left(*data) * self.right(*data)

    def __repr__(self):
        return f"({self.left} * {self.right})"


    @property
    def allows_negative(self) -> bool:
        return self.left.allows_negative or self.right.allows_negative

class Negation(Component):
    """ Component negation """
    def __init__(self, parent: Component):
        self.parent = parent

    def __call__(self, *data):
        return -self.parent(*data)

    def __repr__(self):
        return f"-{self.parent}"

    @property
    def allows_negative(self) -> bool:
        return True

class Sum(Component):
    """ Component addition """
    def __init__(self, left: Component, right: Component):
        self.left = left
        self.right = right

    def __call__(self, *data):
        return self.left(*data) + self.right(*data)

    def __repr__(self):
        return f"({self.left} + {self.right})"

class Constant(Component):
    """ Constant to be used in algebra"""
    def __init__(self, value: int | float | complex):
        self.value = value

    def __call__(self, *data):
        return self.value

    def __repr__(self):
        return repr(self.value)

class SabComponent(Component):
    """ Component of the scattering matrix Sab"""
    def __init__(self, a: int, b: int, name: str):
        self.a = a
        self.b = b
        self.name = name

    def __call__(self,
                 s_ab: np.ndarray,
                 q_vectors: np.ndarray,
                 rlu_to_cart: np.ndarray,
                 multimagnon: np.ndarray | None = None):

        return s_ab[self.a, self.b, :, :]

    def __repr__(self):
        return self.name

class PMComponent(Component):
    """ Base class of P and M matrices and components """
    def __init__(self, a: int, b: int, name: str, z_dir: ArrayLike | None = None):
        self.a = a
        self.b = b
        self.name = name

        self.z_dir_input = z_dir
        self.z_dir = [0,0,1] if z_dir is None else z_dir
        self.z_dir = np.array(self.z_dir, dtype=float).reshape(-1)

        if len(self.z_dir) != 3:
            raise ValueError("Expected z_dir to be length 3")

        self.z_dir /= np.sqrt(np.sum(self.z_dir ** 2))

    def __getitem__(self, item):
        if not isinstance(item, tuple):
            raise TypeError("Expected __getitem__ to have a tuple input")

        if not len(item) == 3:
            raise ValueError("Expected __getitem__ to be called with three components")

        return PMComponent(self.a, self.b, self.name, item)

    def __call__(self,
                 s_ab: np.ndarray,
                 q_vectors: np.ndarray,
                 rlu_to_cart: np.ndarray,
                 multimagnon: np.ndarray | None = None):

        if q_vectors is None:
            # We have been given a z direction
            z = np.array(s_ab).reshape(-1)
            if len(z) != 3:
                raise ValueError("Expected z to be a 3 vector")

            return self.__class__(a=self.a, b=self.b, name=self.name, z_dir=z)

        else:
            # We have been given data to calculate the output for

            # Calculate the coordinate frame

            # x direction is the q vector in cartesian coordinates
            z_dir = self.z_dir

            x_dir = q_vectors @ rlu_to_cart
            x_dir /= np.sqrt(np.sum(x_dir ** 2, axis=1)).reshape(-1, 1)

            y_dir = np.cross(z_dir, x_dir) # Right handed
            y_dir /= np.sqrt(np.sum(y_dir ** 2, axis=1)).reshape(-1, 1)

            z_dir = np.tile(z_dir.reshape(1, 3), (q_vectors.shape[0], 1))

            #
            # Blume–Maleev frame
            #

            from_bm = np.array([x_dir, y_dir, z_dir])

            # Move q index from middle to front for inverse
            from_bm = from_bm.transpose(1, 0, 2)

            to_bm = np.linalg.inv(from_bm)

            # Move q index back to the end again
            from_bm = from_bm.transpose(1,2,0)
            to_bm = to_bm.transpose(1,2,0)

            return self._calculation(s_ab, to_bm, from_bm)

    def _calculation(self, s_ab, to_bm, from_bm):
        """ Do the specific calculation for this class """
        raise NotImplementedError("_calculation not implemented")

    def __repr__(self):
        if self.z_dir_input is None:
            return self.name
        else:
            return self.name + f"[{self.z_dir_input[0]},{self.z_dir_input[1]},{self.z_dir_input[2]}]"

def m_calculation(s_ab, to_bm):
    # The calculation is to_bm s_ab to_bm^T = P S P^T
    # Lets write this in einstein notation:

    #  M_[il,mn] = P_[ij,n] S_[jk,mn] P_[lk,n]
    return np.einsum('ijn, jkmn, lkn -> lim n', to_bm, s_ab, to_bm)


class MAll(PMComponent):
    """ Mab matrix"""
    def __init__(self, z_dir: ArrayLike | None = None, **kwargs):
        super().__init__(a=0, b=0, name="M", z_dir=z_dir) # a and b are not used

    def _calculation(self, s_ab, to_bm, from_bm):
        return m_calculation(s_ab, to_bm)

class MSPerp(PMComponent):
    """ Sperp """
    def __init__(self, z_dir: ArrayLike | None = None, **kwargs):
        super().__init__(a=0, b=0, name="Sperp", z_dir=z_dir) # a and b are not used

    def _calculation(self, s_ab, to_bm, from_bm):
        m = m_calculation(s_ab, to_bm)
        return m[1,1, :, :] + m[2,2, :, :]


class MComponent(PMComponent):
    """ Component of Mab matrix"""
    def _calculation(self, s_ab, to_bm, from_bm):
        m = m_calculation(s_ab, to_bm)

        return m[self.a, self.b, :, :]

class PAll(PMComponent):
    """ Pab matrix"""

    def __init__(self, z_dir: ArrayLike | None = None, **kwargs):
        super().__init__(a=0, b=0, name="P", z_dir=z_dir) # a and b are not used

class PComponentHalf(PMComponent):
    """ Half polarisation components """

    def __init__(self, a: int, name: str, z_dir: ArrayLike | None = None, **kwargs):
        super().__init__(a=a, b=0, name=name, z_dir=z_dir) # b is not used

class PComponentFull(PMComponent):
    """ Full polarisation components """

# List of components

Sxx = SabComponent(0, 0, "Sxx")
Sxy = SabComponent(0, 1, "Sxy")
Sxz = SabComponent(0, 2, "Sxz")
Syx = SabComponent(1, 0, "Syx")
Syy = SabComponent(1, 1, "Syy")
Syz = SabComponent(1, 2, "Syz")
Szx = SabComponent(2, 0, "Szx")
Szy = SabComponent(2, 1, "Szy")
Szz = SabComponent(2, 2, "Szz")

Sperp = MSPerp()

#Pxx = PMComponent(1, 1, "Pxx")

Mab = MAll()
Mxx = MComponent(0,0,"Mxx")
Mxy = MComponent(0,1,"Mxy")
Mxz = MComponent(0,2,"Mxz")
Myx = MComponent(1,0,"Myx")
Myy = MComponent(1,1,"Myy")
Myz = MComponent(1,2,"Myz")
Mzx = MComponent(2,0,"Mzx")
Mzy = MComponent(2,1,"Mzy")
Mzz = MComponent(2,2,"Mzz")


if __name__ == "__main__":
    import numpy as np

    #print(Pxx[1,2,3])

    data = np.array([["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]])

    print(Sxx - Szy)

    print(Sxx(data))
    print(Syz(data))
