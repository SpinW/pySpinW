from typing import overload

import numpy as np
from numpy._typing import ArrayLike


class Component:

    def __add__(self, other: "Component"):
        return Sum(self, other)

    def __neg__(self):
        return Negation(self)

    def __sub__(self, other: "Component | float | int"):
        return self.__add__(other.__neg__())

    def __mul__(self, other: "Component | float | int"):
        if isinstance(other, (float, int)):
            return Product(self, Constant(other))
        return Product(self, other)

    def __rmul__(self, other: "Component | float | int"):
        if isinstance(other, (float, int)):
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
    def __init__(self, left: Component, right: Component):
        self.left = left
        self.right = right

    def __call__(self, *data):
        return self.left(*data) + self.right(*data)

    def __repr__(self):
        return f"({self.left} + {self.right})"

class Constant(Component):
    def __init__(self, value):
        self.value = value

    def __call__(self, *data):
        return self.value

    def __repr__(self):
        return repr(self.value)

class SabComponent(Component):
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

            return PMComponent(self.a, self.b, self.name, z)

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

class MComponent(PMComponent):
    def __init__(self, name: str, z_dir: ArrayLike | None = None):
        super().__init__(a=0, b=0, name=name, z_dir=z_dir) # a and b are not used

    def _calculation(self, s_ab, to_bm, from_bm):
        # The calculation is to_bm s_ab to_bm^T = P S P^T
        # Lets write this in einstein notation:

        # #  M_[il,mn] = P_[ij,n] S_[jk,mn] P_[kl,n]
        # return np.einsum('ijn, jkmn, kln -> ilm n', to_bm, s_ab, to_bm)

        #  M_[il,mn] = P_[ij,n] S_[jk,mn] P_[lk,n]
        return np.einsum('ijn, jkmn, lkn -> ilm n', to_bm, s_ab, to_bm)

class PComponent():
    pass

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

#Pxx = PMComponent(1, 1, "Pxx")

Mab = MComponent("M")


if __name__ == "__main__":
    import numpy as np

    #print(Pxx[1,2,3])

    data = np.array([["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]])

    print(Sxx - Szy)

    print(Sxx(data))
    print(Syz(data))
