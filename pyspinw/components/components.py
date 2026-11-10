from typing import overload

import numpy as np


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

class BottomLevelComponent(Component):
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

class PolarisationComponent(Component):

    def __init__(self, a: int, b: int, name: str, z_dir: np.ndarray | None = None):
        self.a = a
        self.b = b
        self.name = name
        self.z_dir = [0,0,1] if z_dir is None else z_dir
        self.z_dir = np.array(z_dir, dtype=float)
        self.z_dir /= np.sqrt(np.sum(z_dir ** 2))

    @overload
    def __call__(
            self,
            s_ab: np.ndarray,
            q_vectors: np.ndarray,
            rlu_to_cart: np.ndarray,
            multimagnon: np.ndarray | None = None
    ): ...

    @overload
    def __call__(
            self,
            z_dir: np.ndarray
    ): ...


    def __call__(self,
                 s_ab: np.ndarray,
                 q_vectors: np.ndarray | None = None,
                 rlu_to_cart: np.ndarray | None = None,
                 multimagnon: np.ndarray | None = None):

        if (q_vectors is None and rlu_to_cart is not None) or (q_vectors is not None and rlu_to_cart is None):
            raise ValueError("Expected both or neither `rlu_to_cart` and `q_vectors` to be None")

        if multimagnon is not None and q_vectors is None:
            raise ValueError("`multimagnon` should be None if `q_vectors` is None")

        if q_vectors is None:
            # We have been given a z direction
            z = np.array(s_ab).reshape(-1)
            if len(z) != 3:
                raise ValueError("Expected z to be a 3 vector")

            return PolarisationComponent(self.a, self.b, self.name, z)

        else:
            # We have been given data to calculate the output for

            # Calculate the coordinate frame

            # x direction is the q vector in cartesian coordinates
            xdir = q_vectors @ rlu_to_cart
            xdir /= np.sqrt(np.sum(xdir ** 2))




# List of components

Sxx = BottomLevelComponent(0, 0, "Sxx")
Sxy = BottomLevelComponent(0, 1, "Sxy")
Sxz = BottomLevelComponent(0, 2, "Sxz")
Syx = BottomLevelComponent(1, 0, "Syx")
Syy = BottomLevelComponent(1, 1, "Syy")
Syz = BottomLevelComponent(1, 2, "Syz")
Szx = BottomLevelComponent(2, 0, "Szx")
Szy = BottomLevelComponent(2, 1, "Szy")
Szz = BottomLevelComponent(2, 2, "Szz")

#Pxx


if __name__ == "__main__":
    import numpy as np
    data = np.array([["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]])

    print(Sxx - Szy)

    print(Sxx(data))
    print(Syz(data))
