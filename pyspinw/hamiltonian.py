"""Selection of different Hamiltonians"""

from pyspinw._base import Hamiltonian

# pylint: disable=R0903


class GeneralHamiltonian(Hamiltonian):
    """Hamiltonian of the form: H = Σ s.P.s + Q.s

    More explicitly, the equation is ... TODO

    """


class HeisenbergMagnet(GeneralHamiltonian):
    """Heisenberg Ferromagnet Hamiltonian of the form H = - Σ J s.s

    More explicitly, the equation is ... TODO
    """


class BiQuadratic(Hamiltonian):
    """Biquadratic Hamiltonian of the form H = - Σ A (s.s)^2

    More explicitly, the equation is ... TODO
    """
