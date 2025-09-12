"""Instrument specification"""

from typing import Callable
from dataclasses import dataclass

from numpy.typing import NDArray

@dataclass
class Instrument:
    """Defines an instrument

    This class is responsible for instrument specific details, including resolution, energy and q binning, etc.
    """

    energy: NDArray
    resolution: Callable[[float], float] = lambda energy: np.max(self.energy)*0.02

    @staticmethod
    def from_ResINS(instrument_name: str) -> 'Instrument':
        """Instantiate an instrument from ResINS.

        Parameters
        ----------
        instrument_name: str
            The name of the instrument.
        """
        raise NotImplementedError("ResINS compatibility has not yet been implemented.")
