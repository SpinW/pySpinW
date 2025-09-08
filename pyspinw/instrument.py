"""Instrument specification"""

from dataclasses import dataclass

class Instrument:
    """Defines an instrument

    This class is responsible for instrument specific details, including resolution, energy and q binning, etc.
    """

    resolution: float

    @classmethod
    def from_ResINS(instrument_name: str) -> Instrument:
        """Instantiate an instrument from ResINS.

        Parameters
        ----------
        instrument_name: str
            The name of the instrument.
        """
        raise NotImplementedError("ResINS compatibility has not yet been implemented.")
