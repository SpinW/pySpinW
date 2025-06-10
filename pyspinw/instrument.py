"""Instrument specification"""

from dataclasses import dataclass

class Instrument:
    """Defines an instrument

    This class is responsible for instrument specific details, including resolution, energy and q binning, etc.
    """

    resolution: float
