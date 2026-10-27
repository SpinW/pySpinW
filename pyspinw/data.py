"""Data object - placeholder."""

import numpy as np

class Data:
    """Placeholder.

    Parameters
    ----------
    data : Any
        Data payload.

    Attributes
    ----------
    data : Any
        Data payload.
    """

    def __init__(self, data):
        self.data = data

    @property
    def q(self) -> np.ndarray:
        """Q Values."""
        raise NotImplementedError()
