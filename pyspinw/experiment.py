"""Different kinds of experimental setups"""

from pyspinw.sample import Sample
from pyspinw._base import Data

class Experiment:
    """The setup of a neutron experiment."""

    def __init__(self, sample: Sample, instrument: Instrument | None = None):
        self.sample = sample
        self.instrument = instrument

    def calculate(self):
        """Calculate energies for the experiment."""
        raise NotImplementedError

    def fit(self, data: Data):
        """Fit the experimental model to a dataset."""
        raise NotImplementedError
