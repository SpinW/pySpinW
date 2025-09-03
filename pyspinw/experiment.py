"""Different kinds of experimental setups"""

from pyspinw.sample import Sample
from pyspinw._base import Data

import numpy as np

from instrument import Instrument
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext


class Experiment(SPWSerialisable):
    """The setup of a neutron experiment."""

    def __init__(self, sample: Sample, instrument: Instrument | None = None):
        self.sample = sample
        self.instrument = instrument

    def calculate(self, input_q: list[np.ndarray], n_q: int, field: None = None) -> np.ndarray:
        """Calculate energies for the experiment."""
        generated_q = self.sample.generate_q(input_q, n_q, field, self.instrument.resolution)
        energies = self.sample.hamiltonian.energies(generated_q)

        return energies

    def fit(self, data: Data):
        """Fit the experimental model to a dataset."""
        q_values = data.q
        raise NotImplementedError()

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"sample": self.sample.serialise(),
                "instrument": self.instrument.serialise()}

    @staticmethod
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        return Experiment(
            Sample.deserialise(data["sample"]),
            Instrument.deserialise(data["instrument"])
        )