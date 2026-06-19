"""Define experiment-level coordination between samples, instruments, and data.

This module is not fully implemented yet.
"""

from pyspinw.sample import Sample
from pyspinw.data import Data

import numpy as np

from instrument import Instrument
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext


class Experiment(SPWSerialisable):
    """Represent a neutron scattering experiment.

    An experiment combines a physical sample with an optional instrument model.
    The sample provides the magnetic model to evaluate, while the instrument
    provides resolution information for the simulated measurement.
    """

    serialisation_name = "experiment"

    def __init__(self, sample: Sample, instrument: Instrument | None = None):
        """Create an experiment from a sample and optional instrument.

        Parameters
        ----------
        sample : Sample
            Sample model that owns the Hamiltonian and q-point generation logic.
        instrument : Instrument | None, optional
            Instrument model used to provide resolution information. If omitted,
            methods that require instrument resolution need to handle the missing
            instrument before calling this object.
        """
        self.sample = sample
        self.instrument = instrument

    def calculate(self, input_q: list[np.ndarray], n_q: int, field: None = None) -> np.ndarray:
        """Calculate simulated energies for requested momentum inputs.

        Parameters
        ----------
        input_q : list[np.ndarray]
            Momentum-space inputs passed to the sample for q-point generation.
        n_q : int
            Number of q-points to generate from the requested momentum inputs.
        field : None, optional
            External field passed through to the sample q-point generator.

        Returns
        -------
        np.ndarray
            Energies and intensities calculated from the generated q-points.
        """
        generated_q = self.sample.generate_q(input_q, n_q, field, self.instrument.resolution)
        energies = self.sample.hamiltonian.energies_and_intensities(generated_q)

        return energies

    def fit(self, data: Data):
        """Fitting has not been implemented yet."""
        q_values = data.q
        raise NotImplementedError()

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"sample": self.sample.serialise(), "instrument": self.instrument.serialise()}

    @staticmethod
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        return Experiment(Sample.deserialise(data["sample"]), Instrument.deserialise(data["instrument"]))
