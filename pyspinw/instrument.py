"""Instrument specification"""

from pyspinw.measurement import Measurement


class Instrument:
    """ Instrument Model"""

    def resolution(self):
        """ Resolution of the instrument"""

    def required_q_points(self, measurement: Measurement):
        """Q points required for applying the resolution to a given measurement"""

    def create_histogram(self, measurement: Measurement):
        """ Do binning """


class PerfectResolution(Instrument):
    """ Instrument with perfect resolution"""

    def create_histogram(self, measurement: Measurement):
        """ Do Binning"""

