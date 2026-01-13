""" Description of measurements """

class Measurement:
    """ Provides the q and E parameters specific to a given measurement """

    def q_points(self):
        """ Q positions """

    def q_energy_centroids(self):
        """ Centres of the bins in q-Energy space"""



class Slice(Measurement):
    """ A slice in 4D space"""

class OrthogonalSlice(Slice):
    """ A slice in 4D space with energy and q independent """

    def __init__(self, path: Path, energies):
        pass
