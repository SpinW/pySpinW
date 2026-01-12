class Measurement:
    """ Provides the q and E parameters specific to a given measurement """

    def q_points(self):
        pass

    def q_energy_centroids(self):
        pass



class Slice(Measurement):
    pass

class OrthogonalSlice(Slice):
    def __init__(self, path: Path, energies):
        pass