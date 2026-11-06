""" Nearest neighbour searching """

import numpy as np

from pyspinw.site import Site
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances

from pyspinw.util import IncrementalPointHistogram

# Strategy:
#  Pick a radius
#  Fill out the search space until every cell that has a corner inside the radius is added
#  Find radius of all the points within that search space
#
# Stopping condition: we have a layer of cells around the central point where no
#  point is closer than the target radius

relative_neighbours = [
    (i, j, k)
    for i in range(-1, 2)
    for j in range(-1, 2)
    for k in range(-1, 2)
    if i != 0 or j != 0 or k != 0
]

class CellData:
    """ Used to contain the data for each cell in the nearest neighbour search """

    def __init__(self, distance_sq):
        self.distance_sq = distance_sq
        self.min_distance_sq = np.min(distance_sq)

    def __repr__(self):
        s = ", ".join([str(x) for x in self.distance_sq])
        return f"({self.min_distance_sq} | {s})"


def add_triple(a: tuple[int, int, int], b: tuple[int, int, int]):
    """ Pairwise addition of the contents of two triples """
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]

def while_true_iter():
    """ Iterator that keeps iterating forever """
    while True:
        yield None

class NeighbourSearch:
    """ Does a nearest neighbour search, and serves as a data storage class for the output """

    def __init__(self, start_site: Site, search_sites: list[Site], unit_cell: UnitCell, n: int=1,
                 shell_tolerance=tolerances.SAME_SITE_ABS_TOL, max_iters: int | None = None):

        if len(search_sites) < 1:
            raise ValueError("Need at least one search site")

        if not isinstance(n, int):
            raise TypeError("Expected n to be an integer")

        if n < 0:
            raise ValueError("Expected n to be positive")

        self.unit_cell = unit_cell

        self.start_site = start_site
        self.search_sites = search_sites

        self.start_in_search_set = start_site.unique_id in set(site.unique_id for site in search_sites)

        if not self.start_in_search_set:
            n -= 1

        self.start_position = unit_cell.lattice_units_to_cartesian([start_site.ijk])
        self.search_positions = unit_cell.lattice_units_to_cartesian(
                                    np.array([site.ijk for site in search_sites], dtype=float))


        self.relative_positions = self.search_positions - self.start_position

        to_calculate = [(0,0,0)]
        found_n_groups = False

        # We want to have at least one point in the search, whatever happens
        current_distance_sq = unit_cell.main_diagonal_length**2

        search_data: dict[tuple[int,int,int], CellData] = {}

        histogram = IncrementalPointHistogram(tolerance=shell_tolerance)

        # The outer loop sets the target distance
        iteration_counter = while_true_iter() if max_iters is None else range(max_iters)
        for _ in iteration_counter:

            # The inner loop iterates until we have found all the cells larger than a certain amount
            while True:

                # print("Searching:", to_calculate)

                # Check all the distances in the parts
                for lu_offset, cartesian_offset in zip(to_calculate,
                                                       unit_cell.lattice_units_to_cartesian(
                                                        np.array(to_calculate).reshape(-1, 3))):
                    differences = self.relative_positions + cartesian_offset
                    search_data[lu_offset] = CellData(np.sum(differences**2, axis=1))

                    # Add the data to the histogram
                    for index, distance_sq in enumerate(search_data[lu_offset].distance_sq):
                        histogram.add(distance_sq, (lu_offset, index))

                # Find new sites to search
                expanded = set()
                for offset, value in search_data.items():
                    if value.min_distance_sq < current_distance_sq:
                        for relative in relative_neighbours:
                            expanded.add(add_triple(offset, relative))

                # Get just the ones we have not searched
                to_calculate = [offset for offset in expanded if offset not in search_data]

                # When all the cells on the outside have all entries above the set distance, break
                if len(to_calculate) == 0:
                    break

            # # Debugging: search data
            # for key in search_data:
            #     print(key, search_data[key].min_distance_sq >= current_distance_sq, search_data[key])

            # # Debugging: plot
            # import matplotlib.pyplot as plt
            # fig = plt.figure(f"Distance {current_distance_sq}")
            # ax = fig.add_subplot(111, projection="3d")
            # inside = [key for key in search_data if search_data[key].min_distance_sq < current_distance_sq]
            # outside = [key for key in search_data if search_data[key].min_distance_sq >= current_distance_sq]
            #
            # ax.scatter(*zip(*inside), color='b')
            # ax.scatter(*zip(*outside), color='r')
            # plt.show()
            #
            # # Debugging: histogram data
            # for index, distance in enumerate(sorted(histogram._groups.keys())):
            #     print(index, distance, histogram._groups[distance])

            # Calculate the new target distance

            groups = sorted(histogram._groups.keys())
            assert len(groups) > 0, "Should always have groups"

            n_groups = len(groups)

            if n_groups <= n:
                # If we don't have n groups yet, redo the search with the greatest search value,
                # this should always increase as we will automatically add new layers

                # print(f"{n_groups}/{n}: {current_distance_sq}")

                current_distance_sq = groups[-1]

            else:
                # Otherwise, the distance we want is for the n_th group
                current_distance_sq = groups[n]

                # We need to redo the search for this n
                if found_n_groups:
                    break

                found_n_groups = True

        else:
            raise RuntimeError("Reached maximum iterations")

        # Deal with the indexing change based on whether the start site is in the search sites
        if not self.start_in_search_set:
            self.groups = [[]] + histogram.groups()[:n+1]

        else:
            self.groups = histogram.groups()[:n+1]


    def neighbour_sites(self):
        """ Get the group information in terms of sites"""
        output = []
        for group in self.groups:
            output.append([(self.search_sites[index], offset) for offset, index in group])

        return output



if __name__ == "__main__":

    if False:
        sites = [
            Site(0.2, 0.3, 0.4),
            Site(0.5, 0.5, 0.5),
            Site(0.2, 0.4, 0.6),
            Site(0.8, 0.7, 0.6)
        ]

        unit_cell = UnitCell(1,0.2,1)

        NeighbourSearch(sites[1], sites, unit_cell, n=6)

    if True:
        unit_cell = UnitCell(1, 1, 1)
        sites = [Site(0.5, 0.5, 0.5)]
        search_results = NeighbourSearch(sites[0], sites, unit_cell, n=10)

        for n, neighbours in enumerate(search_results.neighbour_sites()):
            print(f"{n}:", neighbours)
