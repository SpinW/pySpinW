""" Load the data from crysFML"""

from importlib import resources
from dataclasses import dataclass
from multiprocessing.managers import Value

import numpy as np
import re


@dataclass
class PointOperation:
    number: int
    name: str
    string_form: str
    matrix: np.ndarray


def parse_point_operation(line: str) -> PointOperation:
    parts = [s for s in line.split(" ") if s != ""]

    try:
        n = int(parts[0])

        point_op_label = parts[1]
        point_op_string = parts[2]
        point_op_matrix = np.array([int(x) for x in parts[3:12]]).reshape(3, 3)

        return PointOperation(
            number=n,
            name=point_op_label,
            string_form=point_op_string,
            matrix=point_op_matrix)

    except ValueError as ex:
        raise ValueError(f"Failed to parse '{line}'") from ex

def split_stringy_line(line: str) -> list[str]:
    """ Split a line containing strings """
    # Lots of ways to do this, this is simple
    # Find all matching a number with whitespace, or a string

    matches = re.findall(r'(?:\d+|"(?:[^"\\]|\\.)*")', line)

    return matches


point_operations = {}
hexagonal_point_operations = {}


with resources.open_text("pyspinw.data", "magnetic_data.txt") as file:


    # Read the square operations
    for i in range(48):
        line = file.readline().strip()
        op = parse_point_operation(line)
        point_operations[op.number] = op

    # Read the hexagonal operations
    for i in range(24):
        line = file.readline().strip()
        op = parse_point_operation(line)
        hexagonal_point_operations[op.number] = op

    for group_id in range(1651):
        print(group_id)

        # First line, names and numbers
        line = file.readline().strip()

        parts = split_stringy_line(line)

        bns_number_part_1 = int(parts[0])
        bns_number_part_2 = int(parts[1])
        bns_number_string = parts[2][1:-1]

        uni_label = parts[3][1:-1]
        bns_label = parts[4][1:-1]

        og_number_part_1 = int(parts[5])
        og_number_part_2 = int(parts[6])
        og_number_part_3 = int(parts[7])

        og_number_string = parts[8][1:-1]
        og_label = parts[9][1:-1]

        # Second line, the type
        line = file.readline().strip()
        group_type = int(line)

        # Now we have group type dependent parsing, type 4 is exceptional
        if group_type == 4:

            # Type 4 parsing
            line = file.readline().strip()
            parts = [s for s in line.split(" ") if s != ""]

            bnsog_point_op = np.array([int(x) for x in parts[:9]]).reshape(3, 3)
            bnsog_origin_num = np.array([int(x) for x in parts[9:12]])
            bnsog_origin_denom = int(parts[12])

        #
        # Next section, operations
        #

        # Operations 1, number of operations
        line = file.readline().strip()
        n_operations = int(line)

        # Operations 2, operations details, multiple operations in one line


        bns_point_op = []
        bns_translation_num = []
        bns_translation_denom = []
        bns_time_inversion = []

        print(bns_label, len(parts), n_operations)

        operation_lines = n_operations // 4
        operation_sizes = [4]*operation_lines + [n_operations % 4]
        operation_sizes = [s for s in operation_sizes if s > 0]

        for line_length in operation_sizes:
            line = file.readline().strip()
            parts = [s for s in line.split(" ") if s != ""]
            for i in range(line_length):
                bns_point_op.append(int(parts[6*i])) # 6n + 0
                bns_translation_num.append([int(parts[6*i + j + 1]) for j in range(3)]) # 6n + 1,2,3
                bns_translation_denom.append(int(parts[6*i+4])) # 6n + 4
                bns_time_inversion.append(int(parts[6*i+5])) # 6n + 5

        # Lattice Vectors 1, number of lattice vectors
        line = file.readline().strip()
        n_lattice_vectors = int(line)

        # Lattice Vectors 2, actual lattice vectors
        line = file.readline().strip()
        parts = [s for s in line.split(" ") if s != ""]

        lattice_vectors_num = []
        lattice_vectors_denom = []

        for i in range(n_lattice_vectors):
            lattice_vectors_num.append([int(parts[4*i + j]) for j in range(3)])
            lattice_vectors_denom.append(int(parts[4*i + 3]))

        # Wyckoff sites 1, number of sites
        line = file.readline().strip()
        n_wyckoff = int(line)

        # Wyckoff loop
        wyckoff_sites = []
        for i in range(n_wyckoff):

            # Wyckoff outer loop
            line = file.readline().strip()
            parts = split_stringy_line(line)

            n_positions = int(parts[0])
            multiplicity = int(parts[1])
            label = parts[2].strip()

            # Wyckoff inner loop
            positions_num = []
            positions_denom = []
            positions_xyz = []
            positions_mag = []
            for j in range(n_positions):
                line = file.readline().strip()
                parts = [s for s in line.split(" ") if s != ""]

                positions_num.append([int(parts[k]) for k in range(3)]) # 10n + 0,1,2
                positions_denom.append(int(parts[3])) # 10n + 3
                positions_xyz.append([int(parts[k + 4]) for k in range(3)]) # 10n + 4,5,6
                positions_mag.append([int(parts[k + 7]) for k in range(3)]) # 10n + 7,8,9

            wyckoff_sites.append((positions_num, positions_denom, positions_xyz, positions_mag))

        # Last part, more type 4 dependent stuff, OG representations

        if group_type == 4:
            # OG line 1, number of sites
            line = file.readline().strip()
            n_og_operations = int(line)

            operation_lines = n_og_operations // 4
            operation_sizes = [4] * operation_lines + [n_og_operations % 4]
            operation_sizes = [s for s in operation_sizes if s > 0]

            # OG operators line

            og_point_op = []
            og_translation_num = []
            og_translation_denom = []
            og_time_inversion = []

            for line_length in operation_sizes:
                line = file.readline().strip()
                parts = [s for s in line.split(" ") if s != ""]

                for i in range(line_length):

                    og_point_op.append(int(parts[6 * i]))  # 6n + 0
                    og_translation_num.append([int(parts[6 * i + j + 1]) for j in range(3)])  # 6n + 1,2,3
                    og_translation_denom.append(int(parts[6 * i + 4]))  # 6n + 4
                    og_time_inversion.append(int(parts[6 * i + 5]))  # 6n + 5

            # OG Lattice Vectors 1, number of lattice vectors
            line = file.readline().strip()
            n_og_lattice_vectors = int(line)

            # OG Lattice Vectors, actual lattice vectors
            line = file.readline().strip()
            parts = [s for s in line.split(" ") if s != ""]

            og_lattice_vectors_num = []
            og_lattice_vectors_denom = []

            for i in range(n_og_lattice_vectors):
                og_lattice_vectors_num.append([int(parts[4 * i + j]) for j in range(3)])
                og_lattice_vectors_denom.append(int(parts[4 * i + 3]))

            # Wyckoff sites 1, number of sites
            line = file.readline().strip()
            n_og_wyckoff = int(line)

            # Wyckoff loop
            og_wyckoff_sites = []
            for i in range(n_og_wyckoff):

                # Wyckoff outer loop
                line = file.readline().strip()
                parts = split_stringy_line(line)

                n_positions = int(parts[0])
                multiplicity = int(parts[1])
                label = parts[2].strip()

                # Wyckoff inner loop
                positions_num = []
                positions_denom = []
                positions_xyz = []
                positions_mag = []
                for j in range(n_positions):
                    line = file.readline().strip()
                    parts = [s for s in line.split(" ") if s != ""]

                    positions_num.append([int(parts[k]) for k in range(3)])  # 10n + 0,1,2
                    positions_denom.append(int(parts[3]))  # 10n + 3
                    positions_xyz.append([int(parts[k + 4]) for k in range(3)])  # 10n + 4,5,6
                    positions_mag.append([int(parts[k + 7]) for k in range(3)])  # 10n + 7,8,9

                og_wyckoff_sites.append((positions_num, positions_denom, positions_xyz, positions_mag))
