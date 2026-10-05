""" Magnetic Spacegroups"""


from pyspinw.serialisation import SPWSerialisable
from pyspinw.symmetry.group import SymmetryGroup, database, ExactMatch
from pyspinw.symmetry.operations import MagneticOperation


class MagneticSpaceGroup(SymmetryGroup, SPWSerialisable):
    """ Representation of a magnetic space group"""

    serialisation_name = "magnetic_group"

    @property
    def space_operations(self):
        """ The SpaceOperations for this group """
        return self.parent_spacegroup.operations


    def __init__(self, operations: list[MagneticOperation]):
        self.operations = operations
        self.parent_spacegroup = self.get_spacegroup(operations)


    def __repr__(self):
        """repr"""
        return f"MagneticSpaceGroup({self.number}, {self.symbol})"

    @staticmethod
    def get_spacegroup(operations):
        """ Get the spacegroup that corresponds to this magnetic group """
        space_operations = [operation.space_operation().text_form for operation in operations]

        match_data = database.spacegroups_with_operations("; ".join(space_operations))

        if isinstance(match_data, ExactMatch):
            return match_data.spacegroup

        else:
            raise ValueError("Failed to find spacegroup exactly matching this magnetic group")

    @staticmethod
    def closure(operations: list[MagneticOperation], centerings: list[MagneticOperation]):
        """ Find the closure of a list of operations, based on the operation/centering convention """
        all_operations = set(centering.and_then(op)
                             for centering in centerings
                              for op in operations)

        last_operations_count = 0
        while last_operations_count < len(all_operations):
            last_operations_count = len(all_operations)

            new_operations = set(
                                op1.and_then(op2)
                                for op1 in all_operations
                                for op2 in all_operations)

            all_operations = all_operations.union(new_operations)

        return MagneticSpaceGroup(list(all_operations))

