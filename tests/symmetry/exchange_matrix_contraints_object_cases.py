""" Various examples of allowed exchange matrices determined independently using Mathematica"""
from pyspinw import UnitCell
from pyspinw.symmetry.operations import SpaceOperation
from pyspinw.symmetry.symmetry_checking import ExchangeMatrixConstraints

ops1 = "{x, -y, z}, {-x, y, z}, {y, x, z}, {-y, x, z}, {-x, -y, z}, {x, y, z}, {y, -x, z}, {-y, -x, z}"
sym1 = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 2]]
asym1 = [0, 0, 0]

ops2 = "{-x, -y, z}, {x, -y, z}, {-x, y, z}, {x, y, z}"
sym2 = [
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 3]
]
asym2 = [0, 0, 0]

ops3 = "{y, x, z}, {x, y, z}"
sym3 = [
    [1, 2, 3],
    [2, 1, 3],
    [3, 3, 4]
]
asym3 = [0,1,-1]

ops4 = "{x, y, z}"
sym4 = [
    [1, 2, 3],
    [2, 4, 5],
    [3, 5, 6]
]
asym4 = [1,2,3]

ops = [ops1, ops2, ops3, ops4]
sym = [sym1, sym2, sym3, sym4]
asym = [asym1, asym2, asym3, asym4]

def string_to_ops(string: str) -> list[SpaceOperation]:
    parts = [part.strip() for part in string.replace("}", "").split(r"{") if part.strip() != ""]
    parts = [part[:-1] if part.endswith(",") else part for part in parts]

    return [SpaceOperation.from_text(part) for part in parts]

cell = UnitCell(1,1,1)

for ops_string in ops:
    ops_list = string_to_ops(ops_string)
    transforms = [op.point_operation_in_cartesian(cell) for op in ops_list]
    constraints: ExchangeMatrixConstraints = ExchangeMatrixConstraints(transforms)
    constraints.print_summary()

    print(constraints.free)
    print(constraints.zero)
