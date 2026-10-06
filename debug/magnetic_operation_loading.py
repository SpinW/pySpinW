from pyspinw.symmetry.magnetic_group import build_group
from pyspinw.symmetry.operations import MagneticOperation

operations_string = "x,y,z,+1; -x+1/2,y+1/2,-z+1/2,+1; -x,-y,-z,+1; x+1/2,-y+1/2,z+1/2,+1"
operations_string = "x,y,z,+1; -x+1/2,y+1/2,-z+1/2,+1; -x,-y,-z,+1"#; x+1/2,-y+1/2,z+1/2,+1"
centering_string = "x,y,z,+1; x+1/2,y+1/2,z,-1"

mag_ops = []
centering_ops = []

for string, output in [(operations_string, mag_ops), (centering_string, centering_ops)]:
    for op_str in string.split(";"):
        op_str = op_str.strip()

        op = MagneticOperation.from_text(op_str)

        print(op)

        output.append(op)

magnetic_group = build_group(mag_ops, centering_ops)
space_group = magnetic_group.get_spacegroup()

print(space_group)