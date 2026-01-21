class AtomDisplayMode:
    SAME = 0
    TYPE_INFORMED = 1

class CouplingDisplayMode:
    SOLID = 0
    WIREFRAME = 1

class DisplayOptions:
    show_moments: bool = True
    show_atoms: bool = True
    show_anisotropies: bool = True
    show_couplings: bool = True
    show_edge_couplings: bool = True

    coupling_mode: CouplingDisplayMode = CouplingDisplayMode.SOLID
    atom_mode: AtomDisplayMode = AtomDisplayMode.SAME

    atom_scaling: float = 1.0
    moment_scaling: float = 1.0
    coupling_scaling: float = 1.0
    anisotropy_scaling: float = 1.0

