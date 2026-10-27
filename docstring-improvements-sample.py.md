# Docstring Overview Improvements Needed: `sample.py`

## `Sample3D.spaghetti_plot` (line 242)
- **Current overview**: `Create a spaghetti diagram with intensity as colorfill overplotted by mode energies.`
- **Issue**: The `scale` parameter (`scale: str = 'linear'`) is accepted in the method signature but is never used in the implementation.
- **Suggested clarification needed**: Clarify whether `scale` is intended to be implemented (e.g., logarithmic intensity scaling) or deprecated.

## `MagneticFieldPowder` (line 1174)
- **Current overview**: `Powder in magnetic field (slightly different parallelisation method to Powder).`
- **Issue**: The class is an unimplemented stub subclass of `Sample1D` with no methods or parameters.
- **Suggested clarification needed**: Clarify whether this class is a placeholder for future implementation and document its intended API once implemented.

## `TwoMagnon` (line 1177)
- **Current overview**: `Two magnon excitations in a single crystal.`
- **Issue**: The class is an unimplemented stub subclass of `Sample3D` with no methods or parameters.
- **Suggested clarification needed**: Clarify whether this class is a placeholder for future implementation and document its intended API once implemented.
