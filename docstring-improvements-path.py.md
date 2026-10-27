# Docstring Overview Improvements Needed: `path.py`

## `Path` (line 11)
- **Current overview**: `Path through q-space.`
- **Issue**: The overview does not clarify what coordinate system `points` are expected in (reciprocal lattice units vs. inverse Angstroms) or how `convert_to_lattice_units_with` transforms them (the inline comment says `# Convert to RLU`, the parameter name says `convert_to_lattice_units_with`, and the implementation calls `unit_cell.lattice_units_to_cartesian(self._points) # inverse transform`).
- **Suggested clarification needed**: Clarify the default coordinate convention of `points` and the exact direction/purpose of `convert_to_lattice_units_with`.
