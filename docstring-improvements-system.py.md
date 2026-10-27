# Docstring Overview Improvements Needed: `system.py`

## `MonoclinicB`, `MonoclinicA`, `MonoclinicC` (lines 294, 346, 398)
- **Current overview**: `Monoclinic lattice system.`
- **Issue**: All three monoclinic lattice system classes have the identical summary `"Monoclinic lattice system."` without distinguishing their unique axis convention (`b`, `a`, or `c`).
- **Suggested clarification needed**: Clarify the unique axis convention in the summary of each class (e.g., unique axis `b`, `a`, or `c`).

## `find_unit_cell_type` (line 710)
- **Current overview**: `Find the lattice system for a unit cell; it should be a length zero or one list, but it might not be.`
- **Issue**: The overview states "it should be a length zero or one list, but it might not be", which is conversational and ambiguous about when zero or multiple lattice systems are returned.
- **Suggested clarification needed**: Explain under what geometric conditions a `UnitCell` matches zero or multiple `LatticeSystem` definitions (e.g., rhombohedral vs. hexagonal/trigonal or degenerate constraints).
