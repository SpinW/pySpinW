# Docstring Overview Improvements Needed: `energy_minimisation.py`

## `OneFixed` (line 86)
- **Current overview**: `Constraint Generator.`
- **Issue**: The overview only states that the class is a constraint generator without explaining the specific constraint pattern it produces (fixing a single site at `index` while leaving all other sites free).
- **Suggested clarification needed**: Clarify that `OneFixed` generates a constraint list where the site at `index` is assigned a `Fixed` constraint and all remaining sites are assigned `Free` constraints.
