# Docstring Overview Improvements Needed: `exchangegroup.py`

## `BiDirectionFilter` (line 97)
- **Current overview**: `Selects vectors in a given direction.`
- **Issue**: Identical summary to `InDirectionFilter`, failing to distinguish that `BiDirectionFilter` accepts vectors both parallel and antiparallel (`|n . v|`) to `direction`.
- **Suggested clarification needed**: Clarify that `BiDirectionFilter` selects vectors collinear (either parallel or antiparallel) with the specified axis.
