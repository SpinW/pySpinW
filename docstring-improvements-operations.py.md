# Docstring Overview Improvements Needed: `operations.py`

## `SpaceOperation.symmorphic` (line 138)
- **Current overview**: `Whether this operation is symmorphic (i.e. has non-integer translation).`
- **Issue**: The parenthetical explanation "(i.e. has non-integer translation)" contradicts both the standard crystallographic definition of a symmorphic operation and the implementation (`self._symmorphic = np.allclose(self.translation, 0.0)`), which returns `True` when the fractional translation is zero (integer lattice translation).
- **Suggested clarification needed**: Confirm whether the parenthetical should read "(i.e. has no non-integer translation)" or "(i.e. has zero fractional translation)".

## `MagneticOperation.from_text` (line 646)
- **Current overview**: `Parse a string form of an expression, e.g. 'x,x-y,z+1/2'.`
- **Issue**: The example `'x,x-y,z+1/2'` only has 3 comma-separated parts, whereas `MagneticOperation.from_text` requires 4 comma-separated parts (including the time reversal factor `+1` or `-1`, e.g. `'x,x-y,z+1/2,-1'`).
- **Suggested clarification needed**: Update the overview example to include the fourth component representing time reversal.
