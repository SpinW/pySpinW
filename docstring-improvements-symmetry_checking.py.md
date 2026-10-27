# Docstring Overview Improvements Needed: `symmetry_checking.py`

## `MatrixConstraints.print_summary` (line 286)
- **Current overview**: `Print out a summary of the exchange properties.`
- **Issue**: `MatrixConstraints` is the base class for both exchange and anisotropy constraints, not just exchange properties.
- **Suggested clarification needed**: Generalize the overview to "Print out a summary of the matrix constraints."

## `AnisotropyMatrixConstraints.text_summary` (line 436)
- **Current overview**: `Textual summary of symmetry allowed exchange properties.`
- **Issue**: Appears to be a copy-paste from `ExchangeMatrixConstraints.text_summary`; this method summarizes symmetry-allowed anisotropy matrices, not exchange properties.
- **Suggested clarification needed**: Update the overview to refer to symmetry-allowed anisotropy properties.
