# Docstring Overview Improvements Needed: `anisotropy.py`

## `Anisotropy` (line 37)
- **Current overview**: `Represent the anisotropy term at a given site.`
- **Issue**: Does not specify the mathematical convention or sign of the single-ion anisotropy term in the spin Hamiltonian (e.g., whether $\mathcal{H}_{\text{anis}} = \sum_i \mathbf{S}_i^T A_i \mathbf{S}_i$ or $-\sum_i \mathbf{S}_i^T A_i \mathbf{S}_i$, and in what coordinate frame `anisotropy_matrix` is defined).
- **Suggested clarification needed**: Document the Hamiltonian sign convention, energy units, and coordinate basis (Cartesian vs. lattice) for the $3 \times 3$ anisotropy matrix $A$.

## `Anisotropy._obeys_symmetry` (line 82)
- **Current overview**: `Main part of the check for whether the anisotropy matrix obeys symmetry.`
- **Issue**: Written as a noun phrase describing its role relative to `obeys_symmetry` rather than explaining the mathematical check performed on the symmetrised matrix $A + A^T$ under Cartesian point operations.
- **Suggested clarification needed**: Clarify that this helper checks whether the symmetrised anisotropy tensor $A + A^T$ is invariant ($R (A + A^T) R^T = A + A^T$) under the supplied point-group operations in Cartesian coordinates.

## `Anisotropy.symmetry_fill` (line 189)
- **Current overview**: `Make multiple copies of this anisotropy so that symmetry is satisfied.`
- **Issue**: Does not explain how duplicate site anisotropies are resolved or how generated terms are specialised; additionally, the parameter `specialisation_rounding_exponent` is accepted in the signature but unused in the method body.
- **Suggested clarification needed**: Clarify how symmetry-equivalent sites are populated (one anisotropy per symmetry-related site) and whether `specialisation_rounding_exponent` should be forwarded to `specialise_anisotropy` or removed.

## `Anisotropy.updated` (line 341)
- **Current overview**: `Return a copy of this anisotropy term with variables replaced.`
- **Issue**: The parameter description for `name` states "If omitted, current name is reused", but the default argument in the signature is `name: str = ""` while the implementation checks `if name is None`, causing an omitted `name` to overwrite the existing name with `""`.
- **Suggested clarification needed**: Confirm whether `name` should default to `None` (matching `AxisMagnitudeAnisotropy.updated`) or whether the docstring should describe the empty-string default behavior.

## `AxisMagnitudeAnisotropy` (line 405)
- **Current overview**: `Represent a uniaxial anisotropy term.`
- **Issue**: Does not state the explicit outer-product relation $A = a \, (\hat{\mathbf{n}} \otimes \hat{\mathbf{n}})$ or which sign of `a` ($a < 0$ vs. $a > 0$) corresponds to easy-axis versus easy-plane anisotropy.
- **Suggested clarification needed**: Specify the formula relating `a` and `direction` to `anisotropy_matrix` and state the sign convention for easy-axis and easy-plane anisotropy.

## `specialise_anisotropy` (line 614)
- **Current overview**: `Find the narrowest anisotropy subclass to fit the exchange.`
- **Issue**: Contains a copy-paste error referring to "the exchange" instead of "the anisotropy".
- **Suggested clarification needed**: Replace "exchange" with "anisotropy" and briefly explain the subclass trial order in `_all_anisotropies`.
