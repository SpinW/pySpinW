# Docstring Overview Improvements Needed: `exchange.py`

## `ExchangeBaseDeserialisation` (line 32)
- **Current overview**: `Hold basic properties of an exchange term.`
- **Issue**: The overview does not explain that this dataclass is an internal deserialisation container returned by `Exchange._base_deserialisation` to bundle deserialised base fields for `Exchange` and its subclasses.
- **Suggested clarification needed**: Clarify that this is an intermediate data structure used during deserialisation of `Exchange` objects rather than a user-facing exchange representation.

## `Exchange.generalise` (line 221)
- **Current overview**: `Convert to most general Exchange class.`
- **Issue**: The summary is terse and does not clarify whether it returns a new `Exchange` instance or modifies the object in place, nor how subclass parameters are represented in the returned general `Exchange`.
- **Suggested clarification needed**: Specify that it returns a new base `Exchange` instance preserving the full 3x3 `exchange_matrix`, sites, offset, name, and metadata.

## `Exchange._obeys_symmetry` (line 377)
- **Current overview**: `Main logic for symmetry checking.`
- **Issue**: Vague overview that describes the role of the method in the call graph rather than the mathematical condition being checked.
- **Suggested clarification needed**: Explain the symmetry invariance criteria for both `identity_operations` (operations preserving site order, $R M R^{-1} = M$) and `inversion_operations` (operations swapping the two sites, $R M R^{-1} = M^T$), and clarify the coordinate frame (lattice vs. Cartesian) used in each check.

## `Exchange.symmetry_copy` (line 468)
- **Current overview**: `Copy this exchange using symmetry operations.`
- **Issue**: Does not explain how a space-group operation mapping `(self.site_1, self.site_2)` to `(site_1, site_2)` transforms the exchange matrix, nor why `cell_offset` is passed manually rather than computed from the symmetry operation (in contrast to `symmetry_fill`).
- **Suggested clarification needed**: Clarify how `cell_offset` is handled and when a user should call `symmetry_copy` vs. `symmetry_fill`.

## `Exchange._specialise` (line 654) and subclass `_specialise` methods (lines 822, 1008, 1160, 1328, 1478, 1657)
- **Current overview**: `Convert this to a specialised exchange type.` / `Create a specialised version of this exchange.`
- **Issue**: Does not state the contract of `_specialise`—namely that it inspects `exchange.exchange_matrix` to test whether it matches the specific algebraic form of the target subclass, returning a new subclass instance if so or `None` if not.
- **Suggested clarification needed**: State the specific matrix structure required by each subclass and note that `None` is returned when the input exchange matrix does not match that structure.
