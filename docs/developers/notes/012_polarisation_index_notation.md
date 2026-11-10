Polarised Cross-Sections — Index (Einstein) Notation
===================================================

This note formalises what `pyspinw/polarisation.py` (`calculate_polarised_intensity`) computes,
using index notation with the summation convention (repeated indices are summed).
A companion note, `013_polarisation_matrix_notation.md`, gives the same calculation in matrix notation.

Data layout and independence
----------------------------

The input `Sab` has shape `(3, 3, n_modes, n_q)`. It is best read as a 2D array, over
modes $`\ell`$ and q-points $`q`$, of $`3\times3`$ complex matrices:

```math
S_{ij}^{(\ell, q)}, \qquad i, j \in \{1,2,3\} \ \text{(Cartesian)}
```

All quantities below are computed **independently for each $`(\ell, q)`$ pair**. The only
thing shared across modes is the Blume–Maleev frame, which depends on $`q`$ alone. From here on
we fix $`(\ell, q)`$ and drop those labels, writing $`S_{ij}`$.

Index conventions:

* $`i, j, k, l, m, n`$ are Cartesian indices (the frame `Sab` is expressed in).
* $`\alpha, \beta \in \{x, y, z\}`$ label the Blume–Maleev axes.
* $`a, b`$ label reciprocal-lattice (r.l.u.) components.
* $`\varepsilon_{ijk}`$ is the Levi-Civita symbol (`LC` in the code), $`\delta_{ij}`$ the Kronecker delta.

Inputs per $`(\ell, q)`$
----------------------

| Symbol | Code | Meaning |
|---|---|---|
| $`S_{ij}`$ | `Sab[:, :, l, q]` | spin–spin correlation matrix |
| $`S_\perp`$ | `Sperp[q, l]` | unpolarised intensity (scalar, supplied by the caller) |
| $`q_a`$ | `q_vectors[q]` | momentum transfer in r.l.u. |
| $`R_{ai}`$ | `rlu_to_cart` | r.l.u. → Cartesian transform |
| $`n_i`$ | `normal` | normal to the horizontal scattering plane |

$`S_\perp`$ is not computed in this module. For consistency with the formulae below it should
equal the perpendicular projection $`S_\perp = \Pi_{ij} S_{ij}`$ (with $`\Pi`$ defined below).

Blume–Maleev frame
------------------

```math
\hat{x}_i = \frac{q_a R_{ai}}{\lVert q_b R_{bj} \rVert}, \qquad
\hat{z}_i = n_i, \qquad
\hat{y}_i = \frac{\varepsilon_{ijk}\, \hat{z}_j \hat{x}_k}{\lVert \varepsilon_{ljk}\, \hat{z}_j \hat{x}_k \rVert}
```

(The code writes $`\hat{y} = -\hat{x}\times\hat{z}`$, which is the same as $`\hat{z}\times\hat{x}`$.)

The frame matrix has the basis vectors as columns:

```math
P_{i\alpha} = e^{(\alpha)}_i, \qquad e^{(x)} = \hat{x},\ e^{(y)} = \hat{y},\ e^{(z)} = \hat{z}
```

and $`P^{-1}_{\alpha i}`$ is its inverse, $`P^{-1}_{\alpha i} P_{i\beta} = \delta_{\alpha\beta}`$.
If $`\hat{z} \perp \hat{x}`$ and $`\lVert\hat z\rVert = 1`$ the frame is orthonormal and right-handed, so
$`P^{-1}_{\alpha i} = P_{i\alpha}`$. The code does not enforce this (see *Implementation notes*).

Auxiliary quantities
--------------------

Symmetric and antisymmetric parts (plain transpose, **not** Hermitian conjugate):

```math
S^{S}_{ij} = \tfrac12 (S_{ij} + S_{ji}), \qquad
S^{A}_{ij} = \tfrac12 (S_{ij} - S_{ji})
```

Chiral (axial) vector of $`S`$:

```math
c_i = \varepsilon_{ijk} S_{jk} = \varepsilon_{ijk} S^{A}_{jk},
\qquad\text{and conversely}\qquad
S^{A}_{jk} = \tfrac12\, \varepsilon_{jkn} c_n
```

The inverse follows from $`\varepsilon_{jkn}\varepsilon_{nab} = \delta_{ja}\delta_{kb} - \delta_{jb}\delta_{ka}`$.
For Hermitian $`S`$ ($`S_{ji} = S_{ij}^*`$), $`c`$ is purely imaginary, so $`i\,c`$ is real.

Projector perpendicular to $`\mathbf{Q}`$:

```math
\Pi_{ij} = \delta_{ij} - \hat{x}_i \hat{x}_j
```

### The chiral term (`mat3`)

The code builds the rank-3 tensor `2*mat1 + LC` and contracts it with $`S^A`$:

```math
h_i = \left(\varepsilon_{ijk} - 2\,\varepsilon_{ijm}\,\hat{x}_m \hat{x}_k\right) S^{A}_{jk}
```

This simplifies. The first term is $`c_i`$. For the second, substitute $`S^A_{jk} = \tfrac12\varepsilon_{jkn}c_n`$
and use $`\varepsilon_{ijm}\varepsilon_{jkn} = \delta_{mk}\delta_{in} - \delta_{mn}\delta_{ik}`$:

```math
-2\,\varepsilon_{ijm}\hat{x}_m\hat{x}_k\,\tfrac12\varepsilon_{jkn}c_n
= -\left(\hat{x}_k\hat{x}_k\, c_i - \hat{x}_i\, \hat{x}_n c_n\right)
= -c_i + \hat{x}_i\,(\hat{x}_n c_n)
```

so

```math
\boxed{\,h_i = \hat{x}_i\,\hat{x}_n c_n\,}
```

The chiral term is the component of the chiral vector along $`\mathbf{Q}`$, pointing along $`\mathbf{Q}`$.

Outputs
-------

### `Mab`: projection of $`S`$ into the Blume–Maleev frame

```math
M_{\alpha\beta} = P^{-1}_{\alpha i}\, S_{ij}\, P^{-1}_{\beta j}
```

In an orthonormal frame this is $`M_{\alpha\beta} = e^{(\alpha)}_i S_{ij} e^{(\beta)}_j`$.

### `Pa`: half polarisation (incident beam polarised along $`\alpha`$)

```math
I_\alpha = S_\perp + i\, P_{i\alpha}\, h_i = S_\perp + i\,(e^{(\alpha)}_i \hat{x}_i)(\hat{x}_n c_n)
```

Note this uses $`P_{i\alpha}`$ (the transpose), not $`P^{-1}_{\alpha i}`$.

### `Pab`: full polarisation

```math
\mathcal{P}_{\alpha\beta}
= 2\,P^{-1}_{\alpha i}\,\Pi_{ij}\,S^{S}_{jk}\,\Pi_{kl}\,P_{l\beta}
\;-\; \delta_{\alpha\beta}\,S_\perp
\;-\; i\,P^{-1}_{\alpha i}\,h_i
```

The last term is independent of $`\beta`$. Mapping to the code:

| Term | Code |
|---|---|
| $`\Pi_{kl}P_{l\beta}`$ | `mat4` (`Pi - outer(...)`) |
| $`\Pi_{ij}S^S_{jk}\Pi_{kl}P_{l\beta}`$ | `Pab1 + Pab2` (`Pab2` supplies the $`-\hat x_i \hat x_j`$ part of the left $`\Pi`$) |
| $`P_{i\beta}S_\perp`$ | `Pab3`, which becomes $`\delta_{\alpha\beta}S_\perp`$ after applying $`P^{-1}`$ |
| $`h_i`$ (broadcast over $`\beta`$) | `Pab4` |
| $`P^{-1}_{\alpha i}(\cdots)`$ | final `einsum('ij,jkm', invP, ...)` |

The first index $`\alpha`$ (the first letter in `Pxy`) is the row. Because the chiral term appears in
row $`x`$ for every column, $`\alpha`$ behaves as the *final* polarisation direction and $`\beta`$ as
the *incident* one. These are polarisation-weighted cross-sections; they are not normalised by
the total cross-section.

Orthonormal frame: explicit components
--------------------------------------

When the frame is orthonormal and right-handed, $`P_{i\alpha}\Pi_{ij}P_{j\beta} = \mathrm{diag}(0,1,1)_{\alpha\beta}`$,
$`S_\perp = M_{yy} + M_{zz}`$ and $`\hat{x}_n c_n = M_{yz} - M_{zy}`$. Writing $`M^{S}_{\alpha\beta} = \tfrac12(M_{\alpha\beta}+M_{\beta\alpha})`$:

```math
I_x = S_\perp + i\,(M_{yz} - M_{zy}), \qquad I_y = I_z = S_\perp
```

```math
\mathcal{P}_{xx} = -S_\perp - i\,(M_{yz}-M_{zy}), \qquad
\mathcal{P}_{xy} = \mathcal{P}_{xz} = -i\,(M_{yz}-M_{zy}), \qquad
\mathcal{P}_{yx} = \mathcal{P}_{zx} = 0
```

```math
\mathcal{P}_{yy} = M_{yy} - M_{zz}, \qquad
\mathcal{P}_{zz} = M_{zz} - M_{yy}, \qquad
\mathcal{P}_{yz} = \mathcal{P}_{zy} = 2M^{S}_{yz}
```

For Hermitian $`S`$: $`i(M_{yz}-M_{zy}) = -2\,\mathrm{Im}\,M_{yz}`$ and $`2M^S_{yz} = 2\,\mathrm{Re}\,M_{yz}`$, so all of these are real.

Combining components
--------------------

A component string such as `'Pxy + 0.33*Pyz'` gives terms $`(w_t, X_t)`$, where $`X_t`$ is one element of
$`M`$, $`I`$ or $`\mathcal{P}`$. The returned intensity, per $`(\ell, q)`$, is

```math
I^{(\ell,q)} = \sum_t w_t\, \left(X_t^{(\ell,q)}\right)^*
```

Note the complex conjugate (`_sum_components` calls `.conj()`). With `components='all'`, the code instead
returns the unconjugated arrays `(Pab, Pa, Mab)`, with shapes `(3, 3, n_modes, n_q)`,
`(3, n_modes, n_q)` and `(3, 3, n_modes, n_q)`.

Implementation notes
--------------------

These came up while deriving the above. They are observations about the current code, not part of the formalism.

1. **Inconsistent use of `rlu_to_cart`.** $`\hat{x}`$ uses `q @ rlu_to_cart` ($`q_a R_{ai}`$), but the
   two-vector plane specification uses `rlu_to_cart @ u` ($`R_{ia}u_a`$). Only one of these can be right
   for a given convention. `spinwave.py` also uses `q @ rlu_to_cart`.
2. **The frame need not be orthonormal.** A single normal vector given in the string is not normalised,
   and $`\hat z`$ is never orthogonalised against $`\hat x`$. In that case $`P^{-1} \neq P^T`$, and the half-polarisation
   formula (which uses $`P^T`$) is inconsistent with the full-polarisation and `Mab` formulae (which use $`P^{-1}`$).
3. **Conjugation differs between modes.** Specific components are conjugated; `'all'` is not.
4. **The full/half mixing check never fires.** `typelen = len(pt[1][0])` is always 1. It was presumably meant to be `len(pt[1])`.
5. **Vector parsing.** The regex in `_str_to_vec` accepts only non-negative integers, so `[0 -1 0]` and `[0.5 0 0]` fail.
6. **$`\mathbf{q}=0`$** gives a zero-length $`\hat x`$ and NaNs.
