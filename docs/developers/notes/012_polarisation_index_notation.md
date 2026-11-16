Polarised Cross-Sections — Index (Einstein) Notation
===================================================

This note formalises what `pyspinw/polarisation.py` (`calculate_polarised_intensity`) computes,
using index notation with the summation convention (repeated indices are summed).
A companion note, `013_polarisation_matrix_notation.md`, gives the same calculation in matrix notation.

Each equation is followed by its geometric meaning. The last section, *Overall interpretation
and compact form*, rewrites the whole calculation as the Blume–Maleev equations for the
perpendicular correlation tensor. If you only want the result, read that section.

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

Geometrically, $`S_{ij}`$ is the correlation $`\langle m_i m_j^*\rangle`$ of the Fourier
component $`\mathbf{m}`$ of the magnetisation for that mode and wavevector. It is a dyadic (a rank-2 tensor), not just an array of numbers.
Conjugation conventions differ, but the structure below does not depend on them.

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
equal $`S_\perp = \Pi_{ij} S_{ij}`$, with $`\Pi`$ defined below.

Blume–Maleev frame
------------------

```math
\hat{x}_i = \frac{q_a R_{ai}}{\lVert q_b R_{bj} \rVert}, \qquad
\hat{z}_i = n_i, \qquad
\hat{y}_i = \frac{\varepsilon_{ijk}\, \hat{z}_j \hat{x}_k}{\lVert \varepsilon_{ljk}\, \hat{z}_j \hat{x}_k \rVert}
```

(The code writes $`\hat{y} = -\hat{x}\times\hat{z}`$, which is the same as $`\hat{z}\times\hat{x}`$.)

*Geometry:* $`\hat x`$ is along $`\mathbf{Q}`$. $`\hat z`$ is vertical, normal to the scattering plane.
$`\hat y`$ lies in the scattering plane, perpendicular to $`\mathbf{Q}`$, and completes a right-handed set.
The $`yz`$ plane is therefore the plane perpendicular to $`\mathbf{Q}`$, the only plane the neutron is sensitive to.

The frame matrix has the basis vectors as columns:

```math
P_{i\alpha} = e^{(\alpha)}_i, \qquad e^{(x)} = \hat{x},\ e^{(y)} = \hat{y},\ e^{(z)} = \hat{z}
```

$`P^{-1}_{\alpha i}`$ is its inverse, $`P^{-1}_{\alpha i} P_{i\beta} = \delta_{\alpha\beta}`$.
The rows of $`P^{-1}`$ form the dual basis: $`P^{-1}_{\alpha i}v_i`$ are the expansion coefficients of $`v`$
in the basis $`e^{(\alpha)}`$. If $`\hat{z} \perp \hat{x}`$ and $`\lVert\hat z\rVert = 1`$, the frame is orthonormal and
$`P^{-1}_{\alpha i} = P_{i\alpha}`$, so the coefficients are simply projections, $`e^{(\alpha)}_i v_i`$.
The code does not enforce this (see *Implementation notes*).

Auxiliary quantities
--------------------

### Projector perpendicular to $`\mathbf{Q}`$

```math
\Pi_{ij} = \delta_{ij} - \hat{x}_i \hat{x}_j
```

*Geometry:* $`\Pi`$ projects orthogonally onto the $`yz`$ plane. Neutrons scatter only from
$`\mathbf{M}_\perp = \Pi\,\mathbf{m}`$, the part of the moment perpendicular to $`\mathbf{Q}`$. The corresponding
**perpendicular correlation tensor** is

```math
T_{ij} = \Pi_{ik}\,S_{kl}\,\Pi_{lj} \;\;\bigl(= \langle M_{\perp i} M_{\perp j}^*\rangle\bigr)
```

The code never forms $`T`$ explicitly, but every polarised output depends on $`S`$ only through $`T`$
(see the final section). In particular $`S_\perp = \Pi_{ij}S_{ij} = T_{ii}`$ is the trace of $`T`$:
the total unpolarised magnetic intensity, $`\langle|\mathbf{M}_\perp|^2\rangle`$.

### Symmetric and antisymmetric parts

Plain transpose, **not** Hermitian conjugate:

```math
S^{S}_{ij} = \tfrac12 (S_{ij} + S_{ji}), \qquad
S^{A}_{ij} = \tfrac12 (S_{ij} - S_{ji})
```

*Geometry:* for Hermitian $`S`$, $`S^S = \mathrm{Re}\,S`$ is the "ellipse" part, the real quadratic form whose
principal axes are the axes of the moment's precession ellipse. $`S^A = i\,\mathrm{Im}\,S`$ holds the handedness of
the precession, and is equivalent to a single vector:

### Chiral (axial) vector

```math
c_i = \varepsilon_{ijk} S_{jk} = \varepsilon_{ijk} S^{A}_{jk},
\qquad\text{and conversely}\qquad
S^{A}_{jk} = \tfrac12\, \varepsilon_{jkn} c_n
```

The inverse follows from $`\varepsilon_{jkn}\varepsilon_{nab} = \delta_{ja}\delta_{kb} - \delta_{jb}\delta_{ka}`$.

*Geometry:* $`c_i = \langle(\mathbf{m}\times\mathbf{m}^*)_i\rangle`$. It is normal to the plane in which the
moment precesses, and its sign gives the sense of rotation. For Hermitian $`S`$, $`c`$ is purely imaginary,
so $`i\,c`$ is real. $`c = 0`$ for collinear (linearly polarised) fluctuations.

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
\boxed{\,h_i = \hat{x}_i\,\hat{x}_n c_n = \varepsilon_{ijk}\,T_{jk}\,}
```

*Geometry:* $`h`$ is the chiral vector of the **perpendicular** moment,
$`\mathbf{h} = \langle\mathbf{M}_\perp\times\mathbf{M}_\perp^*\rangle`$. Both $`\mathbf{M}_\perp`$ and $`\mathbf{M}_\perp^*`$
lie in the $`yz`$ plane, so their cross product must point along $`\mathbf{Q}`$. Only the component of the
chirality along $`\mathbf{Q}`$ survives, which is why polarised neutrons measure only that component.

Outputs
-------

### `Mab`: $`S`$ expressed in the Blume–Maleev frame

```math
M_{\alpha\beta} = P^{-1}_{\alpha i}\, S_{ij}\, P^{-1}_{\beta j}
\qquad\Longleftrightarrow\qquad
S_{ij} = M_{\alpha\beta}\, e^{(\alpha)}_i e^{(\beta)}_j
```

*Geometry:* $`M_{\alpha\beta}`$ are the components of the tensor $`S`$ in the dyad basis
$`e^{(\alpha)}\otimes e^{(\beta)}`$. This is a change of basis only; no projection is applied. In an orthonormal frame,
$`M_{\alpha\beta} = e^{(\alpha)}_i S_{ij} e^{(\beta)}_j`$ is the correlation between the moment components along
$`e^{(\alpha)}`$ and $`e^{(\beta)}`$. Its $`yz`$ block is $`T`$ in that frame.

### `Pa`: half polarisation (incident beam polarised along $`e^{(\alpha)}`$)

```math
I_\alpha = S_\perp + i\, P_{i\alpha}\, h_i = S_\perp + i\,(e^{(\alpha)}_i \hat{x}_i)(\hat{x}_n c_n)
```

Note this uses $`P_{i\alpha}`$ (the transpose), not $`P^{-1}_{\alpha i}`$.

*Geometry:* this is the Blume–Maleev cross-section $`\sigma(\mathbf{p}) = S_\perp + i\,\mathbf{p}\cdot\mathbf{h}`$
evaluated for an incident polarisation $`\mathbf{p} = e^{(\alpha)}`$. Since $`\mathbf{h}\parallel\mathbf{Q}`$,
only the component of $`\mathbf{p}`$ along $`\mathbf{Q}`$ sees the chirality. $`I_y`$ and $`I_z`$ are just the
unpolarised intensity.

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

*Geometry:* since $`\Pi S^S \Pi = T^S`$, column $`\beta`$ is the Cartesian vector

```math
\sigma\mathbf{p}_f\big|_{\mathbf{p}=e^{(\beta)}} = \left(2T^S_{ij} - T_{kk}\,\delta_{ij}\right)e^{(\beta)}_j \;-\; i\,h_i
```

expressed in the Blume–Maleev frame by $`P^{-1}`$. This is the Blume–Maleev scattered polarisation (times
cross-section) for incident polarisation $`e^{(\beta)}`$. So $`\alpha`$ (the first letter in `Pxy`) labels the
*final* polarisation component and $`\beta`$ the *incident* direction. The two parts have distinct meanings:

* $`A_{ij} = 2T^S_{ij} - T_{kk}\delta_{ij}`$ is a linear map that **rotates and depolarises** the incident
  polarisation. It sends $`\hat x \mapsto -S_\perp\hat x`$: polarisation along $`\mathbf{Q}`$
  is always reversed (all magnetic scattering is spin-flip for $`\mathbf{p}\parallel\mathbf{Q}`$). Within the $`yz`$ plane it
  is traceless and symmetric, so it acts as a scaled reflection. For a collinear moment $`\mathbf{M}_\perp = |M_\perp|\hat{\mathbf{u}}`$,
  $`A = |M_\perp|^2(2\hat u\hat u^T - I)`$, a rotation by $`\pi`$ about $`\hat{\mathbf{u}}`$: the polarisation precesses half a turn about $`\mathbf{M}_\perp`$.
* $`-i\,h_i`$ is a **created** polarisation along $`\mathbf{Q}`$, independent of the incident one. Chiral
  fluctuations polarise even an unpolarised beam.

These values are polarisation-weighted cross-sections, not polarisations: they are not divided by $`\sigma`$.

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

*Geometry:* $`M_{xx}`$, $`M_{xy}`$ and $`M_{xz}`$ never appear, because moment components along $`\mathbf{Q}`$ are invisible.
$`\mathcal{P}_{yy} = M_{yy} - M_{zz}`$ is the non-spin-flip minus spin-flip intensity for $`\mathbf{p}\parallel\hat y`$, which gives the
anisotropy of the fluctuations within the $`yz`$ plane. $`\mathcal{P}_{yz}`$ measures the correlation between the $`y`$ and $`z`$
components. The $`x`$ row combines the guaranteed spin flip with the chiral term.

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

Overall interpretation and compact form
---------------------------------------

All of the above checks numerically against the code.

### Everything depends on $`T`$ alone

The polarised outputs need only three objects, all built from the perpendicular correlation tensor
$`T_{ij} = \Pi_{ik}S_{kl}\Pi_{lj}`$:

```math
\underbrace{T_{kk}}_{S_\perp:\ \text{intensity}}, \qquad
\underbrace{T^S_{ij}}_{\text{ellipse of } \mathbf{M}_\perp}, \qquad
\underbrace{h_i = \varepsilon_{ijk}T_{jk}}_{\text{chirality of } \mathbf{M}_\perp,\ \parallel \mathbf{Q}}
```

The intermediate steps in the code (separate $`\Pi`$ products in `Pab1` and `Pab2`, the rank-3 tensor in `mat3`)
are just ways of computing these three objects without forming $`T`$.

### The Blume–Maleev equations

For an arbitrary incident polarisation $`p_j`$ (magnetic scattering only, no nuclear terms):

```math
\sigma(p) = T_{kk} + i\,p_j h_j,
\qquad
\sigma p^f_i(p) = \underbrace{\left(2T^S_{ij} - T_{kk}\delta_{ij}\right)}_{A_{ij}} p_j \;-\; i\,h_i
```

The code evaluates these for $`p = e^{(\beta)}`$ and expresses the results in the Blume–Maleev frame:

```math
I_\beta = \sigma(e^{(\beta)}), \qquad
\mathcal{P}_{\alpha\beta} = P^{-1}_{\alpha i}\,\sigma p^f_i(e^{(\beta)}) = P^{-1}_{\alpha i}\left(A_{ij}P_{j\beta} - i\,h_i\right)
```

The relative sign of the chiral term in $`\sigma`$ and $`\sigma\mathbf{p}^f`$ ($`+i`$ and $`-i`$) matches Blume–Maleev. The absolute sign
depends on whether $`S_{jk}`$ means $`\langle m_j m_k^*\rangle`$ or $`\langle m_j^* m_k\rangle`$, and on the conjugation in `_sum_components`.

### One $`4\times4`$ matrix

Write $`(\sigma, \sigma\mathbf{p}^f)`$ as a 4-vector, indexed $`\mu, \nu \in \{0,1,2,3\}`$ with $`0`$ the intensity
component. The scattering is then a single linear map on $`(1, \mathbf{p})`$:

```math
\begin{pmatrix}\sigma \\ \sigma p^f_i\end{pmatrix}
= K_{\mu\nu}\begin{pmatrix}1 \\ p_j\end{pmatrix},
\qquad
K = \begin{pmatrix} T_{kk} & i\,h_j \\ -i\,h_i & 2T^S_{ij} - T_{kk}\delta_{ij}\end{pmatrix}
```

$`K`$ is Hermitian when $`S`$ is. `Pa` is row 0 of $`K`$ applied to the basis vectors, and `Pab` is rows 1–3.
This form is the natural object to store or test if the code is refactored.

### Orthonormal frame: a $`2\times2`$ problem

In an orthonormal frame, $`T`$ is $`M`$ with its $`x`$ row and column set to zero. Let $`t`$ be the $`yz`$ block of $`M`$,
with entries $`t_{ab}`$ for $`a, b\in\{y,z\}`$. Then:

```math
S_\perp = t_{aa}, \qquad
\hat x\cdot\mathbf{h} = \varepsilon_{ab}\,t_{ab} = t_{yz} - t_{zy}, \qquad
\mathcal{P}_{ab} = 2t^S_{ab} - t_{cc}\,\delta_{ab}
```

and the $`x`$ row and column of $`\mathcal{P}`$ contain only $`-S_\perp`$ (spin flip) and $`-i\,\hat x\cdot\mathbf{h}`$ (chirality).
The whole calculation is therefore: rotate $`S`$ into the frame, take the $`2\times2`$ block perpendicular to $`\mathbf{Q}`$,
and read off its trace, its antisymmetric part and its traceless symmetric part.
