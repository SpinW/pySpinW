Polarised Cross-Sections — Matrix Notation
==========================================

This note formalises what `pyspinw/polarisation.py` (`calculate_polarised_intensity`) computes,
using matrix/vector notation. A companion note, `012_polarisation_index_notation.md`, gives the
same calculation in index notation, with the step-by-step reductions and a list of implementation issues.

Each equation is followed by its geometric meaning. The last section, *Overall interpretation
and compact form*, rewrites the whole calculation as the Blume–Maleev equations for the
perpendicular correlation matrix. If you only want the result, read that section.

Data layout and independence
----------------------------

`Sab` has shape `(3, 3, n_modes, n_q)`: a 2D array, over modes $`\ell`$ and q-points $`q`$, of
complex $`3\times3`$ matrices $`\mathsf{S}^{(\ell,q)}`$. Every formula below acts on one
$`(\ell, q)`$ pair at a time and involves only $`3\times3`$ matrices and 3-vectors. The frame depends on
$`q`$ only and is shared by all modes at that $`q`$. We fix $`(\ell, q)`$ and write $`\mathsf{S}`$.

Geometrically, $`\mathsf{S} = \langle\mathbf{m}\mathbf{m}^\dagger\rangle`$ is the correlation (outer product) of
the Fourier component $`\mathbf{m}`$ of the magnetisation for that mode and wavevector.

Notation: bold lower case for column 3-vectors, sans-serif capitals for $`3\times3`$ matrices,
$`{}^T`$ for the plain transpose (no conjugation), $`\mathbf{1} = (1,1,1)^T`$, $`\mathsf{I}`$ for the identity.

Inputs per $`(\ell, q)`$
----------------------

| Symbol | Code | Meaning |
|---|---|---|
| $`\mathsf{S}`$ | `Sab[:, :, l, q]` | spin–spin correlation matrix |
| $`S_\perp`$ | `Sperp[q, l]` | unpolarised intensity (scalar, supplied by the caller; should equal $`\mathrm{tr}(\Pi\mathsf{S})`$) |
| $`\mathbf{q}`$ | `q_vectors[q]` | momentum transfer, r.l.u. |
| $`\mathsf{R}`$ | `rlu_to_cart` | r.l.u. → Cartesian |
| $`\mathbf{n}`$ | `normal` | normal to the horizontal plane |

Blume–Maleev frame
------------------

```math
\hat{\mathbf{x}} = \frac{\mathsf{R}^T\mathbf{q}}{\lVert\mathsf{R}^T\mathbf{q}\rVert}, \qquad
\hat{\mathbf{z}} = \mathbf{n}, \qquad
\hat{\mathbf{y}} = \frac{\hat{\mathbf{z}}\times\hat{\mathbf{x}}}{\lVert\hat{\mathbf{z}}\times\hat{\mathbf{x}}\rVert}
```

(`q @ rlu_to_cart` is $`\mathsf{R}^T\mathbf{q}`$ as a column vector.)

```math
\mathsf{P} = \begin{pmatrix} \hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \end{pmatrix}
\quad\text{(columns)}
```

*Geometry:* $`\hat{\mathbf{x}}\parallel\mathbf{Q}`$, $`\hat{\mathbf{z}}`$ is vertical, and $`\hat{\mathbf{y}}`$ completes a right-handed
set in the scattering plane. The $`yz`$ plane is the plane perpendicular to $`\mathbf{Q}`$. $`\mathsf{P}`$ maps
Blume–Maleev coordinates to Cartesian ones. $`\mathsf{P}^{-1}`$ maps back, and its rows are the dual basis. For an orthonormal,
right-handed frame, $`\mathsf{P}`$ is a rotation and $`\mathsf{P}^{-1} = \mathsf{P}^T`$.

Auxiliary quantities
--------------------

### Projector perpendicular to $`\mathbf{Q}`$

```math
\Pi = \mathsf{I} - \hat{\mathbf{x}}\hat{\mathbf{x}}^T,
\qquad
\mathsf{T} = \Pi\,\mathsf{S}\,\Pi
```

*Geometry:* $`\Pi`$ is the orthogonal projection onto the $`yz`$ plane, so $`\mathbf{M}_\perp = \Pi\mathbf{m}`$ is the part of the
moment perpendicular to $`\mathbf{Q}`$, the only part the neutron sees. $`\mathsf{T} = \langle\mathbf{M}_\perp\mathbf{M}_\perp^\dagger\rangle`$
is the **perpendicular correlation matrix**. The code never forms it explicitly, but every polarised output depends
on $`\mathsf{S}`$ only through $`\mathsf{T}`$. $`S_\perp = \mathrm{tr}\,\mathsf{T} = \langle|\mathbf{M}_\perp|^2\rangle`$.
In an orthonormal frame, $`\mathsf{P}^T\Pi\mathsf{P} = \mathsf{D} \equiv \mathrm{diag}(0, 1, 1)`$.

### Symmetric and antisymmetric parts

```math
\mathsf{S}^S = \tfrac12(\mathsf{S} + \mathsf{S}^T), \qquad
\mathsf{S}^A = \tfrac12(\mathsf{S} - \mathsf{S}^T)
```

*Geometry:* for Hermitian $`\mathsf{S}`$, $`\mathsf{S}^S = \mathrm{Re}\,\mathsf{S}`$ is a real quadratic form whose principal axes are
the axes of the moment's precession ellipse. $`\mathsf{S}^A = i\,\mathrm{Im}\,\mathsf{S}`$ holds the sense of precession.

### Chiral vector

The chiral vector $`\mathbf{c}`$ is the axial vector of $`\mathsf{S}^A`$:

```math
\mathbf{c} = \begin{pmatrix} S_{yz} - S_{zy} \\ S_{zx} - S_{xz} \\ S_{xy} - S_{yx} \end{pmatrix},
\qquad
\mathsf{S}^A = -\tfrac12[\mathbf{c}]_\times,
\qquad
\mathsf{S}^A\mathbf{v} = \tfrac12\,\mathbf{v}\times\mathbf{c}
```

Here $`[\mathbf{c}]_\times`$ is the cross-product matrix, $`[\mathbf{c}]_\times\mathbf{v} = \mathbf{c}\times\mathbf{v}`$.

*Geometry:* $`\mathbf{c} = \langle\mathbf{m}\times\mathbf{m}^*\rangle`$ is normal to the plane of precession, and its sign is the
handedness. It is purely imaginary for Hermitian $`\mathsf{S}`$, and zero for collinear fluctuations.

### The chiral term (`mat3`)

The code computes

```math
\mathbf{h} = \mathbf{c} - 2\,(\mathsf{S}^A\hat{\mathbf{x}})\times\hat{\mathbf{x}}
```

Using $`\mathsf{S}^A\hat{\mathbf{x}} = \tfrac12\hat{\mathbf{x}}\times\mathbf{c}`$ and
$`(\mathbf{a}\times\mathbf{b})\times\mathbf{d} = \mathbf{b}(\mathbf{a}\cdot\mathbf{d}) - \mathbf{a}(\mathbf{b}\cdot\mathbf{d})`$:

```math
-2(\mathsf{S}^A\hat{\mathbf{x}})\times\hat{\mathbf{x}} = (\mathbf{c}\times\hat{\mathbf{x}})\times\hat{\mathbf{x}}
= \hat{\mathbf{x}}(\hat{\mathbf{x}}^T\mathbf{c}) - \mathbf{c}
```

so

```math
\boxed{\,\mathbf{h} = \hat{\mathbf{x}}\hat{\mathbf{x}}^T\mathbf{c} = \mathrm{axial}(\mathsf{T})\,}
```

*Geometry:* $`\mathbf{h} = \langle\mathbf{M}_\perp\times\mathbf{M}_\perp^*\rangle`$ is the chirality of the perpendicular moment. Both factors lie in
the $`yz`$ plane, so the cross product points along $`\mathbf{Q}`$. Only the component of $`\mathbf{c}`$ along $`\mathbf{Q}`$ is observable.

Outputs
-------

### `Mab`: $`\mathsf{S}`$ in the Blume–Maleev frame

```math
\mathsf{M} = \mathsf{P}^{-1}\,\mathsf{S}\,\mathsf{P}^{-T}
\qquad\Longleftrightarrow\qquad
\mathsf{S} = \mathsf{P}\,\mathsf{M}\,\mathsf{P}^T
\qquad(\mathsf{M} = \mathsf{P}^T\mathsf{S}\mathsf{P}\ \text{if orthonormal})
```

*Geometry:* this is a change of basis only, with no projection. $`\mathsf{M}`$ holds the coefficients of $`\mathsf{S}`$ on the dyads
$`\mathbf{e}_\alpha\mathbf{e}_\beta^T`$. In an orthonormal frame it is $`\mathsf{S}`$ rotated into the frame, and
$`\mathsf{D}\mathsf{M}\mathsf{D}`$ is $`\mathsf{T}`$ in the frame.

### `Pa`: half polarisation, a 3-vector indexed by the incident direction

```math
\mathbf{I} = S_\perp\,\mathbf{1} + i\,\mathsf{P}^T\mathbf{h}
```

*Geometry:* entry $`\alpha`$ is the Blume–Maleev cross-section $`\sigma(\mathbf{p}) = S_\perp + i\,\mathbf{p}\cdot\mathbf{h}`$ for incident
polarisation $`\mathbf{p} = \mathbf{e}_\alpha`$. Because $`\mathbf{h}\parallel\mathbf{Q}`$, only polarisation along $`\mathbf{Q}`$ detects chirality.

### `Pab`: full polarisation, a $`3\times3`$ matrix (row = first letter of `Pab`)

```math
\mathcal{P} = \mathsf{P}^{-1}\left(2\,\Pi\,\mathsf{S}^S\,\Pi\,\mathsf{P} \;-\; S_\perp\,\mathsf{P} \;-\; i\,\mathbf{h}\,\mathbf{1}^T\right)
= 2\,\mathsf{P}^{-1}\Pi\,\mathsf{S}^S\,\Pi\,\mathsf{P} \;-\; S_\perp\,\mathsf{I} \;-\; i\,(\mathsf{P}^{-1}\mathbf{h})\,\mathbf{1}^T
```

The three bracketed terms correspond to `2*(Pab1 + Pab2)`, `Pab3` and `1j*Pab4` in the code.

*Geometry:* since $`\Pi\mathsf{S}^S\Pi = \mathsf{T}^S`$, this is

```math
\mathcal{P} = \mathsf{P}^{-1}\left(\mathsf{A}\,\mathsf{P} - i\,\mathbf{h}\,\mathbf{1}^T\right),
\qquad
\mathsf{A} = 2\,\mathsf{T}^S - (\mathrm{tr}\,\mathsf{T})\,\mathsf{I}
```

Column $`\beta`$ is the Blume–Maleev scattered polarisation (times cross-section), $`\sigma\mathbf{p}_f = \mathsf{A}\mathbf{p} - i\mathbf{h}`$,
for incident $`\mathbf{p} = \mathbf{e}_\beta`$, expressed in frame coordinates. Rows are therefore the final direction and
columns the incident one. The two parts are:

* $`\mathsf{A}`$ is the **rotation/depolarisation** of the incident polarisation. It sends $`\hat{\mathbf{x}} \mapsto -S_\perp\hat{\mathbf{x}}`$: polarisation
  along $`\mathbf{Q}`$ is always flipped. On the $`yz`$ plane it is symmetric and traceless, so it acts as a scaled reflection. For a
  collinear moment $`\mathbf{M}_\perp = |M_\perp|\hat{\mathbf{u}}`$, $`\mathsf{A} = |M_\perp|^2(2\hat{\mathbf{u}}\hat{\mathbf{u}}^T - \mathsf{I})`$, a rotation by $`\pi`$ about $`\hat{\mathbf{u}}`$.
* $`-i\mathbf{h}\mathbf{1}^T`$ is polarisation **created** along $`\mathbf{Q}`$ by chiral fluctuations. It is independent of the incident polarisation,
  so the same vector appears in every column.

The values are not normalised by the total cross-section.

Orthonormal frame: everything in terms of $`\mathsf{M}`$
------------------------------------------------------

With $`\mathsf{P}^{-1} = \mathsf{P}^T`$, $`\mathsf{M}^S = \tfrac12(\mathsf{M}+\mathsf{M}^T)`$,
$`S_\perp = M_{yy}+M_{zz}`$ and $`\chi \equiv \hat{\mathbf{x}}^T\mathbf{c} = M_{yz} - M_{zy}`$:

```math
\mathbf{I} = S_\perp\mathbf{1} + i\chi\,\mathbf{e}_x
= \begin{pmatrix} S_\perp + i\chi \\ S_\perp \\ S_\perp \end{pmatrix}
```

```math
\mathcal{P} = 2\,\mathsf{D}\,\mathsf{M}^S\,\mathsf{D} - S_\perp\mathsf{I} - i\chi\,\mathbf{e}_x\mathbf{1}^T
= \begin{pmatrix}
-S_\perp - i\chi & -i\chi & -i\chi \\
0 & M_{yy} - M_{zz} & 2M^S_{yz} \\
0 & 2M^S_{yz} & M_{zz} - M_{yy}
\end{pmatrix}
```

For Hermitian $`\mathsf{S}`$, $`i\chi = -2\,\mathrm{Im}\,M_{yz}`$ and $`2M^S_{yz} = 2\,\mathrm{Re}\,M_{yz}`$, so every entry is real.

*Geometry:* the $`x`$ row and column of $`\mathsf{M}`$ never appear, because components along $`\mathbf{Q}`$ are invisible. The $`yz`$ block of
$`\mathcal{P}`$ is the traceless part of the in-plane correlations, measuring the anisotropy ($`M_{yy}-M_{zz}`$) and the cross-correlation ($`M_{yz}`$)
of the fluctuations perpendicular to $`\mathbf{Q}`$. The first row combines the guaranteed spin flip with the chiral term.

Combining components
--------------------

For terms $`(w_t, X_t)`$ parsed from the component string, with $`X_t`$ an entry of $`\mathsf{M}`$,
$`\mathbf{I}`$ or $`\mathcal{P}`$:

```math
I^{(\ell,q)} = \sum_t w_t\,\overline{X_t^{(\ell,q)}}
```

Equivalently, for full-polarisation or `M` terms, $`I = \mathrm{tr}(\mathsf{W}^T\overline{\mathcal{X}})`$,
where $`\mathsf{W}`$ is the $`3\times3`$ weight matrix with entries $`w_t`$ at the selected positions.
`components='all'` returns $`(\mathcal{P}, \mathbf{I}, \mathsf{M})`$ stacked over $`(\ell, q)`$, without conjugation.

Batched form
------------

Because each $`(\ell,q)`$ is independent and the frame depends only on $`q`$, all of the above can be
vectorised with the matrix axes leading and $`(\ell, q)`$ trailing. Take `Pinv` with shape `(n_q, 3, 3)`,
`P` with shape `(n_q, 3, 3)`, `Pi` (the projector $`\Pi`$) with shape `(n_q, 3, 3)`, `xhat` with shape `(n_q, 3)`,
`S` with shape `(3, 3, n_modes, n_q)` and `Sperp` with shape `(n_q, n_modes)`:

```python
M     = np.einsum('qai,ijlq,qbj->ablq', Pinv, S, Pinv)
SS    = (S + S.transpose(1, 0, 2, 3)) / 2
c     = np.einsum('ijk,jklq->ilq', LC, S)
chi   = np.einsum('qi,ilq->lq', xhat, c)                     # x̂·c
h     = xhat.T[:, None, :] * chi                             # (3, n_modes, n_q)
I     = Sperp.T[None] + 1j * np.einsum('qia,ilq->alq', P, h)
Pfull = (2 * np.einsum('qai,qij,jklq,qkm,qmb->ablq', Pinv, Pi, SS, Pi, P)
         - np.eye(3)[:, :, None, None] * Sperp.T[None, None]
         - 1j * np.einsum('qai,ilq->alq', Pinv, h)[:, None])
```

This matches the loop implementation's `'all'` output to numerical precision. It was checked with random Hermitian
$`\mathsf{S}`$, an orthonormal frame, and $`S_\perp = \mathrm{tr}(\Pi\mathsf{S})`$.

Overall interpretation and compact form
---------------------------------------

All of the identities in this section were checked numerically against the code.

### Everything depends on $`\mathsf{T} = \Pi\mathsf{S}\Pi`$

The polarised outputs use only three properties of $`\mathsf{T}`$, each with a direct geometric meaning:

| Quantity | Formula | Meaning |
|---|---|---|
| $`S_\perp`$ | $`\mathrm{tr}\,\mathsf{T}`$ | total magnetic intensity, $`\langle\lvert\mathbf{M}_\perp\rvert^2\rangle`$ |
| $`\mathsf{T}^S`$ | $`\tfrac12(\mathsf{T}+\mathsf{T}^T)`$ | precession ellipse of $`\mathbf{M}_\perp`$ (real part, if Hermitian) |
| $`\mathbf{h}`$ | $`\mathrm{axial}(\mathsf{T}) = \hat{\mathbf{x}}\hat{\mathbf{x}}^T\mathbf{c}`$ | chirality of $`\mathbf{M}_\perp`$, always $`\parallel\mathbf{Q}`$ |

The separate $`\Pi`$ products (`Pab1`, `Pab2`) and the rank-3 tensor (`mat3`) are just roundabout ways of
computing these without forming $`\mathsf{T}`$.

### The Blume–Maleev equations

For an arbitrary incident polarisation $`\mathbf{p}`$ (magnetic scattering only):

```math
\sigma(\mathbf{p}) = \mathrm{tr}\,\mathsf{T} + i\,\mathbf{p}\cdot\mathbf{h},
\qquad
\sigma\mathbf{p}_f(\mathbf{p}) = \mathsf{A}\,\mathbf{p} - i\,\mathbf{h},
\qquad
\mathsf{A} = 2\mathsf{T}^S - (\mathrm{tr}\,\mathsf{T})\,\mathsf{I}
```

The outputs are these evaluated on the frame vectors, then written in frame coordinates:

```math
\mathbf{I} = \sigma(\mathbf{e}_\alpha)_{\alpha=x,y,z}, \qquad
\mathcal{P} = \mathsf{P}^{-1}\begin{pmatrix}\sigma\mathbf{p}_f(\hat{\mathbf{x}}) & \sigma\mathbf{p}_f(\hat{\mathbf{y}}) & \sigma\mathbf{p}_f(\hat{\mathbf{z}})\end{pmatrix}, \qquad
\mathsf{M} = \mathsf{P}^{-1}\mathsf{S}\mathsf{P}^{-T}
```

The relative sign of the chiral term ($`+i`$ in $`\sigma`$, $`-i`$ in $`\sigma\mathbf{p}_f`$) matches Blume–Maleev. The absolute sign depends on the
convention for $`\mathsf{S}`$ ($`\langle\mathbf{m}\mathbf{m}^\dagger\rangle`$ or its transpose) and on the conjugation applied in `_sum_components`.

### One $`4\times4`$ matrix

Stack the intensity and polarisation as a 4-vector. The scattering is then one linear map on $`(1, \mathbf{p})`$:

```math
\begin{pmatrix}\sigma \\ \sigma\mathbf{p}_f\end{pmatrix}
= \mathsf{K}\begin{pmatrix}1 \\ \mathbf{p}\end{pmatrix},
\qquad
\mathsf{K} = \begin{pmatrix} \mathrm{tr}\,\mathsf{T} & i\,\mathbf{h}^T \\ -i\,\mathbf{h} & 2\mathsf{T}^S - (\mathrm{tr}\,\mathsf{T})\,\mathsf{I}\end{pmatrix}
```

$`\mathsf{K}`$ is Hermitian when $`\mathsf{S}`$ is. Both polarised outputs are blocks of the $`4\times3`$ matrix
$`\mathrm{diag}(1, \mathsf{P}^{-1})\,\mathsf{K}\,\begin{pmatrix}\mathbf{1}^T \\ \mathsf{P}\end{pmatrix}`$: the top row is `Pa`, and the lower $`3\times3`$ block is `Pab`.
This is the natural single object to compute, store or test if the module is refactored.

### Orthonormal frame: a $`2\times2`$ problem

In an orthonormal frame, $`\mathsf{P}^T\mathsf{T}\mathsf{P} = \mathsf{D}\mathsf{M}\mathsf{D}`$ is $`\mathsf{M}`$ with its $`x`$ row and column zeroed. Let
$`\mathsf{t}`$ be the $`2\times2`$ $`yz`$ block of $`\mathsf{M}`$, and $`\mathsf{J} = \begin{pmatrix}0&1\\-1&0\end{pmatrix}`$. Then

```math
S_\perp = \mathrm{tr}\,\mathsf{t}, \qquad
\hat{\mathbf{x}}\cdot\mathbf{h} = \mathrm{tr}(\mathsf{J}^T\mathsf{t}) = t_{yz} - t_{zy}, \qquad
\mathcal{P}_{yz\text{ block}} = 2\mathsf{t}^S - (\mathrm{tr}\,\mathsf{t})\,\mathsf{I}_2
```

and the $`x`$ row and column of $`\mathcal{P}`$ contain only $`-S_\perp`$ (spin flip) and $`-i\,\hat{\mathbf{x}}\cdot\mathbf{h}`$ (chirality). In words:
rotate $`\mathsf{S}`$ into the Blume–Maleev frame, keep the $`2\times2`$ block perpendicular to $`\mathbf{Q}`$, and split it into its trace
(intensity), its antisymmetric part (chirality) and its traceless symmetric part (rotation of the polarisation).
