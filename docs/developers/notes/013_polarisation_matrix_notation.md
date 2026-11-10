Polarised Cross-Sections — Matrix Notation
==========================================

This note formalises what `pyspinw/polarisation.py` (`calculate_polarised_intensity`) computes,
using matrix/vector notation. A companion note, `012_polarisation_index_notation.md`, gives the
same calculation in index notation, with the step-by-step reductions and a list of implementation issues.

Data layout and independence
----------------------------

`Sab` has shape `(3, 3, n_modes, n_q)`: a 2D array, over modes $`\ell`$ and q-points $`q`$, of
complex $`3\times3`$ matrices $`\mathsf{S}^{(\ell,q)}`$. Every formula below acts on one
$`(\ell, q)`$ pair at a time and involves only $`3\times3`$ matrices and 3-vectors. The frame depends on
$`q`$ only and is shared by all modes at that $`q`$. We fix $`(\ell, q)`$ and write $`\mathsf{S}`$.

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
\quad\text{(columns)}, \qquad
\Pi = \mathsf{I} - \hat{\mathbf{x}}\hat{\mathbf{x}}^T
```

For an orthonormal, right-handed frame, $`\mathsf{P}^{-1} = \mathsf{P}^T`$ and
$`\mathsf{P}^T\Pi\mathsf{P} = \mathsf{D} \equiv \mathrm{diag}(0, 1, 1)`$.

Auxiliary quantities
--------------------

```math
\mathsf{S}^S = \tfrac12(\mathsf{S} + \mathsf{S}^T), \qquad
\mathsf{S}^A = \tfrac12(\mathsf{S} - \mathsf{S}^T)
```

The chiral vector $`\mathbf{c}`$ is the axial vector of $`\mathsf{S}^A`$:

```math
\mathbf{c} = \begin{pmatrix} S_{yz} - S_{zy} \\ S_{zx} - S_{xz} \\ S_{xy} - S_{yx} \end{pmatrix},
\qquad
\mathsf{S}^A = -\tfrac12[\mathbf{c}]_\times,
\qquad
\mathsf{S}^A\mathbf{v} = \tfrac12\,\mathbf{v}\times\mathbf{c}
```

Here $`[\mathbf{c}]_\times`$ is the cross-product matrix, $`[\mathbf{c}]_\times\mathbf{v} = \mathbf{c}\times\mathbf{v}`$.

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
\boxed{\,\mathbf{h} = \hat{\mathbf{x}}\hat{\mathbf{x}}^T\mathbf{c}\,}
```

That is, the chiral vector projected onto $`\mathbf{Q}`$.

Outputs
-------

### `Mab`: $`\mathsf{S}`$ in the Blume–Maleev frame

```math
\mathsf{M} = \mathsf{P}^{-1}\,\mathsf{S}\,\mathsf{P}^{-T}
\qquad(\;=\mathsf{P}^T\mathsf{S}\mathsf{P}\ \text{if orthonormal})
```

### `Pa`: half polarisation, a 3-vector indexed by the incident direction

```math
\mathbf{I} = S_\perp\,\mathbf{1} + i\,\mathsf{P}^T\mathbf{h}
```

### `Pab`: full polarisation, a $`3\times3`$ matrix (row = first letter of `Pab`)

```math
\mathcal{P} = \mathsf{P}^{-1}\left(2\,\Pi\,\mathsf{S}^S\,\Pi\,\mathsf{P} \;-\; S_\perp\,\mathsf{P} \;-\; i\,\mathbf{h}\,\mathbf{1}^T\right)
= 2\,\mathsf{P}^{-1}\Pi\,\mathsf{S}^S\,\Pi\,\mathsf{P} \;-\; S_\perp\,\mathsf{I} \;-\; i\,(\mathsf{P}^{-1}\mathbf{h})\,\mathbf{1}^T
```

The three bracketed terms correspond to `2*(Pab1 + Pab2)`, `Pab3` and `1j*Pab4` in the code.
The outer product $`\mathbf{h}\mathbf{1}^T`$ puts the chiral term in every column. Rows therefore behave as the
final polarisation direction and columns as the incident one. The values are not normalised by the total cross-section.

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

This is a restatement of the formulae above for reference. It has not been checked against the existing loop implementation.
