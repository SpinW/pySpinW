# The Maths behind Linear Spin Wave Theory Calculations

This documentation describes the maths behind linear spin wave theory calculations as implemented in the (Matlab) [SpinW](https://github.com/spinw/spinw) code.
A fuller description can be found in the papers of [S. Petit (2011)](https://doi.org/10.1051/sfn/201112006) and [Toth and Lake (2015)](https://doi.org/10.1088/0953-8984/27/16/166002).
(We shall mostly be using the notation of Petit as this is what is used in the SpinW Matlab code).


## Introduction

Spin waves are collective excitations (normal modes) of a lattice of atomic spins
(atoms with unfilled electronic shells and hence a net magnetic moment) coupled by exchange interactions.
These excitations can be described using a semiclassical approach using an equation of motion in which
the internal magnetic field generated by the moment on one atom causes a torque on its neighbour. [citation needed]
The resulting wave is a precession of the spins about their ordered direction, with a net phase between site.
Another way of thinking about these excitations is as (quantised) quasiparticles in a field theory.
This is "linear spin wave theory" and its quasiparticles are "magnons".


## The Holstein-Primakoff transformation

Consider a single atom with total spin quantum number $S$, which can be in states labelled by $`S_z=-S,-S+1,...,S-1,S`$.
In the ordered state, this spin is in the state with maximum $`S_z=S`$.
A small deviation from this (i.e. initiating a spin wave) will change the state to $`S_z=S-1`$.
This can be described by the spin lowering operator $`\hat{S}^-`$.
Likewise, restoring the ordered state from this deviation is described by the spin raising operator $`\hat{S}^+`$.

Now, the _Holstein-Primakoff_ transformation is a mapping between these lowering and raising operators
to bosonic creation $`\hat{b}^{\dagger}`$ and annihilation $`\hat{b}`$ operators as follows (we drop the hats):


$$S^+ = b \sqrt{2S - b^{\dagger}b}$$
$$S^- = b^{\dagger} \sqrt{2S - b^{\dagger}b}$$
$$S^z = (S - b^{\dagger}b)$$

(We've also neglected the spin quantum $n\hbar$ in this treatment compared to standard texts in accordance with Petit and Toth and Lake).
(Note also that the $`\sqrt{2S}`$ in the above refers to the total spin quantum number $S$, which is also often termed the "spin length".
This is because in some of the literature it is taken to represent the magnitude of the ordered moment rather than a (half integral) quantum number.
As such, SpinW does not restrict $S$ to taking integer or half-integer values.
The main effect of changing $S$ is to scale the magnon energies.

The number of magnons is $`b^{\dagger}b`$ and $`S^z`$ is the projection of the spins along the local ordered moment direction.
We thus see that when no magnons are excited this corresponds to the fully ordered state,
and that as more magnons are excited the spins become canted perpendicular to this direction,
as the ladder (raising/lowering) operators can be related to the $x$ and $y$ spin components by $`S^{\pm} = S^x \pm iS^y`$.

The term in the square root is usually expanded in a Taylor series in practical calculations,
and usually only the first order (linear) term is retained which is equivalent to neglecting the $`b^{\dagger}b`$ term in the square root.
This is strictly only valid when $2S$ is large (to see this, rearrange to get $`\sqrt{\frac{1-b^{\dagger}b}{2S}}`$).
Thus, _linear spin wave theory_ is said to be only valid for large $S$ systems.


## Local spin directions

We see in the above that mapping to the bosonic operators $`b^{\dagger}, b`$ requires the Hamiltonian
to be described in terms of the local spin ordered moment direction, since it describes small deviations from this direction.

We thus define a set of rotation matrices $`R_i^{\alpha}`$ which transforms a spin vector $`\mathbf{S}'_i`$
in the local coordinate system (where $z$ is along the ordered moment direction)
to a vector $`\mathbf{S}_i`$ in a Cartesian coordinate system connected to the crystal lattice
(In SpinW, this Cartesian system is defined by $x \parallel a$, $z \perp  (a, c)$ and $y \perp (x, z$)):

```math
\mathbf{S}_i = R_i \mathbf{S}'_i
```

where $`R_i`$ is a $3 \times 3$ matrix.

Additionally, to more easily map to the operators $S^z$, $S^-$ ($`b^{\dagger}`$) and $S^+$ ($b$) operators above,
we will define the following vectors from the columns of $`R_i`$ for each spin:

```math
\mathbf{z}_i = R_i^1 + i R_i^2
```
```math
\boldsymbol{\eta}_i = R_i^3
```

that is $`z_i`$ is formed from the first and second column of $`R_i`$ whilst $`\eta_i`$ from the third column of $`R_i`$.
This is so that we can express the spin vector (in the local coordinate system) in terms of the bosonic operators as:

```math
\mathbf{S}'_i = \sqrt{\frac{S}{2}}\left( \mathbf{z}_i^* b_i + \mathbf{z}_ib_i^{\dagger} \right) + \boldsymbol{\eta}_i \left( S_i - b_i^{\dagger} b_i \right)
```

where we have made the linear approximation and taken:

$$S^x = \frac{(S^+ + S^-)}{2} = \sqrt{\frac{S}{2}}(b + b^{\dagger}) $$
$$S^y = \frac{(S^+ - S^-)}{2i} = \frac{\sqrt{S}}{i\sqrt{2}}(b - b^{\dagger}) $$
$$S^z = S - b^{\dagger}b $$


## The Hamiltonian

The spin-Hamiltonian used in SpinW is:

```math
\mathcal{H} = \sum_{m \neq n} \sum_{i \neq j} \mathbf{S}^{\intercal}_{im} J_{im, jn} \mathbf{S}_{jn} 
    + \sum_{m} \sum_{i} \mathbf{S}^{\intercal}_{im} A_{im} \mathbf{S}_{im} 
    - \mu_B \mathbf{H}^{\intercal} \sum_{m} \sum_{i} g_i \mathbf{S}_{im} 
```

where the $\intercal$ symbol indicates matrix transpose, and
the first term is the pairwise exchange interaction with exchange "constant" $`J_{im,jn}`$,
the second term is the single-ion anisotropy with anistotropy "constant" $`A_{im}`$ and
the third term is the Zeeman energy when a finite magnetic field $`\mathbf{H}`$ is applied to the sample, and
the g-tensor $`g_i`$ is a measure of the sample anisotropy for each magnetic ion in the unit cell
(note it does not depend on $m$) in an applied field.
The "constants" are called such for historical reasons but they are inputs to our model and can be
treated as variable parameters which may be fitted to data.
$`J_{im,jn}`$, $`A_{im}`$ and $`g_i`$ are all tensors represented by $3 \times 3$ matrices.

Now, if we express $`A_{im} = J_{im, jn}\delta_{ij}\delta_{mn} = J_{im}`$ we can express the first two terms together as:

```math
\begin{array}{rcl}
\mathcal{H} &=& \sum_{m,n} \sum_{i,j} \mathbf{S}^{\intercal}_{im} J_{im,jn} \mathbf{S}_{jn} \\
 & = & \sum_{m,n} \sum_{i, j} \mathbf{S}^{\intercal}(\mathbf{r}_{i}+\mathbf{r}_m)
    J(\mathbf{r}_{i}+\mathbf{r}_m, \mathbf{r}_{j}-\mathbf{r}_n) \mathbf{S}(\mathbf{r}_{j}+\mathbf{r}_n)
\end{array}
```

where labels $i$ and $j$ denote sites within a magnetic (super)-unit cell whilst labels $m$ and $n$ index the unit cells themselves.
Note that the sum $`\sum_{ij}`$ is finite whilst the sum $`\sum_{mn}`$ is infinite.
Now, because SpinW deals exclusively with periodic (crystalline) systems, we can Fourier transform the above equation to replace
the infinite sum over unit cells with another (infinite) sum over momentum space $`\mathbf{q}`$-points.
This is actually a step forward because the "Hamiltonian" at each $`\mathbf{q}`$-point can be diagonalised to obtain wave-like solutions
which are the magnon modes which are measured by inelastic neutron scattering.

The Fourier transform of the spin vector at site $(i,m)$ is:

```math
\mathbf{S}_{i,m}(\mathbf{q}) = \sum_m \left[ \sum_i \mathbf{S} \exp(-i\mathbf{r}_i\cdot\mathbf{q}) \right] \exp(-i\mathbf{r}_m\cdot\mathbf{q})
```

and similarly for the $(j,n)$ indexed terms $`\mathbf{S}_{j,n}`$.
In the next part we will denote the term in square brackets as $`\mathbf{S}_i(\mathbf{q})`$.
In a similar vein, the Fourier transform of the exchange interaction becomes:

```math
J_{ij,mn}(\mathbf{q}) = \sum_m \sum_n
    \left[ \sum_i \sum_j J \exp(-i(\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{q}) \right] \exp(-i(\mathbf{r}_m - \mathbf{r}_n)\cdot\mathbf{q})
```

As with the spins, the terms in the square brackets will be denoted $`J_{ij}(\mathbf{q})`$.
Note that for single-ion anisotropy terms $`A_i = J_{ii}`$ where $i=j$ and $`\mathbf{r}=\mathbf{0}`$ the phase factor in $`J_{ij}(\mathbf{q})`$ is unity.

We can now expressed the Fourier transform of the spin Hamiltonian:

```math
\mathcal{H}(\mathbf{q}) = \sum_{m,n} \left[ \sum_{i,l} \mathbf{S}_i(\mathbf{q}) J_{ij}(\mathbf{q}) \mathbf{S}_j(\mathbf{q}') \right]
    \exp(-i\mathbf{r}_m\cdot\mathbf{q}) \exp(-i\mathbf{r}_n\cdot\mathbf{q}') \exp(-i(\mathbf{r}_m - \mathbf{r}_n)\cdot\mathbf{q})
```

Then we can use the identity $`\sum_{r} \exp(i(\mathbf{q}-\mathbf{q}')\cdot\mathbf{r}) = \delta_{\mathbf{qq}'}`$ to cancel the
$`\mathbf{r}_n`$ terms leaving a single Fourier series in $`\mathbf{r}_m`$.
If we now take the inverse Fourier transform of this, we can replace the sum $`\sum_{mn}`$ by a sum over $\mathbf{q}$, giving:

```math
\mathcal{H} = \sum_{\mathbf{q}} \sum_{i,j} \mathbf{S}_i(\mathbf{q}) J_{ij}(\mathbf{q}) \mathbf{S}_j(\mathbf{q})
```

where the sum over $`\mathbf{q}`$ extends over both positive and negative vectors within the first Brillouin zone.

Expressing the spin vectors in terms of the boson operators using the Holstein-Primakoff transformation described in the last section
(and noting that because of the complex conjugation the Fourier transform of the creation operator $`b^{\dagger}`$ is
$`b^{\dagger}(-\mathbf{q})`$):

```math
\mathcal{H} = \sum_{\mathbf{q}} \sum_{i,j} \left[
    \sqrt{\frac{S_{i}}{2}}\left( \mathbf{z}_{i}^* b_{i}(\mathbf{q}) + \mathbf{z}_{i}b_{i}^{\dagger}(-\mathbf{q}) \right)
                           + \boldsymbol{\eta}_i \left( S_{i} - b_{i}(\mathbf{q}) b_{i}^{\dagger}(-\mathbf{q}) \right)
    \right]^{\intercal} J_{ij}(\mathbf{q}) \left[
    \sqrt{\frac{S_{j}}{2}}\left( \mathbf{z}_{j}^* b_{j}(\mathbf{q}) + \mathbf{z}_{j}b_{j}^{\dagger}(-\mathbf{q}) \right)
                           + \boldsymbol{\eta}_j \left( S_{j} - b_{j}(\mathbf{q}) b_{j}^{\dagger}(-\mathbf{q}) \right)
    \right]
```

where the operators $`b_i`$, $`b^{\dagger}_i`$ and vectors $`\mathbf{z}_i`$ and $`\boldsymbol{\eta}_i`$ relate to each site $i$
in the magnetic unit cell in real-space, and $`J_{ij}`$ is a tensor represented by a $3 \times 3$ matrix.
(Note that the vector $\mathbf{z}$ is complex and we use the asterisk to denote complex conjugation
and the $\intercal$ symbol to indicate vector or matrix transpose).
In addition, in the above summation we explicitly consider terms $`i=j`$ which corresponds to the on-site anisotropy terms
$`\sum_i \mathbf{S}_i A_i \mathbf{S}_i`$ where $`A_i`$ is also a tensor represented by $`3\times 3`$ matrix.
Thus in this treatment, $`J_{ij}`$ is a superset and includes both the exchange and anisotropy terms.

Noting that the boson operators obey the commutation relation

```math
[b_{i}, b_{j}^{\dagger}] = \delta_{i, j}
```

and that the $`\mathbf{z}`$ and $`\boldsymbol{\eta}`$ vectors are perpendicular, we can rewrite the Hamiltonian as a matrix equation
if we define a column vector $\mathbf{X}$ as composed of the boson operators within the unit cell:

```math
\mathbf{X}_m = \left(b_{m,1}, b_{m,2}, \ldots, b_{m,L}, b^{\dagger}_{m,1}, b^{\dagger}_{m,2}, \ldots, b^{\dagger}_{m,L} \right)
```

where $L$ is the number of sites within the unit cell, giving:

```math
\mathcal{H} = \sum_{\mathbf{q}} \mathbf{X}^{*\intercal}(\mathbf{q}) \mathrm{h}(\mathbf{q}) \mathbf{X}(\mathbf{q})
```

where the hermitian matrix $`\mathrm{h}(\mathbf{q})`$ is

```math
\mathrm{h}(\mathbf{q}) = \left[ \begin{array}{cc}
\mathrm{A}(\mathbf{q}) - \mathrm{C} && \mathrm{B}(\mathbf{q}) \\
\mathrm{B}^{*\intercal}(\mathbf{q}) && \mathrm{A}^*(\mathbf{q}) - \mathrm{C}
\end{array} \right]
```

where the $`(i,j)`$ elements are:

```math
\mathrm{A}(\mathbf{q})_{ij} = \frac{\sqrt{S_i S_j}}{2} \mathbf{z}_i^{\intercal} \mathrm{J}_{ij}(\mathbf{q}) \mathbf{z}_j^*
```
```math
\mathrm{B}(\mathbf{q})_{ij} = \frac{\sqrt{S_i S_j}}{2} \mathbf{z}_i^{\intercal} \mathrm{J}_{ij}(\mathbf{q}) \mathbf{z}_j
```
```math
\mathrm{C}_{ij} = \delta_{ij} \sum_{\mathcal{l}} S_{\mathcal{l}} \boldsymbol{\eta}_i^{\intercal} \mathrm{J}_{i\mathcal{l}}(\mathbf{q}=\mathbf{0}) \boldsymbol{\eta}_{\mathcal{l}}
```

where only terms which are quadratic (two-operator) in the boson operators, e.g. $`b^{\dagger} b`$, have been retained.
This is because the expectation value of single-operator terms (e.g. $`b^{\dagger}`$) vanish and we ignore higher order terms in _linear_ spin wave theory.

We can use the same treatment above for the Zeeman term, $`\mu_B \mathbf{H}^{\intercal} \sum_m \sum_i g_i \mathbf{S}_{im}`$.
The Fourier transform yields $`\mu_B \mathbf{H}^{\intercal} \sum_{\mathbf{q}} \sum_i g_i \mathbf{S}_i(\mathbf{q})`$,
which with the Holstein-Primakoff transformation becomes:

```math
\mathcal{H}_{\mathrm{Zeeman}} = - \sum_{\mathbf{q}} \sum_i \mu_B \mathbf{H}^{\intercal} g_i \boldsymbol{\eta}_i b^{\dagger}_i(\mathbf{q}) b_i(\mathbf{q})
```

where again only two-operator terms have been retained.
This yields an additional matrix element

```math
\mathrm{A}_{\mathrm{Zeeman}}(\mathbf{q})_{ij} = - \frac{1}{2} \mu_B \mathbf{H}^{\intercal} \delta_{ij} g_i \boldsymbol{\eta}_i
```

which should be added to the $`\mathrm{A}(\mathbf{q})_{ij}`$ term.

The eigenvalues of $`\mathrm{h}(\mathbf{q})`$ are the magnon energies and its eigenvectors can be used to calculate the neutron cross-section.

We should note here that while the anisotropy terms $`A_i`$ are often said to add a "constant to the diagonal" of the "Hamiltonian",
this is only true when the $`A_i`$ tensors are aligned with the moment direction $`\boldsymbol{\eta}_i`$ such that only the $`C_{ij}`$ terms contribute.
Let us take a cubic ferromagnetic system with one spin in the unit cell which is aligned along $[0,0,1]$.
In this case, $`\mathbf{z} = [1, j, 0]`$ and $`\boldsymbol{\eta} = [0, 0, 1]`$.
If

```math
A_i = J_{ii} = \left( \begin{array}{ccc} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{array} \right)
```

then in the equation above, the $`\mathrm{A}(\mathbf{q})_{ii}`$ and $`\mathrm{B}(\mathbf{q})_{ij}`$ terms will be zero
because $`A_i\mathbf{z}_i = A_i\mathbf{z}_i^* = 0`$. However, if

```math
A_i = J_{ii} = \left( \begin{array}{ccc} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{array} \right)
```

then $`\mathrm{C}(\mathbf{q})_{ii}`$ will be zero and $`\mathrm{A}(\mathbf{q})_{ii}`$ and $`\mathrm{B}(\mathbf{q})_{ij}`$ will both be non-zero
so we will get non-diagonal terms due to the single-ion anisotropy.

Now, often we speak of $`\mathrm{h}(\mathbf{q})`$ as the "Hamiltonian" but strictly the Hamiltonian (the operator which yields the total energy)
is $`\mathcal{H}`$ which is the sum of $`\mathrm{h}(\mathbf{q})`$ over all $`\mathbf{q}`$ in the Brillouin zone.
Nonetheless, in order to calculate the inelastic neutron spectra, the quantity we need to calculate (and diagonalise) is $`\mathrm{h}(\mathbf{q})`$.


### Code for calculating the "Hamiltonian" matrix

Despite the disclaimer above, in the rest of the text we will use "Hamiltonian matrix" to refer to the hermitian matrix $`\mathrm{h}(\mathbf{q})`$.
As can be seen from the equations for the submatrices $\mathrm{A}(\mathbf{q})$, $\mathrm{B}(\mathbf{q})$, and $\mathrm{C}$ above,
the calculation can naturally be divided into a $`\mathbf{q}`$-independent part
(involving the prefactor $`\frac{\sqrt{S_i S_j}}{2}`$ and the $\mathbf{z}$ and $\boldsymbol{\eta}$ vectors)
and a $`\mathbf{q}`$-dependent part (involving the Fourier transform of the exchange interactions $`\mathrm{J}_{ij}(\mathbf{q})`$ where $i,j$ label sites within the unit cell).

The code to calculate the Hamiltonian is entirely in the [spinwave.m](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%2540spinw/spinwave.m) file.
First, the `intmatrix` method of the `spinw` object is called to compute $`J_{ij}`$ in real-space
(in [line 508](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L508)).
Next the magnetic structure vectors $`\mathbf{z}`$ and $`\boldsymbol{\eta}`$ are calculated from the magnetic propagation vector and magnetic basis
(in [lines 555-590](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L555-L590)).

`spinwave.m` uses Matlab vectorisation (using `bsxfun`) in order to speed up calculations, so it needs to expand these vectors to cover the full basis
(there are $L$ spins in the unit cell but the boson operator basis is $2L$ because it includes both creation and annihilation operators)
using `repmat` and then compute the vector transpose (`zedL` and `etaL`) using `permute` in
[lines 608-612](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L608-L612).

It next computes the $\mathbf{q}$-independent parts of the submatrices $`\mathrm{A}`$ (called `AD0` in the code),
$`\mathrm{B}`$ (called `BC0`), and $`\mathrm{C}`$ (called `A20` [upper left submatrix] and `D20` [lower right submatrix] in the code) in
[lines 614-641](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L614-L641).
These lines also include computing the indices $i$ (`atom1`) and $j$ (`atom2`) which are eventually passed to `accumarray` to actually construct the matrix.

[Lines 856-934](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L614-L641)
contain the calculation of the $\mathbf{q}$-dependent part of the Hamiltonian and forming it into a square matrix.
First the phase factor $\exp(i\boldsymbol{\delta}\cdot\mathbf{q})$ is computed in
[line 871](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L871),
then each submatrix $`\mathrm{A}`$ (split into `A1` for upper left and `D1` for lower right pars) and $`\mathrm{B}`$ is multiplied by this phase factor
in [lines 874-876](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L874-L876) using `bsxfun`.

Finally `accumarray` is used to construct the square matrix and the diagonal submatrix $`\mathrm{C}`$ added in
[lines 898-900](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L898-L900).


### Takin implementation

As an alternative to the vectorized Matlab code, there is an independent implementation of the formulism in Toth & Lake in the
[takin](https://arxiv.org/pdf/1903.02632) code which
[uses a loop over the sites of the unit cell](https://github.com/ILLGrenoble/magpie/blob/master/tlibs2/libs/magdyn/hamilton.h#L175-L264)
which is more transparent with respects to the equations in the paper.

<!--

## The spin-spin correlation function

The logical next step (and the one taken by most treatment) is to diagonalise the "Hamiltonian" matrix $`h(\mathbf{q})`$
to find the magnon energies $`E(\mathbf{q})`$ at each $`\mathbf{q}`$ point.
However,

-->

## The Bogoliubov transformation

The next step is to diagonalise the "Hamiltonian" matrix $`h(\mathbf{q})`$ to find the magnon energies (the eigenvalues of $`h(\mathbf{q})`$,
and more importantly to find the linear combination of boson operators $`b_i(\mathbf{q})`$ corresponding to these energies so as to compute the neutron structure factor.
Because the (new) boson operators (also) have to obey the commutation relation we cannot simply use the standard eigenvalue decomposition algorithms.
Instead we need to compute


Instead SpinW uses one of two algorithms:

* That of [Colpa](http://dx.doi.org/10.1016/0378-4371%2878%2990160-7) if the $`h(\mathbf{q})`$ matrix is explicitly treated as Hermitian.
  Note that this algorithm will give an [error](https://github.com/SpinW/spinw/blob/73f604ab4d84084d872d1f5fdf46dbc54e14cdd7/swfiles/%40spinw/spinwave.m#L993-L996)
  if the $`h(\mathbf{q})`$ matrix is not only Hermitian but also positive definite.
* That of [White et al.](https://doi.org/10.1103/PhysRev.139.A450) for a general $`h(\mathbf{q})`$ matrix.
  Note that in this case SpinW can return (unphysical) imaginary magnon energies.
