Calculating the Spin-Spin Correlation
=====================================

The core calculation in pySpinW calculates the spin-spin correlation matrix $S_{ab}$ as follows

$S_{ab} = \left<\psi_{aj}\mid c_{jk} e^{i q (r_j - r_k)} \mid \psi_{bk}\right>$ 

where $c_{jk}$ are the $q$-independent terms (see 001), 
the complex exponential terms are the $q$-dependent phase factors,
and $\psi$ are the eigenvectors produced during the calculation (i.e.
the wavefunctions describing the magnon modes in terms of 
boson creation and annihilation operators). 

This object is this autocorrelation of the changes in the magnetic 
structure that correspond to a magnetic excitation (magnon) from the ground 
state structure to an excited state. 

Processing the Scattering Matrix 
================================

Neutron Scattering
------------------

Neutrons interact with spins via magnetic fields. 
Whilst the spin-spin correlation matrix captures the changes
in spins that correspond to a magnetic excitation, the
fundamental structure of magnetic fields prevents some of
these from happening. Namely, some correlations between 
spin states do not have a magnetic field associated with them,
and so cannot interact with magnetic dipoles, such as the one
carried by a neutron. In other words, there is a magnetic 
selection rule.

The magnetic selection rule follows directly from Maxwell's
equations, specifically Gauss's Law for Magnetism:

$\nabla \cdot \vec{B} = 0$

When fourier transformed, the first order derivatives
become a product, giving

$\vec{q} \cdot \vec{B}(\vec{q}) = 0$

i.e. the component of the magnetic field parallel to $\vec{q}$
is zero.

So, to get the neutron scattering intensity, we need to
remove this component.

Applying the following transform to a quantity removes the component 
in the $\vec{q}$ direction:

$\mathbf{P}(\vec{q}) = \mathbf{I} - \hat{q}\hat{q}^T$

And we can view the magnetic field as this projection applied to
$M$, the magnetisation density.

$\vec B(\vec q)=\mu_0 \mathbf{P}(\vec q)\vec M(\vec q)$

and the spin-spin correlation matrix as the expectation

$S_{ab} = \left< \vec{M}_a \vec{M}_b^* \right>$

So, to transform $S_{ab}$ we need to do the projection
for each rank of the tensor. Call the projected spin-spin
correlation $T_{ab}$.

Switching to Einstein notation now,

$\mathbf T = \mathbf P \mathbf S \mathbf P$

The neutron scattering intensity proportional to the trace of this,
we call this $S_\perp$.

$S_\perp = \mathrm{Tr}(\mathbf T) = \mathrm{Tr}(\mathbf P \mathbf S \mathbf P) = \mathrm{Tr}(\mathbf P \mathbf S)$

The last equality follows because $\mathbf P$ is a projection, 
and thus idempotent. 
When we calculate the trace it is effectively applied twice, which, because of idempotency,
is the same as applying it once. In pySpinW (and the MATLAB version), we calculate $S_\perp$ with
this single application of the projection.

Polarisation
------------

There are two common polarising devices in a polarised neutron scatting experiment, 
a "polariser" which selects the spin direction of the neutrons before scattering
by the sample, and an "analyser", which selects them after scattering.

There are also two types of experiment commonly conducted, both use an analyser,
but one, the "half-polarised" case lacks a polariser, where the "fully-polarised"
case has both.

The polariser and analyser can both be defined by a polarisation direction, and a 
magnitude between zero and one, with zero being unpolarised and one being complete polarisation.
We can therefore specify a polarisation setup in two vectors, the half-polarised case
being a special case of fully-polarised, with the polariser vector equal to zero.

In this section we will work towards an expression for the polarised neutron scattering intensity
in the following form

$S_\text{polarised} = \sum_{a,b} T_{ab} Z_{ab}$

where $T_{ab}$ is the projected spin-spin correlation defined above, and $Z_{ab}$ is a "polarisation tensor"
that we will now derive.

Density Matrix Formalism
------------------------

The mathematics of polarisation is much easier in a density matrix formalism. In this formalism,
the polarisation state is represented as a linear combination of 4 Pauli matrices.

$\mathbf\sigma_0 = \left(\begin{array}{cc} 1 & 0 \\ 0 & 1\end{array}\right) = \mathbf I$

$\mathbf\sigma_1 = \left(\begin{array}{cc} 0 & 1 \\ 1 & 0\end{array}\right) = \mathbf \sigma_x$

$\mathbf\sigma_2 = \left(\begin{array}{cc} 0 & -i \\ i & 0\end{array}\right) = \mathbf \sigma_y$

$\mathbf\sigma_3 = \left(\begin{array}{cc} 1 & 0 \\ 0 & -1\end{array}\right) = \mathbf \sigma_z$

This representation requires that the matrix is Hermitian (note, no imaginary numbers on the diagonal 
and from $\sigma_1$ and $\sigma_2$ off diagonals are necessarily complex conjugates).

We have 

$\mathbf \rho = a_0 \mathbf \sigma_0 + a_1 \mathbf \sigma_0 + a_2 \mathbf \sigma_0 + a_3 \mathbf \sigma_3$  

When representing polarisations, we have a constraint that the matrix representing the polarisation state has
trace of 1. As all the Pauli matrices except $\sigma_0$ have trace 0, this constraint means that

$1 = \mathrm{Tr}(\rho) = a_0 \mathrm{Tr}(\sigma_0) = 2 a_0$

so $a_0 = 1/2$ all the time.

The other coefficients, $a_1$, $a_2$ and $a_3$ are the degree of polarisation in each cartesian direction, x,y and z
respectively.


