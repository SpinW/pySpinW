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
but one, the "half-polarised" case lacks a polariser or analyser, where the "fully-polarised"
case has both.

The polariser and analyser can both be defined by a polarisation direction, and a 
magnitude between zero and one, with zero being unpolarised and one being complete polarisation.
We can therefore specify a polarisation setup in two vectors, the half-polarised case
being a special case of fully-polarised, with the polariser or analyser vector equal to zero.

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

$\mathbf \rho = a_0 \mathbf \sigma_0 + a_1 \mathbf \sigma_0 + a_2 \mathbf \sigma_0 + a_3 \mathbf \sigma_3 = a_0 \mathbf \sigma_0 \cdot \sigma$  

Where $\sigma$ is a vector of $\sigma_x$, $\sigma_y$ and $\sigma_z$.

When representing polarisations, we have a constraint that the trace of the matrix is the intensity, which initially
we will assume is 1. As all the Pauli matrices except $\sigma_0$ have trace 0, this constraint means that

$1 = \mathrm{Tr}(\rho) = a_0 \mathrm{Tr}(\sigma_0) = 2 a_0$

So, the initial state has $a_0 = 1/2$, regardless of polarisation.

The other coefficients, $a_1$, $a_2$ and $a_3$ are the degree of polarisation in each cartesian direction: x, y and z
respectively. These form the vector discussed above describing the polarisation. If the beam
is unpolarised, then we have

$\rho_0 = \frac{1}{2} \mathbf I$

Polarising Elements in Density Matrix Formalism
-----------------------------------------------

Next we'll work through the mathematics of polarising elements, demonstrating how the formalism works a little more.

A polarising element in the density matrix formalism has the form

$\Pi(\vec a) = \frac{1}{2} \left( \mathbf I + \vec{a} \cdot \mathbf \sigma \right)$ 

and we can apply this as matrix to a state to get a new polarisation. For example,
a polariser in $z$ would have a $\Pi$ matrix of

$\Pi(\vec z) = \frac{1}{2} \left(\mathbf I + \sigma_z\right)$

Applying to an unpolarised state would give (note that as the Pauli matrices are all Hermitian,
$\sigma_i^\dagger = \sigma_i$)

$\Pi(\vec z) \rho_0 \Pi(\vec z)^\dagger$ $= \frac{1}{2} \left(\mathbf I + \sigma_z\right) \left(\frac{1}{2} I\right) \frac{1}{2} \left(\mathbf I + \sigma_z\right)$

$= \frac{1}{8} \left(\mathbf I + \sigma_z\right) \left(\mathbf I + \sigma_z\right)$

$= \frac{1}{8} \left( \mathbf I + 2\sigma_z + \sigma_z^2\right)$

The square of the Pauli matrices are the identity, so

$= \frac{1}{8} \left( \mathbf I + 2\sigma_z + \mathbf I\right) = \frac{1}{4}\left( \mathbf I + \sigma_z \right)$

The measured intensity is the trace of this

$\mathrm{Tr}(\Pi(\vec z) \rho_0 \Pi(\vec z)^\dagger) = \frac{1}{2}$

So, complete polarisation removes half the intensity, as expected.

Spin-Spin Correlation in Terms of Density Matrices
==================================================

Technically, we can't write the spin-spin correlation matrix as an operation in the density matrix formalism,
because it is a quadratic term. What we can do, however, is use the components that we called $M$ above, knowing
that they will "square" to the spin-spin correlation, and hope that whatever we calculate turns out to only
include square terms (it will).

We can define a density matrix operator for the vector quantity $M$ as follows:

$F = M\cdot\sigma$

Note, this does not have an identity matrix like the projection vector, this is 
where nuclear scattering would appear. As we're only concerned with magnetic changes,
this nuclear term is absent.

The magnetic selection rule is independent of this, so the maths is the same with or without it,
only we will use $T$ instead of $S$ in the final calculation, and it will factor out.

Combining the Terms
===================

With these matrices established, we can build a chain of components, with $\Pi_p$ and $\Pi_a$ being the 
polarisation operators for the polariser and analyser. $\rho_d$ is density matrix at the detector

$\Sigma = \Pi_a F \Pi_p$

and 

$\rho_d = \Sigma \rho_0 \Sigma^\dagger$

We can expand $\Sigma$

$= \frac{1}{4} \left( \mathbf I + \vec{p_a} \cdot \mathbf \sigma \right) \left( M\cdot\sigma\right) \left( \mathbf I + \vec{p_p} \cdot \mathbf \sigma \right)$

$=\frac{1}{4}\left(\mathbf M\cdot\sigma + (\vec p_a \cdot \sigma) (M\cdot\sigma) + (M\cdot\sigma) (\vec p_p \cdot \sigma) + (\vec p_a \cdot \sigma) (M\cdot\sigma) (\vec p_p \cdot \sigma) \right)$

There's lots of terms here that are products of inner products of $\sigma$, it will be useful to have a rules for
working with them. The most useful is

$(A \cdot \sigma)(B \cdot \sigma) = (A\cdot B) I + i (A \times B) \cdot \sigma$

this can be applied twice to get

$(A \cdot \sigma)(B \cdot \sigma)(C \cdot \sigma) = ((A\cdot B)C + (A\cdot C)B - (B\cdot C)A)\sigma + i ((A\times B)\cdot C)\mathbf I$


so we have

$\Sigma = \frac{1}{4}\left(\mathbf M\cdot\sigma + (\vec p_a\cdot M) \mathbf I + i(\vec p_a\times M) \cdot \sigma + (M \cdot \vec p_p) \mathbf I + i(M \times \vec p_p) \cdot \sigma + ((p_a\cdot M)p_p + (p_a\cdot p_p)M - (M\cdot p_p)p_a)\sigma + i ((p_a\times M)\cdot p_p)\mathbf I\right)$

collecting $\sigma$ and $\mathbf I$ terms, and applying anti-symmetry of the cross product gives

$= \frac{1}{4}\left( (M + ((p_a\cdot M)p_p + (p_a\cdot p_p)M - (M\cdot p_p)p_a) + i (p_a \times M) - i (p_p \times M)) \sigma + (p_a \cdot M + p_p \cdot M + i((p_a \times M)\cdot p_p)) \mathbf I \right)$

using distributivity laws

$= \frac{1}{4}\left( (M + ((p_a\cdot M)p_p + (p_a\cdot p_p)M - (M\cdot p_p)p_a) + i ((p_a - p_p) \times M)) \cdot \sigma + ((p_a + p_p) \cdot M + i((p_a \times M)\cdot p_p)) \mathbf I \right)$

We can introduce some notation to make this simpler, let $U$ be the coefficient of the identity and $V$ be the 
coefficient of $\sigma$, i.e.

$\Sigma = U I + V \cdot \sigma$

then

$U = \frac{1}{4}\left( (M + ((p_a\cdot M)p_p + (p_a\cdot p_p)M - (M\cdot p_p)p_a) + i ((p_a - p_p) \times M)) \right)$

which is a scalar

$V = \frac{1}{4}\left( ((p_a + p_p) \cdot M + i((p_a \times M)\cdot p_p)) \mathbf I \right)$

which is a vector.

When we calculate $\rho_d$, we'll have

$\rho_d = \left(U I + V \cdot \sigma\right) \left(\frac{1}{2} \mathbf I \right) \left(U^* I + V^* \cdot \sigma\right)$

$= \frac{1}{2}(UU^*\mathbf I + UV^*\cdot\sigma + U^*V\cdot\sigma + (V\cdot\sigma)(V^*\cdot\sigma))$

The last term can be simplified

$(V\cdot\sigma)(V^*\cdot\sigma)) = (V\cdot V^*) \mathbf I + i (V \times V^*) \cdot \sigma = |V|^2 \mathbf I + i (V \times V^*) \cdot \sigma$

so

$\rho_d = \frac{1}{2} \left((|U|^2 + |V|^2)\mathbf I + (UV^* + U^*V + i(V\times V^*))\cdot\sigma \right)$










