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
When we calculate the trace it is effectively applied twice, which
is the same as applying it once. In pySpinW (and the MATLAB version), we calculate $S_\perp$ with
this single application of the projection.

Polarisation
------------




