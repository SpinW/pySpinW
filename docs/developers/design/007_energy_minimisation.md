# Classical Energy Minimisation

## Problem

We will consider the Hamiltonian which includes bilinear terms (B), quadratic terms (A), and linear terms (C).
Physically, the bilinear terms are couplings between different sites, the quadratic (single-ion anisotropy) terms
are self-couplings (i.e. same site, but maybe different cell), and the linear terms describe the effect of an
applied  magnetic field.

$H = \sum_{i \neq j} A_{ij} S_i^T S_j + \sum_i S_i^T B_i S_i + \sum_i S_i \cdot C_i $

### Preliminaries

It will be useful to consider this in parts,

$\alpha = \sum_{i \neq j} A_{ij} S_i^T S_j$

$\beta = \sum_i S_i^T B_i S_i$ 

$\gamma = \sum_i S_i \cdot C_i$

so, $H = \alpha + \beta + \gamma$

We will need the derivative of this with respect to the cartesian coordinates of the spins, $S = (S^x, S^y, S^z)$. We have

$\frac{\partial H}{\partial S_k} = \frac{\partial \alpha}{\partial S_k} + \frac{\partial \beta}{\partial S_k} + \frac{\partial \gamma}{\partial S_k}$

we can work through this term by term, firstly, the terms in $\alpha$ are all linear with respect to $S_i$ as $i \neq j$, but their order matters.

$\frac{\partial \alpha}{\partial S_k} = \sum_{ij} (S_i^T A_{ik})^T + A_{kj} S_j = \sum_{ij} A_{ik}^T S_i + A_{kj} S_j$

as $\alpha$ is linear with respect to $S_i$, this is constant with respect to $S_i$.

Next we can look at $\beta$, where we have something very similar, but note that it is linear, rather than constant in $S_k$

$\frac{\partial \beta}{\partial S_k} = \sum_{i} \delta_{ik} (B_i^T + B_i) S_i = \sum_k (B_k^T + B_k) S_k$

where $\delta$ is the Kronecker delta. Finally, the $\gamma$ term is linear, and has the constant value

$\frac{\partial \gamma}{\partial S_k} = \sum_i \delta_{ik} S_i \cdot C_i = C_k$






## Challenges

### Non-unique ground state

The lowest energy state will not be, in general, unique. For example, consider two spins with a ferromagnetic Heisenberg coupling. They will be at an energy minimum when the spins are aligned. 
However, there will be two degrees of freedom with which they can achieve this, it doesn't matter which direction they face.

There's not much we can do about this, however, it will be a bigger problem in stocastic optimisation and it should not affect the spinwave calculation in a way other than changing the coordinates.

This motivates being able to fix the orientation of moments during the optimisation.

### Gimbal Lock

Another challenge is that global parameterisations on the surface of the sphere, or circle, will be discontinuous and potentially singular.

We can avoid this in the optimisation by exploiting the fact that we will only want small rotations. This means we can pick a coordinate system at each
step, work out how to move the moments in that system, then pick a new one in the next step. 



## Coordinate Systems

We want a coordinate system that is defined locally around the current state.

Let $R$ be the a rotation from $z = [0, 0, 1]$ to the current moment, i.e. $S = R z$

### Unconstrained

For unconstrained systems we want to rotate around z first, we can choose a rotation about the x and y axes, resulting in a parameterised rotation of:

$m = R R_x(\alpha) R_y(\beta) z = R [sin(\beta), -sin(\alpha) cos(\beta), cos(\alpha) cos(\beta)]$

To work out the generalised force applied to this, we'll need the derivative of it with respect to $\alpha$ and $\beta$ at the starting point.

$\frac{\partial S}{\partial \alpha}\vert_{\alpha=\beta=0} = R [0,-1,0] $

$\frac{\partial S}{\partial \beta}\vert_{\alpha=\beta=0} = R [1,0,0] $

The generalised force on each site is then just $\frac{dH}{dS_k} \frac{dS_k}{d\alpha_k}$, and similarly for $\beta$. We only need to consider the contributions from couplings/anisotropies/fields that involve the site under consideration.


### Constrained to plane

Note: We will actually constrain rotations to rotations about an axis, so that a moment that is out of plane will not be forced into the plane. Rotations will happen around the axis perpendicular to it. 
With randomisation, however, the moment will be forced into the plane.

When constrained to a plane, things are a bit easier as we don't have to worry about gimbal lock. It will suffice to just rotate about the specified axis.

The matrix for rotation about an axis $v$ by angle $\theta$ can be written as

$R(v, \theta) = (1 - \cos\theta) v \otimes v + cos\theta + (\sin\theta)\epsilon(-v)$

where $\epsilon$ is the "tripple product matrix" such that $X^T \epsilon(Z) Y = (X \times Y) \cdot Z$.
The derivative of this with respect to $\theta$ at $\theta=0$ is simply

$\frac{\partial R}{\partial \theta} \vert_{\theta=0} = \epsilon(-v)$

We therefore have the following for the change of the spin with respect to $\theta$:
$\frac{\partial S}{\partial \theta} = \frac{\partial}{\partial \theta} R S = \frac{\partial R}{\partial \theta} S = \epsilon(-v) S$

# Supercells

Supercells introduce an added complexity. For the supercells in spinW, the spin at a given position in the cell is some function of "partial spins," $T_1 ... T_n$. The spin in cell indexed by $a$

$S = f(T_1 ... T_n)$

Each of these $T_i$ components can be rotated independently, with their magnitudes preserved, as such, each component needs to be parameterised separately.

Other than this, the the minimisation proceeds as before, but with

$\frac{\partial S}{\partial \alpha_k} = \sum_{i=0}^n \frac{\partial f}{\partial T_i} \cdot \frac{\partial T_i}{\partial \alpha_k}  = \frac{\partial f}{\partial T_k} \frac{\partial T_k}{\partial \alpha_k} $

## Specific Supercells

The derivative $\frac{\partial f}{\partial T_k}$ depends on the choice of supercell, in all cases of the supercells implemented in pySpinW it is a linear map 

$\frac{\\partial f}{\partial T_k}: \mathbb{R}^3 \to \mathbb{R}^3$

but sometimes it is a scalar, and sometimes a matrix. For efficiency we can code for either of these options, for uniformity we can multiply scalar values by the identity matrix.

### Trivial supercell

`TrivialSupercell`s are the identity and the derivative is also the identity. They only have one input component

### Transformation supercells

`TransformationSupercell`s are a linear map on one component, where $f(T) = M T$, where $M$ is matrix that depends on the which cell is selected in the supercell. These also only have one input component.

### Summation supercells

`SummationSupercell`s have multiple components, which are defined by vectors $T_i$ with propagation vectors $p_i$ and phases $\varphi_i$, such that

$f(T_1...T_n) = \sum_i T_i exp(2 \pi i p_i \cdot r + \varphi_i)$

Where $r$ is the position of the cell within the unit cell. The derivative can be expressed by a scalar

$\frac{\partial f}{\partial T_k} = exp(2\pi i p_i \cdot r + \varphi_i)$

### Incommensurate supercells

Incommensurate supercells, such as `RotationSupercell` cannot be used in this optimisation.







