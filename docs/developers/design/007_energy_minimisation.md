# Classical Energy Minimisation

## Problem

We will think of the Hamiltonian in form, which is defined in terms of the mathematical form, not the physical meaning.
There are mixed bilinear terms (B), quadratic terms (A), and linear terms (C).
From a physical perspective, the bilinear terms are couplings between different sites, the quadratic terms
are self-couplings (i.e. same site, but maybe different cell) and anisotropies, and the
linear terms describe the effect of the magnetic field.

$H = \sum_{i \neq j} A_{ij} S_i^T S_j + \sum_i S_i^T B_i S_i + \sum_i S_i \cdot C_i $

### Preliminaries

It will be useful to consider this in parts,

$\alpha = \sum_{i \neq j} A_{ij} S_i^T S_j$

$\beta = \sum_i S_i^T B_i S_i$ 

$\gamma = \sum_i S_i \cdot C_i$

so, $H = \alpha + \beta + \gamma$

We will need the derivative of this with respect to the cartesian coordinates of the spins, $S = (S^x, S^y, S^z)$. We have

$\frac{dH}{dS_k} = \frac{d\alpha}{dS_k} + \frac{d\beta}{dS_k} + \frac{d\gamma}{dS_k}$

we can work through this term by term, firstly, the terms in $\alpha$ are all linear with respect to $S_i$ as $i \neq j$, but their order matters.

$\frac{d\alpha}{dS_k} = \sum_{ij} (S_i^T A_{ik})^T + A_{kj} S_j = \sum_{ij} A_{ik}^T S_i + A_{kj} S_j$

as $\alpha$ is linear with respect to $S_i$, this is constant with respect to $S_i$.

Next we can look at $\beta$, where we have something very similar, but note that it is linear, rather than constant in $S_k$

$\frac{d\beta}{dS_k} = \sum_{i} \delta_{ik} (B_i^T + B_i) S_i = \sum_k (B_k^T + B_k) S_k$

where $\delta$ is the Kronecker delta. Finally, the $\gamma$ term is linear, and has the constant value

$\frac{d\gamma}{dS_k} = \sum_i \delta_{ik} S_i \cdot C_i = C_k$






## Challenges

### Non-unique ground state

The lowest energy state will not be, in general, unique. For example, consider two spins with a ferromagnetic Heisenberg coupling. They will be at an energy minimum when the spins are aligned. 
However, there will be two degrees of freedom with which they can achieve this, it doesn't matter which direction they face.

There's not much we can do about this

### Gimbal Lock

Another challenge is that global parameterisations on the surface of the sphere, or circle, will be discontinuous and potentially singular.

We can avoid this in the optimisation by exploiting the fact that we will only want small rotations. This means we can pick a coordinate system at each
step, work out how to move the moments in that system, then pick a new one in the next step. 



## Coordinate Systems

### Unconstrained

### Constrained to plane








