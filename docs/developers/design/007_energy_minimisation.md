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

Next we can look at $\beta$, where we have something very similar

$\frac{d\beta}{dS_k} = \sum_{i} \delta_{ik} (B_i^T + B_i) S_i = \sum_k (B_k^T + B_k) S_k$

where $\delta$ is the Kronecker delta. Finally, the $\gamma$ term is linear, and has the constant value

$\frac{d\gamma}{dS_k} = \sum_i \delta_{ik} S_i \cdot C_i = C_k$






## Challenges

### Non-unique ground state

### Gimbal Lock



## Coordinate Systems

### Unconstrained

### Constrained to plane







