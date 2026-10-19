Higher Order Terms in the Hamiltonian
=====================================

Begin with the linearised spin operator

$\bar{S}_i = \sqrt{\frac{S_i}{2}} \left(z_i^* b_i + z_i b^\dagger \right) + \eta_i\left(S_i - b_i b^\dagger_i\right)$

the overall strategy is then to express the Hamiltonian in the form

$H = \sum_q X^\dagger h X$

where $X$ is a vector of boson operators

$X = \left(b_1, b_2 ... b_N, b^\dagger_1, b^\dagger_2 ... b^\dagger_N\right)$

In the Heisenberg case, we get

$`\bar{S}_i \cdot \bar{S}_j = J \left[\sqrt{\frac{S_i}{2}} \left(z_i^* b_i + z_i b^\dagger \right) + \eta_i\left(S_i - b_i b^\dagger_i\right)\right] \left[\sqrt{\frac{S_i}{2}} \left(z_i^* b_i + z_i b^\dagger \right) + \eta_i\left(S_i - b_i b^\dagger_i\right)\right]`$

From which we get rid of anything that is not quadratic in the boson operators: constant terms are ignored because we only care about energy up to an offset, linear terms are ignored because we are considering perturbations around the ground state (we're at the bottom of an energy well), higher order terms are discarded as an approximation. We're going to square this again, so we cannot yet get rid of the lower order terms

$`= \frac{\sqrt{S_i S_j}}{2} \left(z_i^* b_i + z_i b_i^\dagger\right)^\dagger \left(z_j^* b_j + z_j b_j^\dagger\right)`$

$`+ S_i \sqrt{S_j/2} \eta^\dagger_i \left(z_j^* b_j + z_j b_j^\dagger\right)`$

$`+ S_j \sqrt{S_i/2} \frac{\sqrt{S_i S_j}}{2} \left(z_i^* b_i + z_i b_i^\dagger\right)^\dagger \eta_j`$

$`+ S_i S_j\eta_i^\dagger \eta_j - S_j \eta_i^\dagger b_i^\dagger b_i - S_i \eta_i^\dagger \eta_j b_i^\dagger b_i`$
