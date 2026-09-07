# Anisotropy Ellipsoids for the Viewer 

## Requirements
We want easy-axis matrices, e.g.

$$\left(\begin{array}{ccc} 
-1 & 0 & 0 \\ 
0 & 0 & 0 \\ 
0 & 0 & 0 \\
\end{array}\right)$$

to be long thin ellipses, and easy-plane matrices, e.g.

$$\left(\begin{array}{ccc} 
1 & 0 & 0 \\ 
0 & 0 & 0 \\ 
0 & 0 & 0 \\
\end{array}\right)$$

 to be flat disks.

There's different ways of doing this, but the following seems like the simplest and most stable.

## Adding values

We can add a multiple of the identity to an anisotropy term without affecting the physics.

Let $A$ be an anisotropy matrix, then we have

$H = SAS^T$

Now consider the same with $A$ replaced with $A+kI$, then

$H = S(A + kI)S^T = SAS^T + kSS^T = SAS^T + k|S|^2$

Therefor the energy is changed by a constant term, and so the systems physics remains the same.

## Solution

If we add the smallest eigenvector multiplied by the identity then we assure that the resulting matrix is 
positive-semidefinite, and then we can add a small number to make it slightly bigger for rendering, i.e. 

$$M_\text{rendering} = -(M + (\lambda_\text{min} + \epsilon) I)$$ 

where $\lambda_\text{min}$ is the smallest (most negative) eigenvector of $M$. 

## Scaling

There are various alternatives for scaling, probably the best choice is to find the minimum eigenvector of
$M_\text{rendering}$, which will be $\lambda_\text{max} - \lambda_\text{min}$, then to use that as midpoint of
the slider. This will need to be done across all anisotropies.

## Rendering Matrix

We can use Cholesky decomposition to get a transformation matrix that transforms the sphere model into the
appropriate ellipsoid.