Representation of Superstructures
=================================

There are a number of different ways we can represent magnetic structures exist at a scale beyond the unit cell.

Firstly, we can distinguish between commensurate and incommensurate structures, i.e. those that have a period that
matches up with (can be fully represented with an integer multiple of) the lattice spacing, and those that are not.

Secondly, there is a question of how we represent the moment, for which there are two obvious approaches. 
The first of these two approaches is the usual way of representing them in the literature, where we write
a given moment $m(r)$ in cell $r$ as a sum

$$m(r) = Re \{ \sum_i \varphi_i \exp(2\pi i k \cdot r) \} $$

Lets consider an example from a common use-case: rotation.
If we have a spin that rotates around $(0,0,1)$ with a period of 5 unit cells in the direction (1,0,0) we could write this,
for $x$ th cell the $(1,0,0)$ direction, starting in that direction, as

$$m(x) = Re \{ (1, i, 0) \exp (2\pi i x / 5) \}$$

or we could, if it is preferable, write it in purely real terms as

$$= Re \{ (1, 0, 0) \exp (2\pi i x / 5) \} + Im \{ (0, 1, 0) \exp (2\pi i x / 5) \} = (1,0,0) \cos (2\pi x / 5) + (0,1,0) \sin (2\pi x / 5) $$

There are some details to consider here, specifically, the magnitude of the calculated moment. In our example, we have a moment that is constant

$$|m| = \sqrt{m_x^2 + m_y^2 + m_z^2} = \sqrt{\cos^2(2\pi x / 5) + \sin^2(2\pi x / 5) + 0} = 1$$

However, this is not always the case, it could change in magnitude, for example, if it were the case that $$\phi = (1,0,0)$$, then the moment would be
smaller for $x=1$ than $x=0$, as $|m| = |\cos(2\pi x /5)|$.

We can say then, that rotations are special case of the general setup, and that there are some very restricted conditions. 
In this form, rotations (and reflections, and improper reflections) have the property that

$$ |m| = \text{constant} $$

If we consider a general form for $\phi$:

$$ m(r) = Re \{ (a + bi, c + di, e + fi) \exp(2\pi r \cdot k) \} = (a, c, e) \cos(2\pi r \cdot k) + (b, d, f) \sin(2\pi r \cdot k) $$
$$ |m|^2 = (a^2 + c^2 + e^2) \cos^2(2\pi r \cdot k) + (b^2 + d^2 + f^2) \sin^2(2\pi r \cdot k) + 2(ab + cd + ef) \sin(2\pi r \cdot k) \cos(2\pi r \cdot k) $$

so, to lose the dependence on $r \cdot k$ we require that

$$a^2 + c^2 + e^2 = b^2 + c^2 + d^2 $$

i.e., the magnitude of the real part and imaginary part are equal, in other words
$$|Re\{\phi\}| = |Im\{\phi\}|$$
and then we also have
$$ab + cd + ef = 0$$
which is the requirement the real and imaginary parts are orthogonal, i.e.
$$Re\{\phi\}\cdot Im\{\phi\} = 0$$

which is kind of obvious, I guess. But it is not so easy to assure in practice. 

Multiple $k$ Vectors
--------------------

When we have multiple $\phi_i$ or $k_i$s this specification becomes more difficult.
Adding two or more rotating vectors do not, in general, have a conserved magnitide.

For the specific case of rotations, it makes sense to have a different description.


Alternative Specification
=========================

There is a different way that we could possibly write the $m$ values

$$m(r) = \left(\prod_i M_i^{k \cdot r}\right)m_0$$

where $M$ is a matrix that performs some geometric transformation,
in the case of rotations/reflections it will be an orthogonal matrix.
$m_0$ is a "starting vector". The starting vector makes more sense for
building supercells based on a unit cell.

For example, in the case of the 5 fold rotations we would have

$$m_0 = (1,0,0)$$

and $M$ is the rotation matrix around $z$ by one fifth of a rotation

$$M = \left( \begin{array}{ccc} 
\cos(2\pi/5) & \sin(2\pi/5) & 0 \\
-\sin(2\pi/5) & \cos(2\pi/5) & 0 \\
0 & 0 & 1
\end{array} \right)$$

There are two important things to note about this approach. One is that there
are more parameters in this representation, the other is that the order in which $M$
operations are applied affects the outcome.

If we restrict ourselves to orthogonal matrices, then the degrees of freedom are the same 
as the other representation (i.e. 6 for each $\phi$ or $M$ and 3 for each $k$) except for 
the starting vector introduces an extra three.

Can we convert between the two representations?
===============================================

Yes, however, because of the dependence on $m_0$ this might not be ideal.

