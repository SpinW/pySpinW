# Symmetry Checks

## Transformation of exchanges

A spacegroup object in pySpinW contains a full list of operations, these can be considered ordered pairs of 
transformation matrices and translations

$\mathcal{G} = \{(M_1, T_1), (M_2, T_2), ...\}$

These operations apply to points, but when considering the effect of symmetry on exchanges, we need to think about
how this affects exchange contribution to the overall Hamilton. The exchange terms have a form $ S_i J_{ij} S_j^T $
and are mostly non-spatial, that is to say: these terms are not dependent on the position of the spins in space,
but the spins orientation is dependent on the transformation part of a given symmetry operation.

Let's say that exchange terms are related by a symmetry operation with transformation matrix $M$, and let $i$ and $j$ 
be the indices of the original spins, and $a$ and $b$ be the transformed spins. Then we should find that

$S_a J_{ab} S_b^T = S_i J_{ij} S_j^T $

Adopting the convention that the operations multiply on the left, we know that $S_a = M S_j$, so using 
the fact we have $M^{-1} = M^T$ then we have

$S_i J_{ij} S_j^T = (M^T S_a) J_{ab} (M^T S_b)^T = S_a M^T J_{ab} M S_b$

from which we can see that

$J_{ij} = M^T J_{ab} M$

and so, left and right multiplying by $M$ and $M^T$ respectively

$M J_{ij} M^T = M M^T J_{ab} M^T M = J_{ab} $

## Allowed exchanges

Some sets of exchanges will be forbidden by symmetry, but some exchanges are forbidden even without 
considering them collectively.
To find what these are, we first need to work out which symmetries might potentially leave an exchange unchanged.

As the exchange is defined by a pair of points and an exchange matrix, we can find all the symmetries that 
the point pair satisfies, and then ask which matrices are consistent with these symmetries.

The symmetries, $g \in \mathcal{G}$, of a point pair $(s_1, s_2)$ come in the following two forms

1) direct, where $g(s_1) = s_1$ and $g(s_2) = s_2$
2) swapped, where $g(s_1) = s_2$ and $g(s_2) = s_1$

These have to be dealt with differently, because changing the order of $S_i$ and $S_j$ in the equation
above results in a transposition. 
But they give us what we need to work out the constraints on the exchange matrix. 

### Swapping the sites transposes the exchange matrix

This is very easy to show. It comes from noting that Hamilton should be invariant to the site swapping, that is,
for a bilinear exchange

$S_i J_{ij} S_j^T$

The corresponding swapped $J$, call it $K$, should give the same behaviour with the spins swapped, that is

$S_j K_{ij} S_i^T = S_i J_{ij} S_j^T$

so $K_{ij} = J_{ji}$, or just $K = J^T$

## Getting the constraints on the exchange matrix

We now have a way of listing all the relevant symmetries, and equations that the exchange matrix must obey for all
these symmetries:

* For the unswapped symmetries: $J = M J M^T$ 
* For the swapped symmetries: $J = M J^T M^T$

From these we can build a set of equations on the elements of the matrix, find out which can be changed independently,
which are zero, etc.  The procedure can also be applied to the system as a whole, either to check the symmetry, or provide
a symmetry constrained parametrisation.

### Numerical Strategy

A numerical solution can be obtained by Gaussian elimination / row-reduction.
Essentially, we need to "flatten the equation"
the unswapped constraints are calculated using the following vectorisation identity

$\text{vec}(AXB)=(B^T \otimes A)\text{vec}(X)$

where $\text{vec}$ denotes the vectorisation, i.e. such that

$\text{vec} \left(\begin{array}{ccc} a & b & c \\ d & e & f \\ g & h & i \end{array}\right) = \left(\begin{array}{c} a \\ b \\ c \\ d \\ e \\ f \\ g \\ h \\ i  \end{array}\right)$

So, if $V = \text{vec}(J)$ is the vector of matrix components, then we can translate as

$J \to I V$

$M J M^T \to (M \otimes M) V$

And our symmetry constraints are of the form 

$(M \otimes M - I) V = 0$

The swapped constraints need a "commutation matrix", $K$, defined such that

$\text{vec}(P^T) = K \text{vec}(P)$

With this, for the swapped case we get

$((M \otimes M) K - I) V = 0$

To get a reduced system of equations for out matrix entries we build a $9$-by-$9n$
matrix of $M \otimes M - I$ and $(M \otimes M)K - I$ and row-reduce it.

### Interpreting the row reduction

Once we have our matrix in reduce row echelon form, we can interpret the entries.
The reduced matrix will have non-zero entries in at most the first 9 rows, the leftmost
non-zero value on each row should always be 1 (if not all zeros).

There are three basic cases to consider:
1) *The row has a single non-zero entry*, and as such it corresponds to the equation $x_i = 0$, meaning that
   zero is the only solution
2) *The row has multiple non-zero entries.* This means that the entries are free, but are restricted by an equation of 
   the form $x_i + c_0 x_j + c_1 x_k ... = 0$
3) *The row is all zeros.* This entry provides no constraints at all on the system: $0=0$.

We can also think about which entries are completely free, these will be the ones that have no non-zero
entry in a given column.

### In terms of symmetric and antisymmetric components

This will be enough to give us a simplified set of equations for the system, 
but it is not the preferred form for users. It's more understandable to consider 
the symmetric and antisymmetic parts of the matrix separately.

Let's now write $J$ in terms of matrix entries, with the aim of making a simple linear system in the form $Y = MX$

$J = \left( \begin{array}{ccc} a & b & c \\ b & d & e \\ c & e & f \end{array} \right) + \left( \begin{array}{ccc} 0 & z & -y \\ -z & 0 & x \\ y & -x & 0 \end{array} \right)$

We then want an equation to relate the entries to $V$ (dots used instead of zeros for readability)

$V = \left(\begin{array}{ccccccccc}
1 & \cdot & \cdot & \cdot & \cdot & \cdot & \cdot & \cdot & \cdot \\
\cdot & 1 & \cdot & \cdot & \cdot & \cdot & \cdot & \cdot & 1 \\
\cdot & \cdot & 1 & \cdot & \cdot & \cdot & \cdot & -1 & \cdot \\
\cdot & 1 & \cdot & \cdot & \cdot & \cdot & \cdot & \cdot & -1 \\
\cdot & \cdot & \cdot & 1 & \cdot & \cdot & \cdot & \cdot & \cdot \\
\cdot & \cdot & \cdot & \cdot & 1 & \cdot & 1 & \cdot & \cdot \\
\cdot & \cdot & 1 & \cdot & \cdot & \cdot & \cdot & 1 & \cdot \\
\cdot & \cdot & \cdot & \cdot & 1 & \cdot & -1 & \cdot & \cdot \\
\cdot & \cdot & \cdot & \cdot & \cdot & 1 & \cdot & \cdot & \cdot \\
\end{array}\right) \left(\begin{array}{c} a \\ b \\ c \\ d \\ e \\ f \\ x \\ y \\ z  \end{array}\right)$

To get a system of equations in these terms, we can premultiply before doing the reduction.

## Copying Exchanges

### Uniqueness of symmetry based copying 

*There is not necessarily a single operation that maps the points of an exchange, but the exchange matrix should
change in the same way under all of them as long as it is valid according to that symmetry.*

Proof: For a pair of ordered exchange $(s_1, s_2)$ and $(t_1, t_2)$ we find all the 
operations $g\in\mathcal{G}_x$ in the symmetry group $\mathcal{G}$ that maps them together,
i.e.

$\mathcal{G}_x = \{ g\in\mathcal{G} \mid t_1 = g(s_1), t_2 = g(s_2) \}$

Let's say there is more that there are two operation in $\mathcal{G}_x$: $g$ and $h$,
then $g^{-1}$ and $h^{-1}$ must in the original space group. 
Then also, so must $hg^{-1}$ (as well as all other combinations of $g$, $h$, $g^{-1}$ and $h^{-1}$).

Since $h$ maps the point pair $(s_1, s_2)$ to $(t_1, t_2)$, then $h^{-1}$ must map $(t_1, t_2)$ to $(s_1, s_2)$.
This means that $g h^{-1}$ must map $(t_1, t_2)$ to itself.

We know that the exchange matrix needs to be invariant under all operations that send a point pair to itself.
Remember that is not true of any arbitrary matrix, just those that obey the symmetries of the system.

As $gh^{-1}$ leaves the symmetry constrained exchange matrix unchanged, $(gh^{-1}) h$ and $h$ should do the
same thing to it. Then, because the operations $(gh^{-1}) h$ and $g$ are the same, $h$ and $g$ should do
the same thing to the exchange matrix.