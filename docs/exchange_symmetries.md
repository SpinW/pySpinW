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

$S_i J_{ij} S_j^T $

The corresponding swapped $J$, call it $K$, should give the same behaviour with the spins swapped, that is

$S_j K_{ij} S_i^T = S_i J_{ij} S_j^T$

so $K_{ij} = J_{ji}$, or just $K = J^T$

## Getting the constraints on the exchange matrix

We now have a way of listing all the relevant symmetries, and equations that the exchange matrix must obey for all
these symmetries:

* For the unswapped symmetries: $J = M J M^T$ 
* For the swapped symmetries: $J = M J^T M^T$

From these we can build a set of equations on the elements of the matrix, find out which can be changed independently,
which are zero, etc.

The procedure can also be applied to the system as a whole, either to check the symmetry, or provide
a symmetry constrained parametrisation.

