# Classes for different kinds of magnetic structures

## Context

Diagonalising the bosonic Hamiltonian requires a matrix decomposition that finds $M$ such that

$$MM^\dagger = H$$

The common method is Cholesky decomposition, which requires the matrix to be positive definite. 
This is not always the case in spinW.
In spinW, a small diagonal contribution is added to the input of the Cholesky decomposition to assure the
matrix is positive definite, but this is a hack.

Another alternative, one which does not have the same requirements for positive definiteness, is the LDL decomposition,
in which case we find matrices $L$ and $D$ such that

$$LDL^\dagger = H$$

$D$ is diagonal (and real when $H$ is Hermitian), so we can find $M$ by
$$M = L\sqrt{D}$$

Potentially, this might make the diagonalisation easier too.

LDL is, however, more expensive, here's a plot...

![benchmark_curves.png](supplementary%2F005%2Fbenchmark_curves.png)

Both LDL and Cholesky are approximately linear for large matrices (it's quite likely 
they're actually $n \log n$ or something similar), but Cholesky is about 3 times as fast. The time to
detect failure in Cholesky is more like 5 or 6 times faster than the successful LDL.

It makes sense to try Cholesky first, and if it fails, use LDL.

## Decision

Use Cholesky decomposition and wait for errors, if it errors though LinAlgError, 
run LDL instead.

## Status

Awaiting proper testing

## Advantages

* More rigourous