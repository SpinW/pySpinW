# Chosing a compiled language for the core spin wave computation

## Context

We intend to provide a compiled CPython module to compute the core linear spin wave theory
calculations in parallel to run faster than the native Python implementation.

The are several options for which language to for this purpose:

1. C/C++
2. Rust

Other languages (Fortran, Go) were discounted as they are not as well supported
for building Python modules, and we have less expertise in these languages in the group.
The main advantage of C/C++ is that CPython offers native support for modules written
in C, and the PyBind bindings for C++ is almost as fast.
On the other hand, C/C++ is difficult to write parallel programs in,
and naiive implementations can lead to hard-to-fix bugs.
The architecture of Rust, in contrast, means that many common bugs cannot occur.
Furthermore, Rust was designed with parallel processing in mind and as such
is easier to write such applications in than C++.

## Decision

We will use Rust as the language for the compiled module running the core
linear spin wave theory calculations.

However, pure-Python computation routines will also be retained despite being
slower in order to provide an understandable reference for the computation.

## Status

Accepted

## Advantages

* Rust is faster than C++ and has better parallel processing support
