# Chosing a build system for PySpinW

## Context

We intend to provide a compiled CPython module to compute the core linear spin wave theory
calculations in parallel to run faster than the native Python implementation.
This will require using a build system.
The two most popular systems compatible with the `pyproject.toml` format are:

1. [meson](https://mesonbuild.com/)
2. [cmake](https://cmake.org/) via [scikit-build-core](https://scikit-build-core.readthedocs.io/en/latest/)
3. [cargo](https://doc.rust-lang.org/cargo/) via [maturin](https://www.maturin.rs/)

The advantage of `meson` is that it is itself implemented in Python, has
extensive support for building Python modules, and also natively supports rust.

CMake is a general build system but doesn't natively support build rust projects -
rather there are [extensions](https://github.com/val-ms/cmake-rust-demo)
to wrap rust's `cargo` build system.
CMake does have extensive support for C/C++.

Alternatively we can use a Python build system which directly wraps `cargo` like `maturin`.
However, this would mean we cannot build C/C++ project unless we change the build system.

On the other hand, `meson` directly supports build both rust and C++ modules.

## Decision

We will use `meson` as our build system because it: is implemented in Python and
directly supports building both C/C++ and rust extension modules.

## Status

Accepted

## Advantages

* Implemented in Python
* Directly supports building C/C++ extension modules
* Directly supports building rust extension modules

