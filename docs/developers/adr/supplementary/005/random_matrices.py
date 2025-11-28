import numpy as np

def generate_random_hermitian_matrix(size: int, seed: int | None=None):
    rng = np.random.default_rng(seed)
    m_real = rng.random((size, size))
    m_imag = rng.random((size, size))

    m = m_real + 1j*m_imag

    # might be unnecessary to do both of these, but we need the multiplication
    m = m @ m.conj().T
    m = 0.5*(m + m.conj().T)

    return m

def generate_random_positive_semidefinite_matrix(size: int, seed: int | None=None):
    m = generate_random_hermitian_matrix(size, seed)
    m[0, :] = 0.0
    return m
