import numpy as np
from scipy.linalg import ldl, cholesky
import time

from random_matrices import generate_random_positive_semidefinite_matrix, generate_random_hermitian_matrix

def benchmark(matrix_size: int, sufficient_time=1.0, max_matrices=100):
    matrices = [generate_random_hermitian_matrix(matrix_size) for _ in range(max_matrices)]

    #
    # Cholesky
    #

    total_time = 0
    for i, matrix in enumerate(matrices):
        tic = time.time()

        sqrt = cholesky(matrix)

        toc = time.time()

        total_time += toc - tic

        if total_time > sufficient_time:
            print("Cholesky test completed on sufficient time")
            break

    else:
        print("Chol test completed on max number")

    mean_time_chol = total_time / (i+1)

    #
    # LDL
    #

    total_time = 0
    for i, matrix in enumerate(matrices):
        tic = time.time()

        l, d, perm = ldl(matrix)
        sqrt = l @ np.sqrt(d)

        toc = time.time()

        total_time += toc - tic

        if total_time > sufficient_time:

            print("LDL test completed on sufficient time")
            break

    else:
        print("LDL test completed on max number")

    mean_time_ldl = total_time / (i + 1)

    #
    # Cholesky Error (changes matrix values)
    #

    row = matrix_size // 2

    total_time = 0
    for i, matrix in enumerate(matrices):

        matrix[:, row] = 0.0 # Make matrix not positive definite

        tic = time.time()

        try:
            sqrt = cholesky(matrix)
        except np.linalg.LinAlgError as e:
            pass

        toc = time.time()

        total_time += toc - tic

        if total_time > sufficient_time:
            print("Cholesky error test completed on sufficient time")
            break

    else:
        print("Cholesky error test completed on max number")

    mean_time_chol_error = total_time / (i+1)


    return mean_time_chol, mean_time_ldl, mean_time_chol_error

def run_benchmarks():

    import matplotlib.pyplot as plt

    sizes = [int(x) for x in 2 ** np.linspace(1, 11, 40)]
    chols = []
    ldls = []
    chol_errors = []

    for size in sizes:
        print(size)

        chol_time, ldl_time, err_time = benchmark(size, max_matrices=200, sufficient_time=10.0)

        chols.append(chol_time)
        ldls.append(ldl_time)
        chol_errors.append(err_time)

    plt.loglog(sizes, chols, label="Cholesky")
    plt.loglog(sizes, ldls, label="LDL")
    plt.loglog(sizes, chol_errors, label="Cholesky Positive Definite Error Time")

    plt.xlabel("Matrix size")
    plt.ylabel("Mean time (s)")
    plt.legend()

    plt.show()


if __name__ == "__main__":
    run_benchmarks()
