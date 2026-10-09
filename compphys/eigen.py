"""
Part I, Problem 2 -- eigenvalues of a 5x5 matrix by unshifted QR iteration.

Taken from ``midterm.ipynb`` / ``midterm_prob2a.py``.  Function bodies are
verbatim.  Reorganization changes:
  * the top-level QR loop became ``qr_iteration(matA, max_iteration=100)``
    and returns the (quasi-)triangular matrix plus the iteration count;
  * the top-level block that reads 2x2 blocks off the result became
    ``eigenvalues_from_blocks(matA)``;
  * the (unsuccessful) eigenvector attempt became ``eigenvector_attempt``;
    the right-hand side handed to ``Gauss_Seidel`` is flattened
    (``zeros_vec[:, 0]``) because NumPy 2 no longer accepts a 1-element
    array where a scalar is assigned.  Results are unchanged.
"""

import numpy as np
import cmath

from .nonlinear import Gauss_Seidel


def numSUM(mat):
    """
    Function for counting the number of zeros element in a matrix `mat`
    Tolerance is 1e-4.
    """
    zeros = 0
    for i in range(mat.shape[0]):
        for j in range(mat.shape[0]):
            if abs(mat[i,j]) < 1e-4:
                zeros += 1
    return zeros


def qr_iteration(matA, max_iteration=100):
    matA = np.array(matA, dtype=float)
    sizeA = matA.shape[0] # size matrix, 5

    iteration_counter = 0
    delta_sum_zero = 0
    sum_zero = 0
    Flag = False

    while Flag == False and iteration_counter < max_iteration:
        Q, R = np.linalg.qr(matA) # using built-in function. I hope we're allowed to do this
        matA = R@Q # Flip QR => RQ to for similarity transformation
        iteration_counter += 1
        # below are just stopping criterion, described in the report
        if iteration_counter % sizeA == 0:
            delta_sum_zero = numSUM(matA) - sum_zero
            if delta_sum_zero == 0:
                Flag = True
        sum_zero = numSUM(matA)

    return matA, iteration_counter


#solving an ordinary polynomials of order 2

def characteristics_polynomial(mat):
    a = 1.0
    b = -(mat[0][0]+mat[1][1])
    c = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
    x = (-b + cmath.sqrt(b**2-4*a*c))/2
    return x


def eigenvalues_from_blocks(matA):
    sizeA = matA.shape[0]

    # finding immediate non-zero elements below diagonals
    non_zeros = []
    zeros = []

    for i in range(sizeA-1):
        if abs(matA[i+1][i]) > 1e-3:
            non_zeros.append(i)
    submats = []

    # construction of 2x2 sub-matrices
    for index in non_zeros:
        submats.append(np.array([
            [matA[index][index], matA[index][index+1]],
            [matA[index+1][index], matA[index+1][index+1]]
        ]))

    eigvals = []
    eigvals.append(matA[0][0])

    for mat in submats:
        eigval = characteristics_polynomial(mat)
        eigvals.append(eigval)
        eigvals.append(np.conjugate(eigval))

    return eigvals


# An attempt to find eigenvectors using Gaussian elimination, but unsuccessful

def eigenvector_attempt(matA, eigval):
    matA = np.array(matA, dtype=float)
    sizeA = matA.shape[0]

    A_tilde = matA

    for i in range(sizeA):
        A_tilde[i][i] = matA[i][i] - eigval

    zeros_vec = np.zeros([sizeA, 1])
    eigen_vec1 = Gauss_Seidel(A_tilde, zeros_vec[:, 0], tol=1e-8)
    eigen_vec2 = np.linalg.solve(A_tilde, zeros_vec)

    return eigen_vec1, eigen_vec2, A_tilde
