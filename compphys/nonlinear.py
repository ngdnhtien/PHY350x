"""
Part I, Problem 1 -- a 2x2 system of nonlinear equations.

    cos(xy) + sin(x^4 + y^4) = 1
    x^2 + y^2 + sin(xy)      = 3

Taken from ``midterm.ipynb`` / ``midterm_prob1a.py`` / ``midterm_prob1b.py``.
Function bodies are verbatim.  Reorganization changes:
  * ``plot_zero_contours`` wraps the contour script of ``midterm_prob1a.py``
    (which already used ``np.cos`` / ``np.sin``) in a function.
  * ``jacobian_symbolic`` wraps the sympy block of ``midterm_prob1b.py``.
"""

from sympy import symbols, diff
from sympy import cos as sym_cos, sin as sym_sin

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import random
import copy

cos = np.cos
sin = np.sin
exp = np.exp


def plot_zero_contours(savefig=None):
    x = np.arange(-3.0, 3.0, 0.01)
    y = np.arange(-3.0, 3.0, 0.01)
    X, Y = np.meshgrid(x, y)
    Z1 = np.cos(X*Y) + np.sin(X**4 + Y**4) - 1
    Z2 = X**2 + Y**2 + np.sin(X*Y) - 3
    v = [0, 0]

    fig, ax = plt.subplots()
    CS1 = ax.contour(X, Y, Z1, 0, colors='red')
    CS2 = ax.contour(X, Y, Z2, 0, colors='blue')
    ax.clabel(CS1, inline=True, fontsize=10)
    ax.clabel(CS2, inline=True, fontsize=10)
    ax.set_title(r'Visualizing solutions in $x>-3, y<3$')
    if savefig is not None:
        fig.savefig(savefig, dpi=300)
    return fig, ax


def jacobian_symbolic():
    x, y = symbols('x y')

    f = sym_cos(x*y) + sym_sin(x**4 + y**4) - 1
    g = x**2 + y**2 + sym_sin(x*y) - 3

    f_x = diff(f, x)
    f_y = diff(f, y)
    g_x = diff(g, x)
    g_y = diff(g, y)

    return f_x, f_y, g_x, g_y


# The jacobian
def J(x):
    return np.array([
        [4*x[0]**3*cos(x[0]**4 + x[1]**4) - x[1]*sin(x[0]*x[1]), -x[0]*sin(x[0]*x[1]) + 4*x[1]**3*cos(x[0]**4 + x[1]**4)],
        [2*x[0] + x[1]*cos(x[0]*x[1]), x[0]*cos(x[0]*x[1]) + 2*x[1]]
    ])

# F(x)
def F(x):
    return np.array([cos(x[0]*x[1]) + sin(x[0]**4 + x[1]**4) - 1, x[0]**2 + x[1]**2 + sin(x[0]*x[1]) - 3])

def Gauss_Seidel(A, b, tol):
    """
    Standard Gauss-Seidel method to solve Ax=b
    """

    m = A.shape[0] # A is an n x n matrix
    n = A.shape[1]
    if (m!=n):
        print(r'Matrix $A$ is not square!')
        return

    # initialize x and x_new
    x = np.zeros(n)
    x_n = np.zeros(n)

    # counter for number of iterations
    # norm_x is the conventional Euclidean norm
    # for stopping criterion: ||x-x_n||^2
    iteration_counter = 0
    norm_x = 1

    while (abs(norm_x) > tol) and (iteration_counter < 100):
        for i in range(n):
            x_n[i] = b[i]/A[i,i]
            sum = 0
            # Standard Gauss-Seidel update
            for j in range(n):
                if (j<i): sum+=A[i,j]*x_n[j]
                if (j>i): sum+=A[i,j]*x[j]
            x_n[i] -= sum/A[i,i]
        norm_x = np.linalg.norm(x-x_n, ord=2)
        # update guess x
        for i in range(n):
            x[i]=x_n[i]
        iteration_counter += 1

    return x

def Newton_Raphson(F, J, x, tol):
    """
    Standard Newton-Raphson method to solve a system of
    nonlinear equation. For simplicity we pre-calculated
    the Jacobian matrix J (hence this function is native
    to the above presented system of nonlinear equations)
    """

    Fval = F(x)
    norm_f = np.linalg.norm(Fval, ord=2)
    iteration_counter = 0

    while (abs(norm_f) > tol) and (iteration_counter < 100):
        print(f'{iteration_counter}-th iteration!')
        # Delta = Gauss_Seidel(J(x), -Fval, tol) # My own Gauss-Seidel
        Delta = np.linalg.solve(J(x), -Fval) # Scipy built-in solver
        x = x + Delta # Update x
        Fval = F(x) # Re-compute F
        norm_f = np.linalg.norm(Fval, ord=2) # Euclidean norm of F is used for stopping criterion.
        iteration_counter += 1               # If X is sufficiently close to the exact solution, then F(X) ~ 0
        abs_error = (np.linalg.norm(Delta, ord=2))/(np.linalg.norm(x, ord=2)) # Absolute error
        print(f'Absolute error is {abs_error}!')

    return x, iteration_counter
