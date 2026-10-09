"""
Part III, Topic 5 -- the 1-d heat equation (convection-diffusion with U = 0)
by forward Euler, backward Euler and Crank-Nicolson.

Taken from ``final2.ipynb``.  Function bodies are verbatim; ``u_t0`` is the
initial condition the solvers read from the notebook namespace and now lives
in the same module.
"""

import numpy as np


def u_t0(x):
    return np.sin(np.pi*x)


def heat1d_solver_forward(a, b, dx, tmax, dt, kappa):
    '''
        1-dimensional heat equation solver
        ODE solver: forward Euler

        Args:
            u0: initial condition (for every x, at t=0)
            dx: discretization of space
            dt: discretization of time
            t_max: final time
            kappa: dissipation constant

        Returns:
            [u(t)]: numerical solution
    '''

    x = np.arange(a, b+dx, dx)
    nx = len(x)

    t = np.arange(0, tmax+dt, dt)
    nt = len(t)

    u = np.zeros((nx, nt)) # index 0 = time, index 1 = space

    # boundary condition
    u[0, :] = 0
    u[-1, :] = 0

    # initial condition
    u[:, 0] = u_t0(x)

    gamma = dt/dx**2

    for j in range(1, nt): # loop from 1, since u(t=0) is given, toward the end
        for i in range(1, nx-1): # loop from 2 point to the previous point before last point
            u[i][j] = gamma*u[i+1][j-1] + (1-2*gamma)*u[i][j-1] + gamma*u[i-1][j-1]

    return u


def heat1d_solver_backward(a, b, dx, tmax, dt, kappa):
    '''
        1-dimensional heat equation solver
        ODE solver: backward Euler

        Args:
            u0: initial condition (for every x, at t=0)
            dx: discretization of space
            dt: discretization of time
            t_max: final time
            kappa: dissipation constant

        Returns:
            [u(t)]: numerical solution
    '''

    x = np.arange(a, b+dx, dx)
    nx = len(x)

    t = np.arange(0, tmax+dt, dt)
    nt = len(t)

    u = np.zeros((nx, nt)) # index 0 = time, index 1 = space

    # boundary condition
    u[0, :] = 0
    u[-1, :] = 0

    # initial condition
    u[:, 0] = u_t0(x)

    gamma = dt/dx**2

    matA = np.diag((u.shape[0]-2)*[1+2*gamma], 0) + np.diag((u.shape[0]-3)*[-gamma], -1) + np.diag((u.shape[0]-3)*[-gamma], 1)

    for j in range(1, nt):
        b = u[1:-1, j-1].copy()
        b[0] = b[0] + gamma*u[0, j]
        b[-1] = b[-1] + gamma*u[-1, j]
        sol = np.linalg.solve(matA, b)
        u[1:-1, j] = sol

    return u


def heat_1d_cranknicol(a, b, dx, tmax, dt, kappa):
    '''
        1-dimensional heat equation solver
        ODE solver: Crank-Nicolson method, averaging second order derivative

        Args:
            u0: initial condition (for every x, at t=0)
            dx: discretization of space
            dt: discretization of time
            t_max: final time
            kappa: dissipation constant

        Returns:
            [u(t)]: numerical solution
    '''

    x = np.arange(a, b+dx, dx)
    nx = len(x)

    t = np.arange(0, tmax+dt, dt)
    nt = len(t)

    u = np.zeros((nx, nt)) # index 0 = time, index 1 = space

    # boundary condition
    u[0, :] = 0
    u[-1, :] = 0

    # initial condition
    u[:, 0] = u_t0(x)

    gamma = dt/dx**2

    matA = np.diag([2*(1+gamma)]*(u.shape[0]-2), 0) + np.diag([-gamma]*(u.shape[0]-3), -1) + np.diag([-gamma]*(u.shape[0]-3), 1)
    matB = np.diag([2*(1-gamma)]*(u.shape[0]-2), 0) + np.diag([gamma]*(u.shape[0]-3), -1) + np.diag([gamma]*(u.shape[0]-3), 1)

    for j in range(nt-1):
        b = u[1:-1, j].copy()
        b = np.dot(matB, b)
        b[0] = b[0] + gamma*(u[0, j]+u[0, j+1])
        b[-1] = b[-1] + gamma*(u[-1, j]+u[-1, j+1])
        solution = np.linalg.solve(matA, b)
        u[1:-1,j+1] = solution

    return u
