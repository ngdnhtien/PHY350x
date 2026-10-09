"""
Part II, Chapter 2 -- random number generators for Monte-Carlo simulation:
linear congruential generator, Box-Muller, rejection sampling.

Taken from ``final1.ipynb``.  Function bodies are verbatim (``normal_dist``
was defined inside a plotting cell in the notebook; it now lives here).
"""

import numpy as np


def LCGinator(a, c, m, N):
    '''
        Args:
            a: multiplier
            c: increment
            m: modulus

        Returns:

    '''

    x = np.zeros(N)
    x[0] = int(N*np.random.uniform(0, 1))

    for n in range(1, N):
        x[n] = (a*x[n-1]+c) % m

    x = np.array(x)/max(x)

    return x


def boxmuller_inator(N: int, U1, U2):
    '''
        Args:
            N: how many random numbers do you want?

        Returns:

    '''

    G1 = np.zeros(N)
    G2 = np.zeros(N)

    for n in range(N):
        G1[n] = np.sqrt(-2*np.log(U1[n]))*np.sin(2*np.pi*U2[n])
        G2[n] = np.sqrt(-2*np.log(U1[n]))*np.cos(2*np.pi*U2[n])

    return G1, G2


def normal_dist(x, mu, sigma):
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/2*sigma**2)


def rho(x):
    if x > 0.1:
        return (2/3)*x**(-1/3)
    else:
        return 15


def maxwell_boltzmann(a, v):
    return (np.sqrt(2/np.pi)*v**2)*np.exp(-v**2/(2*a**2))/a**3


def sampler_ginator(a, b, N: int, sample_from: str, Fmax):
    '''
        Args:
            N:
            U: uniform sampling
            G: Gaussian sampling
            Fmax: the maximum value of the target function rho in the range
        Return:

    '''

    if sample_from == 'uniform':
        count = 0
        x_rho = []
        while count < N:
            xrand = np.random.uniform(a, b, N)
            yrand = np.random.uniform(0, Fmax, N)
            for n in range(N):
                if yrand[n] <= maxwell_boltzmann(a=2, v=xrand[n]):
                    count += 1
                    x_rho.append(xrand[n])

        return x_rho
