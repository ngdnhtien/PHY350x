"""
Part II, Chapter 1 -- the Monte-Carlo family: pi estimation and MC integration.

Taken from ``final1.ipynb``.  Function bodies are verbatim.  Reorganization
change: in the notebook the integrators read the integrand ``f`` from the
notebook's global namespace (and ``f`` was redefined before each one); here
``f`` is passed explicitly as the first argument.
"""

import numpy as np


def pi_est(N0: int, style: str):
    '''
        Estimation of the Pi constant π = 3.1415926535...

        Args:
            N0: number of sample points
            style: could be 'hm' or 'sp'

        Return:
            Pi estimation
    '''
    xrand = np.random.uniform(0, 1, N0)
    yrand = np.random.uniform(0, 1, N0)
    if style == 'hm':
        hit = 0
        for i in range(N0):
            if (xrand[i]**2+yrand[i]**2) < 1:
                hit += 1

        return (hit/N0)*4

    if style == 'sp':
        f_mean = np.mean([np.sqrt(1-x**2) for x in xrand])

        return f_mean*4


def MCintegral1D(f, style: str, N: int, xmin, xmax, ymax):
    '''
        Args:

        Returns:
    '''

    xrand = np.random.uniform(xmin, xmax, N)

    if style == 'hm':
        yrand = np.random.uniform(0, ymax, N)
        hit_pos, hit_both = 0, 0
        for n in range(N):
            if 0 < yrand[n] < f(xrand[n]):
                hit_pos += 1
            if 0 < yrand[n] < np.abs(f(xrand[n])):
                hit_both += 1
        int_positive = (hit_pos/N)*ymax*(xmax-xmin)
        int_both = (hit_both/N)*ymax*(xmax-xmin)
        int_negative = int_both - int_positive
        int_final = int_positive - int_negative

        return int_final

    if style == 'sp':
        fmean = 0
        for x in xrand:
            fmean += f(x)
        fmean = fmean/N
        int_final = fmean*(xmax-xmin)

        return int_final


def MCintegral2D(f, style: str, N: int, xmin, xmax, ymin, ymax, zmax):
    '''
        Args
        Returns

    '''

    xrand = np.random.uniform(xmin, xmax, N)
    yrand = np.random.uniform(ymin, ymax, N)
    zrand = np.random.uniform(0, zmax, N)

    if style == 'sm':
        fmean = 0
        for n in range(N):
            fmean += f(xrand[n], yrand[n])
        fmean = fmean/N
        int_final = fmean*(xmax-xmin)*(ymax-ymin)

        return int_final

    if style == 'hm':
        hits = 0
        for n in range(N):
            if zrand[n] < f(xrand[n], yrand[n]):
                hits += 1

        int_final = (hits/N)*(xmax-xmin)*(ymax-ymin)*(zmax-0)

        return int_final

    if style == 'sp':
        fmean = 0
        for n in range(N):
            if (xrand[n]-0.5)**2+yrand[n]**2 - 1 < 0:
                fmean += f(xrand[n], yrand[n])
        fmean = fmean/N
        int_final = fmean*(xmax-xmin)*(ymax-ymin)

        return int_final


def MCintegral3D(f, style: str, N: int, xmin, xmax, ymin, ymax, zmin, zmax):
    '''
        Args
        Returns

    '''

    xrand = np.random.uniform(xmin, xmax, N)
    yrand = np.random.uniform(ymin, ymax, N)
    zrand = np.random.uniform(zmin, zmax, N)

    if style == 'sm':
        fmean = 0
        for n in range(N):
            fmean += f(xrand[n], yrand[n], zrand[n])
        int_final = (fmean/N)*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)

        return int_final

    if style == 'sp':
        fmean = 0
        for n in range(N):
            if (xrand[n]**2+yrand[n]**2 <= 1) and (-1 <= zrand[n] <= 1):
                fmean += f(xrand[n], yrand[n], zrand[n])
        int_final = (fmean/N)*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)

        return int_final
