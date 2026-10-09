"""
Part II, Chapters 1 and 4 -- radioactive decay: one mother / one daughter,
Poisson statistics, the three-nuclei Bateman chain, and a branching chain
(211Rn -> 211At / 207Po -> 211Po / 207Bi -> 207Pb).

Taken from ``final1.ipynb``.  Function bodies are verbatim.  Reorganization
changes:
  * ``np.math.factorial`` (removed in NumPy 2) is read through ``math``;
  * the top-level branching-chain script became ``branching_chain(...)``
    whose arguments are the constants that headed that cell, and which
    returns the time axis and the six populations in the order
    Rn211, At211, Po207, Po211, Bi207, Pb207.
"""

import math
import numpy as np


def decay_func(N0: int, T: float, t: float):
    '''
        Nuclear decay function
    '''
    p = np.log(2)/T

    return N0*np.exp(-p*t)


def decay_mc(N0: int, T: float, multiplier: int):
    '''
        Simulating the decay process (1 mother, 1 daughter)

        Args:
            N0 - Initial population
            T - decay period

        Returns:
            Survival population
    '''
    p = np.log(2)/T # decay probability
    N = [True for i in range(N0)]

    T_end = multiplier*T
    time_simulation = np.arange(1, T_end, 1)

    survival_pop = []
    for t in time_simulation:
        r = np.random.uniform(0, 1, N0)
        N = N&(r>p)
        survival_pop.append(sum(N))

    return time_simulation, survival_pop


def poisson_dist(k, mu):
    rho = (np.exp(-mu)*mu**k)/(math.factorial(k))
    return rho


def three_nuclei_chain(N0: int, lambda_1, lambda_2, end_time):
    '''
        Simulating three nuclei chain decay

    '''

    nucleus_1 = [True for i in range(N0)]
    nucleus_2 = [False for i in range(N0)]
    nucleus_3 = [False for i in range(N0)]

    time_simulation = np.arange(1, end_time, 1)

    pop1 = np.zeros(len(time_simulation))
    pop2 = np.zeros(len(time_simulation))
    pop3 = np.zeros(len(time_simulation))

    for t in time_simulation:
        r1 = np.random.uniform(0, 1, N0)
        nucleus_1_decay = nucleus_1 & (r1 <= lambda_1)
        nucleus_1 = nucleus_1 & (r1 > lambda_1)

        r2 = np.random.uniform(0, 1, N0)
        nucleus_2_decay = nucleus_2 & (r2 <= lambda_2)
        nucleus_2 = (nucleus_2 & (r2 > lambda_2)) | nucleus_1_decay

        nucleus_3 = nucleus_3 | nucleus_2_decay

        pop1[t-1] = sum(nucleus_1)
        pop2[t-1] = sum(nucleus_2)
        pop3[t-1] = sum(nucleus_3)

    return time_simulation, pop1, pop2, pop3


def branching_chain(N0=10000, pRn_At=0.74,
                    lambda_Rn211=np.log(2)/15, lambda_At211=np.log(2)/7.2,
                    lambda_Po207=np.log(2)/5.7, lambda_Po211=0.99999,
                    lambda_Bi207=0.00001, dt=1, tmax=200):
    Rn211 = [True for i in range(N0)]
    At211 = [False for i in range(N0)]
    Po207 = [False for i in range(N0)]
    Po211 = [False for i in range(N0)]
    Bi207 = [False for i in range(N0)]
    Pb207 = [False for i in range(N0)]

    time_simulation = np.arange(1, tmax, dt)

    pop_Rn211 = np.zeros(len(time_simulation))
    pop_At211 = np.zeros(len(time_simulation))
    pop_Po207 = np.zeros(len(time_simulation))
    pop_Po211 = np.zeros(len(time_simulation))
    pop_Bi207 = np.zeros(len(time_simulation))
    pop_Pb207 = np.zeros(len(time_simulation))

    for t in time_simulation:

        r1 = np.random.uniform(0, 1, N0)
        Rn211_decay = Rn211 & (r1 <= lambda_Rn211)
        Rn211 = Rn211 & (r1 > lambda_Rn211)
        pop_Rn211[t-1] = sum(Rn211)

        p_Rn_At_at_t = np.random.uniform(0, 1)
        r2 = np.random.uniform(0, 1, N0)

        if p_Rn_At_at_t < pRn_At:
            At211_decay = At211 & (r2 <= lambda_At211)
            At211 = (At211 & (r2 > lambda_At211)) | Rn211_decay

            Po207_decay = Po207 & (r2 <= lambda_Po207)
            Po207 = (Po207 & (r2 > lambda_Po207))

        else:
            At211_decay = At211 & (r2 <= lambda_At211)
            At211 = (At211 & (r2 > lambda_At211))

            Po207_decay = Po207 & (r2 <= lambda_Po207)
            Po207 = (Po207 & (r2 > lambda_Po207)) | Rn211_decay

        pop_At211[t-1] = sum(At211)
        pop_Po207[t-1] = sum(Po207)

        r3 = np.random.uniform(0, 1, N0)
        Po211_decay = Po211 & (r3 <= lambda_Po211)
        Po211 = (Po211 & (r3 > lambda_Po211)) | At211_decay
        pop_Po211[t-1] = sum(Po211)

        r4 = np.random.uniform(0, 1, N0)
        Bi207_decay = Bi207 & (r4 <= lambda_Bi207)
        Bi207 = (Bi207 & (r4 > lambda_Bi207)) | Po207_decay
        pop_Bi207[t-1] = sum(Bi207)

        Pb207 = Pb207 | Bi207_decay | Po211_decay
        pop_Pb207[t-1] = sum(Pb207)

    return time_simulation, pop_Rn211, pop_At211, pop_Po207, pop_Po211, pop_Bi207, pop_Pb207
