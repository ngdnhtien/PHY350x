"""
Part III, Topic 4 -- molecular simulation: a box of free-flying molecules,
pressure read off the momentum they hand to the walls.

Taken from ``final2.ipynb``.  The body is verbatim.  Reorganization change:
the top-level script became ``simulate_gas(...)`` whose arguments are the
constants that headed that cell; it returns the collision times ``tc`` and
``pressure_after_n_collision``.  The ``print`` inside the loop is kept, as
in the original.
"""

import numpy as np

pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt


def simulate_gas(L=1, Num=200, Temp=300, mass=6.7e-27, kB=1.38e-23, NumCol=100):
    vth = np.sqrt(2*kB*Temp/mass)
    tcollnt = 2*sqrt(2)*L/vth
    delV = 0
    Tn = 0

    x = np.random.uniform(-L, L, Num)
    y = np.random.uniform(-L, L, Num)
    z = np.random.uniform(-L, L, Num)

    phi = np.random.uniform(-pi, pi, Num)
    theta = np.random.uniform(0, pi, Num)

    Vx = vth*sin(theta)*cos(phi)
    Vy = vth*sin(theta)*sin(phi)
    Vz = vth*cos(phi)

    tx = np.zeros(Num)
    ty = np.zeros(Num)
    tz = np.zeros(Num)
    tc = np.zeros(NumCol)

    pressure_after_n_collision = np.zeros(NumCol)

    for n in range(NumCol):
        tcol = tcollnt
        for i in range(Num):
            tx[i] = abs((x[i]-L*np.sign(Vx[i]))/Vx[i])
            if tx[i] < tcol:
                tcol = tx[i].copy()
                mol = i
            ty[i] = abs((y[i]-L*np.sign(Vy[i]))/Vy[i])
            if ty[i] < tcol:
                tcol = ty[i].copy()
                mol = i
            tz[i] = abs((z[i]-L*np.sign(Vz[i]))/Vz[i])
            if tz[i] < tcol:
                tcol = tz[i].copy()
                mol = i

        # found smallest time needed for one collision to happen (regardless of direction)
        Tn += tcol
        tc[n] = Tn

        if tx[mol] <= ty[mol]:
            if tx[mol] <= tz[mol]:
                Vx[mol] = -Vx[mol]
                delV += 2*abs(Vx[mol])
            else:
                Vz[mol] = -Vz[mol]
                delV += 2*abs(Vz[mol])
        else:
            Vy[mol] = -Vy[mol]
            delV += 2*abs(Vy[mol])

        pressure_at_n = mass*delV/tcol
        x += Vx*tcol
        y += Vy*tcol
        z += Vz*tcol

        pressure_after_n_collision[n] = mass*delV/(8*L*Tn)
        print(pressure_at_n)

    return tc, pressure_after_n_collision
