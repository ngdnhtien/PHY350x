"""
compphys -- the routines behind an undergraduate computational-physics course.

The function bodies in this package were lifted verbatim (docstrings and
comments included) out of three Jupyter notebooks written during the course.
Only what was strictly needed to make them importable was changed; every such
change is listed at the top of the module it lives in.

Part I   -- Iterate : nonlinear.py, eigen.py
Part II  -- Sample  : montecarlo.py, random_generators.py, decay.py
Part III -- Evolve  : gas.py, heat.py
"""

from . import nonlinear, eigen, montecarlo, random_generators, decay, gas, heat

__all__ = ["nonlinear", "eigen", "montecarlo", "random_generators", "decay", "gas", "heat"]
