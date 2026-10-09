# misc

Files that are not part of the story as told by `notebooks/` — kept so that
nothing is lost.

`original_files/` holds the six files exactly as they came out of the course,
byte for byte:

| File | What it is | Where it went |
|---|---|---|
| `midterm.ipynb` | Midterm: nonlinear roots + eigenvalues | → `notebooks/01_…`, `compphys/nonlinear.py`, `compphys/eigen.py` |
| `final1.ipynb` | Final, part 1: Monte-Carlo chapters 1–4 | → `notebooks/02_…`, `compphys/montecarlo.py`, `random_generators.py`, `decay.py` |
| `final2.ipynb` | Final, part 2: molecular simulation + heat equation | → `notebooks/03_…`, `compphys/gas.py`, `compphys/heat.py` |
| `midterm_prob1a.py` | Stand-alone script: contour plot of Problem 1 | duplicate of the notebook cells; its `np.cos`/`np.sin` form is the one used in `compphys/nonlinear.py` |
| `midterm_prob1b.py` | Stand-alone script: Jacobian, Gauss–Seidel, Newton–Raphson | duplicate of the notebook cells |
| `midterm_prob2a.py` | Stand-alone script: QR iteration + eigenvalues | duplicate of the notebook cells |

Note that `midterm.ipynb` as originally saved does not run top-to-bottom (its
contour cell uses `cos`/`sin` before `numpy` aliases are defined, and it uses
the `%matplotlib notebook` backend); the reorganized notebook 01 does.
