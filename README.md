# comp-phys

Code from an undergraduate computational-physics course, kept as memorabilia:
it is most probably the last code I wrote by hand, line by line.

Read in order, the exercises tell one story — three ways of teaching a computer to do physics:

| Act | Notebook | What happens |
|---|---|---|
| **Iterate** | [`notebooks/01_iterate_roots_and_eigenvalues.ipynb`](notebooks/01_iterate_roots_and_eigenvalues.ipynb) | Roots of a nonlinear system by Newton–Raphson (with a home-made Gauss–Seidel), eigenvalues of a 5×5 matrix by QR iteration — and an eigenvector attempt that converges to zero. |
| **Sample** | [`notebooks/02_sample_monte_carlo_and_decay_chains.ipynb`](notebooks/02_sample_monte_carlo_and_decay_chains.ipynb) | π, integrals and radioactive decay from random numbers; where random numbers come from (LCG, Box–Muller, rejection sampling); error estimates; Bateman decay chains and a branching chain ²¹¹Rn → … → ²⁰⁷Pb. |
| **Evolve** | [`notebooks/03_evolve_gas_and_heat_equation.ipynb`](notebooks/03_evolve_gas_and_heat_equation.ipynb) | A box of molecules followed collision by collision (whose pressure runs away), and the 1-d heat equation by forward Euler, backward Euler and Crank–Nicolson (where forward Euler reaches 10⁷¹). |

The physics write-up, in the form of a short PRA-style article, is [`report/main.pdf`](report/main.pdf) (source: `report/main.tex`).

## Layout

```
notebooks/   the three acts — the prose is the original, word for word; the notebooks call the package below
compphys/    the function bodies, lifted verbatim from the notebooks into one module per topic
figures/     every plot the course produced, named by content (see figures/README.md)
literature/  the primary reference for each method, with DOIs, and the BibTeX file
report/      the article
misc/        the six original files, untouched, including three stand-alone midterm scripts
```

## Running it

The notebooks already contain their original outputs, so nothing needs to be run to read them. To re-run:

```
pip install -r requirements.txt
jupyter lab notebooks/
```

The notebooks find `compphys/` on their own (no installation step). The first notebook uses `qiskit` only to pretty-print a matrix; without it, it falls back to plain `numpy` output.

## License

MIT — see [`LICENSE`](LICENSE).
