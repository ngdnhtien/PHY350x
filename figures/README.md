# Figures

Every plot produced during the course, exactly as it came out of the notebooks
(the SVG is the notebook's own output; PNG and PDF are conversions of it, PDF
is what the report uses). Nothing was re-plotted or restyled.

Numbering follows the story: `01` = Iterate, `02` = Sample, `03` = Evolve.

| File | What it shows | Produced by |
|---|---|---|
| `01_nonlinear_system_zero_contours` | Zero-level sets of the two equations cos(xy)+sin(x⁴+y⁴)=1 (red) and x²+y²+sin(xy)=3 (blue) on −3<x,y<3; the crossings are the roots. Originally saved as `a.png`. | notebook 01, Problem 1 |
| `02a_decay_mc_vs_analytic` | Survival population of 10⁵ nuclei, T½ = 3.8, Monte-Carlo (blue) vs N₀·2^(−t/T) (orange). | notebook 02, Ch. 1 Topic 2 |
| `02b_decay_counts_poisson` | Histogram of the number of decays per unit time over 10⁵ trials, with the Poisson mass function on top. | notebook 02, Ch. 1 Topic 2 |
| `02c_lcg_uniform_histogram` | 10⁴ numbers from the linear congruential generator a=16807, c=0, m=2³¹−1. | notebook 02, Ch. 2 Topic 1 |
| `02d_box_muller_gaussians` | The two Box–Muller outputs G₁, G₂ built from LCG uniforms, against N(0,1). | notebook 02, Ch. 2 Topic 2 |
| `02e_rejection_target_vs_gaussian_envelope` | The "ugly" target ρ(x)=⅔x^(−1/3) (blue) next to a scaled Gaussian (orange): why rejection sampling needs a cut-off. | notebook 02, Ch. 2, rejection method |
| `02f_rejection_sampling_histogram_fmax15` | First rejection-sampling run (uniform proposal on [10⁻⁵,10], Fmax=15). | notebook 02, Ch. 2, rejection method |
| `02g_maxwell_boltzmann_rejection_sampling` | Rejection-sampled speeds (orange histogram) against the Maxwell–Boltzmann density with a=2 (blue). | notebook 02, Ch. 2, rejection method |
| `02h_bateman_three_nuclei_chain` | Three-isotope Bateman chain, 10⁴ nuclei, T½ = 14.8 and 16.1: pop1 → pop2 → pop3. | notebook 02, Ch. 4 |
| `02i_kr_decay_counts_poisson` | Number of Kr decays in one hour over 2·10⁴ trials, with the Poisson curve. | notebook 02, Ch. 4 |
| `02j_rn211_branching_decay_chain` | Branching chain ²¹¹Rn → ²¹¹At (74 %) / ²⁰⁷Po → ²¹¹Po / ²⁰⁷Bi → ²⁰⁷Pb, 10⁴ nuclei, 200 h. | notebook 02, Ch. 4 |
| `03a_gas_pressure_vs_collisions` | Running pressure estimate of the free-flight gas after each wall collision (it grows linearly: see the report for why). | notebook 03, Topic 4 |
| `03b_heat1d_implicit_explicit_crank_nicolson` | 1-d heat equation, u(x,0)=sin(πx), dx=dt=0.025: backward Euler (left), forward Euler (middle, diverges to 10⁷¹), Crank–Nicolson (right). | notebook 03, Topic 5 |
