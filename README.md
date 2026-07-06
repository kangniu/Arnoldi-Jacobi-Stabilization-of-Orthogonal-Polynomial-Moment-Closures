# Arnoldi-stabilized hyperbolic moment closure numerical cases

This MATLAB package is a first implementation framework for the numerical cases discussed in:

- Morin & McDonald, *Development of globally hyperbolic one-dimensional moment closures based on orthogonal polynomials*, JCP 523 (2025) 113659.
- Vandermonde with Arnoldi / QR-style stabilization idea from the uploaded `vander.pdf` and MATLAB routines.

## What is included

1. `experiment_vandermonde_arnoldi_condition.m`  
   Reproduces a condition-number comparison between monomial Vandermonde matrices and Arnoldi-generated bases.

2. `run_riemann_moment_closures.m`  
   Solves the 1-D Riemann problem from the JCP paper using finite-volume HLL flux.

3. `run_all_cases.m`  
   Runs the condition-number experiment and representative Riemann cases.

4. Closure support:
   - B4 closure with Jacobi/Arnoldi wave-speed computation.
   - B6 closing flux formula included; wave speeds are currently computed by a numerical flux Jacobian.

## Important notes

This is a research-code starting point, not yet a final paper-grade reproduction.

- The JCP paper uses 10,000 cells for the Riemann problem. The default here is smaller (`Nx=600`) for fast testing.
- Direct BGK discrete-velocity reference solution is not yet included.
- Shock-wave steady-state cases are not yet included.
- B6 wave-speed computation should later be replaced by a fully symbolic/Arnoldi Jacobi recurrence implementation.

## Quick start in MATLAB

```matlab
cd ArnoldiMomentClosureCases
run_all_cases
```

For a closer reproduction of the JCP Riemann problem:

```matlab
cfg = default_riemann_config();
cfg.Nx = 10000;
cfg.tEnd = 0.006;
cfg.tau = 1e-3;        % transition regime
cfg.closures = {'B4','B6'};
run_riemann_moment_closures(cfg);
```


## v2 path fix

The output folder is now created by `case_output_dir.m`, using an absolute path
relative to the package directory. This avoids MATLAB `saveas` failures when
the current working directory is not the package directory or when `figures/`
does not already exist.


## v3 stopping criteria and comparison scripts

New stopping criteria in `default_riemann_config.m`:

```matlab
cfg.maxSteps = 20000;
cfg.minDt = 1e-10;
cfg.maxWaveSpeed = 1e6;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
```

The B6 run in v1/v2 can stall because B6 still uses a finite-difference flux
Jacobian. Near non-realizable or near-singular states this numerical Jacobian
can produce enormous wave speeds, forcing the CFL time step to collapse.  This
is exactly why a full Arnoldi/Jacobi implementation for B6/B8 is the next
research task.

New comparison scripts:

```matlab
compare_wave_speeds_B4_original_vs_arnoldi
compare_riemann_B4_original_vs_arnoldi
```

The default `run_all_cases.m` now runs only B4 PDE cases, because B4 has a
clean Arnoldi/Jacobi wave-speed implementation. B6 is kept as an optional
diagnostic with stopping protection.


## v4 next comparison experiments

Run the next comparison suite:

```matlab
run_next_comparison_cases
```

This runs:

1. `experiment_high_order_wave_speed_stability.m`  
   Compares monomial characteristic-polynomial roots with symmetric
   Jacobi eigenvalues for high-order, clustered real spectra.

2. `experiment_quadrature_recovery_vandermonde_vs_arnoldi.m`  
   Recovers quadrature weights from moments and nodes using both a monomial
   Vandermonde system and an Arnoldi-generated basis.

3. `experiment_B6_jacobian_failure_diagnostic.m`  
   Runs the B6 Riemann case with finite-difference flux-Jacobian wave speeds
   and stopping criteria, then plots local wave-speed growth and realizability
   diagnostics.

The first two experiments are designed to isolate the numerical linear
algebra advantage of the Arnoldi/Jacobi realization.  The third experiment
documents why the next implementation stage should replace the B6 numerical
Jacobian by a fully Arnoldi/Jacobi wave-speed calculation.


## v5 corrections

Important correction: the B6 closing flux denominator is now implemented as

```matlab
D = W4 - W3^2 - 1
W6 = -0.5*A/D + 0.5*B/D
```

so that the equilibrium state `(W3,W4,W5)=(0,3,0)` gives `W6=15`.
This is consistent with the B6 construction being chosen to recover
equilibrium.

New B6 prototype files:

```matlab
wave_speeds_B6_jacobi.m
b6_jacobi_coefficients.m
b6_charpoly_coeffs_normalized.m
recurrence_poly_ascending.m
compare_B6_static_original_vs_jacobi.m
```

The B6 Jacobi implementation is a coefficient-matching prototype.  It is
suitable for static tests and diagnostics.  It is not yet optimized for full
PDE runs.

The quadrature recovery experiment has also been corrected.  It now separates:

1. the ill-conditioned monomial system `V' * w = U`,
2. the correct V+A basis formulation `Q' * w = qMom`,
3. an intentionally unsafe post-hoc conversion `qMomBad = C' * U`.

This reflects the central V+A idea: the method achieves high precision by a
basis transformation and by avoiding ill-conditioned linear systems, not by
solving a Vandermonde system and repairing it afterwards.

Run:

```matlab
run_v5_comparison_cases
```


## v6 B6 correction

The B6 closing flux is now derived directly from the recurrence condition

```matlab
b3 = (b1+b2) + 1/2*((a2-a0)^2 + (a2-a1)^2)
```

and solved symbolically for `W6`.  This avoids relying on signs obtained from
visual/PDF parsing of a long formula.  The implemented formula is

```matlab
g  = W3^2 - W4 + 1;
N6 = 3*W3^6 - 10*W3^4*W4 - 2*W3^4 ...
     + 6*W3^3*W5 + 7*W3^2*W4^2 + 8*W3^2*W4 ...
     - W3^2 - 14*W3*W4*W5 - 2*W3*W5 ...
     + 4*W4^3 - 6*W4^2 + 2*W4 + 4*W5^2;
W6 = -0.5*N6/g;
```

The equilibrium check now gives `(W3,W4,W5)=(0,3,0) -> W6=15`, and the
finite-difference Jacobian has real eigenvalues on the diagnostic states.

Run:

```matlab
run_v6_comparison_cases
```


## v7 semi-explicit B6 Jacobi and small-grid PDE comparison

`b6_jacobi_coefficients.m` no longer uses `fminsearch`.  It now recovers
the unknown tail coefficients `a3,a4,a5,b4,b5` semi-explicitly from

```matlab
P6 = A3*P3 - b3*B2*P2
A3 = (lambda-a3)*B2 - b4*(lambda-a5)
B2 = (lambda-a5)*(lambda-a4) - b5
```

The code computes `B2` from the remainder of `P6` modulo `P3`, then obtains
`A3` by polynomial division, and finally recovers the tail Jacobi coefficients
by coefficient comparison.  This implements the requested semi-explicit
recurrence route and avoids nonlinear optimization.

New scripts:

```matlab
run_v7_B6_small_comparison
compare_riemann_B6_small_original_vs_jacobi
run_B6_small_jacobi_only
```

The small-grid B6 comparison uses

```matlab
Nx = 80;
tEnd = 2e-4;
tau = Inf;
```

and compares original finite-difference flux-Jacobian wave speeds with the
V+A/semi-explicit Jacobi wave speeds.


## v8 B6 grid sweep and full-candidate run

New scripts:

```matlab
run_v8_B6_scaling_and_candidate
compare_riemann_B6_grid_sweep
run_B6_jacobi_full_candidate
```

`compare_riemann_B6_grid_sweep` gradually increases grid size and final time:

```matlab
[80,  2e-4]
[160, 2e-4]
[160, 5e-4]
[240, 5e-4]
[240, 1e-3]
```

For each case it compares original finite-difference flux-Jacobian wave speeds
against the V+A/semi-explicit Jacobi wave speeds.

After the sweep passes, a full-candidate Jacobi-only run can be launched by

```matlab
sol = run_B6_jacobi_full_candidate(600,0.006,Inf);
```

or with relaxation regimes, for example

```matlab
sol = run_B6_jacobi_full_candidate(600,0.006,1e-3);
```


## v9 next-stage B6 full benchmark

After the v8 grid/time sweep passes, the recommended next step is a full
B6 V+A/Jacobi regime sweep:

```matlab
T = run_B6_full_regime_sweep(600,0.006);
plot_B6_regime_summary_from_files(600,0.006);
```

This runs the three regimes:

```matlab
tau = Inf;   % free molecular
tau = 1e-3;  % transition
tau = 1e-6;  % continuum
```

and saves per-case MAT files, clean 300-dpi figures, and a summary CSV.

New files:

```matlab
run_v9_B6_full_benchmark_driver
run_B6_full_regime_sweep
run_B6_jacobi_benchmark_case
plot_B6_regime_summary_from_files
compare_B6_original_vs_jacobi_moderate
save_clean_figure
```

For an optional longer original-vs-Jacobi consistency check, run

```matlab
Tcmp = compare_B6_original_vs_jacobi_moderate();
```

`save_clean_figure` uses `exportgraphics` and hides axes toolbar artifacts.


## v10 complete B4/B6/B8 package

The package now contains a complete B4/B6/B8 V+A/Arnoldi--Jacobi workflow:

```matlab
RUN_THIS_FIRST
run_complete_experiment_suite('light')
run_complete_experiment_suite('paper')
run_complete_experiment_suite('full')
```

New in v10:

- `run_complete_experiment_suite.m`: master runner with `light`, `paper`, and `full` modes.
- `run_B4_full_regime_sweep.m`: full B4 V+A/Jacobi sweep.
- `bgk_reference_1d.m`: compact discrete-velocity BGK reference solver.
- `compare_B6_jacobi_with_bgk_reference.m`: compare B6 V+A/Jacobi against BGK.
- `run_B6_vs_BGK_regime_demo.m`: three-regime B6-vs-BGK demo.
- `make_paper_tables_from_results.m`: generate LaTeX tables from CSV summaries.
- `write_latex_table_from_table.m`: simple LaTeX table writer.
- `RUN_THIS_FIRST.m`: command guide.
- `b8_closing_moment_normalized.m`: B8 closing moment generated from the
  hierarchy coefficient rather than from a hand-transcribed formula; it also
  propagates forward sensitivities for the B8 characteristic polynomial.
- `b8_charpoly_coeffs_normalized.m`: B8 characteristic polynomial generated
  from recurrence sensitivities and the normalized-to-raw moment chain rule,
  without centered finite-difference perturbations.
- `b8_jacobi_coefficients.m`: B8 Jacobi-tail recovery using a split
  characteristic polynomial.
- `compare_B8_charpoly_ad_vs_fd.m`: diagnostic comparison between the
  production forward-sensitivity B8 characteristic polynomial and the legacy
  finite-difference polynomial generator.
- `compare_B8_static_original_vs_jacobi.m`,
  `compare_riemann_B8_small_original_vs_jacobi.m`, and
  `compare_riemann_B8_grid_sweep.m`: B8 verification suite.
- `run_original_vs_va_profile_figures.m`: profile plots comparing the original
  finite-difference flux-Jacobian wave-speed implementation with the new
  V+A/Jacobi implementation for the same closure.
- `run_profile_comparison_figures.m`: paper-style profile plots comparing
  B4, B6, B8, and BGK in the free-molecular, transition, and continuum
  regimes.
- `B8_IMPLEMENTATION_NOTE.md`: current B8 implementation status.

Recommended workflow:

```matlab
run_complete_experiment_suite('light')  % quick verification
run_complete_experiment_suite('paper')  % article-level B4/B6 plus B8 verification
```

Optional kinetic reference:

```matlab
run_complete_experiment_suite('full')
```

To regenerate the manuscript profile figures directly:

```matlab
run_original_vs_va_profile_figures()
run_profile_comparison_figures(180,501,0.006)
```

The first command produces implementation-equivalence figures.  In those
plots, the solid curves are the original finite-difference Jacobian method and
the dashed curves are the new V+A/Jacobi method.  The second command produces
model-accuracy figures comparing V+A/Jacobi B4, B6, and B8 against the BGK
reference.

The B8 implementation avoids hand-transcribing the very large B8 formula.  The
closing moment is generated from the recurrence definition of the B8 hierarchy,
the B8 characteristic polynomial is generated by forward sensitivities, and
the Jacobi tail is recovered from a split polynomial identity.  The legacy
finite-difference polynomial routine is kept only as a diagnostic comparator.
