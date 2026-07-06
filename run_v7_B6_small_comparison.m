clear; clc; close all;

fprintf('============================================================\n');
fprintf('V7: small-grid B6 Riemann comparison and semi-explicit Jacobi\n');
fprintf('============================================================\n\n');

fprintf('[1] Verify B6 formula and semi-explicit Jacobi recurrence\n');
verify_B6_formula_equilibrium_and_hyperbolicity;
compare_B6_static_original_vs_jacobi;

fprintf('\n[2] Small-grid B6 Riemann comparison: Nx=80, tEnd=2e-4\n');
compare_riemann_B6_small_original_vs_jacobi;

fprintf('\nV7 B6 small comparison finished. Figures are saved in ./figures\n');
