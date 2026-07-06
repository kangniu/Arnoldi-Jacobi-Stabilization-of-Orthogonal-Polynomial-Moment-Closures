clear; clc; close all;

fprintf('============================================================\n');
fprintf('V6 comparison experiments: corrected V+A basis formulation\n');
fprintf('============================================================\n\n');

fprintf('[1] Corrected quadrature recovery: V+A as basis transformation\n');
experiment_quadrature_recovery_vandermonde_vs_arnoldi;

fprintf('\n[2] Verify corrected B6 formula\n');
verify_B6_formula_equilibrium_and_hyperbolicity;

fprintf('\n[3] B6 static wave-speed comparison with prototype Jacobi implementation\n');
compare_B6_static_original_vs_jacobi;

fprintf('\n[4] Optional corrected B6 diagnostic using finite-difference Jacobian\n');
fprintf('    Run experiment_B6_jacobian_failure_diagnostic manually if needed.\n');

fprintf('\nV6 comparison experiments finished. Figures are saved in ./figures\n');
