clear; clc; close all;

fprintf('============================================================\n');
fprintf('Next comparison experiments for Arnoldi moment closures\n');
fprintf('============================================================\n\n');

fprintf('[1] High-order wave-speed stability test\n');
experiment_high_order_wave_speed_stability;

fprintf('\n[2] Quadrature weight recovery: Vandermonde vs Arnoldi\n');
experiment_quadrature_recovery_vandermonde_vs_arnoldi;

fprintf('\n[3] B6 finite-difference Jacobian diagnostic with stopping criteria\n');
experiment_B6_jacobian_failure_diagnostic;

fprintf('\nNext comparison experiments finished. Figures are saved in ./figures\n');
