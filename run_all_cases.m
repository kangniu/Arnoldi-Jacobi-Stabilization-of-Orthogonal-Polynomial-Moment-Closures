clear; clc; close all;

fprintf('============================================================\n');
fprintf('Arnoldi-stabilized moment closure numerical cases\n');
fprintf('============================================================\n\n');

fprintf('[1] Vandermonde/Arnoldi conditioning experiment\n');
experiment_vandermonde_arnoldi_condition;

fprintf('\n[2] B4 wave-speed comparison: original flux Jacobian vs Arnoldi/Jacobi\n');
compare_wave_speeds_B4_original_vs_arnoldi;

fprintf('\n[3] B4 Riemann comparison: original flux Jacobian vs Arnoldi/Jacobi\n');
compare_riemann_B4_original_vs_arnoldi;

fprintf('\n[4] 1-D Riemann problem: representative B4 cases\n');

cfg = default_riemann_config();
cfg.Nx = 600;
cfg.tEnd = 0.006;
cfg.maxSteps = 20000;
cfg.minDt = 1e-10;
cfg.maxWaveSpeed = 1e6;
cfg.closures = {'B4'};

fprintf('\n  Free molecular regime\n');
cfg.tau = Inf;
cfg.caseName = 'free_molecular_B4';
run_riemann_moment_closures(cfg);

fprintf('\n  Transition regime\n');
cfg.tau = 1e-3;
cfg.caseName = 'transition_B4';
run_riemann_moment_closures(cfg);

fprintf('\n  Continuum regime\n');
cfg.tau = 1e-6;
cfg.caseName = 'continuum_B4';
run_riemann_moment_closures(cfg);

fprintf('\nAll default cases finished. Figures are saved in ./figures\n');
fprintf('\nOptional B6 test with stopping protection:\n');
fprintf('  cfg = default_riemann_config(); cfg.closures={''B6''}; cfg.maxWaveSpeed=1e6; run_riemann_moment_closures(cfg);\n');
