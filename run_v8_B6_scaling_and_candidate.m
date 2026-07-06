clear; clc; close all;

fprintf('============================================================\n');
fprintf('V8: B6 grid sweep and V+A/Jacobi full-candidate run\n');
fprintf('============================================================\n\n');

fprintf('[1] B6 grid/time sweep: original Jacobian vs V+A/Jacobi\n');
T = compare_riemann_B6_grid_sweep;
disp(T);

fprintf('\n[2] Optional full-candidate B6 V+A/Jacobi run\n');
fprintf('    Default command:\n');
fprintf('      sol = run_B6_jacobi_full_candidate(600,0.006,Inf);\n');
fprintf('    To keep runtime moderate, this is not launched automatically.\n');

fprintf('\nV8 sweep finished. Figures and CSV are saved in ./figures\n');
