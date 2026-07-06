clear; clc; close all;

fprintf('============================================================\n');
fprintf('V9: B6 full benchmark driver\n');
fprintf('============================================================\n\n');

fprintf('[1] Recommended immediate run: full B6 V+A/Jacobi regime sweep\n');
fprintf('    This will run tau=Inf, 1e-3, 1e-6 with Nx=600, tEnd=0.006.\n');
fprintf('    Command:\n');
fprintf('      T = run_B6_full_regime_sweep(600,0.006);\n\n');

fprintf('[2] After the sweep finishes, create combined regime plots:\n');
fprintf('      plot_B6_regime_summary_from_files(600,0.006);\n\n');

fprintf('[3] Optional longer original-vs-Jacobi check:\n');
fprintf('      Tcmp = compare_B6_original_vs_jacobi_moderate();\n\n');

fprintf('This driver intentionally does not launch the full benchmark automatically,\n');
fprintf('so that you can start it explicitly and monitor the output.\n');
