%RUN_THIS_FIRST Entry point for the complete code package.
%
% Fast smoke test:
%   run_complete_experiment_suite('light')
%
% Article-level experiments:
%   run_complete_experiment_suite('paper')
%
% Full suite including BGK reference demo:
%   run_complete_experiment_suite('full')
%
% Direct final B6 benchmark:
%   T = run_B6_full_regime_sweep(600,0.006);
%   plot_B6_regime_summary_from_files(600,0.006);
%
% Optional B6-vs-BGK reference:
%   Tbgk = run_B6_vs_BGK_regime_demo(300,501,0.006);

fprintf('Recommended commands:\n');
fprintf('  run_complete_experiment_suite(''light'')\n');
fprintf('  run_complete_experiment_suite(''paper'')\n');
fprintf('  run_complete_experiment_suite(''full'')\n');
