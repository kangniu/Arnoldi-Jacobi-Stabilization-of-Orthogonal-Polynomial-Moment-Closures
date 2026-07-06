function run_complete_experiment_suite(mode)
%RUN_COMPLETE_EXPERIMENT_SUITE Master runner for the B4/B6/B8 V+A code package.
%
% mode:
%   'light' : fast verification suite, suitable for checking installation.
%   'paper' : recommended article-level suite excluding BGK reference.
%   'full'  : includes article suite plus BGK reference demo.
%
% Examples:
%   run_complete_experiment_suite('light')
%   run_complete_experiment_suite('paper')
%   run_complete_experiment_suite('full')

if nargin < 1 || isempty(mode)
    mode = 'light';
end

mode = lower(string(mode));

fprintf('============================================================\n');
fprintf('Complete B4/B6/B8 V+A/Jacobi experiment suite: %s\n',mode);
fprintf('============================================================\n\n');

switch mode
    case "light"
        fprintf('[1] V+A conditioning and quadrature recovery\n');
        experiment_vandermonde_arnoldi_condition;
        experiment_quadrature_recovery_vandermonde_vs_arnoldi;

        fprintf('\n[2] B4 consistency\n');
        compare_wave_speeds_B4_original_vs_arnoldi;
        compare_riemann_B4_original_vs_arnoldi;

        fprintf('\n[3] B6 static and small PDE consistency\n');
        verify_B6_formula_equilibrium_and_hyperbolicity;
        compare_B6_static_original_vs_jacobi;
        compare_riemann_B6_small_original_vs_jacobi;

        fprintf('\n[4] B8 static consistency\n');
        verify_B8_formula_equilibrium_and_hyperbolicity;
        compare_B8_charpoly_ad_vs_fd(10);
        compare_B8_static_original_vs_jacobi;

    case "paper"
        fprintf('[1] V+A conditioning and quadrature recovery\n');
        experiment_vandermonde_arnoldi_condition;
        experiment_quadrature_recovery_vandermonde_vs_arnoldi;

        fprintf('\n[2] High-order root stability diagnostic\n');
        experiment_high_order_wave_speed_stability;

        fprintf('\n[3] B4 original-vs-Jacobi verification\n');
        compare_wave_speeds_B4_original_vs_arnoldi;
        compare_riemann_B4_original_vs_arnoldi;

        fprintf('\n[4] B6 static, small-grid, and grid/time sweep\n');
        verify_B6_formula_equilibrium_and_hyperbolicity;
        compare_B6_static_original_vs_jacobi;
        compare_riemann_B6_small_original_vs_jacobi;
        compare_riemann_B6_grid_sweep;

        fprintf('\n[5] B8 static, small-grid, and grid/time sweep\n');
        verify_B8_formula_equilibrium_and_hyperbolicity;
        compare_B8_charpoly_ad_vs_fd(20);
        compare_B8_static_original_vs_jacobi;
        compare_riemann_B8_small_original_vs_jacobi;
        compare_riemann_B8_grid_sweep;

        fprintf('\n[6] Full B4/B6 V+A/Jacobi regime sweeps\n');
        run_B4_full_regime_sweep(600,0.006);
        run_B6_full_regime_sweep(600,0.006);
        plot_B6_regime_summary_from_files(600,0.006);

        fprintf('\n[7] Original Jacobian versus V+A/Jacobi profile figures\n');
        run_original_vs_va_profile_figures;

        fprintf('\n[8] Paper-style B4/B6/B8/BGK profile figures\n');
        run_profile_comparison_figures(180,501,0.006);

        fprintf('\n[9] Generate LaTeX tables from available CSV files\n');
        make_paper_tables_from_results;

    case "full"
        run_complete_experiment_suite('paper');

        fprintf('\n[10] B6 V+A/Jacobi vs BGK reference demo\n');
        run_B6_vs_BGK_regime_demo(300,501,0.006);

        fprintf('\n[11] Generate LaTeX tables including BGK CSV files\n');
        make_paper_tables_from_results;

    otherwise
        error('Unknown mode: %s. Use light, paper, or full.',mode);
end

fprintf('\nComplete suite finished. Outputs are in ./figures and ./tables.\n');
end
