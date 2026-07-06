function plot_B6_regime_summary_from_files(Nx,tEnd)
%PLOT_B6_REGIME_SUMMARY_FROM_FILES Create combined B6 regime plots.
%
% Requires that run_B6_full_regime_sweep(Nx,tEnd) has already been run.

if nargin < 1 || isempty(Nx), Nx = 600; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end

outdir = case_output_dir();
tags = {'free','transition','continuum'};
labels = {'free molecular','transition','continuum'};

data = cell(numel(tags),1);
for k = 1:numel(tags)
    matFile = fullfile(outdir,sprintf('B6_jacobi_%s_solution_Nx%d_t%.0e.mat',tags{k},Nx,tEnd));
    if ~exist(matFile,'file')
        error('Missing file: %s. Run run_B6_full_regime_sweep first.',matFile);
    end
    data{k} = load(matFile);
end

figure('Name','B6 regime density comparison');
for k = 1:numel(tags)
    plot(data{k}.x,data{k}.rho/data{k}.cfg.left.rho,'LineWidth',1.4, ...
        'DisplayName',labels{k});
    hold on;
end
grid on; xlabel('x'); ylabel('\rho/\rho_L');
title(sprintf('B6 V+A/Jacobi density comparison, Nx=%d, t=%.3e',Nx,tEnd));
legend('Location','best');
save_clean_figure(gcf,fullfile(outdir,sprintf('B6_regime_density_comparison_Nx%d_t%.0e.png',Nx,tEnd)));

figure('Name','B6 regime W3 comparison');
for k = 1:numel(tags)
    plot(data{k}.x,data{k}.W(4,:),'LineWidth',1.4,'DisplayName',labels{k});
    hold on;
end
grid on; xlabel('x'); ylabel('W_3');
title(sprintf('B6 V+A/Jacobi W_3 comparison, Nx=%d, t=%.3e',Nx,tEnd));
legend('Location','best');
save_clean_figure(gcf,fullfile(outdir,sprintf('B6_regime_W3_comparison_Nx%d_t%.0e.png',Nx,tEnd)));

figure('Name','B6 regime theta comparison');
for k = 1:numel(tags)
    plot(data{k}.x,data{k}.theta,'LineWidth',1.4,'DisplayName',labels{k});
    hold on;
end
grid on; xlabel('x'); ylabel('\theta');
title(sprintf('B6 V+A/Jacobi temperature comparison, Nx=%d, t=%.3e',Nx,tEnd));
legend('Location','best');
save_clean_figure(gcf,fullfile(outdir,sprintf('B6_regime_theta_comparison_Nx%d_t%.0e.png',Nx,tEnd)));

end
