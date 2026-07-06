function experiment_B6_jacobian_failure_diagnostic()
%EXPERIMENT_B6_JACOBIAN_FAILURE_DIAGNOSTIC
% Diagnostic run showing why high-order closures should not rely on a
% finite-difference flux-Jacobian wave-speed calculation.
%
% The B6 closing flux is included in the code, but B6 wave speeds are still
% computed by a numerical flux Jacobian.  This diagnostic records the stopped
% state and plots the wave-speed distribution and realizability margin.

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = 600;
cfg.tEnd = 0.006;
cfg.tau = Inf;
cfg.caseName = 'B6_diagnostic';
cfg.closures = {'B6'};
cfg.waveMethodOther = 'jacobian';
cfg.maxSteps = 5000;
cfg.minDt = 1e-9;
cfg.maxWaveSpeed = 1e6;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
cfg.printEvery = 100;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,6);

fprintf('\nB6 finite-difference Jacobian diagnostic run\n');
sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobian',cfg);

U = sol.U;
W = raw_to_normalized(U,5);

W3 = W(4,:);
W4 = W(5,:);
W5 = W(6,:);
margin = W4 - W3.^2 - 1;      % Hamburger realizability margin M_4-like.
denB6 = W3.^2 - W4 + 1;       % denominator used in B6 W6 formula.

smaxCell = zeros(size(x));
for i = 1:numel(x)
    lam = wave_speeds(U(:,i),'B6','jacobian');
    smaxCell(i) = max(abs(lam));
end

fprintf('  status       = %s\n',sol.status);
fprintf('  final time   = %.6e\n',sol.t);
fprintf('  steps        = %d\n',sol.nSteps);
fprintf('  max speed    = %.4e\n',max(smaxCell));
fprintf('  min margin   = %.4e\n',min(margin));
fprintf('  min |denB6|  = %.4e\n',min(abs(denB6)));
fprintf('  max |W3|     = %.4e\n',max(abs(W3)));
fprintf('  max |W4|     = %.4e\n',max(abs(W4)));
fprintf('  max |W5|     = %.4e\n',max(abs(W5)));

figure('Name','B6 diagnostic wave speed');
semilogy(x,smaxCell+eps,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('local max wave speed');
title(['B6 finite-difference Jacobian diagnostic, status=',sol.status]);
drawnow; saveas(gcf,fullfile(outdir,'B6_diagnostic_wave_speed.png'));

figure('Name','B6 diagnostic margin');
plot(x,margin,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('W_4-W_3^2-1');
title('B6 diagnostic: realizability margin');
drawnow; saveas(gcf,fullfile(outdir,'B6_diagnostic_realizability_margin.png'));

figure('Name','B6 diagnostic denominator');
semilogy(x,abs(denB6)+eps,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('|W_3^2-W_4+1|');
title('B6 diagnostic: denominator in B6 closure');
drawnow; saveas(gcf,fullfile(outdir,'B6_diagnostic_denominator.png'));

end
