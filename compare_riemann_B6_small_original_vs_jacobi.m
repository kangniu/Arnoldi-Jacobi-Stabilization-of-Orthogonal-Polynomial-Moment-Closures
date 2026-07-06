function compare_riemann_B6_small_original_vs_jacobi()
%COMPARE_RIEMANN_B6_SMALL_ORIGINAL_VS_JACOBI
% Small-grid B6 Riemann comparison:
%   Nx   = 80
%   tEnd = 2e-4
%
% The two runs use the same finite-volume HLL solver and B6 closing flux.
% Only the wave-speed evaluation is changed:
%   original: finite-difference flux Jacobian eigenvalues,
%   V+A/Jacobi: semi-explicit B6 Jacobi matrix eigenvalues.

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = 80;
cfg.tEnd = 2e-4;
cfg.tau = Inf;
cfg.CFL = 0.35;
cfg.maxSteps = 5000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
cfg.printEvery = 20;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,6);

fprintf('\nSmall-grid B6 Riemann comparison\n');
fprintf('  Nx=%d, tEnd=%.3e, tau=%g\n',cfg.Nx,cfg.tEnd,cfg.tau);

fprintf('\n  Original finite-difference flux-Jacobian wave speeds\n');
tic;
solOrig = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobian',cfg);
timeOrig = toc;

fprintf('\n  V+A / semi-explicit Jacobi wave speeds\n');
tic;
solJac = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobi',cfg);
timeJac = toc;

[rhoO,uO,thetaO] = raw_to_primitive(solOrig.U);
[rhoJ,uJ,thetaJ] = raw_to_primitive(solJac.U);

WO = raw_to_normalized(solOrig.U,5);
WJ = raw_to_normalized(solJac.U,5);

drho = rhoO-rhoJ;
du = uO-uJ;
dtheta = thetaO-thetaJ;
dW3 = WO(4,:)-WJ(4,:);
dW4 = WO(5,:)-WJ(5,:);
dW5 = WO(6,:)-WJ(6,:);

fprintf('\nB6 small-grid original-vs-Jacobi difference\n');
fprintf('  original status: %s, t=%.6e, steps=%d, cpu=%.2fs\n', ...
    solOrig.status,solOrig.t,solOrig.nSteps,timeOrig);
fprintf('  Jacobi   status: %s, t=%.6e, steps=%d, cpu=%.2fs\n', ...
    solJac.status,solJac.t,solJac.nSteps,timeJac);
fprintf('  Linf density difference = %.4e\n',max(abs(drho)));
fprintf('  L1   density difference = %.4e\n',mean(abs(drho)));
fprintf('  Linf velocity difference = %.4e\n',max(abs(du)));
fprintf('  Linf theta difference    = %.4e\n',max(abs(dtheta)));
fprintf('  Linf W3 difference       = %.4e\n',max(abs(dW3)));
fprintf('  Linf W4 difference       = %.4e\n',max(abs(dW4)));
fprintf('  Linf W5 difference       = %.4e\n',max(abs(dW5)));

figure('Name','B6 small Riemann density');
plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
plot(x,rhoJ/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
xlabel('x'); ylabel('\rho / \rho_L');
title('B6 small-grid Riemann: original Jacobian vs V+A/Jacobi');
legend('Location','best');
drawnow; saveas(gcf,fullfile(outdir,'B6_small_riemann_density_original_vs_jacobi.png'));

figure('Name','B6 small Riemann density difference');
plot(x,drho,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('\rho_{orig}-\rho_{Jac}');
title('B6 small-grid density difference');
drawnow; saveas(gcf,fullfile(outdir,'B6_small_riemann_density_difference.png'));

figure('Name','B6 small Riemann W3');
plot(x,WO(4,:),'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
plot(x,WJ(4,:),'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
xlabel('x'); ylabel('W_3');
title('B6 small-grid normalized heat-flux moment W_3');
legend('Location','best');
drawnow; saveas(gcf,fullfile(outdir,'B6_small_riemann_W3_original_vs_jacobi.png'));

end
