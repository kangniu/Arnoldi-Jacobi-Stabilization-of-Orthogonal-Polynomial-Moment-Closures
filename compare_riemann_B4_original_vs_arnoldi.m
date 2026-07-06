function compare_riemann_B4_original_vs_arnoldi()
%COMPARE_RIEMANN_B4_ORIGINAL_VS_ARNOLDI
% Run the same B4 Riemann problem using:
%   1. original flux-Jacobian eigenvalues,
%   2. Arnoldi/Jacobi wave speeds.
%
% For B4 these should agree closely.  This gives a clean verification before
% extending the Arnoldi/Jacobi machinery to B6/B8.

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = 600;
cfg.tEnd = 0.006;
cfg.tau = Inf;
cfg.maxSteps = 20000;
cfg.minDt = 1e-10;
cfg.maxWaveSpeed = 1e6;
cfg.printEvery = 100;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,4);

fprintf('\nB4 Riemann comparison: original Jacobian wave speeds\n');
solOrig = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B4',cfg.CFL,'jacobian',cfg);

fprintf('\nB4 Riemann comparison: Arnoldi/Jacobi wave speeds\n');
solArn = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B4',cfg.CFL,'jacobi',cfg);

[rhoO,uO,thetaO] = raw_to_primitive(solOrig.U);
[rhoA,uA,thetaA] = raw_to_primitive(solArn.U);

drho = rhoO-rhoA;
du = uO-uA;
dtheta = thetaO-thetaA;

fprintf('\nB4 Riemann original-vs-Arnoldi difference\n');
fprintf('  original status: %s, t=%.6e, steps=%d\n',solOrig.status,solOrig.t,solOrig.nSteps);
fprintf('  Arnoldi  status: %s, t=%.6e, steps=%d\n',solArn.status,solArn.t,solArn.nSteps);
fprintf('  Linf density difference = %.4e\n',max(abs(drho)));
fprintf('  L1   density difference = %.4e\n',mean(abs(drho)));
fprintf('  Linf velocity difference = %.4e\n',max(abs(du)));
fprintf('  Linf theta difference    = %.4e\n',max(abs(dtheta)));

figure('Name','B4 Riemann original vs Arnoldi density');
plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
plot(x,rhoA/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','Arnoldi/Jacobi');
xlabel('x'); ylabel('\rho / \rho_L');
title('B4 Riemann problem: original Jacobian vs Arnoldi/Jacobi');
legend('Location','best');
drawnow; saveas(gcf,fullfile(outdir,'B4_riemann_density_original_vs_arnoldi.png'));

figure('Name','B4 Riemann density difference');
plot(x,drho,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('\rho_{original}-\rho_{Arnoldi}');
title('B4 Riemann density difference');
drawnow; saveas(gcf,fullfile(outdir,'B4_riemann_density_difference.png'));

end
