function row = compare_riemann_B8_small_original_vs_jacobi()
%COMPARE_RIEMANN_B8_SMALL_ORIGINAL_VS_JACOBI Small-grid B8 PDE check.

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = 60;
cfg.tEnd = 1e-4;
cfg.tau = Inf;
cfg.CFL = 0.25;
cfg.maxSteps = 20000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.printEvery = 0;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,8);

fprintf('\nSmall-grid B8 Riemann comparison\n');
fprintf('  Nx=%d, tEnd=%.3e, tau=Inf\n',cfg.Nx,cfg.tEnd);

tic;
solOrig = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B8',cfg.CFL,'jacobian',cfg);
cpuOrig = toc;

tic;
solJac = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B8',cfg.CFL,'jacobi',cfg);
cpuJac = toc;

[rhoO,uO,thetaO] = raw_to_primitive(solOrig.U);
[rhoJ,uJ,thetaJ] = raw_to_primitive(solJac.U);
WO = raw_to_normalized(solOrig.U,7);
WJ = raw_to_normalized(solJac.U,7);

LinfRho = max(abs(rhoO-rhoJ));
L1Rho = mean(abs(rhoO-rhoJ));
LinfU = max(abs(uO-uJ));
LinfTheta = max(abs(thetaO-thetaJ));
LinfW3 = max(abs(WO(4,:)-WJ(4,:)));
LinfW7 = max(abs(WO(8,:)-WJ(8,:)));

fprintf('\nB8 small-grid original-vs-Jacobi difference\n');
fprintf('  original status: %s, t=%.6e, steps=%d, cpu=%.2fs\n', ...
    solOrig.status,solOrig.t,solOrig.nSteps,cpuOrig);
fprintf('  Jacobi   status: %s, t=%.6e, steps=%d, cpu=%.2fs\n', ...
    solJac.status,solJac.t,solJac.nSteps,cpuJac);
fprintf('  Linf density difference = %.4e\n',LinfRho);
fprintf('  L1   density difference = %.4e\n',L1Rho);
fprintf('  Linf velocity difference = %.4e\n',LinfU);
fprintf('  Linf theta difference    = %.4e\n',LinfTheta);
fprintf('  Linf W3 difference       = %.4e\n',LinfW3);
fprintf('  Linf W7 difference       = %.4e\n',LinfW7);

figure('Name','B8 small Riemann density original vs Jacobi');
plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
plot(x,rhoJ/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
xlabel('x'); ylabel('\rho/\rho_L');
title('B8 small-grid density comparison');
legend('Location','best');
drawnow; saveas(gcf,fullfile(outdir,'B8_small_riemann_density_original_vs_jacobi.png'));

row = table(cfg.Nx,cfg.tEnd,string(solOrig.status),string(solJac.status), ...
    solOrig.nSteps,solJac.nSteps,solOrig.t,solJac.t,LinfRho,L1Rho,LinfU, ...
    LinfTheta,LinfW3,LinfW7,cpuOrig,cpuJac, ...
    'VariableNames', {'Nx','tEnd','statusOrig','statusJac','stepsOrig','stepsJac', ...
    'timeOrig','timeJac','LinfRho','L1Rho','LinfU','LinfTheta','LinfW3','LinfW7', ...
    'cpuOrig','cpuJac'});
end
