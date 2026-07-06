function row = run_B8_jacobi_benchmark_case(Nx,tEnd,tau,caseTag)
%RUN_B8_JACOBI_BENCHMARK_CASE Run one B8 V+A/Jacobi Riemann benchmark.

if nargin < 4 || isempty(caseTag)
    if isinf(tau)
        caseTag = 'free';
    else
        caseTag = sprintf('tau_%g',tau);
    end
end

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = Nx;
cfg.tEnd = tEnd;
cfg.tau = tau;
cfg.CFL = 0.30;
cfg.maxSteps = 100000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
cfg.printEvery = 100;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,8);

fprintf('\nB8 V+A/Jacobi benchmark case: %s\n',caseTag);
fprintf('  Nx=%d, tEnd=%.6e, tau=%g\n',Nx,tEnd,tau);

tic;
sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B8',cfg.CFL,'jacobi',cfg);
cpu = toc;

[rho,u,theta] = raw_to_primitive(sol.U);
W = raw_to_normalized(sol.U,7);

rhoFile = fullfile(outdir,sprintf('B8_jacobi_%s_density_Nx%d_t%.0e.png',caseTag,Nx,tEnd));
w3File  = fullfile(outdir,sprintf('B8_jacobi_%s_W3_Nx%d_t%.0e.png',caseTag,Nx,tEnd));
momFile = fullfile(outdir,sprintf('B8_jacobi_%s_moments_Nx%d_t%.0e.png',caseTag,Nx,tEnd));
matFile = fullfile(outdir,sprintf('B8_jacobi_%s_solution_Nx%d_t%.0e.mat',caseTag,Nx,tEnd));

figure('Name',['B8 density ',caseTag]);
plot(x,rho/cfg.left.rho,'LineWidth',1.5); grid on;
xlabel('x'); ylabel('\rho/\rho_L');
title(sprintf('B8 V+A/Jacobi density: %s, Nx=%d, t=%.3e',caseTag,Nx,sol.t));
save_clean_figure(gcf,rhoFile);

figure('Name',['B8 W3 ',caseTag]);
plot(x,W(4,:),'LineWidth',1.5); grid on;
xlabel('x'); ylabel('W_3');
title(sprintf('B8 V+A/Jacobi W_3: %s, Nx=%d, t=%.3e',caseTag,Nx,sol.t));
save_clean_figure(gcf,w3File);

figure('Name',['B8 normalized moments ',caseTag]);
plot(x,W(4,:),'LineWidth',1.3,'DisplayName','W_3'); hold on; grid on;
plot(x,W(5,:),'--','LineWidth',1.3,'DisplayName','W_4');
plot(x,W(6,:),':','LineWidth',1.6,'DisplayName','W_5');
plot(x,W(8,:),'-.','LineWidth',1.2,'DisplayName','W_7');
xlabel('x'); ylabel('normalized moments');
title(sprintf('B8 V+A/Jacobi normalized moments: %s',caseTag));
legend('Location','best');
save_clean_figure(gcf,momFile);

save(matFile,'x','sol','rho','u','theta','W','cfg','cpu');

fprintf('  status=%s, steps=%d, final time %.6e, cpu %.2fs\n', ...
    sol.status,sol.nSteps,sol.t,cpu);
fprintf('  min rho %.4e, min theta %.4e\n',min(rho),min(theta));
fprintf('  max |W3| %.4e, max |W7| %.4e\n',max(abs(W(4,:))),max(abs(W(8,:))));

row = table( ...
    string(caseTag), Nx, tEnd, tau, string(sol.status), sol.nSteps, sol.t, cpu, ...
    min(rho), min(theta), max(abs(W(4,:))), max(abs(W(8,:))), ...
    string(rhoFile), string(w3File), string(momFile), string(matFile), ...
    'VariableNames', {'caseTag','Nx','tEnd','tau','status','steps','finalTime','cpu', ...
    'minRho','minTheta','maxAbsW3','maxAbsW7', ...
    'densityFigure','W3Figure','momentsFigure','solutionFile'} );
end
