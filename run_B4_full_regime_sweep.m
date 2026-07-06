function T = run_B4_full_regime_sweep(Nx,tEnd)
%RUN_B4_FULL_REGIME_SWEEP Full B4 V+A/Jacobi regime sweep.

if nargin < 1 || isempty(Nx), Nx = 600; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end

outdir = case_output_dir();

taus = [Inf, 1e-3, 1e-6];
tags = {'free','transition','continuum'};
rows = cell(numel(taus),1);

fprintf('============================================================\n');
fprintf('B4 full-regime V+A/Jacobi sweep\n');
fprintf('Nx=%d, tEnd=%.6e\n',Nx,tEnd);
fprintf('============================================================\n');

for k = 1:numel(taus)
    cfg = default_riemann_config();
    cfg.Nx = Nx;
    cfg.tEnd = tEnd;
    cfg.tau = taus(k);
    cfg.CFL = 0.35;
    cfg.maxSteps = 80000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
    cfg.printEvery = 100;
    cfg.saveFigures = false;

    x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
    U0 = initialize_riemann_state(x,cfg,4);

    fprintf('\nB4 V+A/Jacobi benchmark case: %s\n',tags{k});
    tic;
    sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B4',cfg.CFL,'jacobi',cfg);
    cpu = toc;

    [rho,u,theta] = raw_to_primitive(sol.U);
    W = raw_to_normalized(sol.U,3);

    figFile = fullfile(outdir,sprintf('B4_jacobi_%s_density_Nx%d_t%.0e.png',tags{k},Nx,tEnd));
    matFile = fullfile(outdir,sprintf('B4_jacobi_%s_solution_Nx%d_t%.0e.mat',tags{k},Nx,tEnd));

    figure('Name',['B4 density ',tags{k}]);
    plot(x,rho/cfg.left.rho,'LineWidth',1.5); grid on;
    xlabel('x'); ylabel('\rho/\rho_L');
    title(sprintf('B4 V+A/Jacobi density: %s, Nx=%d, t=%.3e',tags{k},Nx,sol.t));
    save_clean_figure(gcf,figFile);

    save(matFile,'x','sol','rho','u','theta','W','cfg','cpu');

    rows{k} = table(string(tags{k}),Nx,tEnd,taus(k),string(sol.status), ...
        sol.nSteps,sol.t,cpu,min(rho),min(theta),max(abs(W(4,:))), ...
        string(figFile),string(matFile), ...
        'VariableNames', {'caseTag','Nx','tEnd','tau','status','steps', ...
        'finalTime','cpu','minRho','minTheta','maxAbsW3','densityFigure','solutionFile'});
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,sprintf('B4_full_regime_sweep_Nx%d_t%.0e.csv',Nx,tEnd));
writetable(T,csvFile);

fprintf('\nB4 full-regime sweep summary\n');
disp(T(:,{'caseTag','Nx','tEnd','tau','status','steps','finalTime','cpu','minRho','minTheta','maxAbsW3'}));
fprintf('Summary CSV saved to %s\n',csvFile);
end
