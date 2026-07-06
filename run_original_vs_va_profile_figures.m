function T = run_original_vs_va_profile_figures()
%RUN_ORIGINAL_VS_VA_PROFILE_FIGURES Compare original Jacobian and V+A profiles.
%
% The "original" implementation uses finite-difference flux-Jacobian wave
% speeds.  The V+A implementation uses Arnoldi/Jacobi wave speeds while keeping
% the same closure flux and finite-volume solver.

outdir = case_output_dir();

cases = {
    'B4', 600, 0.006, Inf, 4, 'free_full';
    'B4', 600, 0.006, 1e-3, 4, 'transition_full';
    'B4', 600, 0.006, 1e-6, 4, 'continuum_full';
    'B6', 240, 1e-3, Inf, 6, 'free_short';
    'B8', 100, 2e-4, Inf, 8, 'free_short'
};

rows = cell(size(cases,1),1);

for ic = 1:size(cases,1)
    closure = cases{ic,1};
    Nx = cases{ic,2};
    tEnd = cases{ic,3};
    tau = cases{ic,4};
    nMom = cases{ic,5};
    tag = cases{ic,6};
    tagLabel = strrep(tag,'_',' ');

    cfg = default_riemann_config();
    cfg.Nx = Nx;
    cfg.tEnd = tEnd;
    cfg.tau = tau;
    cfg.CFL = 0.35;
    if strcmpi(closure,'B8')
        cfg.CFL = 0.25;
    end
    cfg.maxSteps = 100000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
    cfg.printEvery = 0;
    cfg.saveFigures = false;

    x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
    U0 = initialize_riemann_state(x,cfg,nMom);

    fprintf('\nOriginal-vs-V+A profile: %s, %s, Nx=%d, tEnd=%.3e\n', ...
        closure,tag,Nx,tEnd);

    tic;
    solOrig = solve_moment_1d(U0,x,tEnd,tau,closure,cfg.CFL,'jacobian',cfg);
    cpuOrig = toc;

    tic;
    solVA = solve_moment_1d(U0,x,tEnd,tau,closure,cfg.CFL,'jacobi',cfg);
    cpuVA = toc;

    [rhoO,uO,thetaO] = raw_to_primitive(solOrig.U);
    [rhoV,uV,thetaV] = raw_to_primitive(solVA.U);
    maxMoment = min(size(solOrig.U,1)-1,7);
    WO = raw_to_normalized(solOrig.U,maxMoment);
    WV = raw_to_normalized(solVA.U,maxMoment);

    LinfRho = max(abs(rhoO-rhoV));
    L1Rho = mean(abs(rhoO-rhoV));
    LinfU = max(abs(uO-uV));
    LinfTheta = max(abs(thetaO-thetaV));
    LinfW3 = max(abs(WO(4,:)-WV(4,:)));

    figFile = fullfile(outdir,sprintf('%s_original_vs_va_%s_Nx%d_t%.0e.png',closure,tag,Nx,tEnd));
    figure('Name',[closure,' original vs V+A ',tagLabel],'Position',[100,100,900,720]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','original: finite-difference Jacobian'); hold on; grid on;
    plot(x,rhoV/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','new: V+A/Jacobi');
    ylabel('\rho/\rho_L');
    title(sprintf('%s original vs V+A/Jacobi, %s',closure,tagLabel));
    legend('Location','best');

    nexttile;
    plot(x,uO,'LineWidth',1.4,'DisplayName','original'); hold on; grid on;
    plot(x,uV,'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
    ylabel('u');

    nexttile;
    plot(x,WO(4,:),'LineWidth',1.4,'DisplayName','original'); hold on; grid on;
    plot(x,WV(4,:),'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
    xlabel('x'); ylabel('W_3');

    save_clean_figure(gcf,figFile);

    matFile = fullfile(outdir,sprintf('%s_original_vs_va_%s_Nx%d_t%.0e.mat',closure,tag,Nx,tEnd));
    save(matFile,'closure','tag','cfg','x','solOrig','solVA','rhoO','rhoV','uO','uV', ...
        'thetaO','thetaV','WO','WV','cpuOrig','cpuVA');

    rows{ic} = table(string(closure),string(tag),Nx,tEnd,tau, ...
        string(solOrig.status),string(solVA.status),solOrig.nSteps,solVA.nSteps, ...
        solOrig.t,solVA.t,cpuOrig,cpuVA,LinfRho,L1Rho,LinfU,LinfTheta,LinfW3, ...
        string(figFile),string(matFile), ...
        'VariableNames', {'closure','caseTag','Nx','tEnd','tau','statusOrig','statusVA', ...
        'stepsOrig','stepsVA','timeOrig','timeVA','cpuOrig','cpuVA','LinfRho','L1Rho', ...
        'LinfU','LinfTheta','LinfW3','figureFile','matFile'});

    fprintf('  status original=%s, V+A=%s, Linf rho=%.3e, L1 rho=%.3e\n', ...
        solOrig.status,solVA.status,LinfRho,L1Rho);
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,'original_vs_va_profile_summary.csv');
writetable(T,csvFile);
disp(T(:,{'closure','caseTag','Nx','tEnd','statusOrig','statusVA','LinfRho','L1Rho','cpuOrig','cpuVA'}));
fprintf('Original-vs-V+A profile CSV saved to %s\n',csvFile);
end
