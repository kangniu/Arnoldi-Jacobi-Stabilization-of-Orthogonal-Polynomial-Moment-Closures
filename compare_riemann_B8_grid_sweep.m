function T = compare_riemann_B8_grid_sweep()
%COMPARE_RIEMANN_B8_GRID_SWEEP Conservative B8 original-vs-Jacobi ramp.

outdir = case_output_dir();

cases = [
     40, 5e-5;
     60, 1e-4;
     80, 1e-4;
    100, 2e-4
];

rows = cell(size(cases,1),1);

fprintf('\nB8 grid/time sweep: original Jacobian vs V+A/Jacobi\n');
fprintf('%6s %10s %10s %10s %10s %10s %12s %12s %12s %12s\n', ...
    'Nx','tEnd','stOrig','stJac','nOrig','nJac','LinfRho','L1Rho','cpuOrig','cpuJac');

for ic = 1:size(cases,1)
    cfg = default_riemann_config();
    cfg.Nx = cases(ic,1);
    cfg.tEnd = cases(ic,2);
    cfg.tau = Inf;
    cfg.CFL = 0.25;
    cfg.maxSteps = 20000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
    cfg.printEvery = 0;
    cfg.saveFigures = false;

    x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
    U0 = initialize_riemann_state(x,cfg,8);

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

    rows{ic} = {cfg.Nx,cfg.tEnd,solOrig.status,solJac.status, ...
        solOrig.nSteps,solJac.nSteps,solOrig.t,solJac.t,LinfRho,L1Rho, ...
        LinfU,LinfTheta,LinfW3,LinfW7,cpuOrig,cpuJac};

    fprintf('%6d %10.2e %10s %10s %10d %10d %12.4e %12.4e %12.3f %12.3f\n', ...
        cfg.Nx,cfg.tEnd,short_status(solOrig.status),short_status(solJac.status), ...
        solOrig.nSteps,solJac.nSteps,LinfRho,L1Rho,cpuOrig,cpuJac);

    if ic == size(cases,1)
        figure('Name','B8 grid sweep final density');
        plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
        plot(x,rhoJ/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
        xlabel('x'); ylabel('\rho/\rho_L');
        title(sprintf('B8 comparison: Nx=%d, tEnd=%.1e',cfg.Nx,cfg.tEnd));
        legend('Location','best');
        drawnow; saveas(gcf,fullfile(outdir,'B8_grid_sweep_final_density.png'));
    end
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'Nx','tEnd','statusOrig','statusJac','stepsOrig','stepsJac','timeOrig','timeJac', ...
     'LinfRho','L1Rho','LinfU','LinfTheta','LinfW3','LinfW7','cpuOrig','cpuJac'});

writetable(T,fullfile(outdir,'B8_grid_sweep_original_vs_jacobi.csv'));
end

function s = short_status(status)
switch status
    case 'finished'
        s = 'finished';
    case 'aborted_large_wave_speed'
        s = 'largeS';
    case 'aborted_small_dt'
        s = 'smallDt';
    case 'aborted_max_steps'
        s = 'maxStep';
    otherwise
        s = status;
end
end
