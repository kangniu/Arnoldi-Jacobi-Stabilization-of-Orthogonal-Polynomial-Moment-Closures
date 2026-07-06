function T = compare_riemann_B6_grid_sweep()
%COMPARE_RIEMANN_B6_GRID_SWEEP
% Grid/time sweep for B6 original Jacobian vs V+A/Jacobi wave speeds.
%
% This is the next step after the Nx=80, tEnd=2e-4 verification.
% It gradually increases resolution and final time while keeping stopping
% criteria active.  Results are printed and saved to a CSV file.

outdir = case_output_dir();

% Conservative ramp.  Increase after these pass.
cases = [
     80, 2e-4;
    160, 2e-4;
    160, 5e-4;
    240, 5e-4;
    240, 1e-3
];

nCase = size(cases,1);
rows = cell(nCase,1);

fprintf('\nB6 grid/time sweep: original Jacobian vs V+A/Jacobi\n');
fprintf('%6s %10s %10s %10s %10s %10s %12s %12s %12s %12s\n', ...
    'Nx','tEnd','stOrig','stJac','nOrig','nJac','LinfRho','L1Rho','cpuOrig','cpuJac');

for ic = 1:nCase
    cfg = default_riemann_config();
    cfg.Nx = cases(ic,1);
    cfg.tEnd = cases(ic,2);
    cfg.tau = Inf;
    cfg.CFL = 0.35;
    cfg.maxSteps = 20000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
    cfg.abortOnSmallDt = true;
    cfg.abortOnLargeWaveSpeed = true;
    cfg.printEvery = 0;
    cfg.saveFigures = false;

    x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
    U0 = initialize_riemann_state(x,cfg,6);

    tic;
    solOrig = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobian',cfg);
    cpuOrig = toc;

    tic;
    solJac = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobi',cfg);
    cpuJac = toc;

    [rhoO,uO,thetaO] = raw_to_primitive(solOrig.U);
    [rhoJ,uJ,thetaJ] = raw_to_primitive(solJac.U);

    WO = raw_to_normalized(solOrig.U,5);
    WJ = raw_to_normalized(solJac.U,5);

    LinfRho = max(abs(rhoO-rhoJ));
    L1Rho = mean(abs(rhoO-rhoJ));
    LinfU = max(abs(uO-uJ));
    LinfTheta = max(abs(thetaO-thetaJ));
    LinfW3 = max(abs(WO(4,:)-WJ(4,:)));
    LinfW4 = max(abs(WO(5,:)-WJ(5,:)));
    LinfW5 = max(abs(WO(6,:)-WJ(6,:)));

    rows{ic} = {cfg.Nx,cfg.tEnd,solOrig.status,solJac.status, ...
        solOrig.nSteps,solJac.nSteps,solOrig.t,solJac.t, ...
        LinfRho,L1Rho,LinfU,LinfTheta,LinfW3,LinfW4,LinfW5,cpuOrig,cpuJac};

    fprintf('%6d %10.2e %10s %10s %10d %10d %12.4e %12.4e %12.3f %12.3f\n', ...
        cfg.Nx,cfg.tEnd,short_status(solOrig.status),short_status(solJac.status), ...
        solOrig.nSteps,solJac.nSteps,LinfRho,L1Rho,cpuOrig,cpuJac);

    % Save representative last case figure.
    if ic == nCase
        figure('Name','B6 grid sweep final density');
        plot(x,rhoO/cfg.left.rho,'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
        plot(x,rhoJ/cfg.left.rho,'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
        xlabel('x'); ylabel('\rho / \rho_L');
        title(sprintf('B6 comparison: Nx=%d, tEnd=%.1e',cfg.Nx,cfg.tEnd));
        legend('Location','best');
        drawnow; saveas(gcf,fullfile(outdir,'B6_grid_sweep_final_density.png'));

        figure('Name','B6 grid sweep final W3');
        plot(x,WO(4,:),'LineWidth',1.4,'DisplayName','Original Jacobian'); hold on; grid on;
        plot(x,WJ(4,:),'--','LineWidth',1.4,'DisplayName','V+A/Jacobi');
        xlabel('x'); ylabel('W_3');
        title(sprintf('B6 W_3 comparison: Nx=%d, tEnd=%.1e',cfg.Nx,cfg.tEnd));
        legend('Location','best');
        drawnow; saveas(gcf,fullfile(outdir,'B6_grid_sweep_final_W3.png'));
    end
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'Nx','tEnd','statusOrig','statusJac','stepsOrig','stepsJac','timeOrig','timeJac', ...
     'LinfRho','L1Rho','LinfU','LinfTheta','LinfW3','LinfW4','LinfW5','cpuOrig','cpuJac'});

writetable(T,fullfile(outdir,'B6_grid_sweep_original_vs_jacobi.csv'));

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
