function T = compare_B6_original_vs_jacobi_moderate()
%COMPARE_B6_ORIGINAL_VS_JACOBI_MODERATE
% Moderate original-vs-Jacobi B6 comparison after v8.
%
% This is optional because the original finite-difference Jacobian is slower.
% It checks a longer time than v8 but still keeps runtime reasonable.

cases = [
    300, 2e-3;
    400, 2e-3;
    400, 3e-3
];

outdir = case_output_dir();
rows = cell(size(cases,1),1);

fprintf('\nModerate B6 original-vs-Jacobi comparison\n');
fprintf('%6s %10s %10s %10s %10s %10s %12s %12s %10s %10s\n', ...
    'Nx','tEnd','stOrig','stJac','nOrig','nJac','LinfRho','L1Rho','cpuOrig','cpuJac');

for ic = 1:size(cases,1)
    cfg = default_riemann_config();
    cfg.Nx = cases(ic,1);
    cfg.tEnd = cases(ic,2);
    cfg.tau = Inf;
    cfg.CFL = 0.35;
    cfg.maxSteps = 50000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
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

    rows{ic} = table(cfg.Nx,cfg.tEnd,string(solOrig.status),string(solJac.status), ...
        solOrig.nSteps,solJac.nSteps,solOrig.t,solJac.t, ...
        LinfRho,L1Rho,LinfU,LinfTheta,LinfW3,LinfW4,LinfW5,cpuOrig,cpuJac, ...
        'VariableNames', {'Nx','tEnd','statusOrig','statusJac','stepsOrig','stepsJac', ...
        'timeOrig','timeJac','LinfRho','L1Rho','LinfU','LinfTheta', ...
        'LinfW3','LinfW4','LinfW5','cpuOrig','cpuJac'});

    fprintf('%6d %10.2e %10s %10s %10d %10d %12.4e %12.4e %10.2f %10.2f\n', ...
        cfg.Nx,cfg.tEnd,short_status(solOrig.status),short_status(solJac.status), ...
        solOrig.nSteps,solJac.nSteps,LinfRho,L1Rho,cpuOrig,cpuJac);
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,'B6_moderate_original_vs_jacobi.csv');
writetable(T,csvFile);
fprintf('CSV saved to %s\n',csvFile);

end

function s = short_status(status)
switch char(status)
    case 'finished'
        s = 'finished';
    case 'aborted_large_wave_speed'
        s = 'largeS';
    case 'aborted_small_dt'
        s = 'smallDt';
    case 'aborted_max_steps'
        s = 'maxStep';
    otherwise
        s = char(status);
end
end
