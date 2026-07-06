function cfg = default_riemann_config()
%DEFAULT_RIEMANN_CONFIG Configuration for the JCP one-dimensional Riemann case.

cfg.xmin = -7.0;
cfg.xmax =  7.0;
cfg.Nx = 600;
cfg.tEnd = 0.006;
cfg.CFL = 0.35;
cfg.tau = 1e-3;
cfg.caseName = 'transition';
cfg.closures = {'B4','B6'};
cfg.waveMethodB4 = 'jacobi';
cfg.waveMethodOther = 'jacobian';
cfg.plotEvery = Inf;

% Robust stopping criteria.
cfg.maxSteps = 20000;
cfg.minDt = 1e-10;
cfg.maxWaveSpeed = 1e6;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
cfg.printEvery = 50;
cfg.saveFigures = true;

% Initial conditions from Morin--McDonald JCP Riemann problem.
cfg.left.rho = 4.696;       % kg/m^3
cfg.left.u   = 0.0;         % m/s
cfg.left.p   = 404400.0;    % Pa

cfg.right.rho = 1.408;      % kg/m^3
cfg.right.u   = 0.0;        % m/s
cfg.right.p   = 101100.0;   % Pa

end
