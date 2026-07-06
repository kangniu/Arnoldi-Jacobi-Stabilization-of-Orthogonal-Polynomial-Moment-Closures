function sol = run_B6_small_jacobi_only()
%RUN_B6_SMALL_JACOBI_ONLY Quick B6 small-grid run using V+A/Jacobi wave speeds.
cfg = default_riemann_config();
cfg.Nx = 80;
cfg.tEnd = 2e-4;
cfg.tau = Inf;
cfg.CFL = 0.35;
cfg.maxSteps = 5000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.printEvery = 20;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,6);
sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobi',cfg);
end
