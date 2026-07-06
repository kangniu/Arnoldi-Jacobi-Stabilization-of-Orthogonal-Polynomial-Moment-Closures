function sol = run_B6_jacobi_full_candidate(Nx,tEnd,tau)
%RUN_B6_JACOBI_FULL_CANDIDATE
% Candidate full B6 Riemann run using only V+A/Jacobi wave speeds.
%
% Default:
%   Nx   = 600
%   tEnd = 0.006
%   tau  = Inf
%
% This avoids the finite-difference Jacobian bottleneck and tests whether
% the semi-explicit Jacobi realization can carry a longer B6 simulation.

if nargin < 1 || isempty(Nx), Nx = 600; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end
if nargin < 3 || isempty(tau), tau = Inf; end

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = Nx;
cfg.tEnd = tEnd;
cfg.tau = tau;
cfg.CFL = 0.35;
cfg.maxSteps = 50000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.abortOnSmallDt = true;
cfg.abortOnLargeWaveSpeed = true;
cfg.printEvery = 100;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,6);

fprintf('\nFull-candidate B6 V+A/Jacobi run\n');
fprintf('  Nx=%d, tEnd=%.3e, tau=%g\n',cfg.Nx,cfg.tEnd,cfg.tau);

tic;
sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobi',cfg);
cpu = toc;

[rho,u,theta] = raw_to_primitive(sol.U);
W = raw_to_normalized(sol.U,5);

fprintf('\nB6 V+A/Jacobi full-candidate summary\n');
fprintf('  status: %s\n',sol.status);
fprintf('  final time: %.6e\n',sol.t);
fprintf('  steps: %d\n',sol.nSteps);
fprintf('  cpu: %.2fs\n',cpu);
fprintf('  min rho: %.4e, min theta: %.4e\n',min(rho),min(theta));
fprintf('  max |W3|: %.4e, max |W4|: %.4e, max |W5|: %.4e\n', ...
    max(abs(W(4,:))),max(abs(W(5,:))),max(abs(W(6,:))));

figure('Name','B6 full candidate density');
plot(x,rho/cfg.left.rho,'LineWidth',1.4); grid on;
xlabel('x'); ylabel('\rho / \rho_L');
title(sprintf('B6 V+A/Jacobi candidate: Nx=%d, t=%.3e, tau=%g',Nx,sol.t,tau));
drawnow; save_clean_figure(gcf,fullfile(outdir,sprintf('B6_jacobi_candidate_density_Nx%d_t%.0e.png',Nx,tEnd)));

figure('Name','B6 full candidate W3');
plot(x,W(4,:),'LineWidth',1.4); grid on;
xlabel('x'); ylabel('W_3');
title(sprintf('B6 V+A/Jacobi W_3: Nx=%d, t=%.3e, tau=%g',Nx,sol.t,tau));
drawnow; save_clean_figure(gcf,fullfile(outdir,sprintf('B6_jacobi_candidate_W3_Nx%d_t%.0e.png',Nx,tEnd)));

end
