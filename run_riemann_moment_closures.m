function results = run_riemann_moment_closures(cfg)
%RUN_RIEMANN_MOMENT_CLOSURES Solve the 1-D Riemann problem for several closures.
%
% The moment model is
%     d_t U + d_x F(U) = C(U),
% where F = [U_1, ..., U_n, U_{n+1}(U)]^T.
%
% Supported closures:
%   B4: four moments U_0,...,U_3, closure U_4.
%   B6: six moments U_0,...,U_5, closure U_6.

if nargin < 1
    cfg = default_riemann_config();
end

outdir = case_output_dir();

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
dx = x(2)-x(1);

results = struct();

figure('Name',['Density ',cfg.caseName]); hold on; grid on;
xlabel('x'); ylabel('\rho / \rho_L');
title(['Normalized density: ',strrep(cfg.caseName,'_',' ')]);

figure('Name',['Velocity ',cfg.caseName]); hold on; grid on;
xlabel('x'); ylabel('u');
title(['Velocity: ',strrep(cfg.caseName,'_',' ')]);

figure('Name',['HeatFlux ',cfg.caseName]); hold on; grid on;
xlabel('x'); ylabel('q / (\rho_L \theta_L^{3/2})');
title(['Normalized heat flux: ',strrep(cfg.caseName,'_',' ')]);

for ic = 1:numel(cfg.closures)
    closure = cfg.closures{ic};
    nMom = closure_num_moments(closure);

    U0 = initialize_riemann_state(x,cfg,nMom);

    fprintf('  Closure %-4s: nMom=%d, Nx=%d, tau=%g\n', ...
        closure,nMom,cfg.Nx,cfg.tau);

    if strcmpi(closure,'B4')
        waveMethod = cfg.waveMethodB4;
    else
        waveMethod = cfg.waveMethodOther;
    end

    sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,closure,cfg.CFL,waveMethod,cfg);
    results.(closure) = sol;

    [rho,u,theta] = raw_to_primitive(sol.U);
    W = raw_to_normalized(sol.U,min(5,nMom-1));
    if size(W,1) >= 4
        qnorm = W(4,:);  % normalized third central moment W_3
    else
        qnorm = zeros(size(rho));
    end

    figure(findobj('Name',['Density ',cfg.caseName]));
    plot(x,rho/cfg.left.rho,'LineWidth',1.4,'DisplayName',closure);

    figure(findobj('Name',['Velocity ',cfg.caseName]));
    plot(x,u,'LineWidth',1.4,'DisplayName',closure);

    figure(findobj('Name',['HeatFlux ',cfg.caseName]));
    plot(x,qnorm,'LineWidth',1.4,'DisplayName',closure);

    fprintf('    status=%s: steps=%d, final time %.6e, min rho %.3e, min theta %.3e\n', ...
        sol.status, sol.nSteps, sol.t, min(rho), min(theta));
end

if ~isfield(cfg,'saveFigures') || cfg.saveFigures
    figure(findobj('Name',['Density ',cfg.caseName]));
    legend('Location','best'); drawnow; saveas(gcf,fullfile(outdir,['density_',cfg.caseName,'.png']));

    figure(findobj('Name',['Velocity ',cfg.caseName]));
    legend('Location','best'); drawnow; saveas(gcf,fullfile(outdir,['velocity_',cfg.caseName,'.png']));

    figure(findobj('Name',['HeatFlux ',cfg.caseName]));
    legend('Location','best'); drawnow; saveas(gcf,fullfile(outdir,['heatflux_',cfg.caseName,'.png']));
end

end
