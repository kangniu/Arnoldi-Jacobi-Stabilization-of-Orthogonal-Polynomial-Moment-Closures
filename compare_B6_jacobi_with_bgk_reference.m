function row = compare_B6_jacobi_with_bgk_reference(Nx,tEnd,tau,Nv,caseTag)
%COMPARE_B6_JACOBI_WITH_BGK_REFERENCE Compare B6 V+A/Jacobi with BGK reference.
%
% Default example:
%   row = compare_B6_jacobi_with_bgk_reference(300,0.006,Inf,501,'free');
%
% For faster tests use Nx=200, Nv=301.

if nargin < 1 || isempty(Nx), Nx = 300; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end
if nargin < 3 || isempty(tau), tau = Inf; end
if nargin < 4 || isempty(Nv), Nv = 501; end
if nargin < 5 || isempty(caseTag)
    if isinf(tau), caseTag = 'free'; else, caseTag = sprintf('tau_%g',tau); end
end

outdir = case_output_dir();

cfg = default_riemann_config();
cfg.Nx = Nx;
cfg.tEnd = tEnd;
cfg.tau = tau;
cfg.CFL = 0.35;
cfg.maxSteps = 80000;
cfg.minDt = 1e-12;
cfg.maxWaveSpeed = 1e8;
cfg.printEvery = 100;
cfg.saveFigures = false;

x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';
U0 = initialize_riemann_state(x,cfg,6);

fprintf('\nB6 V+A/Jacobi vs BGK reference: %s\n',caseTag);

tic;
sol = solve_moment_1d(U0,x,cfg.tEnd,cfg.tau,'B6',cfg.CFL,'jacobi',cfg);
cpuMom = toc;

[rhoM,uM,thetaM] = raw_to_primitive(sol.U);
WM = raw_to_normalized(sol.U,5);

ref = bgk_reference_1d(cfg,'Nx',Nx,'Nv',Nv,'tEnd',tEnd,'tau',tau, ...
    'CFL',0.8,'printEvery',500);

rhoB = ref.rho;
uB = ref.u;
thetaB = ref.theta;
WB = ref.W;

errRhoInf = max(abs(rhoM-rhoB));
errRhoRel1 = mean(abs(rhoM-rhoB))/max(mean(abs(rhoB)),eps);
errUInf = max(abs(uM-uB));
errThetaRel1 = mean(abs(thetaM-thetaB))/max(mean(abs(thetaB)),eps);
errW3Inf = max(abs(WM(4,:)-WB(4,:)));

fprintf('\nB6 vs BGK difference, %s\n',caseTag);
fprintf('  moment status=%s, steps=%d, cpu=%.2fs\n',sol.status,sol.nSteps,cpuMom);
fprintf('  BGK    status=%s, steps=%d, cpu=%.2fs\n',ref.status,ref.nSteps,ref.cpu);
fprintf('  Linf rho difference       = %.4e\n',errRhoInf);
fprintf('  relative L1 rho difference= %.4e\n',errRhoRel1);
fprintf('  Linf velocity difference  = %.4e\n',errUInf);
fprintf('  relative L1 theta diff    = %.4e\n',errThetaRel1);
fprintf('  Linf W3 difference        = %.4e\n',errW3Inf);

rhoFile = fullfile(outdir,sprintf('B6_vs_BGK_%s_density_Nx%d_Nv%d.png',caseTag,Nx,Nv));
uFile = fullfile(outdir,sprintf('B6_vs_BGK_%s_velocity_Nx%d_Nv%d.png',caseTag,Nx,Nv));
thetaFile = fullfile(outdir,sprintf('B6_vs_BGK_%s_theta_Nx%d_Nv%d.png',caseTag,Nx,Nv));
matFile = fullfile(outdir,sprintf('B6_vs_BGK_%s_Nx%d_Nv%d.mat',caseTag,Nx,Nv));

figure('Name',['B6 vs BGK density ',caseTag]);
plot(x,rhoM/cfg.left.rho,'LineWidth',1.5,'DisplayName','B6 V+A/Jacobi'); hold on; grid on;
plot(ref.x,rhoB/cfg.left.rho,'--','LineWidth',1.5,'DisplayName','BGK reference');
xlabel('x'); ylabel('\rho/\rho_L');
title(sprintf('B6 V+A/Jacobi vs BGK density: %s',caseTag));
legend('Location','best');
save_clean_figure(gcf,rhoFile);

figure('Name',['B6 vs BGK velocity ',caseTag]);
plot(x,uM,'LineWidth',1.5,'DisplayName','B6 V+A/Jacobi'); hold on; grid on;
plot(ref.x,uB,'--','LineWidth',1.5,'DisplayName','BGK reference');
xlabel('x'); ylabel('u');
title(sprintf('B6 V+A/Jacobi vs BGK velocity: %s',caseTag));
legend('Location','best');
save_clean_figure(gcf,uFile);

figure('Name',['B6 vs BGK theta ',caseTag]);
plot(x,thetaM,'LineWidth',1.5,'DisplayName','B6 V+A/Jacobi'); hold on; grid on;
plot(ref.x,thetaB,'--','LineWidth',1.5,'DisplayName','BGK reference');
xlabel('x'); ylabel('\theta');
title(sprintf('B6 V+A/Jacobi vs BGK temperature: %s',caseTag));
legend('Location','best');
save_clean_figure(gcf,thetaFile);

save(matFile,'x','sol','rhoM','uM','thetaM','WM','ref','cfg','cpuMom');

row = table(string(caseTag),Nx,Nv,tEnd,tau,string(sol.status),string(ref.status), ...
    sol.nSteps,ref.nSteps,cpuMom,ref.cpu,errRhoInf,errRhoRel1,errUInf, ...
    errThetaRel1,errW3Inf,string(rhoFile),string(uFile),string(thetaFile),string(matFile), ...
    'VariableNames', {'caseTag','Nx','Nv','tEnd','tau','statusMoment','statusBGK', ...
    'stepsMoment','stepsBGK','cpuMoment','cpuBGK','LinfRho','RelL1Rho', ...
    'LinfU','RelL1Theta','LinfW3','densityFigure','velocityFigure', ...
    'thetaFigure','matFile'});
end
