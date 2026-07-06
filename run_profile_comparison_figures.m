function T = run_profile_comparison_figures(Nx,Nv,tEnd)
%RUN_PROFILE_COMPARISON_FIGURES Generate paper-style Riemann profile figures.
%
% The figures compare B4, B6, B8, and a discrete-velocity BGK reference for
% the free-molecular, transition, and continuum regimes.

if nargin < 1 || isempty(Nx), Nx = 180; end
if nargin < 2 || isempty(Nv), Nv = 501; end
if nargin < 3 || isempty(tEnd), tEnd = 0.006; end

outdir = case_output_dir();
taus = [Inf, 1e-3, 1e-6];
tags = {'free','transition','continuum'};
labels = {'free molecular','transition','continuum'};

rows = cell(numel(taus),1);
for k = 1:numel(taus)
    cfg = default_riemann_config();
    cfg.Nx = Nx;
    cfg.tEnd = tEnd;
    cfg.tau = taus(k);
    cfg.CFL = 0.30;
    cfg.maxSteps = 100000;
    cfg.minDt = 1e-12;
    cfg.maxWaveSpeed = 1e8;
    cfg.printEvery = 0;
    cfg.saveFigures = false;

    x = linspace(cfg.xmin,cfg.xmax,cfg.Nx)';

    fprintf('\nProfile comparison case: %s\n',tags{k});

    U4 = initialize_riemann_state(x,cfg,4);
    U6 = initialize_riemann_state(x,cfg,6);
    U8 = initialize_riemann_state(x,cfg,8);

    tic; sol4 = solve_moment_1d(U4,x,tEnd,taus(k),'B4',cfg.CFL,'jacobi',cfg); cpu4 = toc;
    tic; sol6 = solve_moment_1d(U6,x,tEnd,taus(k),'B6',cfg.CFL,'jacobi',cfg); cpu6 = toc;
    tic; sol8 = solve_moment_1d(U8,x,tEnd,taus(k),'B8',cfg.CFL,'jacobi',cfg); cpu8 = toc;

    ref = bgk_reference_1d(cfg,'Nx',Nx,'Nv',Nv,'tEnd',tEnd,'tau',taus(k), ...
        'CFL',0.8,'printEvery',0);

    [rho4,u4,theta4] = raw_to_primitive(sol4.U); W4 = raw_to_normalized(sol4.U,3);
    [rho6,u6,theta6] = raw_to_primitive(sol6.U); W6 = raw_to_normalized(sol6.U,5);
    [rho8,u8,theta8] = raw_to_primitive(sol8.U); W8 = raw_to_normalized(sol8.U,7);

    rhoB = ref.rho; uB = ref.u; thetaB = ref.theta; WB = ref.W;

    figFile = fullfile(outdir,sprintf('profile_comparison_%s_Nx%d_Nv%d_t%.0e.png',tags{k},Nx,Nv,tEnd));
    figure('Name',['Profile comparison ',tags{k}],'Position',[100,100,980,760]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(x,rho4/cfg.left.rho,'LineWidth',1.2,'DisplayName','B4'); hold on; grid on;
    plot(x,rho6/cfg.left.rho,'LineWidth',1.2,'DisplayName','B6');
    plot(x,rho8/cfg.left.rho,'LineWidth',1.2,'DisplayName','B8');
    plot(ref.x,rhoB/cfg.left.rho,'k--','LineWidth',1.2,'DisplayName','BGK');
    ylabel('\rho/\rho_L');
    title(sprintf('Riemann profiles, %s regime',labels{k}));
    legend('Location','best');

    nexttile;
    plot(x,u4,'LineWidth',1.2,'DisplayName','B4'); hold on; grid on;
    plot(x,u6,'LineWidth',1.2,'DisplayName','B6');
    plot(x,u8,'LineWidth',1.2,'DisplayName','B8');
    plot(ref.x,uB,'k--','LineWidth',1.2,'DisplayName','BGK');
    ylabel('u');

    nexttile;
    plot(x,W4(4,:),'LineWidth',1.2,'DisplayName','B4'); hold on; grid on;
    plot(x,W6(4,:),'LineWidth',1.2,'DisplayName','B6');
    plot(x,W8(4,:),'LineWidth',1.2,'DisplayName','B8');
    plot(ref.x,WB(4,:),'k--','LineWidth',1.2,'DisplayName','BGK');
    xlabel('x'); ylabel('W_3');

    save_clean_figure(gcf,figFile);

    rows{k} = table(string(tags{k}),Nx,Nv,tEnd,taus(k), ...
        string(sol4.status),string(sol6.status),string(sol8.status),string(ref.status), ...
        sol4.nSteps,sol6.nSteps,sol8.nSteps,ref.nSteps,cpu4,cpu6,cpu8,ref.cpu, ...
        max(abs(rho4-rhoB)),max(abs(rho6-rhoB)),max(abs(rho8-rhoB)), ...
        mean(abs(rho4-rhoB))/max(mean(abs(rhoB)),eps), ...
        mean(abs(rho6-rhoB))/max(mean(abs(rhoB)),eps), ...
        mean(abs(rho8-rhoB))/max(mean(abs(rhoB)),eps), ...
        string(figFile), ...
        'VariableNames', {'caseTag','Nx','Nv','tEnd','tau', ...
        'statusB4','statusB6','statusB8','statusBGK', ...
        'stepsB4','stepsB6','stepsB8','stepsBGK','cpuB4','cpuB6','cpuB8','cpuBGK', ...
        'LinfRhoB4','LinfRhoB6','LinfRhoB8','RelL1RhoB4','RelL1RhoB6','RelL1RhoB8', ...
        'figureFile'});

    save(fullfile(outdir,sprintf('profile_comparison_%s_Nx%d_Nv%d_t%.0e.mat',tags{k},Nx,Nv,tEnd)), ...
        'cfg','x','sol4','sol6','sol8','ref','rho4','rho6','rho8','rhoB', ...
        'u4','u6','u8','uB','theta4','theta6','theta8','thetaB','W4','W6','W8','WB');
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,sprintf('profile_comparison_summary_Nx%d_Nv%d_t%.0e.csv',Nx,Nv,tEnd));
writetable(T,csvFile);
disp(T(:,{'caseTag','statusB4','statusB6','statusB8','statusBGK', ...
    'LinfRhoB4','LinfRhoB6','LinfRhoB8','RelL1RhoB4','RelL1RhoB6','RelL1RhoB8'}));
fprintf('Profile comparison CSV saved to %s\n',csvFile);
end
