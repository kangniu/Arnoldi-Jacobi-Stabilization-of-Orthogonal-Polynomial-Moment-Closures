function T = run_B6_vs_BGK_regime_demo(Nx,Nv,tEnd)
%RUN_B6_VS_BGK_REGIME_DEMO B6 V+A/Jacobi against BGK reference in 3 regimes.
%
% Defaults are moderate so the script can run on a laptop:
%   Nx=300, Nv=501, tEnd=0.006.
%
% Increase Nx/Nv for final figures after verifying runtime.

if nargin < 1 || isempty(Nx), Nx = 300; end
if nargin < 2 || isempty(Nv), Nv = 501; end
if nargin < 3 || isempty(tEnd), tEnd = 0.006; end

taus = [Inf, 1e-3, 1e-6];
tags = {'free','transition','continuum'};

rows = cell(numel(taus),1);
for k = 1:numel(taus)
    rows{k} = compare_B6_jacobi_with_bgk_reference(Nx,tEnd,taus(k),Nv,tags{k});
end

T = vertcat(rows{:});
outdir = case_output_dir();
csvFile = fullfile(outdir,sprintf('B6_vs_BGK_regime_demo_Nx%d_Nv%d_t%.0e.csv',Nx,Nv,tEnd));
writetable(T,csvFile);

fprintf('\nB6 vs BGK regime demo summary\n');
disp(T(:,{'caseTag','Nx','Nv','tEnd','tau','statusMoment','statusBGK', ...
    'cpuMoment','cpuBGK','LinfRho','RelL1Rho','LinfU','RelL1Theta','LinfW3'}));
fprintf('CSV saved to %s\n',csvFile);
end
