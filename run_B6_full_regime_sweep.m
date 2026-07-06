function T = run_B6_full_regime_sweep(Nx,tEnd)
%RUN_B6_FULL_REGIME_SWEEP Full B6 V+A/Jacobi regime sweep.
%
% Default:
%   Nx   = 600
%   tEnd = 0.006
%
% Runs:
%   tau = Inf   : free molecular
%   tau = 1e-3  : transition
%   tau = 1e-6  : continuum
%
% This is the next benchmark after v8 grid/time verification.

if nargin < 1 || isempty(Nx), Nx = 600; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end

outdir = case_output_dir();

taus = [Inf, 1e-3, 1e-6];
tags = {'free','transition','continuum'};

fprintf('============================================================\n');
fprintf('B6 full-regime V+A/Jacobi sweep\n');
fprintf('Nx=%d, tEnd=%.6e\n',Nx,tEnd);
fprintf('============================================================\n');

rows = cell(numel(taus),1);
for k = 1:numel(taus)
    rows{k} = run_B6_jacobi_benchmark_case(Nx,tEnd,taus(k),tags{k});
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,sprintf('B6_full_regime_sweep_Nx%d_t%.0e.csv',Nx,tEnd));
writetable(T,csvFile);

fprintf('\nB6 full-regime sweep summary\n');
disp(T(:,{'caseTag','Nx','tEnd','tau','status','steps','finalTime','cpu', ...
          'minRho','minTheta','minMargin','maxAbsW3','maxAbsW4','maxAbsW5'}));
fprintf('Summary CSV saved to:\n  %s\n',csvFile);

end
