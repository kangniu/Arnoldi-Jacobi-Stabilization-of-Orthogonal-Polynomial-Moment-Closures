function T = run_B8_full_regime_sweep(Nx,tEnd)
%RUN_B8_FULL_REGIME_SWEEP Full B8 V+A/Jacobi regime sweep.

if nargin < 1 || isempty(Nx), Nx = 180; end
if nargin < 2 || isempty(tEnd), tEnd = 0.006; end

outdir = case_output_dir();
taus = [Inf, 1e-3, 1e-6];
tags = {'free','transition','continuum'};

fprintf('============================================================\n');
fprintf('B8 full-regime V+A/Jacobi sweep\n');
fprintf('Nx=%d, tEnd=%.6e\n',Nx,tEnd);
fprintf('============================================================\n');

rows = cell(numel(taus),1);
for k = 1:numel(taus)
    rows{k} = run_B8_jacobi_benchmark_case(Nx,tEnd,taus(k),tags{k});
end

T = vertcat(rows{:});
csvFile = fullfile(outdir,sprintf('B8_full_regime_sweep_Nx%d_t%.0e.csv',Nx,tEnd));
writetable(T,csvFile);

fprintf('\nB8 full-regime sweep summary\n');
disp(T(:,{'caseTag','Nx','tEnd','tau','status','steps','finalTime','cpu', ...
          'minRho','minTheta','maxAbsW3','maxAbsW7'}));
fprintf('Summary CSV saved to:\n  %s\n',csvFile);
end
