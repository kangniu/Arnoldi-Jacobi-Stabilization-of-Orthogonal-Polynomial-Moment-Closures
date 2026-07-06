function make_paper_tables_from_results()
%MAKE_PAPER_TABLES_FROM_RESULTS Convert available CSV summaries to LaTeX.
%
% This utility looks in ./figures and writes tables to ./tables.

outdir = case_output_dir();
tableDir = fullfile(fileparts(outdir),'tables');
if ~exist(tableDir,'dir'), mkdir(tableDir); end

files = {
    'B6_grid_sweep_original_vs_jacobi.csv', ...
    'B8_grid_sweep_original_vs_jacobi.csv', ...
    'B8_charpoly_ad_vs_fd.csv', ...
    'B6_moderate_original_vs_jacobi.csv'
};

for k = 1:numel(files)
    csvFile = fullfile(outdir,files{k});
    if exist(csvFile,'file')
        T = readtable(csvFile);
        texFile = fullfile(tableDir,[erase(files{k},'.csv'),'.tex']);
        write_latex_table_from_table(T,texFile, ...
            strrep(erase(files{k},'.csv'),'_',' '), ...
            ['tab:',erase(files{k},'.csv')]);
    else
        fprintf('Skipping missing CSV: %s\n',csvFile);
    end
end

% Full-regime CSVs have Nx/tEnd in their names.
d = dir(fullfile(outdir,'B*_full_regime_sweep_*.csv'));
for k = 1:numel(d)
    csvFile = fullfile(d(k).folder,d(k).name);
    T = readtable(csvFile);
    texFile = fullfile(tableDir,[erase(d(k).name,'.csv'),'.tex']);
    write_latex_table_from_table(T,texFile, ...
        strrep(erase(d(k).name,'.csv'),'_',' '), ...
        ['tab:',erase(d(k).name,'.csv')]);
end

d = dir(fullfile(outdir,'B6_vs_BGK_*.csv'));
for k = 1:numel(d)
    csvFile = fullfile(d(k).folder,d(k).name);
    T = readtable(csvFile);
    texFile = fullfile(tableDir,[erase(d(k).name,'.csv'),'.tex']);
    write_latex_table_from_table(T,texFile, ...
        strrep(erase(d(k).name,'.csv'),'_',' '), ...
        ['tab:',erase(d(k).name,'.csv')]);
end

fprintf('Paper tables written to %s\n',tableDir);
end
