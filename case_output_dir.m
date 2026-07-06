function outdir = case_output_dir()
%CASE_OUTPUT_DIR Robust output directory relative to this code package.
root = fileparts(mfilename('fullpath'));
outdir = fullfile(root,'figures');
if ~exist(outdir,'dir')
    [ok,msg] = mkdir(outdir);
    if ~ok
        error('Could not create output directory: %s. Message: %s',outdir,msg);
    end
end
end
