function write_latex_table_from_table(T,filename,caption,label)
%WRITE_LATEX_TABLE_FROM_TABLE Minimal LaTeX table writer.
%
% Numeric values are formatted in scientific notation when appropriate.

if nargin < 3, caption = ''; end
if nargin < 4, label = ''; end

[folder,~,~] = fileparts(filename);
if ~isempty(folder) && ~exist(folder,'dir')
    mkdir(folder);
end

fid = fopen(filename,'w');
if fid < 0
    error('Cannot open %s for writing.',filename);
end

fprintf(fid,'%% Auto-generated table\n');
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n');
if ~isempty(caption)
    fprintf(fid,'\\caption{%s}\n',caption);
end
if ~isempty(label)
    fprintf(fid,'\\label{%s}\n',label);
end
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{%s}\n',repmat('c',1,width(T)));
fprintf(fid,'\\toprule\n');

names = T.Properties.VariableNames;
for j = 1:numel(names)
    if j > 1, fprintf(fid,' & '); end
    fprintf(fid,'%s',escape_latex(names{j}));
end
fprintf(fid,'\\\\\n\\midrule\n');

for i = 1:height(T)
    for j = 1:width(T)
        if j > 1, fprintf(fid,' & '); end
        val = T{i,j};
        fprintf(fid,'%s',format_table_value(val));
    end
    fprintf(fid,'\\\\\n');
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table}\n');
fclose(fid);
fprintf('LaTeX table written to %s\n',filename);
end

function s = format_table_value(val)
if iscell(val)
    val = val{1};
end
if isstring(val)
    s = char(val);
elseif ischar(val)
    s = val;
elseif isnumeric(val)
    if isempty(val)
        s = '';
    elseif isscalar(val)
        if isinf(val)
            s = '$\\infty$';
        elseif abs(val) >= 1e4 || (abs(val) > 0 && abs(val) < 1e-3)
            s = sprintf('%.3e',val);
        else
            s = sprintf('%.4g',val);
        end
    else
        s = mat2str(val);
    end
else
    s = char(string(val));
end
s = escape_latex(s);
s = strrep(s,'\$\\infty\$','$\infty$');
end

function s = escape_latex(s)
s = char(s);
s = strrep(s,'_','\_');
s = strrep(s,'%','\%');
end
