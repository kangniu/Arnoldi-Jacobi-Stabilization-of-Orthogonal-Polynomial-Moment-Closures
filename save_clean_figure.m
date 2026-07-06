function save_clean_figure(fig,filename)
%SAVE_CLEAN_FIGURE Export figure without MATLAB axes toolbar artifacts.
%
% Usage:
%   save_clean_figure(gcf,'figures/name.png')

[folder,~,~] = fileparts(filename);
if ~isempty(folder) && ~exist(folder,'dir')
    mkdir(folder);
end

try
    % Hide axes toolbars if present.
    ax = findall(fig,'type','axes');
    for k = 1:numel(ax)
        try
            ax(k).Toolbar.Visible = 'off';
        catch
        end
    end

    exportgraphics(fig,filename,'Resolution',300);
catch
    saveas(fig,filename);
end
end
