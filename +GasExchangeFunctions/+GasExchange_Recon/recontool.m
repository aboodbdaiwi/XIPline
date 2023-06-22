function [varargout] = recontool(varargin)
%% Setup Figure Layout
% Create a figure without toolbar and menubar.
figure_handle = figure('Name','Recon Tool v1.0','NumberTitle','off');
figure_handle.MenuBar = 'none'; % hide menu bar

% Setup resize function
set(figure_handle,'ResizeFcn',@resize);

% Create structure of handles
h = guihandles(figure_handle);

%Create image uipanel
root_panel = uipanel('BorderType','none','Units','pixels','Parent',figure_handle);

% Create Tab group
root_tabgroup = uitabgroup(root_panel);

%% Create ReconModel tab
reconModeltab = uitab(root_tabgroup,'Title','1) Recon Model','Units','pixels');
GasExchange_Recon.ReconModel.ui.addUI(figure_handle,reconModeltab);

% Save structure of handles
h = guidata(figure_handle);

% Nicely give back output
if(nargout > 0)
    varargout{1} = figure_handle;
end

    function resize(source,callbackdata)
        % Get size of parent window
        newSize = get(source,'Position')
        
        set(root_panel,'Position',[1 1 newSize(3) newSize(4)]);
    end
end