function h = addUI(figure_handle, parent)
% get structure of handles
h = guidata(figure_handle);

% Create a dropdown menu for the kernel function
reconModelDropdown = uicontrol(parent, 'Style', 'popup',...
    'String', {'Simple Gridding (default)','Conjugate Gradient'},...
    'Units','pixels', 'Callback', @setReconModel);

% Setup resize function
set(parent,'ResizeFcn',@resize);

% Save structure of handles
guidata(figure_handle);

    function setReconModel(source,callbackdata)
        % Clear any children objects
        delete(get(reconModelDropdown,'Children'));
        clear cggui; 
        
        switch source.Value
            case 1
                % No extra options for simple recon
            case 2
                % load cg GUI
                cggui = load_cg_gui();
            otherwise
                error('There shouldnt be this many elements in the list');
        end
    end

    function resize(source,callbackdata)
        % Get size of parent window
        parentSize = get(parent,'Position');
        
        set(reconModelDropdown,'Position',[1 1 parentSize(3) 80]);
        
    end

    function cggui = load_cg_gui()
        % Get size of parent window
        parentSize = get(parent,'Position');
        
        cggui = uicontrol(parent, 'Style', 'edit',...
            'String', '10', 'Units','pixels','Position',[1 1 parentSize(3) 10]);
    end
end