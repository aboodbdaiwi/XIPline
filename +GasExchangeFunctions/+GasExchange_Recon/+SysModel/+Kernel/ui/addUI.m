function h = addUI(parent)
    % Get size of parent window 
    parentSize = get(parent,'Position');
    
    % Create a dropdown menu for the kernel function
    h = uicontrol(parent, 'Style', 'popup',...
           'String', {'Gaussian (default)','Kaiser-Bessel','Sinc'},...
           'Position', [20 340 100 50],...
           'Callback', @setKernelFunction);
       
    function setKernelFunction(source,callbackdata)
        switch source.Value
            case 1
                % Gaussian
                disp('gaussian');
            case 2
                % Kaiser-Bessel
                disp('kb');
            case 3
                % Sinc
                disp('sinc');
            otherwise 
                error('There shouldnt be this many elements in the list');
        end
    end
end