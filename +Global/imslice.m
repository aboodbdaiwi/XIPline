% IMSLICE - Displays an N-dimmensional dataset one slice at a time.
%
% If the dataset is 2D, the dataset is shown as an image.
% If the dataset is 3D, the volume can be scrolled through any dimmension
% If the dataset is ND, the volume can be scrolled through a slice in any
%      of the first 3 dimensions (assumed to be spatial), then there are
%      sliders to change the other dimmensions.
%
% figure_handle = imslice(data_in)
%
% figure_handle = imslice(data_in, title)
%
% Copyright: 2012 Scott Haile Robertson.
% Website: www.ScottHaileRobertson.com
%
% $Revision: 2.0 $
% $Date: Nov 28, 2012 $
%
% 4/5/2013 Works on N-dimmensions with uneven dimmensions
%
% TODO 1) Install script
%      2) Allow non integer pixel dimmensions
%      3) linking behavior
%      4) Make all text white
%

function [varargout] = imslice(data_in, varargin)

if(nargin > 1)
    title_name = varargin{1};
else
    title_name = 'Image';
end


%% Setup Figure Layout
% Create a figure without toolbar and menubar.
figure_handle = gcf();
%delete(figure_handle);
figure(figure_handle);
set(figure_handle,'defaultTextColor',[1 1 1]);
set(figure_handle,'DefaultAxesColorOrder',[1 1 1]);
set(figure_handle,'defaultLineColor',[1 1 1]);
set(figure_handle,'defaultAxesColor',[1 1 1]);
set(figure_handle,'defaultAxesColor',[1 1 1]);
set(figure_handle,'defaultAxesColor',[1 1 1]);
set(figure_handle,'defaultAxesColor',[1 1 1]);



% Give plot title
set(figure_handle, 'Name',title_name,'NumberTitle','off');

% Create structure of handles
h = guihandles(figure_handle);

%Create image uipanel
% h.im_uipanel = uipanel('BorderType','none','Units','pixels','BackgroundColor','black','Parent',figure_handle);
h.im_uipanel = uipanel('BorderType','none','Units','pixels','Parent',figure_handle);

% Create image axes
h.plot_axes = axes('Parent',h.im_uipanel);

if(isempty(data_in))
    return;
end

%Save the input data
if(isa(data_in,'Volume'))
    h.vol = data_in;
else
    h.vol = Global.Volume(data_in);
end

% Save structure of handles
guidata(figure_handle,h);

% Nicely give back output
if(nargout > 0)
    varargout{1} = figure_handle;
end

%% Check data type
switch h.vol.NDims
    case 2
        % Display image in the axes and get a handle to the image
        h.hImage = imagesc(h.vol.Data, 'Parent',h.plot_axes);
    otherwise
        h = buildNDGUI(h, figure_handle);
end

% Set default plot axes styles
colormap(gray);
axis image;
min_v = min(h.vol.Data(:));
max_v = max(h.vol.Data(:));
if(min_v == max_v)
    if(max_v > 0)
        min_v = 0;
    else
        max_v = 1;
    end
end
set(h.plot_axes,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[],...
    'CLimMode','manual','CLim',[min_v max_v], 'xcolor','green','ycolor','green');

set(figure_handle,'MenuBar','figure');  % Display standard menu bar menus.
% Add Pixel Information tool, specifying image as parent
% h.hpixinfo = impixelinfo(h.im_uipanel, h.hImage);

%% Display tools in toolbar
set(figure_handle,'Toolbar','figure');
h.toolbar = findall(figure_handle,'Type','uitoolbar');

% Get roots for where to find icons
[iconRoot, iconRootMATLAB] = ipticondir;

% Profile tool
uipushtool(h.toolbar,'CData',imread([iconRoot filesep() 'profile.png']),...
    'Separator','on','TooltipString','Generate line profile',...
    'HandleVisibility','off', 'ClickedCallback',@profile_callback);

    function profile_callback(hObject,eventdata)
        improfile(h.hImage);
    end

% Window width/window level tool
uipushtool(h.toolbar,'CData',imread([iconRoot filesep() 'tool_contrast.png']),...
    'Separator','on','TooltipString','Adjust window width/level',...
    'HandleVisibility','off', 'ClickedCallback',@ww_wl_callback);
    function ww_wl_callback(hObject,eventdata)
        imcontrast(h.hImage);
    end

% Save structure of handles
guidata(figure_handle,h);
end

%% Functions
function h = buildNDGUI(h,figure_handle)
%Change default resize
set(figure_handle, 'ResizeFcn',@figResize);

%Create scroll uipanel
h.scroll_uipanel = uipanel('BorderType','none','Units','pixels','Parent',figure_handle);

% Save structure of handles
guidata(figure_handle,h);

%Add sliderbar to show other slices
h.slice_slider = uicontrol(h.scroll_uipanel, 'Style', 'Slider', ...
    'String', 'Image No.', 'Callback', @slider_callback, ...
    'Units', 'pixels');

h.dimension_selection = uicontrol(h.scroll_uipanel, 'Style', 'popup', 'String', ...
    [h.vol.DimNames{1} '|' h.vol.DimNames{2} '|' h.vol.DimNames{3}],...
    'Units','pixels','Callback',@changeDim);
h.dim_val_txt = uicontrol(h.scroll_uipanel, 'Style', 'edit', ...
    'String', '','Callback',@slice_pos_callback,'Units', 'pixels');

% Create a slider bar, dimmension title, and index update string for each
% extra dimmension
if(h.vol.NDims > 3)
    for i=4:(h.vol.NDims)
        h.dim_slider(i-3) = uicontrol(h.scroll_uipanel, 'Style', 'slider', ...
            'String', ['Dim' num2str(i)], 'Callback', @slider_callback, ...
            'Units', 'pixels');
        h.dim_name_txt(i-3) = uicontrol(h.scroll_uipanel, 'Style', 'text', ...
            'String', h.vol.DimNames{i}, ...
            'Units', 'pixels');
        h.dim_val_extra_txt(i-3) = uicontrol(h.scroll_uipanel, 'Style', 'edit', ...
            'String', 'test','Callback',@slice_pos_callback,'Units', 'pixels');
        
        max_v = h.vol.Dims(3);
        step_v = [1 min(round(max_v/4),20)]/(max_v-1);
        
        set(h.dim_slider(i-3),'Max',max_v,'Min',1,'Value',1,'SliderStep',step_v);
    end
end
%Initialize Scroll bar
h.default_scroll_vals = ones(1,3);
updateScrollRange();

% Save structure of handles
guidata(figure_handle,h);

%Initialize Image
updateImage();

% Save structure of handles
guidata(figure_handle,h);
    function slice_pos_callback(hObject,eventdata)
        %h_ = guidata(hObject);
        slice_str = regexp(get(hObject,'String'),'/','split');
        slice_str  = slice_str{1};
        if(hObject == h.dim_val_txt)
            dim = get(h.dimension_selection, 'Value');
            ndims = h.vol.Dims(dim);
            if(~isnumeric(slice_str))
                newSlice = str2num(slice_str);
                if((newSlice > 0)&&(newSlice<=ndims))% Check that its within bounds
                    set(h.slice_slider, 'Value',newSlice);
                    updateImage();
                end
            end
            % Make sure the format is nice
            formatted_string = [num2str(get(h.slice_slider, 'Value')) '/' num2str(ndims)];
            set(hObject,'String',formatted_string);
                                
        elseif(h.vol.NDims > 3)
            for i=1:(h.vol.NDims-3)
                if(h.dim_val_extra_txt(i) == hObject)
                    dim = get(h.dim_name_txt(i), 'Value');
                    ndims = h.vol.Dims(dim+3);
                    
                    if(~isnumeric(slice_str))
                        newSlice = str2num(slice_str);
                        if((newSlice > 0)&&(newSlice<=ndims))% Check that its within bounds
                            set(h.dim_slider(i), 'Value',newSlice)
                            updateImage();
                        end
                    end
                    
                    % Make sure the format is nice
                    formatted_string = [num2str(get(h.dim_slider(i), 'Value')) '/' num2str(ndims)];
                    set(hObject,'String',formatted_string);
                end
            end
        end
    end

    function slider_callback(hObject,eventdata)
        updateImage();
    end

    function changeDim(hObject,eventdata)
        %Change scrollbar range
        updateScrollRange();
        
        updateImage();
    end

    function updateScrollRange()
        dim = get(h.dimension_selection, 'Value');
        maxVal = h.vol.Dims(dim);
        curVal = h.default_scroll_vals(dim);
        stepVal = [1 min(round(maxVal/4),20)]/(maxVal-1);
        set(h.slice_slider,'Max',maxVal,'Min',1,'Value',curVal,'SliderStep',stepVal);
    end

    function updateImage()
        slice = round(get(h.slice_slider, 'Value'));
        dim = get(h.dimension_selection, 'Value');
        
        if(h.vol.NDims > 3)
            main_slice = slice;
            main_dim = dim;
            slice = zeros(1,h.vol.NDims-2);
            dim = zeros(1,h.vol.NDims-2);
            slice(1) = main_slice;
            dim(1) = main_dim;
            
            for i=1:(h.vol.NDims-3)
                slice(i+1) = round(get(h.dim_slider(i), 'Value'));
                dim(i+1) = i+3;
                new_string = [num2str(slice(i+1)) '/' num2str(h.vol.Dims(dim(i+1)))];
                % Only update string if its changed
                if(~strcmp(new_string,get(h.dim_val_extra_txt(i),'String')))
                    set(h.dim_val_extra_txt(i),'String',new_string);
                end
            end
            % Set text only if its changed
            new_main_string = [num2str(main_slice) '/' num2str(h.vol.Dims(main_dim))];
            if(~strcmp(new_main_string,get(h.dim_val_txt,'String')))
                set(h.dim_val_txt,'String',[num2str(main_slice) '/' num2str(h.vol.Dims(main_dim))]);
            end
        else
            % Only change text if its changed
            new_main_string = [num2str(slice) '/' num2str(h.vol.Dims(dim))];
            if(~strcmp(new_main_string,get(h.dim_val_txt,'String')))
                set(h.dim_val_txt,'String',new_main_string);
            end
        end
        h.default_scroll_vals(dim(1)) = slice(1);
        if(isfield(h,'hImage'))
            h.hImage =  Global.showImageSlice(h.vol,slice,dim,h.hImage);
        else
            h.hImage = Global.showImageSlice(h.vol,slice,dim);
        end
    end

    function figResize(hObject,eventdata)
        fig_sz = get(figure_handle, 'Position');
        
        %Build panel layout
        single_scroll_height = 25;
        scroll_panel_height = single_scroll_height*(h.vol.NDims-2);
        set(h.scroll_uipanel, 'Position',[1,1,fig_sz(3),scroll_panel_height]);
        set(h.im_uipanel, 'Position',[1,2+scroll_panel_height,fig_sz(3),fig_sz(4)-scroll_panel_height]);
        
        %Build Slider panel layout
        scroll_sz = get(h.scroll_uipanel, 'Position');
        dimmension_sel_width = 50; %50
        position_width = 60; %60 
        set(h.dim_val_txt, 'Position', [1,scroll_panel_height-single_scroll_height,position_width,single_scroll_height]);
        set(h.slice_slider, 'Position', [2+position_width,scroll_panel_height-single_scroll_height,scroll_sz(3)-position_width-dimmension_sel_width,single_scroll_height]);
        set(h.dimension_selection, 'Position',[2+scroll_sz(3)-dimmension_sel_width,scroll_panel_height-single_scroll_height,dimmension_sel_width,single_scroll_height]);

        %Add slider bars for extra dims
        if(h.vol.NDims > 3)
            for i=1:(h.vol.NDims-3)
                set(h.dim_val_extra_txt(i), 'Position', [1,scroll_panel_height-(i+1)*single_scroll_height,position_width,single_scroll_height]);
                set(h.dim_slider(i), 'Position', [1+position_width,scroll_panel_height-(i+1)*single_scroll_height,scroll_sz(3)-position_width-dimmension_sel_width,single_scroll_height]);
                set(h.dim_name_txt(i), 'Position', [2+scroll_sz(3)-dimmension_sel_width,scroll_panel_height-(i+1)*single_scroll_height,dimmension_sel_width,single_scroll_height]);
            end
        end
    end
set(figure_handle,'position',[560 350 560 421]); % resize the figure for 3 dim or more %%Abood 1/21/2020
end
