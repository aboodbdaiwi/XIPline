function varargout = proton_segmentation(varargin)
% PROTON_SEGMENTATION MATLAB code for proton_segmentation.fig
%      PROTON_SEGMENTATION, by itself, creates a new PROTON_SEGMENTATION or raises the existing
%      singleton*.
%
%      H = PROTON_SEGMENTATION returns the handle to a new PROTON_SEGMENTATION or the handle to
%      the existing singleton*.
%
%      PROTON_SEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROTON_SEGMENTATION.M with the given input arguments.
%
%      PROTON_SEGMENTATION('Property','Value',...) creates a new PROTON_SEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before proton_segmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to proton_segmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help proton_segmentation

% Last Modified by GUIDE v2.5 19-Oct-2017 13:34:53
%
% W. Zha @ July, 2016
%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @proton_segmentation_OpeningFcn, ...
    'gui_OutputFcn',  @proton_segmentation_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before proton_segmentation is made visible.
function proton_segmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to proton_segmentation (see VARARGIN)

% Choose default command line output for proton_segmentation
handles.output = hObject;
handles = initial_settings(hObject,handles);
update_figures(hObject, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes proton_segmentation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = proton_segmentation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in seg_proton.
function seg_proton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_proton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.proton(:)) ~= 0 && sum(handles.proton(:)) ~= 0
    if isfield(handles,'RGseeds')
          [handles.lungsmask,handles.lobemask, ~,handles.trachea,...
        handles.bodymask,handles.RGseeds] = segment_proton_gas(handles.proton,handles.gas,...
        handles.RGseeds,handles.RGthresh,handles.acquisition);
    else
    [handles.lungsmask,handles.lobemask, ~,handles.trachea,...
        handles.bodymask,handles.RGseeds] = segment_proton_gas(handles.proton,handles.gas,...
        [],handles.RGthresh,handles.acquisition);
    end
    
    handles.lungVolume = calculate_lung_volume(handles,handles.lungsmask);
    display(['lungvolume = ' num2str(handles.lungVolume)])
    update_figures(hObject, handles)
end

% --- Executes on button press in add_lung_region.
function add_lung_region_Callback(hObject, eventdata, handles)
% hObject    handle to add_lung_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Select data','Edit proton or gas?','Proton','Gas','Cancel','Proton');
switch button
    case 'Proton'
        figname = 'proton_fig';
    case 'Gas'
        figname = 'gas_fig';
    otherwise
        return;
end
is_recist = get(handles.is_recist,'Value');
if sum(handles.proton(:)) ~= 0
    if is_recist
        figure(1);
        nsl = round(get(handles.sliceSlider,'Value'));
%         imshow(handles.(lower(button))(:,:,nsl),[]);
%         imshow(handles.(lower(button))(:,:,nsl),[0,1]);
% imshow(handles.(lower(button))(:,:,nsl),[-1350,250]);
% imshow(handles.(lower(button))(:,:,nsl),[-700,1100]);
imshow(handles.(lower(button))(:,:,nsl),[handles.rmin, handles.rmax]);

        plot_contours(handles.lungsmask(:,:,nsl));
    else
        axes(handles.(figname));
    end
    % figure(1),
    %     zoom(1.2);
    
    % nsl = round(get(handles.sliceSlider,'Value'));
    
    % imshow(handles.(lower(button))(:,:,nsl),[0,1]);
    % plot_contours(handles.lungsmask(:,:,nsl));
    
    newRegion = roipoly;
    nsl = round(get(handles.sliceSlider,'Value'));
    lungsmask = handles.lungsmask;
    handles.lungsmask_backup = lungsmask;
    
    lungsmask(:,:,nsl) = lungsmask(:,:,nsl) | newRegion;
    lungsmask(lungsmask < 0) = 0;
    handles.lungsmask = lungsmask;
    handles.lungVolume = calculate_lung_volume(handles,handles.lungsmask);
    update_figures(hObject, handles)
end

% --- Executes on button press in livewire_edit.
function livewire_edit_Callback(hObject, eventdata, handles)
% hObject    handle to livewire_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% axes(handles.proton_fig);
button = questdlg('Select data','Livewire on proton or gas?','Proton','Gas','Cancel','Proton');
switch button
    case 'Proton'
        figname = 'proton_fig';
        image_data = handles.proton;
    case 'Gas'
        figname = 'gas_fig';
        image_data = handles.gas;
    otherwise
        return;
end

if sum(image_data(:)) ~= 0
    axes(handles.(figname));
    bwR = livewire(image_data(:,:,handles.nsl));
    %     bwL = livewire(image_data(:,:,handles.nsl));
    handles.lungsmask_backup = handles.lungsmask;
    handles.lungsmask(:,:,handles.nsl) = bwR ;%| bwL;
    handles.lungVolume = calculate_lung_volume(handles,handles.lungsmask);
    update_figures(hObject, handles)
end

% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
nSlices = handles.nSlices;
nsl = get(handles.sliceSlider,'Value');
if nsl <= 0 || nsl > nSlices
    return
elseif nsl < 1
    nsl = 1;
end
set(handles.sliceSlider,'Min',0);
set(handles.sliceSlider,'Max',nSlices);
step1 = 1/nSlices;
step2 = 2/nSlices;
set(handles.sliceSlider, 'SliderStep', [step1 step2]);

handles.nsl = round(nsl);
update_figures(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in save_mask.
function save_mask_Callback(hObject, eventdata, handles)
% hObject    handle to save_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.lungsmask(:)) ~= 0
    lungsmask = handles.lungsmask; %#ok<NASGU>
    handles.lungVolume = calculate_lung_volume(handles,handles.lungsmask);
    dlgtitle ='Save lungsmask to *.mat';
    if isfield(handles,'resultfile_name')
        default_name = strrep(handles.resultfile_name,'results','lungsmask');
    elseif isfield(handles,'gasFileString')
        default_name = char(strcat(handles.gasFileString,'lungsmask'));
    elseif isfield(handles,'lungsmask')
        default_name = char(strcat(handles.protonFileString,'results'));
    else
        default_name = 'lungsmask.mat';
    end
    [handles.maskfile_name,handles.output_path] = uiputfile('*.mat',dlgtitle,default_name);
    if ischar(handles.maskfile_name)&& ischar(handles.output_path)
        output_file = fullfile(handles.output_path,handles.maskfile_name);
        if exist(output_file,'file')
            backup_file = strrep(output_file,'lungsmask','lungsmask_backup');
            movefile(output_file,backup_file ,'f');
            fprintf(' save the existed mask file into %s.\n',backup_file);
            
            
        end
        save(output_file,'lungsmask');
        fprintf(' %s saved.\n',output_file);
    end
end
if sum(handles.LVmask(:))~=0
%     output_path ='K:\Wei\dl\sarp\dot_mat';
    mask_file = sprintf('%s_mask.mat',handles.protonFileString{1}(1:10));
    liverDome = handles.aortamask;
    LVblood = handles.LVmask;
    muscle = handles.Musclemask;
    trachea = handles.lungsmask;
    save(fullfile(output_path,mask_file),'liverDome','LVblood','muscle','trachea');
    fprintf('%s is saved.\n',fullfile(output_path,mask_file));
%     cd('K:\Wei\dl\sarp\dot_mat')
%     save
end

% --- Executes on button press in load_proton.
function load_proton_Callback(hObject, eventdata, handles)
% hObject    handle to load_proton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentPath = pwd;
[handles.proton, handles.protonPath, handles.protonFile,handles.info] = load_image_data(currentPath);
if sum(handles.proton(:)) ~= 0
    %%10/09/2017 UTE
% %     handles.proton = fliplr(imrotate(handles.proton,270));
    if ~isempty(handles.protonFile)
        set(handles.display_proton_file,'String',handles.protonFile);
    end
    if ~isempty(strfind(handles.protonFile,'.nii.gz'))
        handles.protonFileString = strrep(handles.protonFile,'.nii.gz','');
    elseif ~isempty(strfind(handles.protonFile,'.nii'))
        handles.protonsave_FileString = strrep(handles.protonFile,'.nii','');
    else
        handles.protonFileString = handles.protonFile;
    end
     handles.nSlices = size(handles.proton,3);
%     if handles.nSlices ==1
       
        
%         nsl = 1;
%         set(handles.sliceSlider,'Value',nsl);
%         set(handles.sliceSlider,'Max', handles.nSlices);
%         set(handles.sliceSlider,'Min', 0);
%         handles.lungsmask = zeros(size(handles.proton));
%         handles.lobemask = zeros(size(handles.proton));
        
%     elseif handles.nSlices < size(handles.proton,3)
%         handles.nSlices = size(handles.proton,3);
%         handles.lungsmask = zeros(size(handles.proton));
%         handles.lobemask = zeros(size(handles.proton));
        
%     end
    
    
    
    handles.default_size = [size(handles.proton,1),size(handles.proton,2)];
    if ~isempty(handles.info)
        
        handles.voxelSize = handles.info.SliceThickness * prod(handles.info.PixelSpacing)*1e-6;
    else
        handles.voxelSize = [];
    end
    
    is_recist = get(handles.is_recist,'Value');
    is_ute =get(handles.is_ute,'Value');
    if is_ute || is_recist
        num_iter=3;
        delta_t=1/7;
        if is_ute
            kappa = 10;%30;
            handles.proton = fliplr(imrotate(handles.proton,270));
            
        else
            kappa = 10;
        end
        data_flt = zeros(size(handles.proton));
        option=1;
        proton = handles.proton;
        parfor n=1:size(proton,3)
            data_flt(:,:,n)=anisodiff2D(proton(:,:,n),num_iter,delta_t,kappa,option);
        end
        handles.proton = data_flt;
        if is_recist
            handles.proton = scale_raw_ct(handles.proton);
% handles.proton = scale_raw_ct(handles.proton,-700,1500);
%             handles.proton = scale_raw_ct(handles.proton,-160,240);
%              handles.proton = scale_raw_ct(handles.proton,-450,1300);%lung
        end
        
    else
        is_coronal = get(handles.is_coronal,'Value');
        if is_coronal
            intensity_th = locate_histogram_landmark_intensity(handles.proton,.99,'uint16');
            handles.proton(handles.proton>intensity_th)=intensity_th;
        end
                handles.proton = (handles.proton-min(handles.proton(:)))...
                                /(max(handles.proton(:))-min(handles.proton(:)));
    end
    %     handles.voxelSize =[0.82,0.82,2.51]; %handles.info.ReconstructionDiameter/handles.info.Rows*handles.info.ReconstructionDiameter/handles.info.Columns*
    handles.default_size = [size(handles.proton,1),size(handles.proton,2)];
    update_figures(hObject, handles)
    handles.ori_position = handles.proton_fig.Position;
    handles.lungsmask = zeros(size(handles.proton));
        handles.lobemask = zeros(size(handles.proton));
                nsl = 1;
        set(handles.sliceSlider,'Value',nsl);
        set(handles.sliceSlider,'Max', handles.nSlices);
        set(handles.sliceSlider,'Min', 0);
        
    guidata(hObject, handles);
    
end

% --- Executes on button press in load_gas.
function load_gas_Callback(hObject, eventdata, handles)
% hObject    handle to load_gas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentPath = pwd;
[rawgas, handles.gasPath, handles.gasFile] = load_image_data(currentPath);
is_ute = get(handles.is_ute,'Value');

if ~isempty(rawgas)
    if ~isempty(handles.gasFile)
        set(handles.display_gas_file,'String',handles.gasFile);
    end
    if ~isempty(strfind(handles.gasFile,'.nii.gz'))
        handles.gasFileString = strrep(handles.gasFile,'.nii.gz','');
    elseif ~isempty(strfind(handles.gasFile,'.nii'))
        handles.gasFileString = strrep(handles.gasFile,'.nii','');
    else
        handles.gasFileString = handles.gasFile;
    end
    if is_ute
        handles.gas = fliplr(imrotate(double(rawgas),270));
    else
        handles.gas = double(rawgas);
    end
            handles.nSlices = size(handles.gas,3);
        handles.lungsmask = zeros(size(handles.gas));
        handles.lobemask = zeros(size(handles.gas));
%     if handles.nSlices ==1
%         
%         nsl = 1;
%         
%         set(handles.sliceSlider,'Value',nsl);
%         set(handles.sliceSlider,'Max',handles.nSlices);
        
%         handles.lungsmask = zeros(size(handles.gas));
%         handles.lobemask = zeros(size(handles.gas));
%     elseif handles.nSlices < size(handles.gas,3)
%         handles.nSlices = size(handles.gas,3);

        
%     end
    %% Coase estimate of gas mask
    if ~is_ute
        handles.gasmask = estimate_gasmask(handles.gas,handles.poor_SNR);
        
    end
    
    handles.default_size = [size(handles.gas,1),size(handles.gas,2)];
            nsl = 1;
        
        set(handles.sliceSlider,'Value',nsl);
        set(handles.sliceSlider,'Max',handles.nSlices);

    update_figures(hObject, handles)
end

%% -----------------------------------------------
function gasmask = estimate_gasmask(gas,is_poor_SNR)
nClusters = 4;
clustered_gas = kmeans_clustering(gas,nClusters);
if is_poor_SNR
    
    clustered_gas(clustered_gas < 3) = 0;
else
    clustered_gas(clustered_gas < 2) = 0;
end
clustered_gas(clustered_gas > 0) = 1;
[labeled_gas,nGasPieces] = bwlabeln(clustered_gas, 26);
gasmask = zeros(size(gas));
gasmask(clustered_gas == 1) = 1;
for np = 1:nGasPieces
    if sum(labeled_gas(:) == np) <100
        gasmask(labeled_gas == np) = 0;
    end
end

%% -------------------------------------------------



% --- Executes on button press in merge_proton_gas_mask.
function merge_proton_gas_mask_Callback(hObject, eventdata, handles)
% hObject    handle to merge_proton_gas_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.lungsmask(:)) ~= 0 && isfield(handles,'gasmask')
    
    gasmask = handles.gasmask;
    if ~isempty(handles.bodymask)
        gasmask = gasmask.*handles.bodymask;
    end
    if ~isempty(handles.trachea)
        gasmask = gasmask - handles.trachea;
        gasmask(gasmask < 0) = 0;
    end
    handles.lungsmask_backup = handles.lungsmask;
    handles.lungsmask = merge_he_into_proton(handles.lungsmask,gasmask);
%     handles.lungsmask = handles.lungsmask | gasmask;
    update_figures(hObject, handles)
    
end
% --- Executes on button press in remove_lung_region.
function remove_lung_region_Callback(hObject, eventdata, handles)
% hObject    handle to remove_lung_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Select data','Edit proton or gas?','Proton','Gas','Cancel','Proton');
switch button
    case 'Proton'
        figname = 'proton_fig';
    case 'Gas'
        figname = 'gas_fig';
    otherwise
        return;
end
is_recist = get(handles.is_recist,'Value');

if sum(handles.proton(:)) ~= 0
    if is_recist
        figure(1);
        nsl = round(get(handles.sliceSlider,'Value'));
%         imshow(handles.(lower(button))(:,:,nsl),[]);
        imshow(handles.(lower(button))(:,:,nsl),[handles.rmin,handles.rmax]);
        plot_contours(handles.lungsmask(:,:,nsl));
    else
        axes(handles.(figname));
    end
    % figure(1),
    %     zoom(1.2);
    %         zoom on;
    % pause(3);
    % zoom off;
    % nsl = round(get(handles.sliceSlider,'Value'));
    %
    % imshow(handles.(lower(button))(:,:,nsl),[0,1]);
    % plot_contours(handles.lungsmask(:,:,nsl));
    newRegion = roipoly;
    nsl = round(get(handles.sliceSlider,'Value'));
    lungsmask = handles.lungsmask;
    handles.lungsmask_backup = lungsmask;
    
    lungsmask(:,:,nsl) = lungsmask(:,:,nsl) - newRegion;
    lungsmask(lungsmask < 0) = 0;
    handles.lungsmask = lungsmask;
    update_figures(hObject, handles)
end


% --- Executes on button press in undo_edits.
function undo_edits_Callback(hObject, eventdata, handles)
% hObject    handle to undo_edits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'lungsmask_backup')
    handles.lungsmask = handles.lungsmask_backup;
    update_figures(hObject, handles)
    
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
slider_value = get(handles.sliceSlider,'Value');
slider_value = slider_value+ eventdata.VerticalScrollCount;
if slider_value<=0 || slider_value>handles.nSlices
    return
elseif slider_value<1
    slider_value=1;
end
set(handles.sliceSlider,'Value',slider_value);
handles.nsl=slider_value;
update_figures(hObject, handles)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_handles.
function save_handles_Callback(hObject, eventdata, handles)
% hObject    handle to save_handles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if sum(handles.lungsmask(:)) ~= 0
results = struct;

output_fields = segmentation_outputs;
dlgtitle = 'Save handles to *.mat';
if isfield(handles,'maskfile_name')
    default_name = strrep(handles.maskfile_name,'lungsmask','results');
elseif isfield(handles,'gasFileString')
    default_name = char(strcat(handles.gasFileString,'results'));
elseif isfield(handles,'protonFileString')
    default_name = char(strcat(handles.protonFileString,'results'));
else
    default_name = 'results.mat';
end
for nfield = 1:numel(output_fields)
    this_field = output_fields{nfield};
    if isfield(handles,this_field)
        results.(this_field) = handles.(this_field);
    end
end
%     output_path ='K:\Wei\dl\sarp\dot_mat';
%     cd(output_path);
[resultfile_name,output_path] = uiputfile('*.mat',dlgtitle,default_name);

result_file = fullfile(output_path,resultfile_name);
if exist(result_file,'file')
    backup_file = strrep(result_file,'.mat','_backup.mat');
    movefile(result_file,backup_file ,'f');
    fprintf('save existed result file into %s.\n',backup_file);
    
    
end

if output_path ~=0
save(result_file,'results','-v7.3');
fprintf(' %s saved.\n',result_file);
% cd('Y:\density\archive')
end
% end


% --- Executes on button press in load_result.
function load_result_Callback(hObject, eventdata, handles)
% hObject    handle to load_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
is_paired_OE = get(handles.is_OE,'Value');
is_ute = get(handles.is_ute,'Value');
[matfile,matfile_path] = uigetfile('*.mat','Select .mat file');
if matfile~=0
    if is_paired_OE
        load(fullfile(matfile_path, matfile));
        handles.proton = fliplr(imrotate(matched_struct.oxy_hi_v1,270));
        handles.gas = fliplr(imrotate(matched_struct.pse_v1,270));
        handles.lungsmask =  fliplr(imrotate(matched_struct.ventmap_v1,270));
%                 handles.proton = fliplr(imrotate(matched_struct.oxy_hi_v2,270));
%         handles.gas = fliplr(imrotate(matched_struct.pse_v2,270));
%         handles.lungsmask =  fliplr(imrotate(matched_struct.ventmap_v2,270));
        handles.lungsmask(handles.lungsmask>1)=0;
        handles.lobemask = handles.lungsmask;
        handles.nSlices = size(handles.proton,3);
        handles.protonFile = strrep(matfile,'.mat','');
    elseif is_ute
                load(fullfile(matfile_path, matfile));

        handles.proton = reg.oxy_hi;
        handles.gas = permute(reg.oe_reg_cor,[3,2,1]);
        maskfile = strrep(matfile,'oe_results','oxy_hi_lungmask');
        mask_folder = 'K:\Wei\oe\work\lungmasks';
        load(fullfile(mask_folder,maskfile));

        handles.lungsmask = lungsmask(cropped_range.y,cropped_range.x,cropped_range.z);
                handles.lobemask = handles.lungsmask;
        handles.nSlices = size(handles.proton,3);
        handles.protonFile = strrep(matfile,'.mat','');
    else
        proton_seg = load(fullfile(matfile_path, matfile));
        fields = segmentation_outputs;
        for nfield = 1:numel(fields)
            this_field = fields{nfield};
            if isfield(proton_seg.results,this_field)
                handles.(this_field) = proton_seg.results.(this_field);
            end
        end
    end
    clear proton_seg;
    set(handles.sliceSlider,'Value',1);
    set(handles.sliceSlider,'Min',0);
    set(handles.sliceSlider,'Max',handles.nSlices);
    
    handles.nsl = 1;
    if isfield(handles,'gasFile')
        set(handles.display_gas_file,'String',char(handles.gasFile));
    else
        set(handles.display_gas_file,'String',matfile);
    end
    if isfield(handles,'protonFile')
        set(handles.display_proton_file,'String',char(handles.protonFile));
    else
        set(handles.display_proton_file,'String',matfile);
        
    end
    if strcmp(handles.acquisition,'Axial')
        set(handles.is_axial,'Value', 1);
        set(handles.is_coronal,'Value', 0);
    end
    if handles.poor_SNR
        set(handles.is_poor_SNR,'Value', 1);
    end
    if isfield(handles,'enable_recist')
        set(handles.is_recist,'Value',handles.enable_recist);
    end
    if isfield(handles,'enable_ute')
        set(handles.is_ute,'Value',handles.enable_ute);
    end
    if isfield(handles,'RGthresh')
        set(handles.RGthreshold,'String',num2str(handles.RGthresh));
    end
    update_figures(hObject, handles)
    
end

function display_gas_file_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to display_gas_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.gasFile{1} = get(hObject,'String');
set(handles.display_gas_file,'String',handles.gasFile);
% Hints: get(hObject,'String') returns contents of display_gas_file as text
%        str2double(get(hObject,'String')) returns contents of display_gas_file as a double


% --- Executes during object creation, after setting all properties.
function display_gas_file_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to display_gas_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in smooth_mask.
function smooth_mask_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to smooth_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.lungsmask(:)) ~= 0
    handles.lungsmask_backup = handles.lungsmask;
    handles.lungsmask_unsmoothed = handles.lungsmask;
    handles.lungsmask = smooth_lungmask(handles.lungsmask,handles.acquisition);
    update_figures(hObject, handles)
end


function display_proton_file_Callback(hObject, eventdata, handles)
% hObject    handle to display_proton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_proton_file as text
%        str2double(get(hObject,'String')) returns contents of display_proton_file as a double
handles.protonFile{1} =  get(hObject,'String');
set(handles.display_proton_file,'String',handles.protonFile);
update_figures(hObject, handles)

% --- Executes during object creation, after setting all properties.
function display_proton_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_proton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = initial_settings(hObject,handles);
cla(handles.proton_fig,'reset');
cla(handles.gas_fig,'reset');

set(handles.display_proton_file,'String','display proton file');
set(handles.display_gas_file,'String','display gas file');

update_figures(hObject, handles)

%% -----------------------------------------------
function handles = initial_settings(hObject,handles)
output_fields = segmentation_outputs;
output_fields{numel(output_fields)+1} = 'lungsmask_backup';
for nfield = 1:numel(output_fields)
    this_field = output_fields{nfield};
    if isfield(handles,this_field)
        handles = rmfield(handles,this_field);
    end
end
handles.default_size = [256,256];
handles.gas = zeros(handles.default_size);
handles.proton = zeros(handles.default_size);
handles.proton_display = handles.proton;
handles.lungsmask = zeros(handles.default_size);
handles.lobemask = zeros(handles.default_size);
handles.nSlices = 1;
handles.nsl = 0;
set(handles.sliceSlider,'Value',0);
set(handles.sliceSlider,'Min',0);
set(handles.sliceSlider,'Max',1);
set(handles.is_poor_SNR,'Value',0);
% set(handles.is_coronal,'Value',1);
% set(handles.is_axial,'Value',0);
set(handles.is_coronal,'Value',0);
set(handles.is_axial,'Value',1);
handles.current_view = 'axial';
handles.acquisition = 'Coronal';
handles.RGthresh = 0.12;
set(handles.RGthreshold,'String',num2str(handles.RGthresh));
handles.enable_ute = 0;
set(handles.is_ute,'Value',0);
set(handles.is_ute,'Value',0);
% set(handles.is_ute,'Value',1);
handles.enable_recist = 0;
set(handles.is_recist,'Value',0);
handles.enable_OE = 0;
set(handles.is_OE,'Value',0);
handles.poor_SNR = false;



%% -----------------------------------------------
function output_fields = segmentation_outputs
output_fields = {'proton','gas','lungsmask','lobemask','nSlices','gasmask','trachea',...
    'gasPath','protonPath','protonFileString','gasFileString','gasFile','protonFile',...
    'protonFileString','bodymask','maskfile_name','RGseeds','lungsmask_unsmoothed',...
    'RGthresh','poor_SNR','voxelSize','lungVolume','enable_recist','enable_ute',...
    'enable_OE','aortamask','LVmask','Musclemask','current_view','rmin','rmax'};



%% --- Executes on button press in is_poor_SNR.
function is_poor_SNR_Callback(hObject, eventdata, handles)

% hObject    handle to is_poor_SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_poor_SNR
handles.poor_SNR = get(handles.is_poor_SNR,'Value');
if handles.poor_SNR
    handles.RGthresh = 0.03;
    set(handles.RGthreshold,'String',num2str(handles.RGthresh));
else
    handles.RGthresh = 0.08;
    set(handles.RGthreshold,'String',num2str(handles.RGthresh));
    
end
if sum(handles.gas(:))>0
    handles.gasmask = estimate_gasmask(handles.gas,handles.poor_SNR);
end
guidata(hObject,handles);


%% --- Executes on button press in is_coronal.
function is_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to is_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_coronal
is_coronal = get(handles.is_coronal,'Value');
if is_coronal
    handles.acquisition = 'Coronal';
    set(handles.is_axial,'Value',0);
else
    handles.acquisition = 'Axial';
    set(handles.is_axial,'Value',1);
end
guidata(hObject,handles);

%% --- Executes on button press in is_axial.
function is_axial_Callback(hObject, eventdata, handles)
% hObject    handle to is_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_axial
is_axial = get(handles.is_axial,'Value');
if is_axial
    handles.acquisition = 'Axial';
    set(handles.is_coronal,'Value',0);
else
    handles.acquisition = 'Coronal';
    set(handles.is_coronal,'Value',1);
end
guidata(hObject,handles);


%% -----------------------------------------
function RGthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to RGthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RGthreshold as text
%        str2double(get(hObject,'String')) returns contents of RGthreshold as a double
handles.RGthresh = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RGthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RGthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_lung_volume.
function display_lung_volume_Callback(hObject, eventdata, handles)
% hObject    handle to display_lung_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(handles.protonFileString)
is_recist = get(handles.is_recist,'Value');
is_ute = get(handles.is_ute,'Value');
is_paired_OE = get(handles.is_OE,'Value');
if is_recist
    [L,num] = bwlabeln(handles.lungsmask);
    for n=1:num
        this_lesion = zeros(size(L));
        this_lesion(L==n)=1;
        
        lungvolume = calculate_lung_volume(handles,this_lesion);
        intn_range=[handles.rmin,handles.rmax];
        recist = get_RECIST(handles.proton,this_lesion,handles.voxelSize,intn_range);
        if isfield(handles,'voxelSize') && ~isempty(handles.voxelSize)
            fprintf('Lesion@slice#%d volume =%6.2f ml.\n',recist.slice,lungvolume*1000);
            fprintf('RECIST =%6.2f mm.\n',recist.recist);
            fprintf('WHO =%6.2f mm.\n',recist.who);
            
            
        end
    end
elseif is_ute && ~is_paired_OE
    ROIs = {'LVmask','aortamask','Musclemask'};
    for nr = 1:numel(ROIs)
        this_ROI = ROIs{nr};
        if isfield(handles,this_ROI)
            
            pse_values = handles.gas(handles.(this_ROI)==1);
            fprintf('The median PSE of %s is %6.2f.\n',strrep(this_ROI,'mask',''),median(pse_values));
        end
    end
else
    handles.lungVolume = calculate_lung_volume(handles,handles.lungsmask);
    % recist = get_RECIST(handles.proton,handles.lungsmask,handles.voxelSize);
    if isfield(handles,'voxelSize') && ~isempty(handles.voxelSize)
        fprintf('Lung volume =%6.2f L.\n',handles.lungVolume);
        %    fprintf('Lung volume =%6.2f ml.\n',handles.lungVolume*1000);
        %    fprintf('RECIST =%6.2f mm.\n',recist.recist);
        %    fprintf('WHO =%6.2f mm.\n',recist.who);
        
    else
        fprintf('Lung volume = %i voxels.\n',handles.lungVolume);
        
    end
end
guidata(hObject,handles);

%% -------------------------------------------------------
% calculate lung volume
function lungVolume = calculate_lung_volume(handles,lungsmask)
lungVolume = sum(lungsmask(:));
if isfield(handles,'voxelSize') && ~isempty(handles.voxelSize)
    lungVolume = lungVolume*prod(handles.voxelSize)*1e-6;
end



% --- Executes on button press in is_ute.
function is_ute_Callback(hObject, eventdata, handles)
% hObject    handle to is_ute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_ute
is_ute = get(hObject,'Value');
if is_ute
    set(handles.is_recist,'Value',0);
    set(handles.is_OE,'Value',0);
    guidata(hObject,handles);
    
end

%% --- Executes on button press in rendering3D.
function rendering3D_Callback(hObject, eventdata, handles)
% hObject    handle to rendering3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lungsmask = handles.lungsmask;
if sum(lungsmask(:))>0
    is_rendering = get(hObject,'Value');
    if is_rendering ==1
        figure('Name','lung mask 3D rendering','NumberTitle','off');
        rendering3D(lungsmask);
    end
%     guidata(hObject,handles);
end
%% aorta

other_masks = {'aortamask','LVmask','Musclemask'};
if isfield(handles,'aortamask')|| isfield(handles,'LVmask')
color_list = {'y','r','c'};
figure('Name','volume rendering','NumberTitle','off');
hold on;
for n=1:numel(other_masks)
    maskname = other_masks{n};
    
    if isfield(handles,maskname) && sum(handles.(maskname)(:)>0)
        
        patch_plot(handles.(maskname),color_list{n});
    end
    
end
legend('Aorta','LV');
hold off;
set(handles.rendering3D,'Value',0);
end
% guidata(hObject,handles);



%%--- Executes on button press in data_overlay.
function data_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to data_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
is_ute = get(handles.is_ute,'Value');
if is_ute || handles.nSlices>100
    slices = round(linspace(1,handles.nSlices,35));
else
    slices = 1:handles.nSlices;
end
plotting_range = get_plotting_range(handles.lungsmask);

data_struct.proton = handles.proton;
gas_size = size(handles.gas);
if ~isempty(gas_size==1)
    data_struct.he = zeros(size(handles.proton));
else
    data_struct.he = handles.gas;
end
data_struct.lungsmask = handles.lungsmask;
data_struct.lobemask = handles.lungsmask;


button = questdlg('data overlay?','Select data','Mask on Proton','Mask on gas','Gas on proton','Mask on gas');

switch button
    case 'Mask on Proton'
        figure('Name',button,'NumberTitle','off');
        plot_defect_mask(data_struct.proton,data_struct.lungsmask,data_struct,plotting_range,'',slices);
        
    case 'Mask on gas'
        figure('Name',button,'NumberTitle','off');
        
        plot_defect_mask(data_struct.he,data_struct.lungsmask,data_struct,plotting_range,'',slices);
        
        
    case 'Gas on proton'
        figure('Name',button,'NumberTitle','off');
        
        proton_gas_overlay(handles.proton,handles.gas,'',slices);
        
end
set(handles.data_overlay,'Value',0);
guidata(hObject,handles);





% --- Executes on button press in is_recist.
function is_recist_Callback(hObject, eventdata, handles)
% hObject    handle to is_recist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_recist
is_recist = get(hObject,'Value');
if is_recist
    set(handles.is_ute,'Value',0);
    set(handles.is_OE,'Value',0);
    guidata(hObject,handles);
    
    
end


% --- Executes on button press in seeded_region_growing.
function seeded_region_growing_Callback(hObject, eventdata, handles)
% hObject    handle to seeded_region_growing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mid_slice = round(size(handles.proton,3)/2);
set(handles.sliceSlider,'Value',mid_slice);
update_figures(hObject, handles)

button=questdlg('Seed selection','Select a seed','Right Lung','Cancel','Right Lung');
if strcmp(button,'Right Lung')
    [X,Y]=ginput_b(1);
%     handles.RGseeds.right=[Y, X,mid_slice];
        handles.RGseeds.right=[X, Y,mid_slice];

else
    return;
end

button=questdlg('Seed selection','Select a seed','Left Lung','Cancel','Left Lung');
if strcmp(button,'Left Lung')
    [X,Y]=ginput_b(1);
%     handles.RGseeds.left=[Y,X, mid_slice];
    handles.RGseeds.left=[X, Y,mid_slice];

else
    return;
end

seg_proton_Callback(hObject, eventdata, handles);


% --- Executes on button press in convert_to_axial.
function convert_to_axial_Callback(hObject, eventdata, handles)
% hObject    handle to convert_to_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sag_to_ax = [3,2,1];
cor_to_ax = [2,3,1];
backup_fields = {'lungsmask','lobemask','proton','gas','aortamask','LVmask','Musclemask'};
for nf = 1:numel(backup_fields)
    this_field = backup_fields{nf};
    %     ori_field = strrep(this_field,'_ax','');
    if isfield(handles,this_field)
        %            handles.(ori_field) = handles.(this_field);
        %     end
        switch handles.current_view
            case 'coronal'
                handles.(this_field) = fliplr(imrotate(permute(handles.(this_field),cor_to_ax),270));
            case 'sagittal'
                handles.(this_field) = fliplr(imrotate(permute(handles.(this_field),sag_to_ax),270));
        end
    end
end

mid_slice = round(size(handles.proton,3)/2);
set(handles.sliceSlider,'Value',mid_slice);
handles.nSlices = size(handles.proton,3);
set(handles.sliceSlider,'Max',handles.nSlices);

handles.current_view = 'axial';
update_figures(hObject, handles)

% --- Executes on button press in convert_to_coronal.
function convert_to_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to convert_to_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax_to_cor = [3,2,1];
sag_to_cor = [1,3,2];
backup_fields = {'lungsmask','lobemask','proton','gas','aortamask','LVmask','Musclemask'};
for nf = 1:numel(backup_fields)
    this_field = backup_fields{nf};
    if isfield(handles,this_field)
        %            handles.(this_field) = handles.(ori_field);
        %     end
        switch handles.current_view
            case 'axial'
                handles.(this_field) = permute(handles.(this_field),ax_to_cor);
            case 'sagittal'
                handles.(this_field) = permute(handles.(this_field),sag_to_cor);
        end
    end
end



% handles.lungsmask = flip(permute(handles.lungsmask_ax, ax_to_cor),1);
% handles.lobemask = flip(permute(handles.lobemask_ax, ax_to_cor),1);
% handles.proton = permute(handles.proton_ax, ax_to_cor);
% handles.gas = permute(handles.gas_ax, ax_to_cor);
mid_slice = round(size(handles.proton,3)/2);
set(handles.sliceSlider,'Value',mid_slice);
handles.nSlices = size(handles.proton,3);
set(handles.sliceSlider,'Max',handles.nSlices);

handles.current_view = 'coronal';
update_figures(hObject, handles)


% --- Executes on button press in convert_to_sagittal.
function convert_to_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to convert_to_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax_to_sag = [3,1,2];
cor_to_sag = [1,3,2];
backup_fields = {'lungsmask','lobemask','proton','gas','aortamask','LVmask','Musclemask'};
for nf = 1:numel(backup_fields)
    this_field = backup_fields{nf};
    if isfield(handles,this_field)
        %            handles.(this_field) = handles.(ori_field);
        %     end
        switch handles.current_view
            case 'axial'
                handles.(this_field) = permute(handles.(this_field),ax_to_sag);
            case 'coronal'
                handles.(this_field) = permute(handles.(this_field),cor_to_sag);
        end
    end
    
end

% handles.lungsmask = permute(handles.lungsmask_ax, ax_to_sag);
% handles.lobemask = permute(handles.lobemask_ax, ax_to_sag);
% handles.proton = permute(handles.proton_ax, ax_to_sag);
% handles.gas = permute(handles.gas_ax, ax_to_sag);
mid_slice = round(size(handles.proton,3)/2);
set(handles.sliceSlider,'Value',mid_slice);
handles.nSlices = size(handles.proton,3);
set(handles.sliceSlider,'Max',handles.nSlices);

handles.current_view = 'sagittal';
update_figures(hObject, handles)


% --- Executes on button press in get_aorta.
function get_aorta_Callback(hObject, eventdata, handles)
% hObject    handle to get_aorta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'aortamask')
    handles.aortamask = zeros(size(handles.proton));
end
if sum(handles.proton(:)) ~= 0
    figure(1),
    imshow(handles.proton(:,:,handles.nsl),[handles.rmin,handles.rmax])
    bw = roipoly;
    
    if isempty(bw)
        return;
    end
    handles.aortamask(:,:,handles.nsl) = bw;
    
    update_figures(hObject, handles)
end
% --- Executes on button press in get_LV.
function get_LV_Callback(hObject, eventdata, handles)
% hObject    handle to get_LV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'LVmask')
    handles.LVmask = zeros(size(handles.proton));
end
if sum(handles.proton(:)) ~= 0
    figure(1),
    imshow(handles.proton(:,:,handles.nsl),[handles.rmin,handles.rmax])
    bw = roipoly;
    if isempty(bw)
        return;
    end
    
    handles.LVmask(:,:,handles.nsl) = bw;
    
    update_figures(hObject, handles)
end

% --- Executes on button press in propagate_contours.
function propagate_contours_Callback(hObject, eventdata, handles)
% hObject    handle to propagate_contours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button=questdlg('Propagate Contours','Select mask type','aorta','LV','Cancel','Cancel');
if strcmpi(button,'Cancel')
    return;
else
    maskname = sprintf('%smask',button);
end

if isfield(handles,maskname) && sum(handles.(maskname)(:))>0
    
    [hObject,handles] = propagate_contours(hObject,handles,maskname);
    update_figures(hObject, handles)
    
end


% --- Executes on button press in delete_slice_mask.
function delete_slice_mask_Callback(hObject, eventdata, handles)
% hObject    handle to delete_slice_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
this_slice = handles.nsl;
mask_list = {'lungsmask','aortamask','LVmask','Musclemask'};
for n_t = 1:numel(mask_list)
    maskname = mask_list{n_t};
    if isfield(handles,maskname)
        slice_mask = handles.(maskname)(:,:,this_slice);
        if sum(slice_mask(:))>0
            qst_string = sprintf('Delete %s at slice#%d',strrep(maskname,'mask',' mask'),this_slice);
            
            button=questdlg(qst_string,'Delete mask','Yes','No','Yes');
            if strcmp(button,'Yes')
                handles.(maskname)(:,:,this_slice) = zeros(size(slice_mask));
                update_figures(hObject, handles)
                
            end
            
        end
    end
end



% --- Executes on button press in is_OE.
function is_OE_Callback(hObject, eventdata, handles)
% hObject    handle to is_OE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_OE
is_OE = get(hObject,'Value');
if is_OE
    set(handles.is_OE,'Value',1);
    set(handles.is_recist,'Value',0);
    set(handles.is_ute,'Value',0);
    set(handles.is_coronal,'Value',0);
    guidata(hObject,handles);
    
end


% --- Executes on button press in dilate_lungsmask.
function dilate_lungsmask_Callback(hObject, eventdata, handles)
% hObject    handle to dilate_lungsmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'lungsmask') && sum(handles.lungsmask(:))>0
    nSlices = size(handles.lungsmask,3);
    se1=strel('disk',1);
    handles.lungsmask_backup = handles.lungsmask;
    
    for nsl = 1:nSlices
        handles.lungsmask(:,:,nsl) = imdilate(handles.lungsmask(:,:,nsl),se1);
    end
       update_figures(hObject, handles)
        
end


% --- Executes on button press in copy_2Dmask.
function copy_2Dmask_Callback(hObject, eventdata, handles)
% hObject    handle to copy_2Dmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mask_list = {'lungsmask','Musclemask','aortamask','LVmask'};
for nm=1:numel(mask_list)
    mask_name = mask_list{nm};
    if ~isfield(handles,mask_name)
        continue;
    end
    if nm==1
            handles.lungsmask_backup = handles.lungsmask;

    end
   nsl = round(get(handles.sliceSlider,'Value'));
   prev = handles.(mask_name)(:,:,nsl-1);
   post = handles.(mask_name)(:,:,nsl+1);
   mid_sl = round(size(handles.(mask_name),3)/2);
   if sum(prev(:))>0 && sum(post(:))==0
       handles.(mask_name)(:,:,nsl) = handles.(mask_name)(:,:,nsl-1);
   elseif sum(prev(:))==0 && sum(post(:))>0
       handles.(mask_name)(:,:,nsl) = handles.(mask_name)(:,:,nsl+1);
   elseif sum(prev(:))>0 && sum(post(:))>0
        if abs(mid_sl-(nsl-1))>abs(mid_sl-(nsl+1)) 
            handles.(mask_name)(:,:,nsl) = handles.(mask_name)(:,:,nsl+1);
        else
            handles.(mask_name)(:,:,nsl) = handles.(mask_name)(:,:,nsl-1);
        end
   end
end
    update_figures(hObject, handles)


% --- Executes on button press in add_muscle.
function add_muscle_Callback(hObject, eventdata, handles)
% hObject    handle to add_muscle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'Musclemask')
    handles.Musclemask = zeros(size(handles.proton));
end

if sum(handles.proton(:)) ~= 0
    figure(1),
    imshow(handles.proton(:,:,handles.nsl),[handles.rmin,handles.rmax])
    bw = roipoly;
    if isempty(bw)
        return;
    end
    
    handles.Musclemask(:,:,handles.nsl) =  handles.Musclemask(:,:,handles.nsl) +bw;
    
    update_figures(hObject, handles)
end
