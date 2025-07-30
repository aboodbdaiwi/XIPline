function [Ventilation] = VDP_calculationASB(Ventilation,Proton,MainInput)
% this code is based on the ventilation analysis using the adaptive k-means method from Iowa site
%
%
%  added by: Abdullah S. Bdaiwi 2024/14/11
maskarray = double(Ventilation.LungMask);
maskarray(Ventilation.AirwayMask == 1) = 0;
maskarray(Ventilation.VesselMask == 1) = 0;

% defect analysis
temp_handles.proton_used = 'SSFSE'; % hard code for now
temp_handles.proton = Proton.Image;

% Apply dilation
if size(maskarray, 1) <= 128
    se = strel('disk', 2);
elseif size(maskarray, 1) <= 192
    se = strel('disk', 4);
else
    se = strel('disk', 6);    
end
dilated_mask = zeros(size(maskarray));
for i =1:size(maskarray,3)
    dilated_mask(:,:,i) = double(imdilate(maskarray(:,:,i), se) > 0);
end

temp_handles.he = Ventilation.Image; %xe.img;
temp_handles.lungsmask = dilated_mask; % maskarray
temp_handles.vessel_mask = Ventilation.VesselMask;
temp_handles.lobemask = maskarray;
temp_handles.ventmaskbinary = maskarray;
%         temp_handles.kmeans_threshold=1;
%         temp_handles.num_clusters=4;
%         temp_handles.low_intn_perc_cutoff=0.1;
temp_handles.hist_mu = 385;
temp_handles.frangi_thresh = 0.25;%0.15;
temp_handles.frangi_sigmas = 1;%1.5;%
temp_handles.pid = 'XeCTC';
temp_handles.vessel_removal = false;
%         temp_handles.proton_trasform=true;
temp_handles.vessel_thresh = 80;
temp_handles.view = MainInput.SliceOrientation;
temp_handles.helium_voxelsize = prod(MainInput.PixelSpacing)*MainInput.SliceThickness*1e-6;
[defect_mask,~,~,vessel_mask,transformed_proton,low_intn_perc,ventilation_mask] = VentilationFunctions.Akmeans.ventilation_defect_by_kmeans(temp_handles);
% imslice(ventilation_mask)
ventilation_mask = ventilation_mask.*maskarray;
defect_mask = defect_mask.*maskarray;
volume_struct= VentilationFunctions.Akmeans.calculate_percent_defect(defect_mask,temp_handles);

ventP_raw_cor = VentilationFunctions.Akmeans.classify_ventilation_levles(ventilation_mask,temp_handles.lungsmask);
Ventilation.Akmeans_VDP = ventP_raw_cor(1);
Ventilation.Akmeans_LVP = ventP_raw_cor(2);
Ventilation.Akmeans_MVP = ventP_raw_cor(3);
Ventilation.Akmeans_HVP = ventP_raw_cor(3);
Ventilation.Akmeans_ventP_raw_cor = ventP_raw_cor;

fprintf("Adaptive k-means VDP = %f", Ventilation.Akmeans_VDP);
segmentation = ventilation_mask;
Ventilation.Aksegmentation = ventilation_mask;
Ventilation.Akmeans_defect_mask = defect_mask;
%% %Check for non-zero slices along the third dimension

% check images size, 2D or 3D
if size(Ventilation.Image,3) > 40
    Image_3D = 1;
    nrow = ceil(sqrt(size(Ventilation.Image,3)));
    % Initialize an empty vector to store indices of slices with no zeros
    slices_with_no_zeros = zeros(1,size(maskarray,3));    
    % Loop through each slice
    for slice_idx = 1:size(maskarray,3)
        % Extract the current slice
        current_slice = maskarray(:, :, slice_idx);
        
        % Check if all elements in the current slice are non-zero
        if sum(current_slice(:)) > 0
            % Append the slice index to the vector if all elements are non-zero
            slices_with_no_zeros(slice_idx) = 1;
        end
    end
    nonZeroSlices = find(slices_with_no_zeros);
    nonZeroSlice_lng = nonZeroSlices(end) - nonZeroSlices(1);
    nonZeroSlice_space = floor(nonZeroSlice_lng/16);
    slice_indices = nonZeroSlices(1):nonZeroSlice_space:nonZeroSlices(end);
    nonzero_seg = segmentation(:, :, slice_indices);
    nonzero_msk = maskarray(:, :, slice_indices);
    re_seg = reshape(nonzero_seg, size(nonzero_seg, 1), []);
    msk2d = reshape(nonzero_msk,size(nonzero_msk, 1), []);
else
    Image_3D = 0;
    nonzero_slices = any(any(maskarray, 1), 2);
    % %Extract non-zero slices and reshpe to 2D
    nonzero_seg = segmentation(:, :, nonzero_slices);
    nonzero_msk = maskarray(:, :, nonzero_slices);
    re_seg = reshape(nonzero_seg, size(nonzero_seg, 1), []);
    msk2d = reshape(nonzero_msk,size(nonzero_msk, 1), []);
end 

%% plot the clustered colormap
re_seg(msk2d == 0) = NaN;
% re_seg = (re_seg + 1).*msk2d;
% cmap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 0 0 1];
cmap = [1 0 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0];
% % %plot and see the colormap
outputpath = [Ventilation.parentPath, '\VDP_Analysis'];
mkdir(outputpath);
cd(outputpath);

figure;
montage_kmeans = imshow(re_seg, cmap, 'DisplayRange', []);
set(montage_kmeans, 'AlphaData', ~isnan(re_seg))
axis on;
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', 'k')
colorbar(gca, 'Ticks', [1.5,2.5,3.5,4.5,5.5], 'TickLabels', [1, 2,3,4,5], ...
    'FontSize', 10, 'FontWeight', 'bold');
% Save the colormap image
set(gcf,'PaperPositionMode','auto')
% print('kmeans_VDPmap','-dpng','-r300')
saveas(gca,'Adaptive_kmeans_VDPmap.png');
close all;
%% write tiff and read back VDPVent maps
clc
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving VentDefectmap Tiff...')


% Create an output image initialized to zeros (same size as xenon_image)
VDPmap = zeros(size(segmentation));
class_1_inside_lung = (segmentation == 1) & (maskarray == 1);
class_2to4_inside_lung = (segmentation >= 2) & (maskarray == 1);
VDPmap(class_1_inside_lung) = 2;
VDPmap(class_2to4_inside_lung) = 1;

% imslice(VDPmap); colormap(cmap)

cmap = [0 0 0;  0 1 0; 1 0 0];
for slice = 1:size(VDPmap,3)
    [~,~] = Global.imoverlay(squeeze(abs(Proton.ProtonRegistered(:,:,slice))),squeeze(VDPmap(:,:,slice).*2),[1,4],[0,0.8*max(Proton.ProtonRegistered(:))],cmap,1,gca);
    colormap(gca,cmap) 
%     imshow(VDPmap(:,:,slice).*2, cmap, 'DisplayRange', [0 4]);
    Xdata = getframe(gcf);
    X = Xdata.cdata; 
    if (slice == 1)
        imwrite(X,[outputpath,'\AkmeansDefectmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\AkmeansDefectmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving kmeans Defectmap Tiff Completed.')
close all;
cd(outputpath)
% read tiff
tiff_info = imfinfo('AkmeansDefectmap.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
AkmeansDefectmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('AkmeansDefectmap.tif', ii);
    AkmeansDefectmap(:,:,:,ii) = temp_tiff;
end
Ventilation.AkmeansDefectmap = AkmeansDefectmap;

%% write tiff and read back VDPVent maps
clc
tiff2 = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff2);ax2 = axes('Parent',tiff2);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving VentDefectmap2 Tiff...')

VDPmap= defect_mask;

% imslice(VDPmap)
cus_colormap = zeros(100,3);
% colorgradient = 0:0.01:0.99;
cus_colormap(:,1) = 1;
cus_colormap(:,2) = 0;
cus_colormap(:,3) = 0;

for slice = 1:size(VDPmap,3)
    [~,~] = Global.imoverlay(squeeze(abs(Ventilation.Image(:,:,slice))),squeeze(VDPmap(:,:,slice)),[1,100],[0,max(Ventilation.Image(:))],cus_colormap,1,gca);
    colormap(gca,cus_colormap) 
%     imshow(VDPmap(:,:,slice).*2, cmap, 'DisplayRange', [0 4]);
    Xdata = getframe(gcf);
    X = Xdata.cdata; 
    if (slice == 1)
        imwrite(X,[outputpath,'\AkmeansDefectmap2.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\AkmeansDefectmap2.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving Akmeans Defectmap 2 Tiff Completed.')
close all;
cd(outputpath)
% read tiff
tiff_info2 = imfinfo('AkmeansDefectmap2.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
AkmeansDefectmap2 = uint8(zeros(tiff_info2(1).Height ,tiff_info2(1).Width ,3,length(tiff_info2)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info2, 1)
    temp_tiff = imread('AkmeansDefectmap2.tif', ii);
    AkmeansDefectmap2(:,:,:,ii) = temp_tiff;
end
Ventilation.AkmeansDefectmap2 = AkmeansDefectmap2;
%% Output mask boundaries  with ventilation overlays:
% Create a folder for the ventilation images:
outputPath = outputpath;
cd(outputPath);
cus_colormap = zeros(100,3);

cus_colormap(:,1) = 1; % Red channel
cus_colormap(:,2) = 0; % Green channel
cus_colormap(:,3) = 1; % Blue channel

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')

% Normalize the intensity of the original image to fall between [0,1].
MR = Ventilation.Image;
MR2 = MR / max(MR,[], 'all');

for slice=1:size(maskarray,3) %repeat for rest of slices
    maskboundaries = bwboundaries(maskarray(:,:,slice));
    % Check if bwboundaries returned an empty 0x1 cell
    if isempty(maskboundaries) && isequal(size(maskboundaries), [0, 1])
        maskboundaries = zeros(size(MR2(:,:,slice)));
    else
        % Convert boundary points to a binary mask
        maskImg = false(size(maskarray(:,:,slice)));
        for k = 1:length(maskboundaries)
            boundary = maskboundaries{k};
            maskImg(sub2ind(size(maskImg), boundary(:,1), boundary(:,2))) = true;
        end
        % Thicken boundary using morphological dilation
        se = strel('disk', 1);  % Use radius=1 or larger for thicker lines
        thickMask = imdilate(maskImg, se);
        maskboundaries = double(thickMask);
        % maskboundaries = double(maskImg);
    end

    [~,~] = Global.imoverlay(squeeze(abs(MR2(:,:,slice))),squeeze(maskboundaries),[1,100],[0,max(MR2(:))],cus_colormap,1,gca);
    colormap(gca,cus_colormap)     
    Xdata = getframe(gcf);
    X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[outputPath,'\maskboundaries.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputPath,'\maskboundaries.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('maskboundariesTiff Completed.')
close all;
% read tiff
cd(outputPath)
tiff_info = imfinfo('maskboundaries.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
Mask_Vent_boundaries = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 1 : size(tiff_info, 1)
    temp_tiff = imread('maskboundaries.tif', ii);
    Mask_Vent_boundaries(:,:,:,ii) = temp_tiff;
end
Mask_Vent_boundaries = permute(Mask_Vent_boundaries,[1 2 4 3]);
Ventilation.Mask_Vent_boundaries = Mask_Vent_boundaries;
%% Output mask boundaries  with proton overlays:
% Create a folder for the ventilation images:

outputPath = outputpath;
cd(outputPath);
cus_colormap = zeros(100,3);

cus_colormap(:,1) = 1; % Red channel
cus_colormap(:,2) = 0; % Green channel
cus_colormap(:,3) = 1; % Blue channel

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')

% Normalize the intensity of the original image to fall between [0,1].
HImage = Proton.Image;
MR2 = HImage / max(HImage,[], 'all');

for slice=1:size(maskarray,3) %repeat for rest of slices
    maskboundaries = bwboundaries(Ventilation.LungMask(:,:,slice));
    % Check if bwboundaries returned an empty 0x1 cell
    if isempty(maskboundaries) && isequal(size(maskboundaries), [0, 1])
        maskboundaries = zeros(size(MR2(:,:,slice)));
    else
        % Convert boundary points to a binary mask
        maskImg = false(size(Ventilation.LungMask(:,:,slice)));
        for k = 1:length(maskboundaries)
            boundary = maskboundaries{k};
            maskImg(sub2ind(size(maskImg), boundary(:,1), boundary(:,2))) = true;
        end
        % Thicken boundary using morphological dilation
        se = strel('disk', 1);  % Use radius=1 or larger for thicker lines
        thickMask = imdilate(maskImg, se);
        maskboundaries = double(thickMask);
        % maskboundaries = double(maskImg);
    end

    [~,~] = Global.imoverlay(squeeze(abs(HImage(:,:,slice))),squeeze(maskboundaries),[1,100],[0,max(HImage(:))],cus_colormap,1,gca);
    colormap(gca,cus_colormap)     
    Xdata = getframe(gcf);
    X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[outputPath,'\protonmaskboundaries.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputPath,'\protonmaskboundaries.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('protonmaskboundaries Tiff Completed.')
close all;
% read tiff
cd(outputPath)
tiff_info = imfinfo('protonmaskboundaries.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
Mask_Proton_boundaries = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 1 : size(tiff_info, 1)
    temp_tiff = imread('protonmaskboundaries.tif', ii);
    Mask_Proton_boundaries(:,:,:,ii) = temp_tiff;
end
Mask_Proton_boundaries = permute(Mask_Proton_boundaries,[1 2 4 3]);
Ventilation.Mask_Proton_boundaries = Mask_Proton_boundaries;
%% save ppt 
% %how many slices to incluse in the final powerpoint file
sl_include=zeros(1,size(maskarray,3));
for sl=1:size(maskarray,3)
    choose_slice=maskarray(:,:,sl);       
    logica_test2=sum(choose_slice(:));
    if logica_test2>1
         sl_include(sl)=sl;
    else
         sl_include(sl)=0;
    end 
end
sl_include( :, ~any(sl_include,1) ) = [];  %columns
sl_1=sl_include(1);
sl_end=sl_include(end);
NumSliceView = length(sl_include);
if NumSliceView > 16
    NumSliceView = 16;
end
scalefactor = 16/NumSliceView;
if Image_3D ~= 1
    slice_indices = sl_1:1:sl_end;
end
ventimage = Ventilation.Image(:,:,slice_indices);
ventimage = ventimage/max(ventimage(:));
AVentMontage = figure('Name','Vent Image');set(AVentMontage,'WindowState','minimized');
montage(reshape(ventimage,[size(ventimage,1), size(ventimage,2), 1, size(ventimage,3)]),...
    'Size',[1 size(ventimage,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
AVentMontagePosition=get(gcf,'position'); 

% re_seg = (re_seg + 1).*msk2d;
KmeansVDPMontage = figure('Name','Vent Image'); set(KmeansVDPMontage,'WindowState','minimized');
% montage(re_seg,'DisplayRange',[]); colormap(cmap);
kmeansDefectslices = AkmeansDefectmap(:,:,:,slice_indices);
montage(kmeansDefectslices,'Size',[1 size(kmeansDefectslices,4)])
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
AKmeansVDPMontagePosition=get(gcf,'position'); 

AKmeansVDPMontage2 = figure('Name','Vent Image');set(AKmeansVDPMontage2,'WindowState','minimized');
kmeansDefectslices = AkmeansDefectmap2(:,:,:,slice_indices);
montage(kmeansDefectslices,'Size',[1 size(kmeansDefectslices,4)])
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
AKmeansVDPMontagePosition2=get(gcf,'position'); 

parentPath = Ventilation.parentPath;
cd(parentPath)
%Start new presentation
isOpen  = Global.exportToPPTX(); 
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
% Generate filename with today's date
today = datetime('today');
date_str = datestr(today, 'yyyymmdd');
ppt_file_name = ['Ventilation_Analysis_', date_str, '.pptx'];
ReportTitle = ['Ventilation_Analysis_', date_str];
% Create or open the presentation
if isfile(ppt_file_name)
    disp('File existed')
    Global.exportToPPTX('open', ppt_file_name);
    Global.exportToPPTX('switchslide', 1);
else            
    Global.exportToPPTX('new', 'Dimensions', [16 9], ...
        'Title', ReportTitle, ...
        'Author', 'CPIR @ CCHMC');
end
%Add slides
Global.exportToPPTX('addslide'); % angles
Global.exportToPPTX('addtext','AK-means VDP','Position',[6 0 7 1],'Color','b','FontSize',25);

Global.exportToPPTX('addpicture',AVentMontage,'Position',...
    [0 0.5 NumSliceView*scalefactor scalefactor*NumSliceView*(AVentMontagePosition(4)/AVentMontagePosition(3))]);

Global.exportToPPTX('addpicture',KmeansVDPMontage,'Position',...
    [0 2 NumSliceView*scalefactor scalefactor*NumSliceView*(AKmeansVDPMontagePosition(4)/AKmeansVDPMontagePosition(3))]);
Global.exportToPPTX('addtext',['AK-mean VDP = ',num2str(round(Ventilation.Akmeans_VDP,2)),'%'],'Position',[0 2 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');

Global.exportToPPTX('addpicture',AKmeansVDPMontage2,'Position',...
    [0 3.5 NumSliceView*scalefactor scalefactor*NumSliceView*(AKmeansVDPMontagePosition(4)/AKmeansVDPMontagePosition(3))]);
Global.exportToPPTX('addtext',['AK-mean VDP = ',num2str(round(Ventilation.Akmeans_VDP,2)),'%'],'Position',[0 3.5 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');

Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');
            Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
            Global.exportToPPTX('close');
            fprintf('PowerPoint file has been saved\n');               
close all;


%% Output mask images with ventilation overlays:
% Create a folder for the ventilation images:

cd(outputpath);
cus_colormap = zeros(100,3);
% colorgradient = 0:0.01:0.99;
cus_colormap(:,1) = 1; % Red channel
cus_colormap(:,2) = 0; % Green channel
cus_colormap(:,3) = 1; % Blue channel

% Vent_img = Ventilation.Image;
% Vent_img = Vent_img./max(Vent_img(:));
% Vent_image_colored = ones(size(Ventilation.Image,1),size(Ventilation.Image,2),size(Ventilation.Image,3), 3);
% ch1 = 1;
% ch2 = 0;
% ch3 = 1;
% Vent_image_colored(:,:,:,1) = Vent_img .* ch1;
% Vent_image_colored(:,:,:,2) = Vent_img .* ch2;
% Vent_image_colored(:,:,:,3) = Vent_img .* ch3;
% orthosliceViewer((Vent_image_colored));

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')

% Normalize the intensity of the original image to fall between [0,1].
MR = Ventilation.Image;
MR2 = MR / max(MR,[], 'all');
scaledImage = zeros(size(MR2));
scaledImage2 = zeros(size(MR2));
for m = 1:size(MR2,3)
    scaledImage(:,:,m) = MR2(:,:,m);
    scaledImage2(:,:,m) = scaledImage(:,:,m);
%     outputFileName = sprintf('Img_%d.png',m);
%     fullFileName = fullfile(outputPath, outputFileName);
%     imwrite(scaledImage2(:,:,m),fullFileName,'png');
end

% %Mask images
% MR4 = uint16(Ventilation.Image);
% 
% for m2 = 1:size(MR4,3)
%     scaledImage4(:,:,m2) = imadjust(MR4(:,:,m2));
%     scaledImage5(:,:,m2) = adapthisteq(scaledImage4(:,:,m2));
% end
% scaledImage5 = double(scaledImage5);
% scaledImage5 = scaledImage5 - min(scaledImage5(:));
% scaledImage5 = scaledImage5 ./ max(scaledImage5(:));

for slice=1:size(maskarray,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(scaledImage2(:,:,slice))),squeeze(maskarray(:,:,slice)),[1,100],[0,max(scaledImage2(:))],cus_colormap,0.4,gca);
    colormap(gca,cus_colormap)   
%     X = print('-RGBImage',['-r',num2str(size(VentBinMap2,2)/2)]);%2 inches
%      hImage = get( gca, 'Children' ); 
%      X = hImage.CData;   
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[outputpath,'\Mask_Vent_Reg.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\Mask_Vent_Reg.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('Mask_Vent_Reg Tiff Completed.')
close all;
% read tiff
cd(outputpath)
tiff_info = imfinfo('Mask_Vent_Reg.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
Mask_Vent_Reg = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 1 : size(tiff_info, 1)
    temp_tiff = imread('Mask_Vent_Reg.tif', ii);
    Mask_Vent_Reg(:,:,:,ii) = temp_tiff;
end
Mask_Vent_Reg = permute(Mask_Vent_Reg,[1 2 4 3]);
Ventilation.Mask_Vent_Reg = Mask_Vent_Reg;
%% compute VDP per Slice
defect_mask = Ventilation.Akmeans_defect_mask;
defect_mask(defect_mask > 0 & defect_mask < 2) = 1;
defect_mask(defect_mask > 1) = 0;
defectMap = double(maskarray) + defect_mask; 
[vdp_per_slice_local, vdp_per_slice_global] = VentilationFunctions.computeVDPperSlice(defectMap);
Ventilation.AKmeans_vdp_per_slice_local = vdp_per_slice_local;
Ventilation.AKmeans_vdp_per_slice_global = vdp_per_slice_global;
%% Histogram

% Extract voxel values within lung
XeVals = Ventilation.Image(maskarray > 0);
ClusterVals = Ventilation.Aksegmentation(maskarray > 0);

% Define bin edges and cluster colors
edges = linspace(0, prctile(XeVals,99), 200);  % or adjust based on your data range
binCenters = (edges(1:end-1) + edges(2:end)) / 2;

% Define colors (adjust to match your aesthetic)
clusterColors = [
    0.9 0.2 0.2;   % Cluster 1 (red)
    0.95 0.6 0.1;  % Cluster 2 (orange)
    0.4 0.8 0.4;   % Cluster 3 (green)
    0.2 0.2 0.8    % Cluster 5 (blue)
];

% Create figure
figure; hold on;

% Plot histogram for each cluster
for c = 1:4
    vals_c = XeVals(ClusterVals == c);
    hc = histogram(vals_c, edges, ...
        'FaceColor', clusterColors(c,:), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0);
end

% Plot overall density curve (optional)
XeVals = XeVals(XeVals > 0);  % remove zeros and negatives
[f, xi] = ksdensity(XeVals, 'Support', 'positive');
plot(xi, f * numel(XeVals) * (edges(2)-edges(1)), 'k--', 'LineWidth', 1.5);

% Labeling
xlabel('Ventilation Value', 'FontSize', 14);
ylabel('Voxel Count', 'FontSize', 14);
title('Adaptive K-means Clustered Histogram', 'FontSize', 14);
xlim([0 prctile(XeVals,99)]);
set(gca,'FontSize',14)
% Optional: Add vertical cluster boundaries
% (e.g., using means or centroids if available)
box on;
set(gcf,'PaperPosition',[0 0 6 3.5]);
cd(outputpath);
% print([foldername + 'Ventilation_Histogram'],'-dpng','-r300');
saveas(gca,'Adaptive_Kmeans_Histogram.png');
close all;

%%
save_data=[parentPath,'\','Ventilation_Analysis','.mat'];
save(save_data); 


end
