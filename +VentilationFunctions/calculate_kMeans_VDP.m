function [Ventilation] = calculate_kMeans_VDP(Ventilation,Proton,MainInput)

%% MATLAB script to perform K-means VDP calculations
% Inputs:
% MR = ventilation image data
% maskarray = corresponding masks for MR
%
% Following Kirby et al., Academic Radiology, 2012.
%
% Author: Riaz Hussain, PhD
img = Ventilation.UncorrectedImage;
img = (img - min(img(:)))./(max(img(:) - min(img(:))));
maskarray = Ventilation.LungMask;
img = img.*maskarray;
%% %% Load and re-orient data
% % Select nifti image
% [img_name,path] = uigetfile({'*.nii.gz'},'Select Ventillation Nifti file');
% cd(path)
% MR = niftiread(img_name);
% %%
% MR_norm = MR / max(MR,[], 'all');
% I90 = imrotate(MR_norm, 90);
% MR2 = flip(I90, 1);
% montage(MR2);
% % Load mask
% maskarr = niftiread("img_ventilation_mask.nii.gz");
% I90_msk = imrotate(maskarr, 90);
% maskarray = flip(I90_msk, 1);

%% Run 2 rounds of kmeans following Kirby et al., 2012.
% Number of clusters
numClusters = 4;
% Initialize allMask
allMask = ones(size(img));
% Kmeans Round 1
imageVector1 = img(allMask == 1);

try
    clusterCenters1 = linspace(0.2, 0.8, 4)' * (max(imageVector1) - min(imageVector1)) + min(imageVector1);
%     clusterCenters1 = [0.4, 0.55, 0.7, 0.85]' * (max(imageVector1) - min(imageVector1)) + min(imageVector1);
catch
    error('Kmeans: Round 1 failed.');
end
kmeansRound1 = kmeans(imageVector1, numClusters, 'Start', clusterCenters1, ...
                      'MaxIter',1000);
% Kmeans Round 2
imageVector2 = imageVector1(kmeansRound1 == 1);
clusterCenters2 = linspace(0.2, 0.8, 4)' * (max(imageVector2) - min(imageVector2)) + min(imageVector2);
% clusterCenters2 = [0.4, 0.55, 0.7, 0.85]' * (max(imageVector2) - min(imageVector2)) + min(imageVector2);
try
    kmeansRound2 = kmeans(imageVector2, numClusters, 'Start', clusterCenters2, ...
                          'MaxIter',1000);
catch
    error('Kmeans: Round 2 failed.');
end

%% Adjust clusters and make total 5 clusters
kmeansRound2(kmeansRound2 == 2) = 1;
kmeansRound2(kmeansRound2 == 3 | kmeansRound2 == 4) = 2;
kmeansRound1 = kmeansRound1 + 1;
kmeansRound1(kmeansRound1 == 2) = kmeansRound2;
% Initialize arrays to store min and max values for each cluster
minValues = zeros(1, numClusters+1);
maxValues = zeros(1, numClusters+1);
% Calculate min and max values for each cluster
for i = 1:numClusters+1
    clusterIndices = (kmeansRound1 == i);
    minValues(i) = min(imageVector1(clusterIndices));
    maxValues(i) = max(imageVector1(clusterIndices));
end
% Display the results
for i = 1:numClusters+1
    fprintf('Cluster %d: Min = %f, Max = %f\n', i, minValues(i), maxValues(i));
end
% %Create segmentation
segmentation = zeros(size(img));
segmentation(allMask == 1) = kmeansRound1;
segmentation = segmentation .* double(maskarray);
tot_vox=sum(segmentation > 0, 'all');
cluster5=(sum(segmentation == 5, 'all')/tot_vox)*100;
cluster4=(sum(segmentation == 4, 'all')/tot_vox)*100;
cluster3=(sum(segmentation == 3, 'all')/tot_vox)*100;
cluster2_hypo = (sum(segmentation == 2, 'all')/tot_vox)*100;
cluster1_vdp = (sum(segmentation == 1, 'all')/tot_vox)*100;
VDP = cluster1_vdp + cluster2_hypo;
Ventilation.Kmeans_segmentation = segmentation;
Ventilation.KmeansVDP = VDP;
fprintf("VDP = %f", VDP);
fprintf("Cluster 3 = %f, Cluster 4 = %f, Cluster 5 = %f\n", cluster3, cluster4, cluster5);

%% %Check for non-zero slices along the third dimension

% check images size, 2D or 3D
if size(Ventilation.Image,3) > 30
    Image_3D = 1;
    nrow = ceil(sqrt(size(Ventilation.Image,3)));
    % Initialize an empty vector to store indices of slices with no zeros
    slices_with_no_zeros = zeros(1,size(Ventilation.LungMask,3));    
    % Loop through each slice
    for slice_idx = 1:size(Ventilation.LungMask,3)
        % Extract the current slice
        current_slice = Ventilation.LungMask(:, :, slice_idx);
        
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
cmap = [1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0];
% % %plot and see the colormap
outputpath = [Ventilation.outputpath, '\VDP Analysis'];
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
saveas(gca,'kmeans_VDPmap.png');
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
VDPmap = segmentation;
VDPmap(VDPmap>0 & VDPmap<3) = 2;
VDPmap(VDPmap>2) = 1;
% imslice(VDPmap); colormap(cmap)

cmap = [0 0 0;  0 1 0; 1 0 0];
for slice = 1:size(VDPmap,3)
    [~,~] = Global.imoverlay(squeeze(abs(Proton.ProtonRegistered(:,:,slice))),squeeze(VDPmap(:,:,slice).*2),[1,4],[0,0.8*max(Proton.ProtonRegistered(:))],cmap,1,gca);
    colormap(gca,cmap) 
%     imshow(VDPmap(:,:,slice).*2, cmap, 'DisplayRange', [0 4]);
    Xdata = getframe(gcf);
    X = Xdata.cdata; 
    if (slice == 1)
        imwrite(X,[outputpath,'\kmeansDefectmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\kmeansDefectmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving kmeans Defectmap Tiff Completed.')
close all;

% read tiff
tiff_info = imfinfo('kmeansDefectmap.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
kmeansDefectmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('kmeansDefectmap.tif', ii);
    kmeansDefectmap(:,:,:,ii) = temp_tiff;
end
Ventilation.kmeansDefectmap = kmeansDefectmap;
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
ventimage = Ventilation.UncorrectedImage(:,:,slice_indices);
ventimage = ventimage/max(ventimage(:));
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(reshape(ventimage,[size(ventimage,1), size(ventimage,2), 1, size(ventimage,3)]),...
    'Size',[1 size(ventimage,3)],'DisplayRange',[0 0.8]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
VentMontagePosition=get(gcf,'position'); 

% re_seg = (re_seg + 1).*msk2d;
cmap = [0 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0];
re_seg(isnan(re_seg)) = 0;
KmeansVDPMontage = figure('Name','Vent Image');set(KmeansVDPMontage,'WindowState','minimized');
montage(re_seg,'DisplayRange',[]); colormap(cmap);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
KmeansVDPMontagePosition=get(gcf,'position'); 

parentPath = Ventilation.outputpath;
cd(parentPath)
ReportTitle = 'Ventilation_Analysis';
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
ppt_file_name = 'Ventilation_Analysis.pptx';
if isfile(ppt_file_name)
    disp('file existed')
    Global.exportToPPTX('open',ppt_file_name);
    Global.exportToPPTX('switchslide',1);
else            
    Global.exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
end 
%Add slides
Global.exportToPPTX('addslide'); % angles
Global.exportToPPTX('addtext','K-means VDP','Position',[6 0 7 1],'Color','b','FontSize',25);

Global.exportToPPTX('addpicture',VentMontage,'Position',...
    [0 0.5 NumSliceView*scalefactor scalefactor*NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);

Global.exportToPPTX('addpicture',KmeansVDPMontage,'Position',...
    [0 2 NumSliceView*scalefactor scalefactor*NumSliceView*(KmeansVDPMontagePosition(4)/KmeansVDPMontagePosition(3))]);
Global.exportToPPTX('addtext',['K-mean VDP = ',num2str(round(VDP,2)),'%'],'Position',[0 2 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');

Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');
            Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
            Global.exportToPPTX('close');
            fprintf('PowerPoint file has been saved\n');               

close all;
