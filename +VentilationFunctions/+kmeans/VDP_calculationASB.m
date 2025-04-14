function [Ventilation] = VDP_calculationASB(Ventilation,Proton,MainInput)

% load the deformed lung/lobe mask. 
maskarray = Ventilation.LungMask;
% maskarray = double(Ventilation.LungMask + Ventilation.VesselMask);
% maskarray(maskarray > 1) = 0;
% maskarray = double(maskarray);
Ventilation.SliceOrientation = MainInput.SliceOrientation;

% k-means clustering using our previous code.
[Ventilation] = VentilationFunctions.kmeans.clustering_3DASB(Ventilation);
clusters_3D = Ventilation.clusters_3D;

% Xemaskd = Ventilation.Image.*Ventilation.LungMask;
% [heSlicesClusters,heSlicesMask,clusterMin] = VentilationFunctions.CallClustering(Xemaskd);
% clusters_3D = heSlicesClusters;
% Ventilation.clusters_3D = clusters_3D;
% imslice(clusters_3D)

% get the number of He image clusters. Usually, one integral corresponds to one cluster.
clusters_3D_labels = unique(clusters_3D);
clusters_3D_label_length = length (clusters_3D_labels);

% get the number of lobes. Usually, one lobel corresponds to one lobe.
lobemask_labels = unique(maskarray);
lobemask_labels = lobemask_labels(lobemask_labels>0);
lobemask_label_length = length (lobemask_labels);

% create a 2D array to save the regional VDP.
% row represents lobes in the order of left upper---left lower---right upper---right
% middle---right lower 
lobe_wise_VDP = zeros(lobemask_label_length,clusters_3D_label_length + 1);

% interate through each lobe% open and close file by file ID.
for i = 1 : lobemask_label_length
    
    % find the current lobe i.
    lobe = (maskarray==lobemask_labels(i));
    lobe_wise_VDP(i,1) = lobemask_labels(i);
    % interate through each He cluster
    for j = 1 : clusters_3D_label_length
        % find the current cluster.
        cluster = (clusters_3D==clusters_3D_labels(j));
        % find cluster j in the current lobe i. for example, we want to
        % see cluster 0 in lobe 1.=, similarly, we can see cluster 1 in lobe 1
        cluster_i_in_lobe_j = lobe.*cluster;
        % calculate the region of the current cluster, in the units of
        % voxel.
        lobe_wise_VDP(i,j+1) = sum(cluster_i_in_lobe_j(:));
    end
    % normalize the region of current cluster to the current lobe for
    % regional VDP calculation.
    lobe_wise_VDP(i,2:clusters_3D_label_length + 1) = lobe_wise_VDP(i,2:clusters_3D_label_length + 1)/sum(lobe(:))*100;
end
    
% whole lung VDP calculation, we actually calculate the ratio of each
% He image cluster over the whole lung.
wholelung_VDP = zeros (clusters_3D_label_length,1);
% get the whole lung volume.
total_lung_volume = (maskarray~=0);

for i = 1 : clusters_3D_label_length 
        cluster_within_lung = (clusters_3D==clusters_3D_labels(i)).*total_lung_volume;
        wholelung_VDP(i) = sum(cluster_within_lung(:))/sum(total_lung_volume(:))*100;
end

Ventilation.wholelung_VDP = wholelung_VDP;
Ventilation.KmeansVDP = wholelung_VDP(1);
fprintf("k-means VDP = %f", wholelung_VDP(1));
segmentation = clusters_3D;


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
outputpath = [Ventilation.parentPath, '\VDP Analysis'];
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
set(gcf,'PaperPositionMode','auto');  
% print('kmeans_VDPmap','-dpng','-r300')
saveas(gca,'kmeans_VDPmap.png');
close all;
%% write tiff and read back VDPVent maps
clc
% tiff = figure('Visible', 'off');  % Create the figure completely invisible
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized','Visible', 'off');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving VentDefectmap Tiff...')


% Create an output image initialized to zeros (same size as xenon_image)
VDPmap = zeros(size(segmentation));
class_0_inside_lung = (segmentation == 0) & (maskarray == 1);
class_1to4_inside_lung = (segmentation >= 1) & (maskarray == 1);
VDPmap(class_0_inside_lung) = 2;
VDPmap(class_1to4_inside_lung) = 1;
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

VDPmap(class_0_inside_lung) = 1;
VDPmap(class_1to4_inside_lung) = 0;
% imslice(VDPmap)
cus_colormap = zeros(100,3);
% colorgradient = 0:0.01:0.99;
cus_colormap(:,1) = 1;
cus_colormap(:,2) = 0;
cus_colormap(:,3) = 1;

for slice = 1:size(VDPmap,3)
    [~,~] = Global.imoverlay(squeeze(abs(Ventilation.Image(:,:,slice))),squeeze(VDPmap(:,:,slice)),[1,100],[0,max(Ventilation.Image(:))],cus_colormap,1,gca);
    colormap(gca,cus_colormap) 
%     imshow(VDPmap(:,:,slice).*2, cmap, 'DisplayRange', [0 4]);
    Xdata = getframe(gcf);
    X = Xdata.cdata; 
    if (slice == 1)
        imwrite(X,[outputpath,'\kmeansDefectmap2.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\kmeansDefectmap2.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving kmeans Defectmap 2 Tiff Completed.')
close all;
cd(outputpath)
% read tiff
tiff_info2 = imfinfo('kmeansDefectmap2.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
kmeansDefectmap2 = uint8(zeros(tiff_info2(1).Height ,tiff_info2(1).Width ,3,length(tiff_info2)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info2, 1)
    temp_tiff = imread('kmeansDefectmap2.tif', ii);
    kmeansDefectmap2(:,:,:,ii) = temp_tiff;
end
Ventilation.kmeansDefectmap2 = kmeansDefectmap2;
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
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(reshape(ventimage,[size(ventimage,1), size(ventimage,2), 1, size(ventimage,3)]),...
    'Size',[1 size(ventimage,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
VentMontagePosition=get(gcf,'position'); 

% re_seg = (re_seg + 1).*msk2d;
KmeansVDPMontage = figure('Name','Vent Image');set(KmeansVDPMontage,'WindowState','minimized');
% montage(re_seg,'DisplayRange',[]); colormap(cmap);
kmeansDefectslices = kmeansDefectmap(:,:,:,slice_indices);
montage(kmeansDefectslices,'Size',[1 size(kmeansDefectslices,4)])
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
KmeansVDPMontagePosition=get(gcf,'position'); 

KmeansVDPMontage2 = figure('Name','Vent Image');set(KmeansVDPMontage2,'WindowState','minimized');
kmeansDefectslices = kmeansDefectmap2(:,:,:,slice_indices);
montage(kmeansDefectslices,'Size',[1 size(kmeansDefectslices,4)])
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
KmeansVDPMontagePosition2=get(gcf,'position'); 

parentPath = Ventilation.parentPath;
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
Global.exportToPPTX('addtext',['K-mean VDP = ',num2str(round(Ventilation.KmeansVDP,2)),'%'],'Position',[0 2 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');

Global.exportToPPTX('addpicture',KmeansVDPMontage2,'Position',...
    [0 3.5 NumSliceView*scalefactor scalefactor*NumSliceView*(KmeansVDPMontagePosition(4)/KmeansVDPMontagePosition(3))]);
Global.exportToPPTX('addtext',['K-mean VDP = ',num2str(round(Ventilation.KmeansVDP,2)),'%'],'Position',[0 3.5 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');

Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');
            Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
            Global.exportToPPTX('close');
            fprintf('PowerPoint file has been saved\n');               
close all;
segmentation = (segmentation + 1).*maskarray;
Ventilation.Kmeans_segmentation = segmentation;


save_data=[parentPath,'\','Ventilation_Analysis','.mat'];
save(save_data); 


end
