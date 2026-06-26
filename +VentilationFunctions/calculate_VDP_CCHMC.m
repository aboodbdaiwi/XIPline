function [Ventilation] = calculate_VDP_CCHMC(Ventilation,Proton,MainInput)
%% MATLAB script to perform VDP calculation
%
% Authors: Joseph Plummer
%         Abdullah Bdaiwi
% Date: 05/09/2021.

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% VDP Calculation Code:
MR = Ventilation.Image;

maskarray = double(Ventilation.LungMask);
maskarray(Ventilation.AirwayMask == 1) = 0;
maskarray(Ventilation.VesselMask == 1) = 0;

complete = Ventilation.CompleteThresh;
% complete = Ventilation.NoiseCompleteThresh;
incomplete = Ventilation.IncompleteThresh;
hyper = Ventilation.HyperventilatedThresh;
medfilter = Ventilation.MedianFilter;
DataPath = Ventilation.parentPath;

% Ensure trailing file separator
if ~endsWith(DataPath, filesep)
    DataPath = [DataPath filesep];
end

Overall_SNR = Ventilation.SNR_lung;
Ventilation.Threshold.Thresholds = [complete,incomplete,hyper];

cd(DataPath)
close all; 
warning('off','all') % Suppress all the tiff warnings
% Define output folder names based on whether the images are corrected or not:
% if N4_bias_analysis == 'no'
%     foldername = ["Threshold VDP Analysis\"];
% elseif N4_bias_analysis == 'yes'
%     foldername = ["N4 Bias Threshold VDP Analysis\"];
% else
%     disp('User must declare whether or they are analyzing N4 corrected images. Check settings/input_params(), and make sure N4_bias_analysis = yes or no.')
% end
cd(DataPath)
foldername = 'VDP_Analysis';
[~, lastFolder] = fileparts(DataPath);
if ~strcmp(lastFolder, foldername)
    parentPath = fullfile(DataPath, foldername);
else
    parentPath = DataPath;
end
if ~exist(parentPath, 'dir')
    mkdir(parentPath);
end
outputPath = parentPath;
cd(outputPath);

% Ensure only useful signal is operated on:
maskarray = double(maskarray);
MR = double(MR);
NormMR = MR.*(maskarray>0);

% Redefine thresholds to fit within 0-1 limits:
if complete || incomplete || hyper > 1
    complete = complete/100;
    incomplete = incomplete/100;
    hyper = hyper/100;
end

% Run VDP algorithm depending on whether the median filer is applied.
switch medfilter
    case "yes"
        % Calculate the defect ventilation percentage with median filter
        Incomplete_defect = VentilationFunctions.medFilter(double(NormMR<(mean(NormMR(NormMR>0))...
            *incomplete)*(maskarray>0)));
        Complete_defect = VentilationFunctions.medFilter(double(NormMR<(mean(NormMR(NormMR>0))...
            *complete)*(maskarray>0)));
        Hyper_defect = VentilationFunctions.medFilter(double(NormMR>(mean(NormMR(NormMR>0))*hyper)...
            *(maskarray>0)));
        defectArray = Incomplete_defect+Complete_defect+4*Hyper_defect;
    case "no"
        % Calculate the defect ventilation percentage without median filter
        Incomplete_defect = NormMR<(mean(NormMR(NormMR>0))...
            *incomplete)*(maskarray>0);
        Complete_defect = NormMR<(mean(NormMR(NormMR>0))...
            *complete)*(maskarray>0);
        Hyper_defect = NormMR>(mean(NormMR(NormMR>0))*2)...
            *(maskarray>0);
        defectArray = Incomplete_defect+Complete_defect+4*Hyper_defect;
    otherwise 
        fprintf('There is an error with median filter. Enter a value of 1 (yes) or 0 (no) in the input_params().\n')
end

% Calculate bin defect signal for histogram
d_Incomplete = NormMR.*(defectArray == 1); %Incomplete
d_Complete = NormMR.*(defectArray == 2);%Complete
d_Hyper = NormMR.*(defectArray == 4);%Hyper
d_Normal = (NormMR.*(defectArray == 0));%Normal

% Calculating percent of defects
VDP=100*(sum(defectArray==1,'all') + sum(defectArray==2,'all')) / sum(maskarray>0,'all');
Incomplete = 100*sum(defectArray == 1,'all')/sum(maskarray>0,'all');
Complete = 100*sum(defectArray == 2,'all')/sum(maskarray>0,'all');
Hyper = 100*sum(defectArray == 4,'all')/sum(maskarray>0,'all');
Normal = 100*(sum(maskarray>0,'all') - sum(defectArray == 1,'all') - sum(defectArray == 2,'all') ...
         - sum(defectArray == 4,'all'))/sum(maskarray>0,'all');

% Find mean and standard deviation of defects
md_Normal = mean(d_Normal(d_Normal>0));
sd_Normal = std(d_Normal(d_Normal>0));
md_Incomplete = mean(d_Incomplete(d_Incomplete>0));
sd_Incomplete = std(d_Incomplete(d_Incomplete>0));
md_Complete = mean(d_Complete(d_Complete>0));
sd_Complete = std(d_Complete(d_Complete>0));
md_Hyper = mean(d_Hyper(d_Hyper>0));
sd_Hyper = std(d_Hyper(d_Hyper>0));

% Equalize histogram bin counts for each defect
acnts = 50;
a = hist(d_Normal(d_Normal>0),acnts);
bcnts = length(a)*(max(d_Incomplete(d_Incomplete>0))- ...
    min(d_Incomplete(d_Incomplete>0)))/(max(d_Normal(d_Normal>0))-min(d_Normal(d_Normal>0)));
ccnts = length(a)*(max(d_Complete(d_Complete>0))- ...
    min(d_Complete(d_Complete>0)))/(max(d_Normal(d_Normal>0))-min(d_Normal(d_Normal>0)));
dcnts = length(a)*(max(d_Hyper(d_Hyper>0))- ...
    min(d_Hyper(d_Hyper>0)))/(max(d_Normal(d_Normal>0))-min(d_Normal(d_Normal>0)));
if isempty(dcnts) % 91-95 added due to lower Hypervent values/#of pixels causing error
    dcnts = 1;
elseif ((dcnts < 1) && (length(find(d_Hyper)) < 5))
    dcnts = 1;
end 
if isempty(ccnts)
    ccnts = 1;
elseif ((dcnts < 1) && (length(find(d_Complete)) < 5))
    ccnts = 1;
end 
% Create histogram figure of ventilation defects and print/save figure
figure('position',[350 350 750 350]);
hist(d_Normal(d_Normal>0),acnts);
hold on
hist(d_Incomplete(d_Incomplete>0),bcnts);
try
    hist(d_Complete(d_Complete>0),ccnts);
catch
    hist(d_Complete(d_Complete>0));
end
hist(d_Hyper(d_Hyper>0),dcnts)

% Image settings:
h = findobj(gcf,'Type','patch');
set(h(1),'FaceColor','b','EdgeColor','k','facealpha',0.5);
set(h(2),'FaceColor','r','EdgeColor','k','facealpha',0.5);
set(h(3),'FaceColor','y','EdgeColor','k','facealpha',0.5);
set(h(4),'FaceColor','w','EdgeColor','k','facealpha',0.5);
legend1 = sprintf('Normal: %0.2f±%0.2f (%0.1f%%)',md_Normal, sd_Normal, Normal);
legend2 = sprintf('Incomplete: %0.2f±%0.2f (%0.1f%%)',md_Incomplete,sd_Incomplete,Incomplete);
legend3 = sprintf('Complete: %0.2f±%0.2f (%0.1f%%)',md_Complete,sd_Complete, Complete);
legend4 = sprintf('Hyper: %0.2f±%0.2f (%0.1f%%)',md_Hyper,sd_Hyper, Hyper);
legend({legend1 legend2 legend3 legend4});
% title1 = sprintf('Ventilation Histogram');
title2 = sprintf('\nVentilation Defect Percentage %0.1f%%',VDP);
titlemain = [title2,'%'];
titlemain = sprintf(titlemain);
title({titlemain},'Fontweight','bold','FontSize',12);
set(gca,'FontSize',14)

set(gcf,'PaperPosition',[0 0 7.5 3.5]);

% print([foldername + 'Ventilation_Histogram'],'-dpng','-r300');
saveas(gca,'TH_Histogram.png');
%% Output mask images with proton overlays:

tiffFile = fullfile(outputPath,'MaskRegistered.tif');

forceTiffRegen = false;
if exist('MainInput','var') && ...
        isfield(MainInput,'ForceRegenerateReportTiffs') && ...
        MainInput.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end
if isfield(Ventilation,'ForceRegenerateReportTiffs') && ...
        Ventilation.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end

if forceTiffRegen
    if isfile(tiffFile)
        delete(tiffFile);
        fprintf('Deleted stale TIFF: %s\n', tiffFile);
    end

    if isfield(Ventilation,'Mask_Proton_Reg')
        Ventilation = rmfield(Ventilation,'Mask_Proton_Reg');
        fprintf('Removed stale Ventilation.Mask_Proton_Reg\n');
    end
end

if forceTiffRegen || ~isfile(tiffFile) || ~isfield(Ventilation,'Mask_Proton_Reg')
    disp('Saving MaskRegistered Tiff...')

    HImage = Proton.ProtonRegistered;

    frameH = size(HImage,1);
    frameW = size(HImage,2);
    nSlices = size(maskarray,3);

    MaskRegistered = uint8(zeros(frameH,frameW,3,nSlices));

    for slice = 1:nSlices
        lungMaskSlice = Ventilation.LungMask(:,:,slice);

        X = Global.robustOverlayRGB( ...
            abs(HImage(:,:,slice)), ...
            lungMaskSlice, ...
            [0,0.5*max(HImage(:))], ...
            [0 1 0], ...
            0.6);

        if size(X,1) ~= frameH || size(X,2) ~= frameW
            X = imresize(X,[frameH frameW]);
        end

        MaskRegistered(:,:,:,slice) = X;

        if slice == 1
            imwrite(X,tiffFile,'tif', ...
                'Description','Package Version: 1; Cohort: test');
        else
            imwrite(X,tiffFile,'tif','WriteMode','append', ...
                'Description','Package Version: 1; Cohort: test');
        end
    end

    MaskRegistered = permute(MaskRegistered,[1 2 4 3]);
    Ventilation.Mask_Proton_Reg = MaskRegistered;

    disp('MaskRegistered Tiff Completed.')
else
    disp('Using existing MaskRegistered Tiff / Ventilation.Mask_Proton_Reg')
end

MaskRegistered = Ventilation.Mask_Proton_Reg;
% S = orthosliceViewer((MaskRegistered)); %colormap(SixBinMap);
%% Output mask images with ventilation overlays:

tiffFile = fullfile(outputPath,'Mask_Vent_Reg.tif');

forceTiffRegen = false;
if exist('MainInput','var') && ...
        isfield(MainInput,'ForceRegenerateReportTiffs') && ...
        MainInput.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end
if isfield(Ventilation,'ForceRegenerateReportTiffs') && ...
        Ventilation.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end

if forceTiffRegen
    if isfile(tiffFile)
        delete(tiffFile);
        fprintf('Deleted stale TIFF: %s\n', tiffFile);
    end

    if isfield(Ventilation,'Mask_Vent_Reg')
        Ventilation = rmfield(Ventilation,'Mask_Vent_Reg');
        fprintf('Removed stale Ventilation.Mask_Vent_Reg\n');
    end
end

if forceTiffRegen || ~isfile(tiffFile) || ~isfield(Ventilation,'Mask_Vent_Reg')
    disp('Saving Mask_Vent_Reg Tiff...')

    MR = Ventilation.Image;
    scaledImage2 = MR/max(MR,[],'all');

    frameH = size(scaledImage2,1);
    frameW = size(scaledImage2,2);
    nSlices = size(maskarray,3);

    Mask_Vent_Reg = uint8(zeros(frameH,frameW,3,nSlices));

    for slice = 1:nSlices
        lungMaskSlice = Ventilation.LungMask(:,:,slice);

        X = Global.robustOverlayRGB( ...
            abs(scaledImage2(:,:,slice)), ...
            lungMaskSlice, ...
            [0,max(scaledImage2(:))], ...
            [1 0 1], ...
            0.4);

        if size(X,1) ~= frameH || size(X,2) ~= frameW
            X = imresize(X,[frameH frameW]);
        end

        Mask_Vent_Reg(:,:,:,slice) = X;

        if slice == 1
            imwrite(X,tiffFile,'tif', ...
                'Description','Package Version: 1; Cohort: test');
        else
            imwrite(X,tiffFile,'tif','WriteMode','append', ...
                'Description','Package Version: 1; Cohort: test');
        end
    end

    Mask_Vent_Reg = permute(Mask_Vent_Reg,[1 2 4 3]);
    Ventilation.Mask_Vent_Reg = Mask_Vent_Reg;

    disp('Mask_Vent_Reg Tiff Completed.')
else
    disp('Using existing Mask_Vent_Reg Tiff / Ventilation.Mask_Vent_Reg')
end

Mask_Vent_Reg = Ventilation.Mask_Vent_Reg;
% orthosliceViewer((Mask_Vent_Reg));
% orthosliceViewer(Mask_Vent_Reg);
% Global.imslice(Mask_Vent_Reg);


%% Output mask boundaries with ventilation overlays:

cd(outputPath);
tiffFile = fullfile(outputPath,'maskboundaries.tif');

forceTiffRegen = false;
if exist('MainInput','var') && ...
        isfield(MainInput,'ForceRegenerateReportTiffs') && ...
        MainInput.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end
if isfield(Ventilation,'ForceRegenerateReportTiffs') && ...
        Ventilation.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end

if forceTiffRegen
    if isfile(tiffFile)
        delete(tiffFile);
        fprintf('Deleted stale TIFF: %s\n', tiffFile);
    end

    if isfield(Ventilation,'Mask_Vent_boundaries')
        Ventilation = rmfield(Ventilation,'Mask_Vent_boundaries');
        fprintf('Removed stale Ventilation.Mask_Vent_boundaries\n');
    end
end

if forceTiffRegen || ~isfile(tiffFile) || ~isfield(Ventilation,'Mask_Vent_boundaries')
    disp('Saving maskboundaries Tiff...')

    MR = Ventilation.Image;
    MR2 = MR/max(MR,[],'all');

    frameH = size(MR2,1);
    frameW = size(MR2,2);
    nSlices = size(maskarray,3);

    Mask_Vent_boundaries = uint8(zeros(frameH,frameW,3,nSlices));

    for slice = 1:nSlices
        lungMaskSlice = maskarray(:,:,slice);
        maskboundaries = bwboundaries(lungMaskSlice);

        if isempty(maskboundaries)
            boundaryMask = false(size(lungMaskSlice));
        else
            boundaryMask = false(size(lungMaskSlice));
            for k = 1:length(maskboundaries)
                boundary = maskboundaries{k};
                boundaryMask(sub2ind(size(boundaryMask),boundary(:,1),boundary(:,2))) = true;
            end
            boundaryMask = imdilate(boundaryMask,strel('disk',1));
        end

        X = Global.robustOverlayRGB( ...
            abs(MR2(:,:,slice)), ...
            boundaryMask, ...
            [0,max(MR2(:))], ...
            [1 0 1], ...
            1);

        if size(X,1) ~= frameH || size(X,2) ~= frameW
            X = imresize(X,[frameH frameW]);
        end

        Mask_Vent_boundaries(:,:,:,slice) = X;

        if slice == 1
            imwrite(X,tiffFile,'tif', ...
                'Description','Package Version: 1; Cohort: test');
        else
            imwrite(X,tiffFile,'tif','WriteMode','append', ...
                'Description','Package Version: 1; Cohort: test');
        end
    end

    Mask_Vent_boundaries = permute(Mask_Vent_boundaries,[1 2 4 3]);
    Ventilation.Mask_Vent_boundaries = Mask_Vent_boundaries;

    disp('maskboundaries Tiff Completed.')
else
    disp('Using existing maskboundaries Tiff / Ventilation.Mask_Vent_boundaries')
end

Mask_Vent_boundaries = Ventilation.Mask_Vent_boundaries;


%% Output mask boundaries with proton overlays:

tiffFile = fullfile(outputPath,'protonmaskboundaries.tif');

forceTiffRegen = false;
if exist('MainInput','var') && ...
        isfield(MainInput,'ForceRegenerateReportTiffs') && ...
        MainInput.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end
if isfield(Ventilation,'ForceRegenerateReportTiffs') && ...
        Ventilation.ForceRegenerateReportTiffs
    forceTiffRegen = true;
end

if forceTiffRegen
    if isfile(tiffFile)
        delete(tiffFile);
        fprintf('Deleted stale TIFF: %s\n', tiffFile);
    end

    if isfield(Ventilation,'Mask_Proton_boundaries')
        Ventilation = rmfield(Ventilation,'Mask_Proton_boundaries');
        fprintf('Removed stale Ventilation.Mask_Proton_boundaries\n');
    end
end

if forceTiffRegen || ~isfile(tiffFile) || ~isfield(Ventilation,'Mask_Proton_boundaries')
    disp('Saving protonmaskboundaries Tiff...')

    % Important: use registered proton, not unregistered Proton.Image
    HImage = Proton.ProtonRegistered;

    frameH = size(HImage,1);
    frameW = size(HImage,2);
    nSlices = size(maskarray,3);

    Mask_Proton_boundaries = uint8(zeros(frameH,frameW,3,nSlices));

    for slice = 1:nSlices
        lungMaskSlice = Ventilation.LungMask(:,:,slice);
        maskboundaries = bwboundaries(lungMaskSlice);

        if isempty(maskboundaries)
            boundaryMask = false(size(lungMaskSlice));
        else
            boundaryMask = false(size(lungMaskSlice));
            for k = 1:length(maskboundaries)
                boundary = maskboundaries{k};
                boundaryMask(sub2ind(size(boundaryMask),boundary(:,1),boundary(:,2))) = true;
            end
            boundaryMask = imdilate(boundaryMask,strel('disk',1));
        end

        X = Global.robustOverlayRGB( ...
            abs(HImage(:,:,slice)), ...
            boundaryMask, ...
            [0,max(HImage(:))], ...
            [1 0 1], ...
            1);

        if size(X,1) ~= frameH || size(X,2) ~= frameW
            X = imresize(X,[frameH frameW]);
        end

        Mask_Proton_boundaries(:,:,:,slice) = X;

        if slice == 1
            imwrite(X,tiffFile,'tif', ...
                'Description','Package Version: 1; Cohort: test');
        else
            imwrite(X,tiffFile,'tif','WriteMode','append', ...
                'Description','Package Version: 1; Cohort: test');
        end
    end

    Mask_Proton_boundaries = permute(Mask_Proton_boundaries,[1 2 4 3]);
    Ventilation.Mask_Proton_boundaries = Mask_Proton_boundaries;

    disp('protonmaskboundaries Tiff Completed.')
else
    disp('Using existing protonmaskboundaries Tiff / Ventilation.Mask_Proton_boundaries')
end

Mask_Proton_boundaries = Ventilation.Mask_Proton_boundaries;

%% Output ventilation images with defect overlays:
% figure; Global.imslice(MR)
% figure; Global.imslice(scaledImage)
% Normalize the intensity of the original image to fall between [0,1].

% Now add the defect information to the images:
d_Incomplete_rgb = (defectArray == 1)*255;
d_Complete_rgb = (defectArray == 2)*255;
d_Hyper_rgb = (defectArray == 4)*255;
d_Normal_rgb = (defectArray == 0)*255;

% Define colormaps:
cm_incomplete = [0 0 0
    1 1 0]; % yellow
% cm_incomplete = [0 0 0
%     1 0 0];
array4dincomplete = zeros(size(d_Incomplete_rgb,1), size(d_Incomplete_rgb,2),3,...
    size(d_Incomplete_rgb,3));
for slice = 1:size(d_Incomplete_rgb,3)
    array4dincomplete(:,:,:,slice) =  ind2rgb(d_Incomplete_rgb(:,:,slice),cm_incomplete);
end
cm_complete=[0 0 0
    1 0 0]; % red
% cm_complete=[0 0 0
%     1 0 0];
array4dcomplete = zeros(size(d_Complete_rgb,1), size(d_Complete_rgb,2),3,...
    size(d_Complete_rgb,3));
for slice = 1:size(d_Complete_rgb,3)
    array4dcomplete(:,:,:,slice) =  ind2rgb(d_Complete_rgb(:,:,slice),cm_complete);
end
cm_hyper=[0 0 0
    0 0 1]; % change to none for now but you can switch back to blue
array4dhyper = zeros(size(d_Hyper_rgb,1), size(d_Hyper_rgb,2),3,...
    size(d_Hyper_rgb,3));
for slice = 1:size(d_Hyper_rgb,3)
    array4dhyper(:,:,:,slice) =  ind2rgb(d_Hyper_rgb(:,:,slice),cm_hyper);
end
cm_normal=[0 0 0
    1 1 1];
array4dnormal = zeros(size(d_Normal_rgb,1), size(d_Normal_rgb,2),3, size(d_Normal_rgb,3));
for slice = 1:size(d_Normal_rgb,3)
    array4dnormal(:,:,:,slice) =  ind2rgb(d_Normal_rgb(:,:,slice),cm_normal);
end

array4dMRrgb = zeros(size(scaledImage2,1), size(scaledImage2,2), 3,...
    size(scaledImage2,3));
array4dMRrgb(:,:,1,:) = scaledImage2;
array4dMRrgb(:,:,2,:) = scaledImage2;
array4dMRrgb(:,:,3,:) = scaledImage2;

% Turn 3D array of masks into 4D array to eliminate the background from the
% normal ventilation defect array mask
array4dMaskrgb = zeros(size(maskarray,1), size(maskarray,2), 3,...
    size(maskarray,3));
array4dMaskrgb(:,:,1,:) = maskarray;
array4dMaskrgb(:,:,2,:) = maskarray;
array4dMaskrgb(:,:,3,:) = maskarray;

% Eliminate background from 'normal' ventilation defect array mask and add
% all defect array color masks into one array
array4dMRrgb2=  array4dnormal.*array4dMRrgb;
defectArray_rgb2 = array4dincomplete+array4dcomplete+4*array4dhyper;
midslice = round((size(defectArray_rgb2,4)/2));
array4dMRrgb3 = defectArray_rgb2+array4dMRrgb2;
array4dMRrgb4 = uint8(255*array4dMRrgb3);

% save scaledImage2 images
% MR = uint16(Ventilation.Image);
% for m = 1:size(MR,3)
%     scaledImage(:,:,m) = imadjust(MR(:,:,m));
%     scaledImage2(:,:,m) = adapthisteq(scaledImage(:,:,m));
% end
% scaledImage2 = double(scaledImage2);
% scaledImage2 = scaledImage2 - min(scaledImage2(:));
% scaledImage2 = scaledImage2 ./ max(scaledImage2(:));

% save Vent images
% niftiwrite(abs(fliplr(rot90(scaledImage2,3))),[outputPath,'\ScaledVentilationImage.nii'],'Compressed',true);
% MR = double(Ventilation.Image);
% niftiwrite(abs(fliplr(rot90(MR,3))),[outputPath,'\VentilationImage.nii'],'Compressed',true);

% save mask images
% niftiwrite(abs(fliplr(rot90(maskarray,-1))),[outputPath,'\MaskImage.nii'],'Compressed',true);
close all;
%% %% write tiff and read back VDPVent maps
tiffFile = fullfile(outputPath,'VentDefectmap.tif');
if isfile(tiffFile)
    delete(tiffFile);
end

% Fixed output size based on the actual image size
frameH = size(array4dMRrgb4,1);
frameW = size(array4dMRrgb4,2);

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off', ...
    'Resize','off','Visible','off', ...
    'Units','pixels', ...
    'Position',[100 100 frameW frameH], ...
    'Color','black');
ax1 = axes('Parent',tiff, ...
    'Units','pixels', ...
    'Position',[1 1 frameW frameH], ...
    'Visible','off');
disp('Saving VentDefectmap Tiff...')

for slice = 1:size(array4dMRrgb4,4)
    cla(ax1,'reset');
    image(ax1,array4dMRrgb4(:,:,:,slice));
    axis(ax1,'image');
    axis(ax1,'off');
    set(ax1,'Units','pixels','Position',[1 1 frameW frameH]);
    set(ax1,'LooseInset',[0 0 0 0]);

    drawnow;
    % Capture axes only, not the full figure
    X = getframe(ax1).cdata;
    % Force every saved page to be identical size
    if size(X,1) ~= frameH || size(X,2) ~= frameW
        X = imresize(X,[frameH frameW]);
    end
    X = uint8(X);
    if slice == 1
        imwrite(X,tiffFile,'tif', ...
            'Description',strcat('Package Version: ','1','; Cohort: ','test'));
    else
        imwrite(X,tiffFile,'tif','WriteMode','append', ...
            'Description',strcat('Package Version: ','1','; Cohort: ','test'));
    end
end
close(tiff)

disp('Saving VentDefectmap Tiff Completed.')
tiff_info = imfinfo(tiffFile);
if numel(tiff_info) ~= size(array4dMRrgb4,4)
    warning('Expected %d slices, but saved TIFF contains %d slices.', ...
        size(array4dMRrgb4,4), numel(tiff_info));
end
VentDefectmap = uint8(zeros(frameH,frameW,3,numel(tiff_info)));
for ii = 1:numel(tiff_info)
    temp_tiff = imread(tiffFile,ii);
    if size(temp_tiff,1) ~= frameH || size(temp_tiff,2) ~= frameW
        temp_tiff = imresize(temp_tiff,[frameH frameW]);
    end
    if size(temp_tiff,3) == 1
        temp_tiff = repmat(temp_tiff,[1 1 3]);
    end
    VentDefectmap(:,:,:,ii) = uint8(temp_tiff);
end
cd(DataPath)
% S = orthosliceViewer(permute (VentDefectmap, [1,2,4,3])); %colormap(SixBinMap);
%% Generate montages
% check images size, 2D or 3D
MR = Ventilation.Image;
MR2 = MR / max(MR,[], 'all');
scaledImage = MR2;
scaledImage2 = MR2;
if size(Ventilation.Image,3) > 40
    Image_3D = 1;
    nrow = ceil(sqrt(size(scaledImage,3)));
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
else
    Image_3D = 0;
end 
% Xe images:
subplot(3,1,1)
if Image_3D == 1
HPXeMontage = montage(reshape(scaledImage,[size(scaledImage,1),...
    size(scaledImage,2), 1, size(scaledImage,3)]),'Size',...
    [nrow nrow],'DisplayRange',[min(scaledImage(:)) max(scaledImage(:))]);
else
HPXeMontage = montage(reshape(scaledImage,[size(scaledImage,1),...
    size(scaledImage,2), 1, size(scaledImage,3)]),'Size',...
    [1 size(scaledImage,3)],'DisplayRange',[min(scaledImage(:)) max(scaledImage(:))]);
end

title('Xenon Montage')
savedMontage=get(HPXeMontage,'CData');
savedMontage2 = imadjust(savedMontage);
outputFileName2 = 'HPXeMontage.png';
fullFileName2 = char(fullfile(outputPath, outputFileName2));
imwrite(savedMontage2,fullFileName2,'png');

% Mask images:
subplot(3,1,2)
MaskRegistered = permute(MaskRegistered,[1 2 4 3]);
if Image_3D == 1
MaskMontage=montage(reshape(MaskRegistered,[size(MaskRegistered,1),...
    size(MaskRegistered,2), size(MaskRegistered,3), size(MaskRegistered,4)]),...
    'Size',[nrow nrow],'DisplayRange',[]);
else
MaskMontage=montage(reshape(MaskRegistered,[size(MaskRegistered,1),...
    size(MaskRegistered,2), size(MaskRegistered,3), size(MaskRegistered,4)]),...
    'Size',[1 size(MaskRegistered,4)],'DisplayRange',[]);
end
title('Mask Array')
savedMaskMontage = get(MaskMontage,'CData');
% savedMaskMontage2 = imadjust(savedMaskMontage);
outputFileName2 = 'MaskMontage.png';
fullFileName2 = fullfile(outputPath, outputFileName2);
imwrite(savedMaskMontage,fullFileName2,'png');

% Defect images:
subplot(3,1,3)
if Image_3D == 1
DefectArrayMontage=montage(reshape(array4dMRrgb4,[size(array4dMRrgb4,1),...
    size(array4dMRrgb4,2),size(array4dMRrgb4,3),size(array4dMRrgb4,4)]),...
    'Size',[nrow nrow],'DisplayRange',[]);
else
DefectArrayMontage=montage(reshape(array4dMRrgb4,[size(array4dMRrgb4,1),...
    size(array4dMRrgb4,2),size(array4dMRrgb4,3),size(array4dMRrgb4,4)]),...
    'Size',[1 size(array4dMRrgb4,4)],'DisplayRange',[]);
end
title('Defect Array')
savedDefectArrayMontage=get(DefectArrayMontage,'CData');
outputFileName2= 'DefectArrayMontage.png';
fullFileName2 = fullfile(outputPath, outputFileName2);
imwrite(savedDefectArrayMontage,fullFileName2,'png');

% Proton Images
subplot(3,1,3)
Proton.ProtonRegistered = Proton.ProtonRegistered/max(Proton.ProtonRegistered(:));
if Image_3D == 1
ProtonImages=montage(reshape(Proton.ProtonRegistered,[size(Proton.ProtonRegistered,1),...
    size(Proton.ProtonRegistered,2),size(Proton.ProtonRegistered,3)]),...
    'Size',[nrow nrow],'DisplayRange',[0 1]);
else
ProtonImages=montage(reshape(Proton.ProtonRegistered,[size(Proton.ProtonRegistered,1),...
    size(Proton.ProtonRegistered,2),size(Proton.ProtonRegistered,3)]),...
    'Size',[1 size(Proton.ProtonRegistered,3)],'DisplayRange',[0 1]);
end
title('Proton Montage')
savedProtonImages=get(ProtonImages,'CData');
outputFileName2= 'ProtonMontage.png';
fullFileName2 = fullfile(outputPath, outputFileName2);
imwrite(savedProtonImages,fullFileName2,'png');

% Xe_Vent Images
subplot(4,1,4)
if isfield(Ventilation, 'Mask_Vent_Reg')
    Mask_Vent_Reg = Ventilation.Mask_Vent_Reg;
end
Mask_Vent_Reg = permute(Mask_Vent_Reg,[1 2 4 3]);
if Image_3D == 1
Mask_Vent_RegMontage=montage(reshape(Mask_Vent_Reg,[size(Mask_Vent_Reg,1),...
    size(Mask_Vent_Reg,2), size(Mask_Vent_Reg,3), size(Mask_Vent_Reg,4)]),...
    'Size',[nrow nrow],'DisplayRange',[]);
else
Mask_Vent_RegMontage=montage(reshape(Mask_Vent_Reg,[size(Mask_Vent_Reg,1),...
    size(Mask_Vent_Reg,2), size(Mask_Vent_Reg,3), size(Mask_Vent_Reg,4)]),...
    'Size',[1 size(Mask_Vent_Reg,4)],'DisplayRange',[]);
end
title('Mask Array')
savedMask_Vent_RegMontage = get(Mask_Vent_RegMontage,'CData');
% savedMaskMontage2 = imadjust(savedMaskMontage);
outputFileName2 = 'MaskMontage.png';
fullFileName2 = fullfile(outputPath, outputFileName2);
imwrite(savedMask_Vent_RegMontage,fullFileName2,'png');


%% save result in a powerpoint file // Abdullah 1/27/2020
%% 
% cd(outputPath)
close all;
% %how many slices to incluse in the final powerpoint file
sl_include=zeros(1,size(maskarray,3));
for sl=1:size(maskarray,3)
    choose_slice = maskarray(:,:,sl);       
    if sum(choose_slice(:)) > 1
         sl_include(sl) = sl;
    else
         sl_include(sl) = 0;
    end 
end
sl_include( :, ~any(sl_include,1) ) = [];  %columns
sl_1 = sl_include(1);
sl_end = sl_include(end);
NumSliceView = length(sl_include);
if NumSliceView > 16
    NumSliceView = 16;
end
if Image_3D ~= 1
    slice_indices = sl_1:1:sl_end;
end
% if you want the un-segemted slices not to show up on the powepoint, run
% the above code and add this line (,'Indices', sl_1:sl_end) to the end of
% each Montage. 
%%% get rid of the gray frame in each montage
close all;
zeroImage = zeros(size(scaledImage2));  %scaledImage2
zeroMontage = figure('Name','zero Image');set(zeroMontage,'WindowState','minimized');
montage(reshape(zeroImage,[size(zeroImage,1), size(zeroImage,2), 1, size(zeroImage,3)]),...
    'Size',[1 size(zeroImage,3)]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

VentscaledImage = scaledImage2(:,:,slice_indices);
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(reshape(VentscaledImage,[size(VentscaledImage,1), size(VentscaledImage,2), 1, size(VentscaledImage,3)]),...
    'Size',[1 size(VentscaledImage,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
VentMontagePosition=get(gcf,'position'); 

ProtonRegisteredMontage = figure('Name','Proton Image');set(ProtonRegisteredMontage,'WindowState','minimized');
if strcmp(MainInput.NoProtonImage,'yes') == 1 
    ProtonImage = Proton.ProtonRegisteredColored(:,:,slice_indices);
    montage(reshape(ProtonImage,[size(ProtonImage,1), size(ProtonImage,2), 1, size(ProtonImage,3)]),...
    'Size',[1 size(ProtonImage,3)],'DisplayRange',[0 1]);
else 
    ProtonImage = permute(Proton.ProtonRegisteredColored, [1 2 4 3]);
    ProtonImage = ProtonImage(:,:,:,slice_indices);
    montage(reshape(ProtonImage,[size(ProtonImage,1), size(ProtonImage,2),size(ProtonImage,3),...
    size(ProtonImage,4)]),'Size',[1 size(ProtonImage,4)],'DisplayRange',[]);
end
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
ProtonRegisteredPosition = get(gcf,'position'); 

maskarrayimg = MaskRegistered(:,:,:,slice_indices);
MaskMontage = figure('Name','Mask');set(MaskMontage,'WindowState','minimized');
montage(reshape(maskarrayimg,[size(maskarrayimg,1), size(maskarrayimg,2), size(maskarrayimg,3), size(maskarrayimg,4)]),...
    'Size',[1 size(maskarrayimg,4)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
MaskMontagePosition = get(gcf,'position'); 

DefectArray = array4dMRrgb4(:,:,:,slice_indices);
DefectArrayMontage = figure('Name','Defects');set(DefectArrayMontage,'WindowState','minimized');
montage(reshape(DefectArray,[size(DefectArray,1), size(DefectArray,2),size(DefectArray,3),...
    size(DefectArray,4)]),'Size',[1 size(DefectArray,4)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
DefectArrayMontagePosition = get(gcf,'position'); 

protonarrayimg = Proton.ProtonRegistered(:,:,slice_indices);
ProtonImages = figure('Name','Mask');set(ProtonImages,'WindowState','minimized');
if sum(protonarrayimg(:)) > 0
    montage(reshape(protonarrayimg,[size(protonarrayimg,1), size(protonarrayimg,2), size(protonarrayimg,3), size(protonarrayimg,4)]),...
        'Size',[1 size(protonarrayimg,3)],'DisplayRange',[0 max(Proton.ProtonRegistered(:))/2]);
else
    montage(reshape(protonarrayimg,[size(protonarrayimg,1), size(protonarrayimg,2), size(protonarrayimg,3), size(protonarrayimg,4)]),...
        'Size',[1 size(protonarrayimg,3)]);
end
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
ProtonMontagePosition = get(gcf,'position');

Mask_Vent_Regarrayimg = Mask_Vent_Reg(:,:,:,slice_indices);
Mask_Vent_RegMontage = figure('Name','Mask');set(Mask_Vent_RegMontage,'WindowState','minimized');
montage(reshape(Mask_Vent_Regarrayimg,[size(Mask_Vent_Regarrayimg,1), size(Mask_Vent_Regarrayimg,2), size(Mask_Vent_Regarrayimg,3), size(Mask_Vent_Regarrayimg,4)]),...
    'Size',[1 size(Mask_Vent_Regarrayimg,4)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
Mask_Vent_RegMontagePosition = get(gcf,'position'); 

% save ppt 
cd(DataPath)
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

    % % Add human-readable date to the slide (e.g., "June 12, 2025")
    % formatted_date = datestr(today, 'mmmm dd, yyyy');
    % Global.exportToPPTX('addtext', ['Date: ', formatted_date], ...
    %     'Position', [12, 0.2, 4, 1], 'FontSize', 12, ...
    %     'FontWeight', 'normal', 'HorizontalAlignment', 'right');
end

scalefactor = 16/NumSliceView; 
%Add slides
Global.exportToPPTX('addslide'); % Image/mask/VDP
% Global.exportToPPTX('addtext',sprintf(foldername),'Position',[5 0 7 1],'Color','b','FontSize',25);
% Title = ['| VDP Analysis |', MainInput.PatientInfo];
Title = ['VDP Analysis'];
Global.exportToPPTX('addtext',Title ,'Position',[4.2 -0.1 11.5 1],'Color','b','FontSize',25);
Global.exportToPPTX('addpicture',VentMontage,'Position',...
    [0 0.4 NumSliceView*scalefactor scalefactor*NumSliceView*(MaskMontagePosition(4)/MaskMontagePosition(3))]);
Global.exportToPPTX('addpicture',ProtonImages,'Position',...
    [0 1.6 NumSliceView*scalefactor scalefactor*NumSliceView*(ProtonMontagePosition(4)/ProtonMontagePosition(3))]);
Global.exportToPPTX('addpicture',ProtonRegisteredMontage,'Position',...
    [0 3 NumSliceView*scalefactor scalefactor*NumSliceView*(ProtonRegisteredPosition(4)/ProtonRegisteredPosition(3))]);
Global.exportToPPTX('addpicture',Mask_Vent_RegMontage,'Position',...
    [0 4.2 NumSliceView*scalefactor scalefactor*NumSliceView*(Mask_Vent_RegMontagePosition(4)/Mask_Vent_RegMontagePosition(3))]);
Global.exportToPPTX('addpicture',MaskMontage,'Position',...
    [0 5.5 NumSliceView*scalefactor scalefactor*NumSliceView*(MaskMontagePosition(4)/MaskMontagePosition(3))]);
% Global.exportToPPTX('addpicture',DefectArrayMontage,'Position',...
%     [0 4 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);
Global.exportToPPTX('addpicture',DefectArrayMontage,'Position',...
    [0 7 NumSliceView*scalefactor scalefactor*NumSliceView*(DefectArrayMontagePosition(4)/DefectArrayMontagePosition(3))]);
% Global.exportToPPTX('addpicture',ProtonImages,'Position',...
%     [0 5.2 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);

VDP_hist = imread(fullfile(outputPath ,'TH_Histogram.png'));    
%   Global.exportToPPTX('addpicture',VDP_hist,'Position',[0 5 8.5 8.5*(size(VDP_hist,1)/size(VDP_hist,2))]);
Global.exportToPPTX('addtext',['129Xe-Vent Images (SNR=',num2str(Overall_SNR),')'],'Position',[0 0.4 4.2 0.4],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Anatomical Proton','Position',[0 1.6 3 0.4],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Registered Proton Images','Position',[0 3 4 0.4],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Xenon-Mask Overlay ','Position',[0 4.2 3 0.4],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Proton-Mask Overlay','Position',[0 5.5 3 0.4],'Color','r','FontSize',20,'BackgroundColor','k');
vdp_title = sprintf('VDP = %0.1f%%',VDP);
% Global.exportToPPTX('addtext',vdp_title,'Position',[0 3.6 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext',vdp_title,'Position',[0 7 2 0.4],'Color','r','FontSize',20,'BackgroundColor','k');


Global.exportToPPTX('save',fullfile(DataPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');             

close all;
%% compute VDP per Slice
dmap = defectArray;
dmap(dmap == 0) = 3;
dmap = dmap.* maskarray;
array_temp = dmap; % Create a temporary copy
dmap(array_temp == 1) = 2; % Change 1s to 2s
dmap(array_temp == 2) = 1; % Change 2s to 1s
Ventilation.Threshold.defectArray = dmap;

defect_mask = Ventilation.Threshold.defectArray;
defect_mask(defect_mask > 2) = 0;
defect_mask = double(defect_mask > 0);
defectMap = double(Ventilation.LungMask) + defect_mask;  
[vdp_per_slice_local, vdp_per_slice_global] = VentilationFunctions.computeVDPperSlice(defectMap);
Ventilation.Threshold.vdp_per_slice_local = vdp_per_slice_local;
Ventilation.Threshold.vdp_per_slice_global = vdp_per_slice_global;

%% Save workspace data:
legend1 = sprintf('%0.2f±%0.2f (%0.1f%%)',md_Normal, sd_Normal, Normal);
legend2 = sprintf('%0.2f±%0.2f (%0.1f%%)',md_Incomplete,sd_Incomplete,Incomplete);
legend3 = sprintf('%0.2f±%0.2f (%0.1f%%)',md_Complete,sd_Complete, Complete);
legend4 = sprintf('%0.2f±%0.2f (%0.1f%%)',md_Hyper,sd_Hyper, Hyper);

Ventilation.Threshold.VDP = VDP;
Ventilation.Threshold.VentscaledImage = VentscaledImage;
Ventilation.Threshold.DefectArray = DefectArray;

Ventilation.Threshold.VDP_hist = VDP_hist;
Ventilation.Threshold.VentDefectmap = VentDefectmap;
Ventilation.Threshold.THBins = [Complete,Incomplete,Normal,Hyper];
Ventilation.Threshold.legend1 = legend1;
Ventilation.Threshold.legend2 = legend2;
Ventilation.Threshold.legend3 = legend3;
Ventilation.Threshold.legend4 = legend4;

save_title = fullfile(outputPath ,'VDPThresholdAnalysis.mat');
save(save_title);

end