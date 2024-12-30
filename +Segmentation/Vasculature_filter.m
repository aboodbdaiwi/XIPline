
function [Ventilation] = Vasculature_filter(Proton, Ventilation, MainInput)
clc;

warning('off', 'all');
% Load the 3D lung MRI image (assuming the image is already pre-processed)
switch MainInput.vesselImageMode
    case 'xenon'
        options.FrangiScaleRange = [1 1]; % Scale range for sigma (in voxels)
        options.FrangiScaleRatio = 1; % Step size between scales
        options.FrangiBeta = .2; % Controls sensitivity to deviations from tubular shape
        options.BlackWhite = true; % Set to false because vessels are dark in your images
        options.FrangiAlpha = 0.5;
        options.frangi_thresh = MainInput.frangi_thresh;
        options.vessel_thresh = 5;
        % options.FrangiC = 1500;
        image3D = Ventilation.Image;
    case 'proton' % not optimized by Abood 
        options.FrangiScaleRange = [1 1]; % Scale range for sigma (in voxels)
        options.FrangiScaleRatio = 2; % Step size between scales
        options.FrangiBeta = .2; % Controls sensitivity to deviations from tubular shape
        options.BlackWhite = false; 
        options.FrangiAlpha = 0.5;
        options.frangi_thresh = MainInput.frangi_thresh;
        options.vessel_thresh = 5;
        % options.FrangiC = 1500;
        image3D = Proton.Image;    
end
disp(options);
image3D = (image3D - min(image3D(:)))/(max(image3D(:)) - min(image3D(:)));

% image3D = VentilationFunctions.medFilter(image3D);

% Convert image to double for processing
image3D = double(image3D);
lungsmask = Ventilation.LungMask;
[nRows,nCols,nSlices] = size(image3D);

se=strel('disk',2);
se4=strel('disk',4);
se5=strel('disk',10);
se6=strel('disk',6);

if strcmp(MainInput.SliceOrientation,'coronal')
    outter_boundary = VentilationFunctions.Akmeans.split_lung_boundaries_coronal(lungsmask,1:nSlices);
else
    mask_struct.lungsmask = lungsmask;
    outter_boundary = VentilationFunctions.Akmeans.split_lung_boundaries(mask_struct,1:nSlices);
end

% Apply the Frangi filter
vessel_stack = zeros(size(image3D));
VesselMask = zeros(size(image3D));
shrunk_lungmask = zeros(size(image3D));
dilated_lung = zeros(size(image3D));

for i = 1:size(VesselMask, 3)
    [vesselnessSlice, ~] = Segmentation.frangi_filterv2.FrangiFilter2D(image3D(:,:,i), options);
    vessel_stack(:,:,i) = vesselnessSlice; % filtered_at_sigmas;
    shrunk_lungmask(:,:,i) = imerode(imdilate(lungsmask(:,:,i),se4),se6);
    dilated_lung(:,:,i) = imdilate(edge(outter_boundary(:,:,i)),se5);
        
    VesselMask(:,:,i) = double(vesselnessSlice);
end
% figure; imslice(dilated_lung, 'vesselMask1'); % Display the original 3D image

VesselMask = (VesselMask - min(VesselMask(:)))/(max(VesselMask(:)) - min(VesselMask(:)));
lungMaskLogical = logical(Ventilation.LungMask(:));
NFactor = prctile(VesselMask(lungMaskLogical), 99.0);
        VesselMask = VesselMask/NFactor;
        VesselMask(VesselMask > 1) = 1;
% vesselMask1 = vesselMask1.*NFactor;
% vesselMask1(vesselMask1 > NFactor) = 0;
% switch MainInput.vesselImageMode
%     case 'xenon'
%         shrunk_mask = shrunk_lungmask;
%         VesselMask = double(VesselMask.*lungsmask);
%     case 'proton' % not optimized by Abood 
%         radius = 5; % Shrink by 5 voxels
%         semask = strel('sphere', radius); % 3D spherical structuring element
%         shrunk_mask = imerode(lungsmask, semask);       
% end
VesselMask = double(VesselMask.*lungsmask);
VesselMask = (VesselMask - min(VesselMask(:)))/(max(VesselMask(:)) - min(VesselMask(:)));
VesselMask = double(imbinarize(VesselMask,options.frangi_thresh)); % "global" or "adaptive". bright  dark
% vesselMask = vesselMask.*shrunk_lungmask;
VesselMask = VesselMask - dilated_lung;
VesselMask(VesselMask<0)=0;
VesselMask = double(bwareaopen(VesselMask,options.vessel_thresh));
Ventilation.VesselMask = VesselMask;
Ventilation.vessel_stack = vessel_stack;

disp('Save NIFTI')
niftiwrite((fliplr(rot90(VesselMask,-1))),[MainInput.XeDataLocation,'\vessel_mask.nii'],'Compressed',true);
niftiwrite((fliplr(rot90(vessel_stack,-1))),[MainInput.XeDataLocation,'\vessel_stack.nii'],'Compressed',true);

save_data=[MainInput.XeDataLocation,'\','VesselMask','.mat'];
save(save_data);

disp('vesselMask completed')

end