
function [Proton,MainInput,GasExchange] = GeneralRegisterProton_to_Xenon(...
    Proton,...
    XeImage,...
    MainInput,GasExchange)
%   Inputs:
%          
%   Outputs:
%      LungMask
%    'translation' | 'rigid' | 'similarity' | 'affine'
%   Example: 
%
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% match number of slices 

moving1 = double(Proton.Image);
if strcmp(MainInput.AnalysisType,'Ventilation')
    try
        fixed1 = XeImage.Image;
    catch
        fixed1 = XeImage;
    end
elseif strcmp(MainInput.AnalysisType,'GasExchange')
    try
        fixed1 = XeImage.VentImage;
    catch
        fixed1 = XeImage;
    end
end
%     figure; Global.imslice(fixed1);

HImage = Proton.Image;
nSlice = size(fixed1,3);
if strcmp(MainInput.AnalysisType,'GasExchange') && size(HImage,1) <= 128
    HImage = imresize3(HImage, size(fixed1));   
    Proton.ProtonImageHR = imresize3(Proton.ProtonImageHR, size(fixed1));  
    moving1 = HImage;
    ProtonVoxelInfo = MainInput.ProtonVoxelInfo;
    XeVoxelInfo = MainInput.XeVoxelInfo;
    if any(size(moving1) == size(fixed1))
        ProtonVoxelInfo.PixelSize1 = XeVoxelInfo.PixelSize1;
        ProtonVoxelInfo.PixelSize2 = XeVoxelInfo.PixelSize2;
        ProtonVoxelInfo.SliceThickness = XeVoxelInfo.SliceThickness;
    end
elseif strcmp(MainInput.AnalysisType,'GasExchange') && size(HImage,1) > 130
    % Get the size of the arrays
    size1 = size(fixed1);
    size2 = size(HImage);
    
    % Calculate the starting and ending indices for cropping
    startIdx = (size2 - size1) / 2 + 1;
    endIdx = startIdx + size1 - 1;
    
    % Crop the center portion of array2
    HImage = HImage(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3));
    Proton.ProtonImageHR = Proton.ProtonImageHR(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3));
    % Verify the size of the cropped array
    disp(size(HImage));

    moving1 = HImage;
    ProtonVoxelInfo = MainInput.ProtonVoxelInfo;
    XeVoxelInfo = MainInput.XeVoxelInfo;
    if any(size(moving1) == size(fixed1))
        ProtonVoxelInfo.PixelSize1 = XeVoxelInfo.PixelSize1;
        ProtonVoxelInfo.PixelSize2 = XeVoxelInfo.PixelSize2;
        ProtonVoxelInfo.SliceThickness = XeVoxelInfo.SliceThickness;
    end
else
    ProtonVoxelInfo = MainInput.ProtonVoxelInfo;
    XeVoxelInfo = MainInput.XeVoxelInfo;
end
RegistrationType = MainInput.RegistrationType;
if strcmp(MainInput.AnalysisType,'GasExchange')
    DataLocation = fullfile(MainInput.XeDataLocation,'Gas Exchange Analysis');
else
    DataLocation = fullfile(MainInput.XeDataLocation,'Ventilation Analysis');
end
%     figure; Global.imslice(moving1);

cd(DataLocation)
H_RecMatrix = Proton.H_RecMatrix;
% if strcmp(MainInput.AnalysisType,'Ventilation') == 1 
%     mkdir([MainInput.XeDataLocation '\Ventilation Analysis']);
%     DataLocation = [MainInput.XeDataLocation,'\Ventilation Analysfullfis'];
%     MainInput.XeDataLocation = DataLocation;
% elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
%     mkdir([MainInput.XeDataLocation '\GasExchange Analysis']);
%     DataLocation = [MainInput.XeDataLocation,'\GasExchange Analysis'];
%     MainInput.XeDataLocation = DataLocation;
% end

%% 
disp('performing registration, please wait...');
% match number of slices 
if (MainInput.SliceSelection == 1) && any(size(moving1) ~= size(fixed1))
    moving1 = squeeze(moving1(:,:,MainInput.Hstart:MainInput.Hend));
    fixed1 = squeeze(fixed1(:,:,MainInput.Xestart:MainInput.Xeend));
    % maxsize = max(size(fixed1,1), size(fixed1,2));
    % if size(fixed1, 1) ~= size(fixed1, 2) 
    %     resizefixed1 = zeros(maxsize, maxsize, size(fixed1,3));
    %     for i = 1:size(fixed1, 3)
    %         resizefixed1(:,:,i) = imresize(squeeze(fixed1(:,:,i)),[maxsize, maxsize]);
    %     end
    %     fixed1 = resizefixed1;
    % end    
    % if size(moving1, 1) ~= size(fixed1, 1) || size(moving1, 2) ~= size(fixed1, 2)
    %     resizemoving1 = zeros(size(fixed1,1), size(fixed1,2), size(moving1,3));
    %     for i = 1:size(moving1, 3)
    %         resizemoving1(:,:,i) = imresize(squeeze(moving1(:,:,i)),[maxsize, maxsize]);
    %     end
    %     moving1 = resizemoving1;
    % end
    if size(moving1,3) < size(fixed1,3)
        Im1 = fixed1;
        Im2 = moving1;
        [moving1] = Global.match_n_slices(Im1,Im2);        
    elseif size(moving1,3) > size(fixed1,3)
        Im1 = moving1;
        Im2 = fixed1;
        [fixed1] = Global.match_n_slices(Im1,Im2);
    end
    nSlice = size(Im1,3);
end


% % Get the size of the original array
% [Xerows, Xecols, Xeslices] = size(fixed1);
% [Hrows, Hcols, ~] = size(moving1);
% % Check if the third dimension is less than 16
% if Xeslices < 16
%     % Create a new 3D array of size 100x100x16 filled with zeros
%     newfixed1 = zeros(Xerows, Xecols, 16);
%     newmoving1 = zeros(Hrows, Hcols, 16);
%     % Calculate the starting index to place the original array in the middle
%     startSlice = floor((16 - Xeslices) / 2) + 1;
%     endSlice = startSlice + Xeslices - 1;
% 
%     % Place the original array in the middle of the new array
%     newfixed1(:,:,startSlice:endSlice) = fixed1;
%     newmoving1(:,:,startSlice:endSlice) = moving1;
%     fixed = newfixed1;
%     moving = newmoving1;
% else
%     moving = moving1;
%     fixed = fixed1;
%     disp('The third dimension is already 16 or greater. No resizing needed.');
% end

% check image dimenssion if less than 16
if size(moving1,3) < 16  
    moving = zeros(size(moving1,1),size(moving1,2),16);
    moving(:,:,1:size(moving1,3)) = abs(moving1);
else
    moving = moving1;
end
if size(fixed1,3) < 16
    fixed = zeros(size(fixed1,1),size(fixed1,2),16);
    fixed(:,:,1:size(fixed1,3)) = abs(fixed1); 
else
    fixed = fixed1;
end
% figure; imslice(moving) 
% figure; imslice(fixed) 

%% create a lung mask
LungMask = double(Segmentation.SegmentLungthresh(fixed,1,1));
se = strel('sphere', 3);  % Structuring element (sphere with radius 1)

% Perform erosion on the 3D binary mask
LungMask = imerode(LungMask, se);
% figure; imslice(LungMask)
%% Register Multimodal 3-D Medical Images.
% if strcmp(MainInput.AnalysisType,'Ventilation')
%     movingVolume = moving;
% elseif strcmp(MainInput.AnalysisType,'GasExchange')
%     movingVolume = moving(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix);
% end
movingVolume = moving;
fixedVolume = fixed;

% helperVolumeRegistration(fixedVolume,movingVolume);
% centerFixed = size(fixedVolume)/2;
% centerMoving = size(movingVolume)/2;
% figure
% imshowpair(movingVolume(:,:,centerMoving(3)), fixedVolume(:,:,centerFixed(3)));
% title('Unregistered Axial Slice')
[optimizer,metric] = imregconfig('multimodal');
Rfixed  = imref3d(size(fixedVolume),XeVoxelInfo.PixelSize2,XeVoxelInfo.PixelSize1,XeVoxelInfo.SliceThickness);
Rmoving = imref3d(size(movingVolume),ProtonVoxelInfo.PixelSize2,ProtonVoxelInfo.PixelSize1,ProtonVoxelInfo.SliceThickness);
% Rmoving.XWorldLimits
% Rmoving.PixelExtentInWorldX
% Rmoving.ImageExtentInWorldX
optimizer.InitialRadius = 0.00005;
% movingRegisteredVolume = imregister(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
% figure
% imshowpair(movingRegisteredVolume(:,:,centerFixed(3)), fixedVolume(:,:,centerFixed(3)));
% title('Axial Slice of Registered Volume')
% helperVolumeRegistration(fixedVolume,movingRegisteredVolume);
geomtform = imregtform(movingVolume,Rmoving, fixedVolume,Rfixed, RegistrationType, optimizer, metric);
tform = geomtform;
% geomtform.T
% centerXWorld = mean(Rmoving.XWorldLimits);
% centerYWorld = mean(Rmoving.YWorldLimits);
% centerZWorld = mean(Rmoving.ZWorldLimits);
% [xWorld,yWorld,zWorld] = transformPointsForward(geomtform,centerXWorld,centerYWorld,centerZWorld);
% [r,c,p] = worldToSubscript(Rfixed,xWorld,yWorld,zWorld);
ProtonRegistered = imwarp(movingVolume,Rmoving,geomtform,'bicubic','OutputView',Rfixed);
if strcmp(MainInput.AnalysisType,'GasExchange')
    ProtonHRRegistered = imwarp(Proton.ProtonImageHR,Rmoving,geomtform,'bicubic','OutputView',Rfixed);
end
% figure; imslice(ProtonRegistered)
%% stor data
for slice =1:size(fixedVolume,3)
    A = ProtonRegistered(:,:,slice);
    B = fixedVolume(:,:,slice);
    ProtonRegisteredColored(:,:,:,slice) = imfuse(A,B,'falsecolor','ColorChannels','green-magenta'); %'red-cyan'|'green-magenta'
end
ProtonRegistered = ProtonRegistered(:,:,1:nSlice);
ProtonRegisteredColored = ProtonRegisteredColored(:,:,:,1:nSlice);
Proton.ProtonRegistered = ProtonRegistered;
Proton.Rmoving = Rmoving;
Proton.Rfixed = Rfixed;
Proton.optimizer = optimizer;
Proton.geomtform = geomtform;
Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
% figure; orthosliceViewer(Proton.ProtonRegisteredColored)
%% view

if strcmp(MainInput.AnalysisType,'GasExchange')
    %Determine Slices to Plot
    NumPlotSlices = 7;
    if size(fixed,3) < 7
        NumPlotSlices = size(fixed,3);
    end
    %Coronal
    Co_MaskVoxels = squeeze(sum(LungMask,[1,2]));
    Co_StartIndex = find(Co_MaskVoxels,1,'first');
    Co_EndIndex = find(Co_MaskVoxels,1,'last');
    Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
    Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
    Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
    %Axial
    Ax_MaskVoxels = squeeze(sum(LungMask,[2,3]));
    Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
    Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
    Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
    Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
    Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
    fixedMasked = fixed.*LungMask;
    %View
    Registrationfig = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registrationfig,'WindowState','minimized');
    set(Registrationfig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
    tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(abs(moving(:,:,Slices_Co(slice)))/max(abs(moving(:))),abs(fixed(:,:,Slices_Co(slice)))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        if slice == round(NumPlotSlices/2)
            title('Before Registration','FontSize',24)
        end
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(fliplr(rot90(squeeze(abs(moving(Slices_Ax(slice),:,:))),-1))/max(abs(moving(:))),fliplr(rot90(squeeze(abs(fixed(Slices_Ax(slice),:,:))),-1))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(abs(ProtonRegistered(:,:,Slices_Co(slice)))/max(abs(ProtonRegistered(:))),abs(fixed(:,:,Slices_Co(slice)))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        if slice == round(NumPlotSlices/2)
            title('After Registration','FontSize',24)
        end
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(fliplr(rot90(squeeze(abs(ProtonRegistered(Slices_Ax(slice),:,:))),-1))/max(abs(ProtonRegistered(:))),fliplr(rot90(squeeze(abs(fixed(Slices_Ax(slice),:,:))),-1))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
    end
    savefig('Registerationfig.fig')
    close(gcf)
    
    
    Proton.LungMask = LungMask;
    Proton.tform = tform;
    Proton.Slices_Co = Slices_Co; 
    Proton.Slices_Ax = Slices_Ax; 
    Proton.optimizer = optimizer;
    Proton.ProtonRegistered = ProtonRegistered;
    if strcmp(MainInput.AnalysisType,'GasExchange')
        Proton.ProtonHRRegistered = ProtonHRRegistered;
    end
    Proton.ProtonMaskRegistred = LungMask;
    Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
    GasExchange.LungMask = LungMask;
    GasExchange.ProtonRegistered = ProtonRegistered;
    GasExchange.ProtonMaskRegistred = LungMask;
    GasExchange.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);

end
%% save images
%Tiffs (movingRegisteredVolume)

Global.write_imshowpair(ProtonRegistered,fixed(:,:,1:nSlice),DataLocation)
% figure; orthosliceViewer(Proton.ProtonRegisteredColored)

end
