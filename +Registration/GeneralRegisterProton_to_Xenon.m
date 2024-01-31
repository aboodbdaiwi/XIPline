
function [Proton,MainInput] = GeneralRegisterProton_to_Xenon(...
    Proton,...
    XeImage,...
    MainInput)
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
ProtonVoxelInfo = MainInput.ProtonVoxelInfo;
XeVoxelInfo = MainInput.XeVoxelInfo;
RegistrationType = MainInput.RegistrationType;
DataLocation = MainInput.XeDataLocation;


% if strcmp(MainInput.AnalysisType,'Ventilation') == 1 
%     mkdir([MainInput.XeDataLocation '\Ventilation Analysis']);
%     DataLocation = [MainInput.XeDataLocation,'\Ventilation Analysis'];
%     MainInput.XeDataLocation = DataLocation;
% elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
%     mkdir([MainInput.XeDataLocation '\GasExchange Analysis']);
%     DataLocation = [MainInput.XeDataLocation,'\GasExchange Analysis'];
%     MainInput.XeDataLocation = DataLocation;
% end

moving1 = double(Proton.Image);
fixed1 = XeImage.Image;
nSlice = size(fixed1,3);
%% 

% match number of slices 
if MainInput.SliceSelection == 1
    moving1 = squeeze(moving1(:,:,MainInput.Hstart:MainInput.Hend));
    fixed1 = squeeze(fixed1(:,:,MainInput.Xestart:MainInput.Xeend));
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

% check image dimenssion if less than 16
if size(moving1,3) < 16  
    moving = zeros(size(moving1,1),size(moving1,2),16);
    moving (:,:,1:size(moving1,3)) = abs(moving1);
else
    moving = moving1;
end
if size(fixed1,3) < 16
    fixed = zeros(size(fixed1,1),size(fixed1,2),16);
    fixed (:,:,1:size(fixed1,3)) = abs(fixed1); 
else
    fixed = fixed1;
end
% figure; imslice(moving) 
% figure; imslice(fixed) 
%% Register Multimodal 3-D Medical Images.

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
% geomtform.T
% centerXWorld = mean(Rmoving.XWorldLimits);
% centerYWorld = mean(Rmoving.YWorldLimits);
% centerZWorld = mean(Rmoving.ZWorldLimits);
% [xWorld,yWorld,zWorld] = transformPointsForward(geomtform,centerXWorld,centerYWorld,centerZWorld);
% [r,c,p] = worldToSubscript(Rfixed,xWorld,yWorld,zWorld);
ProtonRegistered = imwarp(movingVolume,Rmoving,geomtform,'bicubic','OutputView',Rfixed);

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
%% view
% figure
% x = imshowpair(movingRegisteredVolume(:,:,centerFixed(3)), fixedVolume(:,:,centerFixed(3)));
% title('Axial Slice of Registered Volume')
% disp('done')

%% save images
%Tiffs (movingRegisteredVolume)

Global.write_imshowpair(ProtonRegistered,XeImage.Image,DataLocation)
% figure; orthosliceViewer(Proton.ProtonRegisteredColored)

end
