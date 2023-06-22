
function [movingRegisteredVolume,ProtonMaskRegistered] = RegisterProton_to_Xenon(ProtonImage,XeImage,ProtonVoxelInfo,XeVoxelInfo)
%   Inputs:

%          
%   Outputs:
%      LungMask
%
%   Example: 

%
% 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/
%% intilize 
ProtonMaskRegistered = [];

%% check image dimenssion if less than 16 
moving = ProtonImage;
fixed = XeImage;
if size(ProtonImage,3) < 16  
    moving = zeros(size(ProtonImage,1),size(ProtonImage,2),16);
    moving (:,:,1:14) = abs(ProtonImage);
end
if size(XeImage,3) < 16
    fixed = zeros(size(XeImage,1),size(XeImage,2),16);
    fixed (:,:,1:14) = abs(XeImage); 
end
%% Register Multimodal 3-D Medical Images.

movingVolume = moving;
fixedVolume = fixed;

% helperVolumeRegistration(fixedVolume,movingVolume);

centerFixed = size(fixedVolume)/2;
centerMoving = size(movingVolume)/2;
% figure
imshowpair(movingVolume(:,:,centerMoving(3)), fixedVolume(:,:,centerFixed(3)));
% title('Unregistered Axial Slice')

[optimizer,metric] = imregconfig('multimodal');
Rfixed  = imref3d(size(fixedVolume),XeVoxelInfo.PixelSize2,XeVoxelInfo.PixelSize1,XeVoxelInfo.SliceThickness);
Rmoving = imref3d(size(movingVolume),ProtonVoxelInfo.PixelSize2,ProtonVoxelInfo.PixelSize1,ProtonVoxelInfo.SliceThickness);

% Rmoving.XWorldLimits
% Rmoving.PixelExtentInWorldX
% Rmoving.ImageExtentInWorldX
optimizer.InitialRadius = 0.00005;
movingRegisteredVolume = imregister(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
% figure
% imshowpair(movingRegisteredVolume(:,:,centerFixed(3)), fixedVolume(:,:,centerFixed(3)));
% title('Axial Slice of Registered Volume')
helperVolumeRegistration(fixedVolume,movingRegisteredVolume);

geomtform = imregtform(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
% geomtform.T

% centerXWorld = mean(Rmoving.XWorldLimits);
% centerYWorld = mean(Rmoving.YWorldLimits);
% centerZWorld = mean(Rmoving.ZWorldLimits);
% [xWorld,yWorld,zWorld] = transformPointsForward(geomtform,centerXWorld,centerYWorld,centerZWorld);
% [r,c,p] = worldToSubscript(Rfixed,xWorld,yWorld,zWorld);

movingRegisteredVolume = imwarp(movingVolume,Rmoving,geomtform,'bicubic','OutputView',Rfixed);

%view
figure
x = imshowpair(movingRegisteredVolume(:,:,centerFixed(3)), fixedVolume(:,:,centerFixed(3)));
title('Axial Slice of Registered Volume')
disp('done')

end
