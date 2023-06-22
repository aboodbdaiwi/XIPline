function [LungMask,ProtonSignal,LungSignal,LungStd,NoiseMask] = GasExchange_MakeProtonMask(ProtonImage,NewMask)
%   Syntax:  [LungMask] = GasExchange_MakeProtonMask(ProtonImage,NewMask)
%
%   Inputs:
%      ProtonImage - 3D array
%      NewMask -  1 forces new masks to be made, 0 imports masks if present
%   Outputs:
%      LungMask
%
%   Example: 
%   LoadData_Proton_GasExchange_Philips_Sin('C:\Users\mwillmering\Documents\Subject1')
%
% 
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/
%% Make Proton Mask
if(NewMask==1)
    disp('Making Lung Mask...')
    %Initial mask
    InitMask = ~imbinarize(ProtonImage/(prctile(abs(ProtonImage(:)),99.9)));%initial otsu thresholding; select low intensity
    ProtonSignal = mean(ProtonImage(~InitMask));%determine mean proton signal for thresholding
    %Make mask of noise
    NoiseMask = InitMask;
    NoiseMask(ProtonImage>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung
    temp = imclearborder(NoiseMask);%remove noise
    NoiseMask = NoiseMask - temp;%leave only noise
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = imdilate(NoiseMask,strel('sphere',7));%expand to avoid edges
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = logical(NoiseMask);
    %Initial Lung Mask
    LungInitMask = InitMask - NoiseMask;
    LungInitMask(LungInitMask<0) = 0;
    LungInitMask = logical(LungInitMask);
    LungStd =  std(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
    LungInitMask = imerode(LungInitMask,strel('sphere',5));%only for-certain lung
    LungSignal =  mean(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
    LungInitMask = zeros(size(LungInitMask));
    LungInitMask(ProtonImage<(LungSignal+(ProtonSignal-LungSignal)/2)) = 1;
    LungInitMask = LungInitMask - NoiseMask;
    LungInitMask(LungInitMask<0) = 0;
    Hregions = regionprops3(logical(LungInitMask), 'all'); %get information on remaining regions
    Hregions = sortrows(Hregions,1,'descend'); %order by size
    LungInitMask = zeros(size(LungInitMask)); %start final mask
    for ii=1:size(Hregions,1)
        if (Hregions.Volume(ii) > 0.33*Hregions.Volume(1)) %Select only largest regions
            LungInitMask(Hregions.VoxelIdxList{ii,1}) = 1;
        else
            break
        end
    end
    %Final touches of Lung Mask
    LungMask = LungInitMask;
    LungMask(ProtonImage>(LungSignal+5*LungStd)) = 0;%remove clear proton
    Hregions = regionprops3(logical(LungInitMask), 'all'); %get information on remaining regions
    Hregions = sortrows(Hregions,1,'descend'); %order by size
    LungMask = zeros(size(LungMask)); %start final mask
    for ii=1:size(Hregions,1)
        if (Hregions.Volume(ii) > 0.33*Hregions.Volume(1)) %Select only largest regions
            LungMask(Hregions.VoxelIdxList{ii,1}) = 1;
        else
            break
        end
    end
    LungMask = imerode(logical(LungMask),strel('sphere',2));%avoid edges
    LungMask = imclose(logical(LungMask),strel('sphere',4));%avoid edges
else
    disp('Importing Lung Mask...')
    %Calculate Proton and Lung signals
    InitMask = ~imbinarize(ProtonImage/(prctile(abs(ProtonImage(:)),99.9)));%initial otsu thresholding; select low intensity
    ProtonSignal = mean(ProtonImage(~InitMask));%determine mean proton signal for thresholding
    NoiseMask = InitMask;
    NoiseMask(ProtonImage>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung
    temp = imclearborder(NoiseMask);%remove noise
    NoiseMask = NoiseMask - temp;%leave only noise
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = imdilate(NoiseMask,strel('sphere',7));%expand to avoid edges
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = logical(NoiseMask);
    LungInitMask = InitMask - NoiseMask;
    LungInitMask(LungInitMask<0) = 0;
    LungInitMask = logical(LungInitMask);
    LungStd =  std(ProtonImage(LungInitMask));
    LungInitMask = imerode(LungInitMask,strel('sphere',5));%only for-certain lung
    LungSignal =  mean(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
    %Import Mask
    LungMask = logical(fliplr(rot90(niftiread([DataLocation,'\ProtonMask.nii.gz']),-1)));
end

if LungMask == zeros(size(LungMask))
    %if no mask, create mask of ones
    success = 0;%mark as unsuccessful
    LungMask = logical(ones(size(LungMask)));
end



end