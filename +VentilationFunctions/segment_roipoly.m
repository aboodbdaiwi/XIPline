function [maskarray] = segment_roipoly(MR, MR2)
%% Mask segmentation function using the roipoly tool from MATLAB
% The goal of this function is to segment the images using the roipoly
% function included with MATLAB. This roipoly function is fairly old, and
% MATLAB 2020b released an upgraded version that will operate faster.
%
% Instructions:
% Run: [maskarray] = segment_roipoly(MR, MR2)
%
% Inputs:
% MR = image data for multiple slices.
% MR2 = normalized image data between 0 and 1.
%
% Outputs:
% maskarray = arrays containing the binarized mask for each slice.
%
% Author: Joseph Plummer.
% Date: 05/05/2021.
%
% NOTE: This script utilizes MATLAB's roipoly function to draw masks.
% The same mask can be drawn on ITK snap using the paint brush tool,
% and saving the mask as a 'segmentation' nifti file. All values in
% mask must be set to 0 or 1.
%
% Questions? Contact Joseph Plummer.
fprintf( 'Segment lung parenchyma.\n' );
for p = 1:size(MR2,3)
    scaledImage(:,:,p) = imadjust(MR2(:,:,p));
    scaledImage2(:,:,p) = adapthisteq(scaledImage(:,:,p));
    imshow(scaledImage2(:,:,p),[]);
    HPmask(:,:,p) = roipoly(scaledImage2(:,:,p));
    choice = questdlg('Poly ROI acceptable?', ...
        'Save Poly ROI', ...
        'Yes','No','Yes');
    while ~isequal(choice,'Yes')
        switch choice
            case 'Yes'
                disp(' Outline Next Lung')
                break
            case 'No'
                disp(' Redraw ROI');
                HPmask(:,:,p) = roipoly(scaledImage2(:,:,p));
                break
        end
    end
    
    HPmask2(:,:,p) = roipoly(scaledImage2(:,:,p));
    choice = questdlg('Poly ROI acceptable?', ...
        'Save Poly ROI', ...
        'Yes','No','Yes');
    while ~isequal(choice,'Yes')
        switch choice
            case 'Yes'
                disp(' Outline Next Lung')
                break
            case 'No'
                disp(' Redraw ROI')
                HPmask2(:,:,p) = roipoly(scaledImage2(:,:,p));
                break
        end
    end
    
    mask(:,:,p) = imadd(HPmask(:,:,p),HPmask2(:,:,p));
    maskarray(:,:,p) = Functions.medFilter(mask(:,:,p));
    %imshow(maskarray(:,:,p));
end