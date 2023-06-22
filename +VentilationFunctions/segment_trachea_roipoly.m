function [maskarraytrachea] = segment_trachea_roipoly(MR, maskarray)
%% Trachea segmentation function using the roipoly tool from MATLAB
% The goal of this function is to segment the trachea using the roipoly
% function included with MATLAB. This roipoly function is fairly old, and
% MATLAB 2020b released an upgraded version that will operate faster.
%
% Instructions:
% Run: [maskarray] = segment_roipoly(MR, MR2, foldername)
%
% Inputs:
% MR = image data for multiple slices.
% maskarray = binarized array of masks for each slice.
%
% Outputs:
% maskarraytrachea = arrays containing the binarized trachea mask for each slice.
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
%% 1) Normalize the image data to make it easier to segment:
fprintf( 'Segment trachea.\n' );
figure('Name','Segment trachea')
% Normalize the intensity of the original image to fall between [0,1].
MR2 = MR - min(MR(:));
MR2 = MR / max(MR(:));

%% 2) Segment trachea to use for the linear binning algorithms:
% Note: This entire section can be edited to the new 2020b MATLAB
% segmentation, as opposed to this 'older' roipoly method. 
for p = 1:size(MR2,3)
    scaledImagetrachea(:,:,p) = imadjust(MR2(:,:,p));
    scaledImage2trachea(:,:,p) = adapthisteq(scaledImagetrachea(:,:,p));
    imshow(scaledImage2trachea(:,:,p),[]);
    HPmasktrachea(:,:,p) = roipoly(scaledImage2trachea(:,:,p));
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
                              HPmasktrachea(:,:,p) = roipoly(scaledImage2trachea(:,:,p));
                             break
                      end
                  end
    masktrachea(:,:,p) = (HPmasktrachea(:,:,p));
    
    % Apply median filter to the masked trachea:
    masktrachea2(:,:,p) = VentilationFunctions.medFilter(masktrachea(:,:,p));
    imshow(masktrachea2(:,:,p));
end

% Combine trachea mask to maskarray for lungs:
maskarraytrachea = imadd(maskarray,masktrachea2);

end