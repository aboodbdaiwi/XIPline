function [SNR_slice, Overall_SNR] = calculate_SNR(MR, maskarray, airwaymask)
%% MATLAB script to perform SNR calculation
% This code uses the ventilation images (can be N4 corrected or
% original), and the corresponding masks, to generate an array of SNR
% values for each slice. Advanced features of this code can be edited, such
% as the structural element size for the mask dilations.
%
% Instructions:
% Run: [SNR] = Functions.calculate_SNR(MR, maskarray, parentPath)
%
% Inputs:
% MR = ventilation image data
% maskarray = corresponding masks for MR
%
% Outputs:
% SNR = array of SNR values for each image slice
% Overall_SNR = mean SNR value for entire masked imageset
%
%
% Author: Joseph Plummer.
% Date: 05/09/2021.

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% Calculate SNR
% Create a structuring element to dilate the Xe mask array.
% The structuring element size '28' can be changed as needed to remove
% partial volume pixels and the trachea.
SE = strel('square', 28);
maskarray = double(maskarray);
MR = double(MR);
% Dilate the mask arrary according to the structuring element:
maskarray_dilated = imdilate(maskarray, SE);

% Take the complement of the dilated mask array and create a background
% noise array:
% if ~exist(airwaymask, 'var')
if sum(airwaymask(:)) == 0
    airwaymask = Segmentation.SegmentLungthresh(MR,1,1);
end
backgroundmask = double(imcomplement(maskarray_dilated).*~airwaymask);
background1 = (MR).*(backgroundmask);

% % Label 0s in the background noise array as NaNs
% background1(background1 == 0) = NaN;
% 
% % Label 0s in the Xe vent array as NaNs
% MR2 = MR.*(maskarray>0);
% MR2(MR2 == 0) = NaN;


% Calculate the column size needed to reshape the background noise
% array for accurate StdDev calculation:
%x = size(background1,1) * size(background1,2);

% Calculate the SNR for each slice:
SNR_slice = zeros(1, size(MR,3));
signal_n_avg = zeros(size(SNR_slice));
std_noise = zeros(size(SNR_slice));
mean_noise = zeros(size(SNR_slice));
for n = 1:size(MR,3)
        %Calculating the siganl vector
        mask = maskarray(:,:,n);
        nois_mask = backgroundmask(:,:,n);

        signal_vec = MR(:,:,n).*mask;
        signal_vec = reshape(signal_vec,size(MR,1)*size(MR,2),1);
        signal_vec = nonzeros(signal_vec);
        %signal_vec(mask == 0) = [];
        signal_n_avg(n) = mean(signal_vec);
        signal_n_avg(isnan(signal_n_avg))=0;

        %Calculating the noise vector
        noise_vec= MR(:,:,n).*nois_mask; 
        noise_vec = reshape(noise_vec,size(MR,1)*size(MR,2),1);
        noise_vec = nonzeros(noise_vec);
        %noise_vec(nois_mask == 0) = [];        
        mean_noise(n) = mean(noise_vec);         %the mean of the noise               
        std_noise(n) = std(noise_vec);       %the standard deviation of the noise
        mean_noise(isnan(mean_noise))=0;
        std_noise(isnan(std_noise))=0;

        SNR_slice(n) = round(signal_n_avg(n) / std_noise(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
end
SNR_slice(isnan(SNR_slice))=0;
SNR_slice(isinf(SNR_slice)) = 0;
SNR_slice = SNR_slice(SNR_slice ~= 0);
Overall_SNR = mean(SNR_slice);
disp(['Mean SNR for all masked ventilation slices = ',num2str(Overall_SNR)])

% Generate a text file with the SNR of each slice
C = num2cell(SNR_slice);
filePh = fopen('SNR.txt','w');
fprintf(filePh,'%s\n',C{:});
fclose(filePh);

end