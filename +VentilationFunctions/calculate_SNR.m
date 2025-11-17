function [Ventilation] = calculate_SNR(Ventilation)
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
parentPath = Ventilation.parentPath;
SE = strel('square', 28);
MR = double(Ventilation.Image);
MR = (MR - min(MR(:))) / (max(MR(:)) - min(MR(:)));
maskarray = double(Ventilation.LungMask);

MRvv = double(MR);
NormMR = MRvv.*(maskarray>0);
incomplete = 0.5;
defectmask = VentilationFunctions.medFilter(double(NormMR<(mean(NormMR(NormMR>0))*incomplete)*(maskarray>0)));
Ventilation.SNRdefectmask = defectmask;
% maskarray = double(Ventilation.LungMask + Ventilation.VesselMask);
% maskarray(maskarray > 1) = 0;
% maskarray = double(maskarray);
% figure; imslice(~defectmask)
% Dilate the mask arrary according to the structuring element:
maskarray_dilated = imdilate(maskarray, SE);

% Take the complement of the dilated mask array and create a background
% noise array:
% if ~exist(airwaymask, 'var')
% if sum(airwaymask(:)) == 0
%     airwaymask = Segmentation.SegmentLungthresh(MR,1,1);
% end
% force airway mask
airwaymask = Segmentation.SegmentLungthresh(MR,1,0.3);
% airwaymask = double(MR > 2000);

zerofillings = double(MR == 0);
backgroundmask = double(imcomplement(maskarray_dilated).*~airwaymask.*~zerofillings);
Noise99percentile = prctile(MR(backgroundmask == 1),99.0);

background1 = (MR).*(backgroundmask);
% imslice(backgroundmask)
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
SNRvv_slice = zeros(1, size(MR,3));
signal_n_avg = zeros(size(SNR_slice));
signalvv_avg = zeros(size(SNR_slice));
std_noise = zeros(size(SNR_slice));
mean_noise = zeros(size(SNR_slice));
for n = 1:size(MR,3)
        %Calculating the siganl vector
        mask = maskarray(:,:,n);
        nois_mask = backgroundmask(:,:,n);
        dfctmmask = double(~defectmask(:,:,n));

        signal_vec = MR(:,:,n).*mask;
        signal_vec = reshape(signal_vec,size(MR,1)*size(MR,2),1);
        signal_vec = nonzeros(signal_vec);
        %signal_vec(mask == 0) = [];
        signal_n_avg(n) = mean(signal_vec);
        signal_n_avg(isnan(signal_n_avg))=0;

        signal_vec = MR(:,:,n).*mask.*dfctmmask;
        signal_vec = reshape(signal_vec,size(MR,1)*size(MR,2),1);
        signal_vec = nonzeros(signal_vec);
        %signal_vec(mask == 0) = [];
        signalvv_avg(n) = mean(signal_vec);
        signalvv_avg(isnan(signalvv_avg))=0;

        %Calculating the noise vector
        noise_vec= MR(:,:,n).*nois_mask; 
        noise_vec = reshape(noise_vec,size(MR,1)*size(MR,2),1);
        noise_vec = nonzeros(noise_vec);
        %noise_vec(nois_mask == 0) = [];        
        mean_noise(n) = mean(noise_vec);         %the mean of the noise               
        std_noise(n) = std(noise_vec);       %the standard deviation of the noise
        mean_noise(isnan(mean_noise))=0;
        std_noise(isnan(std_noise))=0;

        SNR_slice(n) = round((signal_n_avg(n) ) / std_noise(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
        SNRvv_slice(n) = round((signalvv_avg(n)) / std_noise(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
        
        % SNR_slice(n) = round((signal_n_avg(n) - mean_noise(n)) / std_noise(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
        % SNRvv_slice(n) = round((signalvv_avg(n) - mean_noise(n)) / std_noise(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
end
SNR_slice(isnan(SNR_slice))=0;
SNR_slice(isinf(SNR_slice)) = 0;
SNR_slice(SNR_slice < 0) = 0;
% SNR_slice = SNR_slice(SNR_slice ~= 0);

SNRvv_slice(isnan(SNRvv_slice))=0;
SNRvv_slice(isinf(SNRvv_slice)) = 0;
SNRvv_slice(SNRvv_slice < 0) = 0;
% SNRvv_slice = SNRvv_slice(SNRvv_slice ~= 0);

% Remove zeros and NaNs before mean/std calculation
vals_signal = MR(maskarray == 1);
vals_signal = vals_signal(vals_signal ~= 0 & ~isnan(vals_signal));

vals_noise = MR(backgroundmask == 1);
vals_noise = vals_noise(vals_noise ~= 0 & ~isnan(vals_noise));

overall_mean = mean(vals_signal, "all");
overall_std  = std(vals_noise, 0, "all");
SNR_lung = round((overall_mean/overall_std), 2)*sqrt(2 - (pi/2));

maskarray2 = maskarray.*~defectmask;
overallvv_mean = mean(MR(maskarray2 == 1), "all", 'omitnan');
SNR_vv = round((overallvv_mean/overall_std), 2)*sqrt(2 - (pi/2));

disp(['Mean SNR for all masked lung slices = ',num2str(SNR_lung)])
disp(['Mean SNR for all masked ventilated slices = ',num2str(SNR_vv)])

cd(parentPath);
% Generate a text file with the SNR of each slice
C = num2cell(SNR_slice);
filePh = fopen('SNRlung.txt','w');
fprintf(filePh,'%s\n',C{:});
fclose(filePh);

C = num2cell(SNRvv_slice);
filePh = fopen('SNRvv.txt','w');
fprintf(filePh,'%s\n',C{:});
fclose(filePh);

Ventilation.SNR_slice = SNR_slice; 
Ventilation.SNRvv_slice = SNRvv_slice;
Ventilation.SNR_lung = SNR_lung; 
Ventilation.SNR_vv = SNR_vv;
Ventilation.NoiseCompleteThresh = round(overall_mean/Noise99percentile);


end