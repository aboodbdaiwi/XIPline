
function [Proton] = LoadData_Proton_GasExchange_Philips_Sin(ProtonDataLocation,PixelShift,Institute,Proton)
%   Will process the resulting proton mri data
%   for V3 and XeCTC versions of the data
%   Requires a folder with the Dissolved_Xe raw/lab/sin files,
%   Dissolved_Xe_Mask raw/lab/sin files, and two (and only two!) sets of 
%   data/list files for the dissolved and proton images
%   If no folder is provided, user will be prompted to select the folder
%
%   Syntax:  LoadData_Gas_GasExchange_Philips_Sin(DataLocation)
%
%   Inputs:
%      ProtonDataLocation - string containing the path of the folder
%
%   Outputs:
%      ProtonImage
%
%   Example: 
%   LoadData_Proton_GasExchange_Philips_Sin('C:\Users\mwillmering\Documents\Subject1')
%
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Co-Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% Determine files to use
if strcmp(Institute,'XeCTC') == 1    
    ScanVersion = 'Xe_CTC';
%     HSinFile = dir([ProtonDataLocation,'\*CTC_Anatomical.sin']);
elseif strcmp(Institute,'CCHMC') == 1                 
    ScanVersion = 'CCHMC_V3';
%     HSinFile = dir([ProtonDataLocation,'\*Dissolved_Xe_Mask_20191008.sin']);
end
HSinFile = dir([ProtonDataLocation,'\*.sin']);

% HSinFile = dir([ProtonDataLocation,'\*Dissolved_Xe_Mask_20191008.sin']);
% ScanVersion = 'CCHMC_V3';
% if size(HSinFile,1) < 1 %didn't find that sequence version
%     HSinFile = dir([DataLocation,'\*CTC_Anatomical.sin']);
%     ScanVersion = 'Xe_CTC';
% end
DataFiles = dir([ProtonDataLocation,'\*.data']);
file_name = DataFiles.name;
HDataFile = DataFiles(1);
disp('Looking For Necessary Files Completed.')

cd(ProtonDataLocation)
mkdir([ProtonDataLocation '\Gas Exchange Analysis']);
outputpath = [ProtonDataLocation '\Gas Exchange Analysis'];
cd(outputpath)
Proton.outputfolder = outputpath;
%% Import Acquisition Information
disp('Importing Acquisition Information...')
HSin = GasExchangeFunctions.loadSIN([HSinFile.folder,'\',HSinFile.name]);

H_AcqMatrix = HSin.scan_resolutions.vals(1);
if strcmp(ScanVersion,'Xe_CTC') %recon 2x interpolated
    H_RecMatrix = 2*H_AcqMatrix;
else %CCHMC V3, recon 2x interpolated
    H_RecMatrix = H_AcqMatrix; %already acquired at 2x resolution
end
disp('Importing Acquisition Information Completed.')
%% Import Data
disp('Importing Data...')
%Load in FIDs
[HData,HInfo] = GasExchangeFunctions.loadLISTDATA([HDataFile.folder,'\',HDataFile.name]);
HData = squeeze(HData);

if strcmp(ScanVersion,'Xe_CTC')
    %Proton k-space
    HKSpace = HData;

    %Get sizes of dimensions
    H_nsamp = size(HKSpace,1);
    H_nprof = size(HKSpace,2);
    H_coils = size(HKSpace,3);
    H_interleaves = 1;
else
    %Proton k-space
    HKSpaceInit = HData;

    %Get sizes of dimensions
    H_nsamp = size(HKSpaceInit,1);
    H_nprof = size(HKSpaceInit,2);
    H_interleaves = size(HKSpaceInit,3);
    H_coils = size(HKSpaceInit,4);
    
    % Put KSpaces in order and logical format
    HOrder = GasExchangeFunctions.GetProjectionOrder(HInfo);
    HKSpace = zeros(H_nsamp,H_nprof*H_interleaves,H_coils);
    for i=1:H_coils
        HKSpace(:,:,i) = GasExchangeFunctions.SortUTEData_AcqOrder(GasExchangeFunctions.reshapeUTEData(HKSpaceInit(:,:,:,i)),[H_nsamp, H_nprof, H_interleaves], HOrder);
    end
end
disp('Importing Data Completed.')

%% Calculate Trajectories
disp('Calculating Trajectories...')
if strcmp(ScanVersion,'Xe_CTC')
    HTraj = GasExchangeFunctions.philipsradialcoords(1.25,2,[HSinFile.folder,'\',HSinFile.name]); %1.25us delay, Haltoned Spiral
    HTraj = permute(HTraj,[4 3 2 1]);
else
    HTraj = GasExchangeFunctions.philipsradialcoords(1.65,1,[HSinFile.folder,'\',HSinFile.name]); %1.65us delay, GM
    HTraj = permute(HTraj,[4 3 2 1]);
    HTraj = GasExchangeFunctions.SortUTETraj_AcqOrder(GasExchangeFunctions.reshapeUTETraj(HTraj), [H_nsamp, H_nprof, H_interleaves], HOrder);%reshape and order
end
disp('Calculating Trajectories Completed.')

%% Account for Exhale/Motion in Proton Image
%Functions to fit
for coil = 1:size(HKSpace,3)
    HPhase(:,coil) = smooth(angle(HKSpace(1,:,coil)),15,'rloess');
end
HPhaseDynamics = rad2deg(mean(unwrap(HPhase,1),2));
HPulses = (1:length(HPhaseDynamics))';

%Fits to data; use only 1st third to avoid exhale/movement - assuming first 5 seconds are uncorrupted
HDynMean = mean(HPhaseDynamics(1000:round(length(HPulses)/3)));
HDynFit = HDynMean*ones(size(HPhaseDynamics(1000:round(length(HPulses)/3))));
HDynRMSE = sqrt(immse(HDynFit,HPhaseDynamics(1000:round(length(HPulses)/3))));
if (3*HDynRMSE<5) 
    HDynRMSE = 5/3; %If small phase varaition, use acceptance window of 5deg
end

%smooth to remove small excursions
HPhaseDynamics_Smoothed = medfilt1(HPhaseDynamics,15,'truncate');

%Determine Proton Data Cutoffs
HDynFit = HDynMean*ones(size(HPulses));%calculate function at all x points
keep = abs(HDynFit - HPhaseDynamics_Smoothed) < 3*HDynRMSE;%label as outlier if >3 rmse away
keep = medfilt1(single(keep),15,'truncate');%smooth to remove small excursions
keep(1:round(length(HPulses)/3)) = 1; %keep all from fit
motion_ind = find(keep==0,1,'first');%find first non keep index

%H Motion
HKSpace_Steady = HKSpace;
HTraj_Steady = HTraj;
if ~isempty(motion_ind)%only remove data if motion
    HKSpace_Steady(:,motion_ind:end,:) = [];
    HTraj_Steady(:,:,motion_ind:end) = [];
else
    motion_ind = length(HPulses);
end

%% Determine FOV Shift
% %not using planned shifts as techs don't seem to do it consistently
% disp('Determining Image offset...')
% PixelShiftinit = [0; 0; 0]; %start with no offset
% converged = 0; %start not converged
% iter = 0;
% while (converged == 0) %run until no change in shift
%     
%     tempIm = GasExchange.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,1),HTraj_Steady/2,PixelShiftinit); %make image; 2x Resolution
% %     tempIm = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace,XeTraj/2,PixelShiftinit); %make image; 2x Resolution
%     [tempIm, ~] = GasExchange.Vent_BiasCorrection_Coarse(tempIm);
%     tempLab = bwlabeln(imbinarize(abs(tempIm)/(0.5*max(abs(tempIm(:)))))); %make mask label
%     center = regionprops3(tempLab); %get centroids
%     center = sortrows(center,1,'descend'); %order by size
%     if (size(center,1) > 1 && table2array(center(2,1)) > 0.15 * table2array(center(1,1)))        
%         CentVox(1) = (min(center.BoundingBox(:,1))+max(center.BoundingBox(:,1)+center.BoundingBox(:,4)))/2;
%         CentVox(2) = (min(center.BoundingBox(:,2))+max(center.BoundingBox(:,2)+center.BoundingBox(:,5)))/2;
%         CentVox(3) = (min(center.BoundingBox(:,3))+max(center.BoundingBox(:,3)+center.BoundingBox(:,6)))/2;
%     else
%         CentVox(1) = center.BoundingBox(1,1)+center.BoundingBox(1,4)/2;
%         CentVox(2) = center.BoundingBox(1,2)+center.BoundingBox(1,5)/2;
%         CentVox(3) = center.BoundingBox(1,3)+center.BoundingBox(1,6)/2;
%     end
%     PixelShift = 0.5*H_RecMatrix - CentVox'; %get pixel shift
%     PixelShift = PixelShift([3, 2, 1]); %put in same order
%     PixelShift(2:3) = -PixelShift(2:3); %correct for reference
%     PixelShift = PixelShift + PixelShiftinit; %sum previous shifts
%     iter = iter + 1;
%     if (sum(abs(PixelShift - PixelShiftinit) > 1) == 0) %if all offsets changed less than 1 pixel, end
%         converged = 1;
%     elseif (iter == 5) %if reached max iterations
%         converged = 1;
%         if (sum(abs(PixelShift - PixelShiftinit) < 3) ~= 3) %if all offsets changed less than 3 pixels, keep, else set to 0
%             PixelShift = [0; 0; 0]; %move back to no offset
%         end
%         disp('Did not converge on image offset...')
%     else %continue loop
%         PixelShiftinit = PixelShift; %set init for next loop
%     end
% end
% disp('Determining Image Offset Completed.')

% PixelShift = [8.5; 1.5; 2]; %start with no offset
%% Proton Image
%Reconstruct Image
disp('Reconstructing UTE Image...')
ProtonImage = zeros(H_RecMatrix*3,H_RecMatrix*3,H_RecMatrix*3,H_coils);%3x overgridding
for i = 1:H_coils
    if strcmp(ScanVersion,'Xe_CTC') %recon 2x interpolated
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon_CTC(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady/2,PixelShift);
    else
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady,PixelShift);
    end
end
ProtonImage = GasExchangeFunctions.SOS_Coil_Combine(ProtonImage);%Combine Coils
disp('Reconstructing UTE Image Competed.')
ProtonImageHR = ProtonImage;
%Correct Bias
disp('Correcting Proton Bias...')
ProtonImage = medfilt3(ProtonImage, [7 7 7]);%remove undersampling artifacts
ProtonImage = imgaussfilt3(ProtonImage, 0.5);%reduce noise
[ProtonImage, ~] = GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImage);

% Determine Levels
ProtonMax = prctile(abs(ProtonImage(:)),99.99);

disp('Proton Bias Corrected.')  
cd(outputpath)
%View
ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
montage(abs(ProtonImageHR(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix)),'DisplayRange',[0 ProtonMax])%unregistered for these
savefig('ProtonMontage.fig')
close(gcf)

Proton.Image = double(ProtonImage);
Proton.ProtonImageHR = double(ProtonImageHR);
Proton.filename = file_name;
Proton.folder = ProtonDataLocation;
Proton.H_RecMatrix = H_RecMatrix;
Proton.ProtonMax = ProtonMax;

end
