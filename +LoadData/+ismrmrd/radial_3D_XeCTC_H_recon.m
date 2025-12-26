
function [Proton] = radial_3D_XeCTC_H_recon(MainInput, GasExchange, Proton)
%% Working with an existing ISMRMRD data set
   
%   Co-Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Institute = MainInput.Institute;
% Institute = 'XeCTC';
% MainInput.HFileName = '20220414_100621_DukeIPF_Mask_Gas_Exchange.h5';
% MainInput.XeDataLocation = 'D:\OneDrive - cchmc\Lab\Xe_App\testtingData\Gasexchange\20220414 ILD-HC-001\ismrmrd';
cd(MainInput.HDataLocation)
if exist(MainInput.HFileName, 'file')
    dset = LoadData.ismrmrd.Dataset(MainInput.HFileName, 'dataset');
else
    error(['File ' MainInput.XeFileName ' does not exist.  Please generate it.'])
end
ProtonDataLocation = MainInput.HDataLocation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr = LoadData.ismrmrd.xml.deserialize(dset.readxml);
%% Read all the data
clc
D = dset.readAcquisition();
%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
if strcmp(MainInput.Scanner,'Siemens')
    enc_Ny = size(D.data,2);
end
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;


%% Determine files to use
if strcmp(Institute,'XeCTC') == 1    
    ScanVersion = 'XeCTC';
%     HSinFile = dir([ProtonDataLocation,'\*CTC_Anatomical.sin']);
elseif strcmp(Institute,'CCHMC') == 1                 
    ScanVersion = 'CCHMC_V3';
%     HSinFile = dir([ProtonDataLocation,'\*Dissolved_Xe_Mask_20191008.sin']);
end

cd(ProtonDataLocation)
if ~isfield(GasExchange, 'outputpath') || isempty(GasExchange.outputpath)
    mkdir([GasDataLocation '\GasExchange_Analysis']);
    outputpath = [GasDataLocation '\GasExchange_Analysis'];
    GasExchange.outputpath = outputpath;
else
    outputpath = GasExchange.outputpath;
end
Proton.outputfolder = outputpath;
%% Import Acquisition Information

H_AcqMatrix = enc_Nx;
H_RecMatrix = 2*H_AcqMatrix;

%% Import Data
disp('Importing Data...')
%Load in FIDs

cellSize = cell2mat(D.data(1,1));
HData = zeros(size(cellSize,1),enc_Ny,size(cellSize,2));
for i = 1:enc_Ny
    HData(:,i,:) = cell2mat(D.data(1,i));
end

%Proton k-space
HKSpace = HData;

%Get sizes of dimensions
H_nsamp = size(HKSpace,1);
H_nprof = size(HKSpace,2);
H_coils = size(HKSpace,3);
H_interleaves = 1;
    
disp('Importing Data Completed.')

%% Calculate Trajectories
disp('Calculating Trajectories...')

cellSize = cell2mat(D.traj(1,1));
HTraj = zeros(size(cellSize,1),size(cellSize,2),enc_Ny);
for i = 1:enc_Ny
    HTraj(:,:,i) = cell2mat(D.traj(1,i));
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
try
    PixelShift = GasExchange.PixelShift; % [1.5;4;-4];
catch
    PixelShift = [0,0,0];
end
disp('Reconstructing UTE Image...')
ProtonImage = zeros(H_RecMatrix*3,H_RecMatrix*3,H_RecMatrix*3,H_coils);%3x overgridding
for i = 1:H_coils
    if strcmp(ScanVersion,'XeCTC') %recon 2x interpolated
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon_CTC(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady/2,PixelShift);
    else
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady,PixelShift);
    end
end
ProtonImage = GasExchangeFunctions.SOS_Coil_Combine(ProtonImage);%Combine Coils
disp('Reconstructing UTE Image Competed.')
if strcmp(MainInput.Scanner,'Siemens')
    ProtonImage = permute(ProtonImage, [2,1,3]);
    ProtonImage = flip(ProtonImage,1);
    ProtonImage = flip(ProtonImage,2);
end
ProtonImageHR = ProtonImage;

%Correct Bias
disp('Correcting Proton Bias...')
ProtonImageLR = medfilt3(ProtonImage, [7 7 7]);%remove undersampling artifacts
ProtonImageLR = imgaussfilt3(ProtonImageLR, 0.5);%reduce noise
[ProtonImageLR, ~] = GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageLR);
[ProtonImageHR, ~] = GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageHR);

% Determine Levels
ProtonMax = prctile(abs(ProtonImageLR(:)),99.99);
figure; imslice(ProtonImage)

disp('Proton Bias Corrected.')  
cd(outputpath)
%View
ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
montage(abs(ProtonImage(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix)),'DisplayRange',[0 ProtonMax])%unregistered for these
savefig('ProtonMontage.fig')
close(gcf)

Proton.Image = double(ProtonImageHR);
Proton.ProtonImageLR = double(ProtonImageLR);
Proton.ProtonImageHR = double(ProtonImageHR);
Proton.filename = MainInput.HFileName;
Proton.folder = ProtonDataLocation;
Proton.H_RecMatrix = H_RecMatrix;
Proton.ProtonMax = ProtonMax;

end