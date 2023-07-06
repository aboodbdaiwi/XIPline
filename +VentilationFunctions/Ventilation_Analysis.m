
function [Ventilation] = Ventilation_Analysis (Ventilation,Proton,MainInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Matlab script for processing VDP Analysis  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% 1) Apply the settings:
% We shall use a custom 'input_params' function to define the initial
% settings:
close all;
f = waitbar(0,'Processing ventilation analysis...');
pause(.1)

settings = VentilationFunctions.input_params(...
1,... % Median filter
0,... % RF correction
1,... % Save data
1,... % Calculate SNR
Ventilation.N4Analysis,... % N4 analysis
Ventilation.IncompleteThresh,... % Incomplete threshold
Ventilation.CompleteThresh,... % Complete threshold
Ventilation.HyperventilatedThresh); % Hyperventilated threshold
% Usage: settings = Functions.input_params(medfilter, RFcorrection,
% savedata, calculateSNR, N4, incomplete, complete, hyper)

MR = Ventilation.Image;
mkdir([MainInput.XeDataLocation '\Ventilation Analysis']);
parentPath = [MainInput.XeDataLocation '\Ventilation Analysis\'];
%parentPath = [MainInput.XeDataLocation,'\'];
FileNames = Ventilation.filename;
maskarray = double(Ventilation.LungMask);
try
    airwaymask = double(Ventilation.AirwayMask);
catch
    airwaymask = zeros(size(Ventilation.LungMask));
end
ventmean = Ventilation.LB_RefMean; % Defined by data collected up to Jan 2021.
ventstd =  Ventilation.LB_RefSD; % Defined by data collected up to Jan 2021.


%% 5) Apply RF correction to images:
% Generally, we choose not to perform RF correction on images if we run N4
% bias correction later on.
switch settings.RF_correction
    case "yes"
        disp('RF correction performed on images.')
        MR = VentilationFunctions.rfCorrection(MR,maskarray);
    case "no"
        disp('RF correction has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end

%% 6) Apply N4 Bias Correction:
% NOTE: There are opportunities to make this function more robust, such as
% making sure it is functional for 2D images as well as 3D (it currently is
% optimized for 3D imagesets). Similarly, we could also provide it with
% additional settings inputs, but I chose to keep it simple and only
% provide the optimal settings as chosen when I was the research assistant.
% -- written by Joey Plummer, 05/10/2021

switch settings.N4_bias_analysis
    case "yes"    
        waitbar(.10,f,'Performing N4 bias analysis...');
        pause(.1)
        disp('Running N4 bias field correction...')
        [N4, Bias] = VentilationFunctions.N4_bias_correction(MR, maskarray, parentPath);
        Ventilation.Bias = Bias;       
        % Output Nifti file of N4 corrected images to save for later use.
        [~,Name,ext] = fileparts(FileNames);
        if iscell(ext) == 0 % Deal with multiple dicom vs single dicom/nifti upload.
            ext = ext;
            Name = Name;
        else
            ext = ext{1};
            Name = Name{1};
        end
        
        if ext == ".gz"
            [~,Name,ext] = fileparts(Name); % Deal with double extension (.nii.gz)
        else
        end
        
        % Apply rotations and flip to account for niftiwrite rotation.
        N4_2 = imrotate(N4,90); 
        N4_2 = flipdim(N4_2,1);
        
        % Save images as Nifti.
        NameN4 = Name + "N4";
%         niftiwrite(abs(N4_2),[parentPath,char(NameN4)]);
        niftiwrite(abs(N4_2),[parentPath + "Ventilation_ImagesN4"]); % Or do this until we have a naming convention
        MR = N4;
    case "no"
        disp('N4 bias field correction has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end

%% 7) Calculate SNR:
switch settings.calculate_SNR
    case "yes"
        waitbar(.20,f,'Calculating SNR...');
        pause(.1)        
        disp('Calculating SNR...')
        [SNR_slice, Overall_SNR] = VentilationFunctions.calculate_SNR(MR, maskarray);
        Ventilation.SNR_slice = SNR_slice;
        Ventilation.Overall_SNR = Overall_SNR;
    case "no"
        disp('SNR calculation has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end
close all;
%% Proton image
if strcmp(MainInput.NoProtonImage,'yes') == 1 
        Proton.ProtonRegistered = zeros(size(MR));
        Proton.ProtonRegisteredColored = zeros(size(MR));
end

%% 8) Calculate VDP:
% Note: VDP can be calculated with or without the median filter performed
% on the mask regions. This dependency will be contained within the
% calculate_VDP function.
if strcmp(Ventilation.ThreshAnalysis,'yes') == 1   % 'yes'; || 'no'
    f = waitbar(.20,'Calculating SNR...');
    waitbar(.50,f,'Performing Thershold Analysis ...');
    pause(.1)    
    medfilter = settings.Median_Filter;
    complete = settings.Complete_threshold;
    incomplete = settings.Incomplete_threshold;
    hyper = settings.Hyper_threshold;
    N4_bias_analysis = settings.N4_bias_analysis;

    % Call the calculate_VDP() function. Make sure you call the correct data
    % depending on N4 bias correction (N4) vs original (MR):
%     if settings.N4_bias_analysis == "yes"
%          [VDP,VentscaledImage,DefectArray,VDP_hist,VentDefectmap] = VentilationFunctions.calculate_VDP(N4, maskarray, complete, incomplete, hyper, medfilter, N4_bias_analysis, parentPath,Overall_SNR,Proton,MainInput);
%     elseif settings.N4_bias_analysis == "no"
%          [VDP,VentscaledImage,DefectArray,VDP_hist,VentDefectmap] = VentilationFunctions.calculate_VDP(MR, maskarray, complete, incomplete, hyper, medfilter, N4_bias_analysis, parentPath,Overall_SNR,Proton,MainInput);
%     else
%         disp('Neither original or N4 bias corrected images were declared. Declare these in the input_params() function.')
%      end

    if strcmp(MainInput.Institute,'XeCTC') == 1 || strcmp(MainInput.Institute,'CCHMC') == 1
        [Ventilation] = VentilationFunctions.calculate_VDP_CCHMC(MR, maskarray, complete, incomplete, hyper, medfilter, N4_bias_analysis, parentPath,Overall_SNR,Ventilation,Proton,MainInput);    
    else    
        [Ventilation] = VentilationFunctions.calculate_VDP(MR, maskarray, complete, incomplete, hyper, medfilter, N4_bias_analysis, parentPath,Overall_SNR,Ventilation,Proton,MainInput);    
    end
    % Display the final VDP:
    disp(['The overall ventilation defect percentage is: ',num2str(Ventilation.VDP),'%'])
end 
close all;
%%  Perform linear binning ventilation analysis:
if strcmp(Ventilation.LB_Analysis,'yes') == 1 
    f = waitbar(.50,'Performing Thershold Analysis ...');
    waitbar(.75,f,'Performing Linear Binning Analysis ...');
    pause(.1)     
    % Note: This section is highly dependent on whether N4 bias correction is
    % applied. N4 bias correction is applied in this case INSIDE the
    % calculate_LB_VDP() function. It uses the same optimization parameters as
    % the N4 correction function used earlier in this script.
    
    % 9) Segment the trachea for linear binning analysis:
%     maskarraytrachea2 = imrotate(maskarraytrachea,90); % niftiwrite likes to rotate and flip, account for this
%     maskarraytrachea2 = flipdim(maskarraytrachea2,1);

    maskarraytrachea = maskarray;
    maskarraytrachea(airwaymask == 1)=0;    
    maskarraytrachea2 = maskarraytrachea;
    
    niftiwrite(maskarraytrachea2,[parentPath + "Lung_and_Trachea_mask"]);
    Ventilation.maskarraytrachea = maskarraytrachea;

    % 10) Perform linear binning ventilation analysis:
    N4_bias_analysis = settings.N4_bias_analysis;

%     if settings.N4_bias_analysis == "yes"
% %         ventmean = 0.6786; % Defined by data collected up to Jan 2021.
% %         ventstd =  0.1395; % Defined by data collected up to Jan 2021.
%         % Important! We don't use the N4 corrected images that we had earlier.
%         % We perform N4 correction using the new trachea + lung mask inside
%         % the calculate_LB_VDP() function. The N4 function uses the same
%         % optimization strategy as before.
%         [LB_VDP, LB_mean, LB_std,VentscaledImage,Ventcolormap,LB_Venthist] = VentilationFunctions.calculate_LB_VDP(MR, maskarray, maskarraytrachea, ventmean, ventstd, parentPath, N4_bias_analysis, Overall_SNR, Proton,MainInput);
% 
%     elseif settings.N4_bias_analysis == "no"
% %         ventmean = 0.5421; % Default (non-N4 corrected, based on data up to Jan 2021)
% %         ventstd = 0.1509; % Default (non-N4 corrected, based on data up to Jan 2021)
%         [LB_VDP, LB_mean, LB_std,VentscaledImage,Ventcolormap,LB_Venthist] = VentilationFunctions.calculate_LB_VDP(MR, maskarray, maskarraytrachea, ventmean, ventstd, parentPath, N4_bias_analysis, Overall_SNR, Proton,MainInput);
%     else
%         disp('Neither original or N4 bias corrected images were declared. Declare these in the input_params() function.')
%     end
    
    [Ventilation] = VentilationFunctions.calculate_LB_VDP(MR, maskarray, maskarraytrachea, ventmean, ventstd, parentPath, N4_bias_analysis, Overall_SNR,Ventilation,Proton,MainInput);    

end
close all;
%% GLRLM analysis
if Ventilation.GLRLM_Analysis == "yes" % 'yes'; || 'no'
    f = waitbar(.75,'Performing Linear Binning Analysis ...');
    waitbar(.90,f,'Performing Texture Analysis ...');
    pause(.1)     
    cd(parentPath)
    if settings.N4_bias_analysis == "yes"
       VentImage = N4.*maskarray;        
    else
       VentImage = MR.*maskarray; 
    end
    [Ventilation] = VentilationFunctions.GLRLM_Analysis(Ventilation,VentImage,parentPath);
end 
close all;
%% %% save maps in mat file
save_data=[parentPath,'\','Ventilation_Analysis','.mat'];
save(save_data);  

disp('ventilation analysis completed');
close all;
end %end of function



% Please add updates here
% 
% 
% 
% 

