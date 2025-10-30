function CCHMC_Db_Diff_Pipeline(MainInput)

warning('off', 'all'); % Turn off all warnings 'off' | 'on'

SubjectID = MainInput.SubjectID;
Age = MainInput.Age;
Sex = MainInput.Sex;
Disease = MainInput.Disease;
ScanDate = MainInput.ScanDate;
ReconType = MainInput.ReconType;

diff_file = MainInput.diff_file;

analysisFolder = MainInput.analysisfolder;

MainInput.OutputPath = analysisFolder;

xedatapath = diff_file;

Ventilation = '';
GasExchange = '';
MainInput.AnalysisType = 'Diffusion';
MainInput.Institute = 'CCHMC'; 

% get segmentation 
MainInput.SE = 1;
MainInput.thresholdlevel = 0.6;

%%

Outputs = [];
Outputs.SUBJECT_ID = SubjectID;
Outputs.Age = Age;
Outputs.Sex = Sex;
Outputs.Disease = Disease;
Outputs.ScanDate = ScanDate;
Outputs.SW_VER = MainInput.ScannerSoftware;
Outputs.DIFF_FILEPATH_NEW = diff_file(57:end);
Outputs.analysispath = analysisFolder(57:end);
Outputs.maindirectory = diff_file(1:56);
Outputs.ReconType = ReconType;
Outputs.MaskPath = analysisFolder(57:end); 
Outputs.Note = MainInput.Note;
Outputs.sernum = MainInput.sernum;

Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
Outputs.AnalysisDate = str2double(datestr(datetime('today'), 'yyyymmdd'));
AnalysisDate = Outputs.AnalysisDate;

%%
analysisSubfolder = MainInput.diff_analysis_folder;
MainInput.analysissesFolder = analysisSubfolder;

% Define file extensions to delete
fileExtensions = {'*.pdf', '*.ppt', '*.pptx'};
% Loop through each file type and delete
for i = 1:numel(fileExtensions)
    files = dir(fullfile(analysisSubfolder, fileExtensions{i}));
    for k = 1:numel(files)
        filePath = fullfile(analysisSubfolder, files(k).name);
        delete(filePath);
        fprintf('Deleted: %s\n', filePath);
    end
end

[path,filename] = fileparts(xedatapath);
XeFullPath = xedatapath;
XeDataLocation = path;
[~,XeFileName,xe_ext] = fileparts(xedatapath);

MainInput.XeFullPath = XeFullPath;
MainInput.XeDataLocation = XeDataLocation;
MainInput.XeFileName = filename;
MainInput.XeDataext = xe_ext;
cd(MainInput.XeDataLocation)

% Outputs.XeFullPath = XeFullPath;
Outputs.XeDataLocation = XeDataLocation(57:end);
Outputs.XeFileName = XeFileName;
Outputs.XeDataext = xe_ext;


%% Load Data
MainInput.NoProtonImage = 1;
MainInput.denoiseXe = 'no';
[Ventilation, Diffusion, GasExchange, Proton, MainInput] = LoadData.LoadReadData(MainInput);
Diffusion.Image = double(Diffusion.Image);

Diffusion.outputpath = analysisSubfolder;
Outputs.AnalysisCode_hash = MainInput.AnalysisCode_hash; 

%% Reorient Diffusion Images if necessary



%%
cd(MainInput.XeDataLocation)

% % % diary Log.txt
MainInput.SegmentationMethod = 'Auto'; % 'Threshold' || 'Manual' || 'Auto'
MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'


MainInput.AIScript = 'Python';
MainInput.thresholdlevel = 1; % 'threshold' 
MainInput.SE = 1;

MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
% MainInput.SliceOrientation = SliceOrientation; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
[Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

if exist('Diffusion.AirwayMask', 'var')
    % skip
else
    Diffusion.AirwayMask = zeros(size(Diffusion.LungMask));
end

%%
%MainInput.PatientAge = MainInput.Age;
Diffusion.ADCFittingType = 'Log Weighted Linear';
Diffusion.ADCLB_Analysis = 'yes';
Diffusion.ADCLB_RefMean = 0.0002*MainInput.Age+0.029; % mean equetion for healthy co. 
Diffusion.ADCLB_RefSD = 5e-5*MainInput.Age+0.0121; 
%%

%%%%%%%%%%%%Diffusion = DiffusionFunctions.Diffusion_Analysis(Diffusion,MainInput);
WinBUGSPath = mfilename('fullpath');
idcs = strfind(WinBUGSPath,filesep);%determine location of file separators
WinBUGSPath = [WinBUGSPath(1:idcs(end)-1),filesep];%remove file
WinBUGSPath = [WinBUGSPath,'\WinBUGS14'];

% WinBUGSPath = 'P:\OneDrive - cchmc\Lab\WinBUGS14';

Diffusion.WinBUGSPath = WinBUGSPath;  
Diffusion.MA_WinBUGSPath = WinBUGSPath;
try
    PatientAge = eval(MainInput.PatientAge);
catch
    PatientAge = MainInput.Age;
end
% ADC_Analysis = Diffusion.ADC_Analysis;
ADC_Analysis = 'yes';
ADCFittingType = Diffusion.ADCFittingType;

bvalues = MainInput.b_values;
ADCLB_Analysis = Diffusion.ADCLB_Analysis;
Diffusion.MorphometryAnalysis = 'yes';
MorphometryAnalysis = Diffusion.MorphometryAnalysis;
% ADCAnalysisType = Diffusion.ADCAnalysisType;  % human | animals
ADCAnalysisType = 'human';
diffimg = Diffusion.Image;
lung_mask = Diffusion.LungMask;
try
    airway_mask = Diffusion.AirwayMask;
catch
    airway_mask = zeros(size(lung_mask));
end


% create a final mask
final_mask = lung_mask;
final_mask(airway_mask==1)=0;

% flip the binary mask to create the noise mask(1 to 0, 0 to 1)
noise_mask = ((~lung_mask)+(~airway_mask))-1;
noise_mask=noise_mask>0;

%% get rid of any slice that has no signal (un-segmented)
slices = size(diffimg,3);
Nbvalues = length(bvalues);
slice_checker=zeros(1,slices);
for check_slice=1:slices
    sl=lung_mask(:,:,check_slice);
    if sum(sl(:))>1
        slice_checker(check_slice)=1;
    elseif sum(sl(:))== 0
        slice_checker(check_slice) = 0;
    end
end
find_slice_checker = find(slice_checker);
Ndiffimg = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker), Nbvalues);
Nlung_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nairway_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nfinal_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nnoise_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
for find_sl=1:length(find_slice_checker)
    choose_sl=find_slice_checker(find_sl);
    Ndiffimg(:,:,find_sl,:)=diffimg(:,:,choose_sl,:);
    Nlung_mask(:,:,find_sl)=lung_mask(:,:,choose_sl);
    Nairway_mask(:,:,find_sl)=airway_mask(:,:,choose_sl);
    Nfinal_mask(:,:,find_sl)=final_mask(:,:,choose_sl);
    Nnoise_mask(:,:,find_sl)=noise_mask(:,:,choose_sl);
end
% figure; imslice(Ndiffimg)
% figure; imslice(Nfinal_mask)
diffimg = Ndiffimg;
lung_mask = Nlung_mask;
airway_mask = Nairway_mask;
final_mask = Nfinal_mask;
noise_mask = Nnoise_mask;
slices = length(find_slice_checker);

Diffusion.lung_mask = lung_mask;
Diffusion.airway_mask = airway_mask;
Diffusion.final_mask = final_mask;
Diffusion.noise_mask = noise_mask;

%% create path
Nbvalues = length(bvalues);
outputpath = [MainInput.diff_analysis_folder '\'];
cd(outputpath)

%% Save iamges and maskes
% diffusion images
niftiwrite(abs(fliplr(rot90(diffimg,-1))),[outputpath,'\DiffusionImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\DiffusionImage.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(diffimg,-1))),[outputpath,'\DiffusionImage.nii'],info,'Compressed',true);
% lung mask
niftiwrite(abs(fliplr(rot90(lung_mask,-1))),[outputpath,'\LungMask.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\LungMask.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(lung_mask,-1))),[outputpath,'\LungMask.nii'],info,'Compressed',true);
% airway mask
niftiwrite(abs(fliplr(rot90(airway_mask,-1))),[outputpath,'\AirwayMask.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\AirwayMask.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(airway_mask,-1))),[outputpath,'\AirwayMask.nii'],info,'Compressed',true);

%% ADC Analysis
% waitbar(.10,f,'Calculating ADC');
% pause(1)
if strcmp(ADC_Analysis, 'yes') == 1
    [ADCmap,ADCcoloredmap,SNR_table,meanADC,stdADC,ADC_hist] = DiffusionFunctions.ADC_Analysis(...
        diffimg,...
        lung_mask,...
        airway_mask,...
        bvalues,...
        ADCFittingType,...
        ADCAnalysisType,...
        WinBUGSPath,...
        outputpath,...
        PatientAge);

    % store result    
    Diffusion.ADCmap = ADCmap;
    Diffusion.ADCcoloredmap = ADCcoloredmap;
    Diffusion.SNR_table = SNR_table;
    Diffusion.meanADC = meanADC;
    Diffusion.stdADC = stdADC;
    Diffusion.ADC_hist = ADC_hist;
end

%% ADC linear binning analysis
% waitbar(.20,f,'Processing Linear Binning Analysis');
% pause(1)
if strcmp(ADCLB_Analysis, 'yes') == 1 && strcmp(ADC_Analysis, 'yes') == 1
    [colorDiffmap,DiffMean,DiffStd,LBADC_hist,BinTable] = DiffusionFunctions.LinearBinningADC(...
        Diffusion.Image,...         % diffusion image
        ADCmap,...                  % ADC amps
        final_mask,...              % lung mask without airways
        PatientAge,...              % Patient age for reference
        outputpath,...              % Results location 
        Diffusion.ADCLB_RefMean,... % mean equetion for healthy co. 
        Diffusion.ADCLB_RefSD);     % std equetion for healthy co. 

    % store result    
    Diffusion.LB_colorDiffmap = colorDiffmap;
    Diffusion.LBADCMean = DiffMean;
    Diffusion.LBADCStd = DiffStd;
    Diffusion.LBADC_hist = LBADC_hist;
    Diffusion.LB_BinTable = BinTable;
    Diffusion.DiffLow1Percent = BinTable{3,2};
    Diffusion.DiffLow2Percent = BinTable{3,3};
    Diffusion.DiffNormal1Percent = BinTable{3,4};
    Diffusion.DiffNormal2Percent = BinTable{3,5};
    Diffusion.DiffHigh1Percent = BinTable{3,6};
    Diffusion.DiffHigh2Percent = BinTable{3,7};

end
%% Morphometry Analysis
if strcmp(MorphometryAnalysis, 'yes') == 1 && Nbvalues > 3
% CM Morphometry analysis
%     waitbar(.30,f,'Processing CM Morphometry Analysis ');
%     pause(1)
    if strcmp(Diffusion.CMMorphometry, 'yes') == 1
        cd(outputpath)
%         [R_map,h_map,r_map,Lm_map,SVR_map,Na_map,DT_map,Dan_map,DL_map,Dmean_map,So_map,MSR_map,Residual_map,LungCMMorphometrySummary]=...
%             DiffusionFunctions.CM_Lung_Morphometry(...
%             diffimg,...
%             final_mask,...
%             noise_mask,...
%             Diffusion.MorphometryAnalysisType,...
%             bvalues,...
%             Diffusion.Do ,...
%             Diffusion.Delta,...
%             outputpath,...
%             WinBUGSPath);
%         
%         % store result
%         Diffusion.DT_map = DT_map;
%         Diffusion.Dan_map = Dan_map;
%         Diffusion.DL_map = DL_map;
%         Diffusion.Dmean_map = Dmean_map;
%         Diffusion.R_map = R_map;
%         Diffusion.h_map = h_map;
%         Diffusion.r_map = r_map;
%         Diffusion.Lm_map = Lm_map;
%         Diffusion.SVR_map = SVR_map;
%         Diffusion.Na_map = Na_map;
%         Diffusion.So_map = So_map;
%         Diffusion.MSR_map = MSR_map;
%         Diffusion.Residual_map = Residual_map;
%         Diffusion.DT_mean = LungCMMorphometrySummary{1,2};
%         Diffusion.DL_mean = LungCMMorphometrySummary{1,3};
%         Diffusion.Dan_mean = LungCMMorphometrySummary{1,4};
%         Diffusion.Dmean_mean = LungCMMorphometrySummary{1,5};
%         Diffusion.R_mean = LungCMMorphometrySummary{1,6};
%         Diffusion.h_mean = LungCMMorphometrySummary{1,7};
%         Diffusion.r_mean = LungCMMorphometrySummary{1,8};
%         Diffusion.Lm_mean = LungCMMorphometrySummary{1,9};
%         Diffusion.SVR_mean = LungCMMorphometrySummary{1,10};
%         Diffusion.Na_mean = LungCMMorphometrySummary{1,11};
%         Diffusion.DT_std = LungCMMorphometrySummary{2,2};
%         Diffusion.DL_std = LungCMMorphometrySummary{2,3};
%         Diffusion.Dan_std = LungCMMorphometrySummary{2,4};
%         Diffusion.Dmean_std = LungCMMorphometrySummary{2,5};
%         Diffusion.R_std = LungCMMorphometrySummary{2,6};
%         Diffusion.h_std = LungCMMorphometrySummary{2,7};
%         Diffusion.r_std = LungCMMorphometrySummary{2,8};
%         Diffusion.Lm_std = LungCMMorphometrySummary{2,9};
%         Diffusion.SVR_std = LungCMMorphometrySummary{2,10};
%         Diffusion.Na_std = LungCMMorphometrySummary{2,11};
         
[R_map,h_map,r_map,Lm_map,SVR_map,Na_map,So_map,LungCMMorphometrySummary]=...
            DiffusionFunctions.CM_MorphometryFit(...
            diffimg,...
            final_mask,...
            noise_mask,...
            Diffusion.MorphometryAnalysisType,...
            bvalues,...
            Diffusion.Do ,...
            Diffusion.Delta,...
            outputpath,...
            WinBUGSPath);
        
        % store result
        Diffusion.R_map = R_map;
        Diffusion.h_map = h_map;
        Diffusion.r_map = r_map;
        Diffusion.Lm_map = Lm_map;
        Diffusion.SVR_map = SVR_map;
        Diffusion.Na_map = Na_map;
        Diffusion.So_map = So_map;
        
        Diffusion.R_mean = LungCMMorphometrySummary{1,2};
        Diffusion.h_mean = LungCMMorphometrySummary{1,3};
        Diffusion.r_mean = LungCMMorphometrySummary{1,4};
        Diffusion.Lm_mean = LungCMMorphometrySummary{1,5};
        Diffusion.SVR_mean = LungCMMorphometrySummary{1,6};
        Diffusion.Na_mean = LungCMMorphometrySummary{1,7};
        
        Diffusion.R_std = LungCMMorphometrySummary{2,2};
        Diffusion.h_std = LungCMMorphometrySummary{2,3};
        Diffusion.r_std = LungCMMorphometrySummary{2,4};
        Diffusion.Lm_std = LungCMMorphometrySummary{2,5};
        Diffusion.SVR_std = LungCMMorphometrySummary{2,6};
        Diffusion.Na_std = LungCMMorphometrySummary{2,7};
        
    end
    %% SEM Morphometry analysis
%     waitbar(.80,f,'Processing SEM Morphometry Analysis ');
%     pause(1)
    if strcmp(Diffusion.SEMMorphometry, 'yes') == 1
        cd(outputpath)
%         [DDC_map,alpha_map,SEMSo_map,LmD_map,LungSEMMorphometrySummary]= DiffusionFunctions.SEM_Morphometry...
%             (diffimg,...
%             final_mask,...
%             noise_mask,...
%             bvalues,...
%             Diffusion.Do ,...
%             Diffusion.Delta,...
%             outputpath,...
%             WinBUGSPath);
        [DDC_map,alpha_map,SEMSo_map,LmD_map,LungSEMMorphometrySummary]= DiffusionFunctions.SEM_MorphometryFit...
            (diffimg,...
            final_mask,...
            noise_mask,...
            bvalues,...
            Diffusion.Do ,...
            Diffusion.Delta,...
            outputpath,...
            WinBUGSPath);
        % store result
        Diffusion.DDC_map = DDC_map;
        Diffusion.alpha_map = alpha_map;
        Diffusion.SEMSo_map = SEMSo_map;
        Diffusion.LmD_map = LmD_map;
        Diffusion.DDC_mean = LungSEMMorphometrySummary{1,2};
        Diffusion.alpha_mean = LungSEMMorphometrySummary{1,3};
        Diffusion.LmD_mean = LungSEMMorphometrySummary{1,4};
        Diffusion.DDC_std = LungSEMMorphometrySummary{2,2};
        Diffusion.alpha_std = LungSEMMorphometrySummary{2,3};
        Diffusion.LmD_std = LungSEMMorphometrySummary{2,4};        
    
    end
elseif strcmp(MorphometryAnalysis, 'yes') == 1 && Nbvalues < 3
    disp('Morphometry cannot be performed with only 2 b values')
end
close all;
%% save maps in mat file
save_data=[outputpath,'Diffusion_Analysis.mat'];
save(save_data);   

% waitbar(1,f,'Finishing');
% pause(1)
% close(f)
end


% Please add updates here
% 
% 
% 
% 
