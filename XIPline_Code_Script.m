%This script is intended for debugging, developing, and testing new tools before implementing them in the main interface.

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 
%% calibration 
[filename, path] = uigetfile('*.*','Select xenon data file');
XeFullPath = [path,filename];
XeDataLocation = path(1:end-1);
[~,~,xe_ext] = fileparts(XeFullPath);

MainInput.XeFullPath = XeFullPath;
MainInput.XeDataLocation = XeDataLocation;
MainInput.XeFileName = filename;
MainInput.XeDataext = xe_ext;
cd(MainInput.XeDataLocation)

% [GasExResults, CalResults] = Calibration.XeCTC_Calibration(MainInput);
[GasExResults, CalResults] = Calibration.XeCTC_Calibration_MRD(MainInput); 
% the rest of the calculations is in the main interface 

%% Image analysis
clc; clear; close all; % start with fresh variables, etc.
MainInput.Institute = '';
MainInput.Scanner = '';
MainInput.ScannerSoftware = '';
MainInput.SequenceType = '';
MainInput.XeDataLocation = '';
MainInput.XeDataType = '';
MainInput.NoProtonImage = '';
MainInput.HDataLocation = '';
MainInput.HDataType = '';
MainInput.AnalysisType = '';
MainInput.SubjectAge = [];
MainInput.RegistrationType = '';

% add patinet info
MainInput.PatientInfo = '';

% 1) choose the type of analysis
MainInput.AnalysisType = 'Ventilation';  % 'Ventilation', 'Diffusion', 'GasExchange'

% 2) Do you have protom images? 
MainInput.NoProtonImage = 1;  % 1: There is no proton images  % 0: There is  proton images   

MainInput.Institute = 'CCHMC';  % 'CCHMC', 'XeCTC', 'Duke'
MainInput.Scanner = 'Philips'; % Siemens, Philips
MainInput.ScannerSoftware = '5.9.0'; % '5.3.1', '5.6.1','5.9.0'
MainInput.SequenceType = '2D GRE'; % '2D GRE', '3D Radial'

% diary Log.txt
[filename, path] = uigetfile('*.*','Select xenon data file');
XeFullPath = [path,filename];
XeDataLocation = path(1:end-1);
[~,Xe_name,xe_ext] = fileparts(XeFullPath);

MainInput.XeFullPath = XeFullPath;
MainInput.XeDataLocation = XeDataLocation;
MainInput.XeFileName = filename;
MainInput.Xe_name = Xe_name;
MainInput.XeDataext = xe_ext;
cd(MainInput.XeDataLocation)

% select calibration
if strcmp(MainInput.Institute,'XeCTC') && (strcmp(MainInput.XeDataext,'.dat') || strcmp(MainInput.XeDataext,'.p'))
    [filename, path] = uigetfile('*.*','Select calibration data file');
    CalFullPath = [path,filename];
    CalDataLocation = path(1:end-1);
    [~,Cal_name,Cal_ext] = fileparts(CalFullPath);

    MainInput.CalFullPath = CalFullPath;
    MainInput.CalDataLocation = CalDataLocation;
    MainInput.CalFileName = filename;
    MainInput.Cal_name = Cal_name;    
    MainInput.CalDataext = Cal_ext;
end


% select proton
if MainInput.NoProtonImage == 0
    [filename, path] = uigetfile('*.*','Select proton data file');
    HFullPath = [path,filename];
    HDataLocation = path(1:end-1);
    [~,H_name,H_ext] = fileparts(HFullPath);

    MainInput.HFullPath = HFullPath;
    MainInput.HDataLocation = HDataLocation;
    MainInput.HFileName = filename;
    MainInput.H_name = H_name;    
    MainInput.HDataext = H_ext;
end

[Ventilation, Diffusion, GasExchange, Proton] = LoadData.LoadReadData(MainInput);

if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
    figure; Global.imslice(Ventilation.Image) 
elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
    figure; Global.imslice(Diffusion.Image)
elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
    figure; Global.imslice(GasExchange.VentImage)
end
if MainInput.NoProtonImage == 0
    figure; Global.imslice(Proton.Image)
end
% Ventilation.Image = flipdim(Ventilation.Image,3); 
% Ventilation.Image = permute(Ventilation.Image,[3 2 1]); % 
% Rotated = imrotate(Diffusion.Image,-90);
      
%% Registration 

% Move from 1st to 3rd dim and flip Z (optional)
% Prot_Test = permute(Ventilation.Image,[3 2 1]); %
% Test3 = flip(Prot_Test,1);
% Ventilation.Image = flipdim(Test3,3);
% figure; imslice(Ventilation.Image)

% Flip Z of Vent and Proton (optional)
% Ventilation.Image = flipdim(Ventilation.Image,3);
% Proton.Image = flipdim(Proton.Image,3); 
% Ventilation.LungMask = flipdim(Ventilation.LungMask,3);

% preform registration 
clc
cd(MainInput.XeDataLocation)
%Proton.Image = P_image; % Use only for interpulation
% 
% diary LogFile_LoadingData
MainInput.RegistrationType = 'affine'; % 'translation' | 'rigid' | 'similarity' | 'affine'
% enable this code in case number of slices is different between xe and H
MainInput.SliceSelection = 0;
if MainInput.SliceSelection == 1
    MainInput.Xestart = 1;
    MainInput.Xeend = 13;
    MainInput.Hstart = 1;
    MainInput.Hend = 10;
end

if strcmp(MainInput.AnalysisType, 'Ventilation')
    Xesize = size(Ventilation.Image);
    Hsize = size(Proton.Image);
    sizeRatio = Hsize./Xesize;
elseif strcmp(MainInput.AnalysisType, 'GasExchange')
    Xesize = size(GasExchange.VentImage);
    Hsize = size(Proton.Image);
    sizeRatio = Hsize./Xesize;
end

useRatio = 1;
if useRatio == 1
    MainInput.XeVoxelInfo.PixelSize1 = sizeRatio(1);
    MainInput.XeVoxelInfo.PixelSize2 = sizeRatio(2);
    MainInput.XeVoxelInfo.SliceThickness = sizeRatio(3);
else % manual
    MainInput.XeVoxelInfo.PixelSize1 = 1.13; %Y
    MainInput.XeVoxelInfo.PixelSize2 = 1.13; %X
    MainInput.XeVoxelInfo.SliceThickness = 1;
end

MainInput.ProtonVoxelInfo.PixelSize1 = 1;
MainInput.ProtonVoxelInfo.PixelSize2 = 1;
MainInput.ProtonVoxelInfo.SliceThickness = 1;

[Proton] = Registration.PerformRegistration(Proton,Ventilation,GasExchange,MainInput);
viewing_img = Proton.ProtonRegisteredColored;

S = orthosliceViewer(viewing_img,...      
                               'DisplayRange', [min(viewing_img(:))  max(viewing_img(:))]);

%% Load Mask (.nii or .gz)
% clc
% [filename, path] = uigetfile('*.*','Select proton data file');
% cd(path)
% Mask = LoadData.load_nii(filename);
% A = Mask.img;
% A = permute(A,[1 3 2]); 
% A = double(squeeze(A));
% A = imrotate(A,-90);
% A = flip(A,2);                   
% Mask = A;  
% figure; Global.imslice(Mask) 

%% segmenting 
clc
cd(MainInput.XeDataLocation)

% diary Log.txt
MainInput.SegmentationMethod = 'Threshold'; % 'Threshold' || 'Manual' || 'Auto'
MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'

MainInput.thresholdlevel = 1; % 'threshold' 
MainInput.SE = 1;

MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
MainInput.SliceOrientation = 'coronal'; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
[Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

figure; Global.imslice(Ventilation.LungMask)

%% Ventilation analysis
clc
close all;
cd(MainInput.XeDataLocation)

% diary LogFile_LoadingData
Ventilation.N4Analysis = 0;
Ventilation.IncompleteThresh = 60;
Ventilation.RFCorrect = 0;
Ventilation.CompleteThresh = 15;
Ventilation.HyperventilatedThresh = 200;
Ventilation.ThreshAnalysis = 'yes'; % 'yes'; || 'no'
Ventilation.LB_Analysis = 'no'; % 'yes'; || 'no'
Ventilation.Kmeans = 'no';  % 'yes'; || 'no'
Ventilation.DDI2D = 'yes';  % 'yes'; || 'no'
Ventilation.DDI3D = 'no';  % 'yes'; || 'no'
Ventilation.DDI3Dx = 3;
Ventilation.DDI3Dy = 3;
Ventilation.DDI3Dz = 3;
Ventilation.DDIDefectMap = 'Threshold';  % 'Threshold'; || 'Linear Binning' || 'Kmeans'
Ventilation.GLRLM_Analysis = 'no'; % 'yes'; || 'no'

if Ventilation.N4Analysis == 1
    Ventilation.LB_RefMean = 0.6786; % Defined by data collected up to Jan 2021.
    Ventilation.LB_RefSD =  0.1395; % Defined by data collected up to Jan 2021.
elseif Ventilation.N4Analysis == 0
    Ventilation.LB_RefMean = 0.5421; % Default (non-N4 corrected, based on data up to Jan 2021)
    Ventilation.LB_RefSD = 0.1509; % Default (non-N4 corrected, based on data up to Jan 2021)
end

[Ventilation] = VentilationFunctions.Ventilation_Analysis(Ventilation, Proton, MainInput);

%%  Diffusion analysis
clc
cd(MainInput.XeDataLocation)
% diary Log.txt
MainInput.PatientAge = '17';
Diffusion.ADCFittingType = 'Log Linear'; % 'Log Weighted Linear' | 'Log Linear' | 'Non-Linear' | 'Bayesian'
Diffusion.ADCAnalysisType = 'human'; % human | animals;  % human | animals
Diffusion.bvalues = '[0, 6.25, 12.5, 18.75, 25]';
% Diffusion.bvalues = '[0, 7.5, 15]';
% Diffusion.bvalues = '[0, 10, 20, 30]'; % for T-Scanner

Diffusion.ADCLB_Analysis = 'yes'; % 'yes'; || 'no'
Diffusion.ADCLB_RefMean = '0.0002*age+0.029'; 
Diffusion.ADCLB_RefSD = '5e-5*age+0.0121'; 

Diffusion.MorphometryAnalysis = 'yes';  % yes || no
Diffusion.MorphometryAnalysisType = 'human'; % human | animals
Diffusion.MA_WinBUGSPath = 'C:\Users\BAS8FL\Desktop\WinBUGS14';
Diffusion.CMMorphometry = 'yes';  % yes || no
Diffusion.SEMMorphometry = 'yes';  % yes || no
Diffusion.Do = 0.14; % cm2/s
Diffusion.Delta = 3.5; % ms

[Diffusion] = DiffusionFunctions.Diffusion_Analysis(Diffusion,MainInput);

% diary off
%% Gas Exchange analysis
clc;
cd(MainInput.XeDataLocation)
% diary Log.txt

MainInput.ImportHealthyCohort = 0; % 0: import CCHMC healthy Ref., 1: import .mat H. Ref., 2: enter values manually
if MainInput.ImportHealthyCohort == 0
    MainInput.HealthyReferenceType = 'Default';
elseif MainInput.ImportHealthyCohort == 1
    GasExchange.ImportHealthyCohort = 'yes';
    [filename, path] = uigetfile('*.mat*','Import Healthy Cohort');
    GasExchange.HealthyCohortFullPath = [path,filename];
    GasExchange.HealthyCohortDataLocation = path(1:end-1);
    [~,~,HealthyCohort_ext] = fileparts(GasExchange.HealthyCohortFullPath);
    GasExchange.HealthyCohortFileName = filename;
    GasExchange.HealthyCohort_ext = HealthyCohort_ext;
    MainInput.HealthyReferenceType = 'Import';
elseif MainInput.ImportHealthyCohort == 2
    GasExchange.VentHealthyMean = 0.51; GasExchange.VentHealthyStd = 0.19;%6
    GasExchange.DissolvedHealthyMean = 0.0075; GasExchange.DissolvedHealthyStd = 0.00125;%6
    GasExchange.BarrierHealthyMean = 0.0049; GasExchange.BarrierHealthyStd = 0.0015;%6
    GasExchange.RBCHealthyMean = 0.0026; GasExchange.RBCtHealthyStd = 0.0010;%6
    GasExchange.RBCBarrHealthyMean = 0.53; GasExchange.RBCBarrHealthyStd = 0.18;%6
    GasExchange.RBCOscHealthyMean = 8.9596; GasExchange.RBCOscHealthyStd = 10.5608;%6
    MainInput.HealthyReferenceType = 'Manual';
end

[GasExchange] = GasExchangeFunctions.GasExchange_Analysis(GasExchange,Proton,MainInput);

disp('Analysis done')
% diary off

%% patient report
clc
MainInput.studyID = '740H';
MainInput.studyType = 'CF';
MainInput.patientID = '001';
MainInput.scanDate = '7/6/2023';
MainInput.xeDoseVolume = '1000';
MainInput.PatientName = 'test test';
MainInput.MRMnumber = '1201201';
MainInput.sex = 'M';
MainInput.age = '25';
MainInput.height = '180';
MainInput.weight = '200';
MainInput.summaryofFindings = 'all good';
MainInput.notes = 'all good';
MainInput.dataAnalyst = 'ASB';
MainInput.processingDate = '7/6/2023';

Global.PatientReport(MainInput)


%% 









