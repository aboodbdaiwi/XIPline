
clc; clear; close all; % start with fresh variables, etc.

% commit test 2 

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
MainInput.PatientInfo = 'CF Perfusion|IRC668-005.V2|11y/o,M|04/01/2022|';



% 1) choose the type of analysis
% MainInput.AnalysisType = 'Ventilation';
% MainInput.AnalysisType = 'Diffusion';
MainInput.AnalysisType = 'GasExchange';

% 2) Do you have protom images? 
% MainInput.NoProtonImage = 'yes';  % There is no proton images 
MainInput.NoProtonImage = 'no';    % There is  proton images 

MainInput.Institute = 'CCHMC'; 
% MainInput.Institute = 'XeCTC'; 
% MainInput.Institute = 'Duke';
MainInput.Scanner = 'Philips';
% MainInput.ScannerSoftware = '5.3.1'; % R-scanner
MainInput.ScannerSoftware = '5.9.0'; % R-scanner
% MainInput.ScannerSoftware = '5.6.1'; % T1-scanner
% MainInput.SequenceType = '2D GRE';
MainInput.SequenceType = '3D Radial';


% diary Log.txt
[filename, path] = uigetfile('*.*','Select xenon data file');
XeFullPath = [path,filename];
XeDataLocation = path(1:end-1);
[~,~,xe_ext] = fileparts(XeFullPath);

MainInput.XeFullPath = XeFullPath;
MainInput.XeDataLocation = XeDataLocation;
MainInput.XeFileName = filename;
MainInput.XeDataext = xe_ext;
cd(MainInput.XeDataLocation)

% select proton
if strcmp(MainInput.NoProtonImage, 'no') == 1
    [filename, path] = uigetfile('*.*','Select proton data file');
    HFullPath = [path,filename];
    HDataLocation = path(1:end-1);
    [~,~,H_ext] = fileparts(HFullPath);

    MainInput.HFullPath = HFullPath;
    MainInput.HDataLocation = HDataLocation;
    MainInput.HFileName = filename;
    MainInput.HDataext = H_ext;
end


[Ventilation, Diffusion, GasExchange, Proton] = LoadData.LoadReadData(MainInput);

if strcmp(MainInput.AnalysisType, 'Ventilation') == 1
    if (strcmp(xe_ext, '.nii') == 0) && (strcmp(xe_ext, '.gz') == 0)
        disp('nii images saved')
        MR = double(Ventilation.Image);
        niftiwrite(abs(fliplr(rot90(MR,3))),[MainInput.XeDataLocation,'\VentilationImage.nii'],'Compressed',true);
        % niftiwrite(abs(MR,3),[MainInput.XeDataLocation,'\VentilationImage.nii'],'Compressed',true);
    end 
end

if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
    figure; Global.imslice(Ventilation.Image) 
elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
    figure; Global.imslice(Diffusion.Image)
elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
    figure; Global.imslice(GasExchange.VentImage)
end
% Ventilation.Image = double(Ventilation.Image)
% Proton.Image = double(Proton.Image)

% Ventilation.Image = imrotate(Ventilation.Image,-180);
% figure; Global.imslice(Ventilation.Image(:,:,1,1))
figure; Global.imslice(Proton.Image)
% figure; Global.imslice(A) 
% figure; Global.imslice(Proton.ProtonRegistered)
% orthosliceViewer(Proton.ProtonRegisteredColored)
% DiffusionDataImage = DiffusionData.Image;
% Ventilation.Image = flipdim(Ventilation.Image,3); 
% Ventilation.Image = permute(Ventilation.Image,[3 2 1]); % 

% Rotated = imrotate(Diffusion.Image,-90);
% Ventilation.Image = flipdim(Ventilation.Image,1); % Original
% figure; Global.imslice(Rotated)

% Diffusion.Image = double(Rotated);
% figure; Global.imslice(Diffusion.Image)

% figure; Global.imslice(Ventilation.Image)
% 
% Proton.Image = Proton.Image(:,:,16:42);
% figure; Global.imslice(T2)

% diary off
% viewing_img = Ventilation.Image;
% S = orthosliceViewer(viewing_img,...      
%                                'DisplayRange', [min(viewing_img(:))  max(viewing_img(:))]);
                      
%% Registration 
clc
cd(MainInput.XeDataLocation)
Vent_slice_indx = 3:12;
Prot_slice_indx = 3:12;
Prot_Test = squeeze(Proton.Image(:,:,Prot_slice_indx)); % numerator
Vent_Test = squeeze(Ventilation.Image(:,:,Vent_slice_indx)); % denomenator
Ratio = size(Prot_Test, 3) ./ size(Vent_Test, 3);

%% -- Have # of slices of Test2 equal Image3D
scale_factor = 0; % if orginal ratio claculation doesn't work, use this as a fudge factor to get the slices you want

Img = squeeze(Prot_Test(:,:,:,1));
Image_size = size(Img); %size of the original image
largest_dim=max(Image_size); %find the max dim
[rows,cols,slices] = size(Img);
% find steps for the meshgrid
stepX=rows/largest_dim; 
stepY=cols/largest_dim;
stepZ = Ratio + scale_factor; % ratio between              
[X,Y,Z] = meshgrid(1:Image_size(2), 1:Image_size(1), 1:Image_size(3));
[X2,Y2,Z2] = meshgrid(1:stepX:rows, 1:stepY:cols, 1:stepZ:slices);
Image_3D = interp3(X, Y, Z, Img, X2, Y2, Z2, 'linear', 0);

P_image = zeros(size(Prot_Test,1),size(Prot_Test,2),size(Ventilation.Image,3));
P_image(:,:,Vent_slice_indx) = Image_3D;
    

figure; Global.imslice(P_image);

%% Move from 1st to 3rd dim and flip Z (optional)
Prot_Test = permute(Ventilation.Image,[3 2 1]); %
Test3 = flip(Prot_Test,1);
Ventilation.Image = flipdim(Test3,3);
figure; imslice(Ventilation.Image)
%% Flip Z of Vent and Proton (optional)
Ventilation.Image = flipdim(Ventilation.Image,3);
Proton.Image = flipdim(Proton.Image,3); 
% Ventilation.LungMask = flipdim(Ventilation.LungMask,3);
%% preform registration 
clc
cd(MainInput.XeDataLocation)

Proton.Image = P_image; % Use only for interpulation
% 
% diary LogFile_LoadingData
MainInput.RegistrationType = 'affine'; % 'translation' | 'rigid' | 'similarity' | 'affine'
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

MainInput.XeVoxelInfo.PixelSize1 = sizeRatio(1);
MainInput.XeVoxelInfo.PixelSize2 = sizeRatio(2);
MainInput.XeVoxelInfo.SliceThickness = sizeRatio(3);

% MainInput.XeVoxelInfo.PixelSize1 = 1.13; %Y
% MainInput.XeVoxelInfo.PixelSize2 = 1.13; %X
% MainInput.XeVoxelInfo.SliceThickness = 1;

MainInput.ProtonVoxelInfo.PixelSize1 = 1;
MainInput.ProtonVoxelInfo.PixelSize2 = 1;
MainInput.ProtonVoxelInfo.SliceThickness = 1;
% 

[Proton] = Registration.PerformRegistration(Proton,Ventilation,GasExchange,MainInput);
viewing_img = Proton.ProtonRegisteredColored;

% viewing_img = Proton.ProtonRegistered;
% figure; Global.imslice(A) 
% figure; Global.imslice(MR)
% S = orthosliceViewer(Proton.ProtonRegisteredColored, 'Colormap', eval('jet'));
% viewing_img = double(abs(viewing_img));
S = orthosliceViewer(viewing_img,...      
                               'DisplayRange', [min(viewing_img(:))  max(viewing_img(:))]);
% diary off
% S = orthosliceViewer(Proton.ProtonRegisteredColored);

% % Get the dimensions of the original image
% original_image = double(LungMask);
% [height, width, depth] = size(original_image);
% 
% % Calculate the cropping boundaries to obtain a centered crop of size 112x112x112
% crop_start = [(height - 112) / 2 + 1, (width - 112) / 2 + 1, (depth - 112) / 2 + 1];
% crop_end = crop_start + [111, 111, 111];
% 
% % Crop the image
% cropped_image = original_image(crop_start(1):crop_end(1), crop_start(2):crop_end(2), crop_start(3):crop_end(3));
% montage(cropped_image)
%% Flip Mask for Coronal Slices (optional)
% Read in Image:
msgboxw('Please, Select Nifti Image File (Mask or Image)',14);
[parentFile,parentPath] = uigetfile('*.nii.gz', 'Select Nifti file');

% Create a storage folder:
c_parentPath = parentPath;
mkdir([c_parentPath 'Coronal to Axial Images']);
sub_parentPath = [c_parentPath,'\Coronal to Axial Images\'];

% Read in the Nifti images
Before = load_nii([parentPath parentFile]); % Original
Before = Before.img;
Before = double(squeeze(Before));

niftiwrite(abs(Before),[sub_parentPath,'CoronalImage_before.nii']);%Image

% Dimension shifted image set:
After = permute(Before,[3 1 2]); % 
Rotated = imrotate(After,-90);
Flipped = flipdim(Rotated,1); % Original
% Normalize the intensity of the original image to fall between [0,1].
% After = After-min(After(:));
% After = After/max(After(:));
niftiwrite(abs(Flipped),[sub_parentPath,'AxialImage_after.nii']);%Image

% Conversion Complete Statement
fprintf( 'Conversion Complete\n' );
msgboxw('Conversion complete - See "Coronal to Axial Images" folder.',14);
%% Load Mask
clc
[filename, path] = uigetfile('*.*','Select proton data file');
cd(path)
% load([path,filename])
% Mask = load_nii('C:\Users\BAS8FL\Documents\DUSTIN DOCS\Image Analysis\IRC186H-411\20211215\1101_CPIR_Ventilation_021121514024081114\Mask 2\IRC186H-411_Corrected_Vent_20211215_mask.nii');
% Mask = double(labels);

Mask = LoadData.load_nii(filename);
% Mask = LoadData.load_untouch_nii(filename);
A = Mask.img;
A = permute(A,[1 3 2]); 
A = double(squeeze(A));
A = imrotate(A,-90);
A = flip(A,2);                   
Mask = A;  
figure; Global.imslice(A) 
figure; Global.imslice(Mask) 
%% 

Ventilation.LungMask = Mask;
Ventilation.AirwayMask = zeros(size(Mask));

Ventilation.LungMask = imrotate(Ventilation.LungMask,-180);

Ventilation.LungMask = permute(Ventilation.LungMask,[3 2 1]); % 
% Rotated = imrotate(After,-90);
Ventilation.LungMask = flipdim(Ventilation.LungMask,1); % Original



figure; Global.imslice(Ventilation.LungMask) 
% figure; Global.imslice(A)

%% segmenting 
clc
cd(MainInput.XeDataLocation)

% diary Log.txt
MainInput.SegmentationMethod = 'Threshold'; % 'Threshold' || 'Manual' || 'Auto'
MainInput.SegmentAnatomy = 'Airway'; % 'Airway'; || 'parenchyma'
MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'

MainInput.thresholdlevel = 1; % 'threshold' 
MainInput.SE = 1; % 'threshold' 

MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
MainInput.SliceOrientation = 'coronal'; % 'transversal' || 'isotropic'
[Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

% [Proton,Ventilation,Diffusion, GasExchange] = Segmentation.loadSegmentedMask(Proton,Ventilation,Diffusion,GasExchange,MainInput);
% Ventilation.AirwayMask = zeros(size(Ventilation.LungMask));

figure; Global.imslice(Ventilation.LungMask)
% figure; Global.imslice(Images)
% figure; Global.imslice(auto_segmented_mask) 
% figure; Global.imslice(Diffusion.LungMask)
% volumeSegmenter(Ventilation.Image) 
% figure; Global.imslice(NormMR)
% diary off 

% [preprocessed_data, MainInput] = Segmentation.preprocess_images_for_auto_segmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);
% cd(MainInput.AutoSegmentPath)
% if strcmp(MainInput.NoProtonImage, 'no') == 1
%     system('predict_mask_Vent_w_H.exe')
% else 
%     system('predict_mask_Vent_wout_H.exe')           
% end
% 
% cd(MainInput.XeDataLocation)

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

MainInput.PatientAge = 11;
Diffusion.ADCFittingType = 'Log Linear'; % 'Log Weighted Linear' | 'Log Linear' | 'Non-Linear' | 'Bayesian'
Diffusion.ADCAnalysisType = 'human'; % human | animals;  % human | animals
Diffusion.bvalues = '[0, 6.25, 12.5, 18.75, 25]';
% Diffusion.bvalues = '[0, 7.5, 15]';
% Diffusion.bvalues = '[0, 10, 20, 30]'; % for T-Scanner

Diffusion.ADCLB_Analysis = 'no'; % 'yes'; || 'no'
Diffusion.ADCLB_RefMean = '0.0002*age+0.029'; 
Diffusion.ADCLB_RefSD = '5e-5*age+0.0121'; 

Diffusion.MorphometryAnalysis = 'yes';  % yes || no
Diffusion.MorphometryAnalysisType = 'human'; % human | animals
Diffusion.MA_WinBUGSPath = 'C:\Users\BAS8FL\Desktop\WinBUGS14';
Diffusion.CMMorphometry = 'yes';  % yes || no
Diffusion.SEMMorphometry = 'no';  % yes || no
Diffusion.Do = 0.14; % cm2/s
Diffusion.Delta = 3.5; % ms


[Diffusion] = DiffusionFunctions.Diffusion_Analysis(Diffusion,MainInput);

% diary off
%% Gas Exchange analysis
clc;
cd(MainInput.XeDataLocation)
% diary Log.txt

MainInput.ImportHealthyCohort = 1;
if MainInput.ImportHealthyCohort == 1
    GasExchange.ImportHealthyCohort = 'yes';
    [filename, path] = uigetfile('*.mat*','Import Healthy Cohort');
    GasExchange.HealthyCohortFullPath = [path,filename];
    GasExchange.HealthyCohortDataLocation = path(1:end-1);
    [~,~,HealthyCohort_ext] = fileparts(GasExchange.HealthyCohortFullPath);
    GasExchange.HealthyCohortFileName = filename;
    GasExchange.HealthyCohort_ext = HealthyCohort_ext;
else
    GasExchange.VentHealthyMean = 0.51; GasExchange.VentHealthyStd = 0.19;%6
    GasExchange.DissolvedHealthyMean = 0.0075; GasExchange.DissolvedHealthyStd = 0.00125;%6
    GasExchange.BarrierHealthyMean = 0.0049; GasExchange.BarrierHealthyStd = 0.0015;%6
    GasExchange.RBCHealthyMean = 0.0026; GasExchange.RBCtHealthyStd = 0.0010;%6
    GasExchange.RBCBarrHealthyMean = 0.53; GasExchange.RBCBarrHealthyStd = 0.18;%6
    GasExchange.RBCOscHealthyMean = 8.9596; GasExchange.RBCOscHealthyStd = 10.5608;%6
end



[GasExchange] = GasExchangeFunctions.GasExchange_Analysis(GasExchange,Proton,MainInput);

disp('Analysis done')
% diary off



%% 
clc
[Images, MainInput] = Segmentation.preprocess_images_for_auto_segmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);
cd(MainInput.AutoSegmentPath)
%system('predict_mask_Vent_w_H.exe')
cd(MainInput.XeDataLocation)

figure; Global.imslice(squeeze(Images))


 %% 
clc                     
ppt = actxserver('PowerPoint.Application');


presentation = ppt.Presentations.Open('D:\OneDrive - cchmc\Lab\Xe_App\testtingData\diffdata\Subject1\Diffusion Analysis\Diffusion_Analysis.pptx');
outputPath = 'D:\OneDrive - cchmc\Lab\Xe_App\testtingData\diffdata\Subject1\Diffusion Analysis\\Diffusion_Analysis.pdf';
presentation.SaveAs(outputPath, 32);

presentation.Close();
ppt.Quit();
delete(ppt);

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









