


%% Make Proton Mask (trained ML model method)
clc
disp('Making Lung Mask (trained ML model method)...')

%Write Proton image to external file for Python-based ML segmentation

if isfile('Mat2Py_IMLvars.csv') == 0  %double-checking if already exists
    VentImage =  double(VentImage);
%     ProtonImage_112 =  ProtonImage(H_RecMatrix+1:2*H_RecMatrix,H_RecMatrix+1:2*H_RecMatrix,H_RecMatrix+1:2*H_RecMatrix); %taking only the central 112x112x112 points (currently only capable of modeling this)
    writematrix(ProtonImage_112,'Mat2Py_IMLvars.csv') 
end
% figure; Global.imslice(ProtonImage_112)
%% 

clc
pathToMod = fileparts(which('AutoSegmentation3D.py'));
if count(py.sys.path,pathToMod)==0
    insert(py.sys.path,int32(0),pathToMod)
    disp('Uploaded Segmentation3DModule')
end

lungmaskname = 'AutoMask.nii.gz';
MODELLOC = "D:\Github\HP_Xe_Analysis_App\+Segmentation\AutoSegment_3DGasExchange_Xe_200e.hdf5";
Model_ML = "AutoSegment_3DGasExchange_Xe_200e.hdf5";
threshold = 0.2;
command = strcat('AutoSegmentation3D.py', {' '}, 'VentImage.mat', {' '},MODELLOC,{' '},num2str(threshold),{' '},lungmaskname);
system(command);

% AutoLungMask = system(strcat('AutoSegmentation3D.py', {' '}, 'Mat2Py_IMLvars.csv', {' '},MODELLOC,{''}, Model_ML,{' '},num2str(threshold),{' '},lungmaskname));
%% 
clc
pythonScript = 'AutoSegmentation3D.py';
InputImg = 'VentImage.mat'; %'Mat2Py_IMLvars.csv'; %'VentImage.mat';
MODELLOC = "D:\Github\HP_Xe_Analysis_App\+Segmentation\AutoSegment_3DGasExchange_Xe_200e.hdf5";
threshold = 0.2;
OutputName = 'AutoMask.nii.gz';
command = join(['python', pythonScript,'',InputImg,'',MODELLOC,'',num2str(threshold),'',OutputName]);

% system(command);

%% 
clc
% C:\Users\bda5ik\AppData\Local\Programs\Python\Python311
pyenv('Version','C:\Users\bda5ik\AppData\Local\anaconda3\envs\asbML\python.exe')


%% 
Mat2Py_preprocessing = zeros(1, 112, 112, 112, 1);
img = VentImage./max(VentImage(:));
Mat2Py_preprocessing(1,:,:,:,1) = img;
% Images = (permute(Images, [1,4, 2, 3]));
% Images = Images(:,:,:,:,1);
size(Mat2Py_preprocessing)

Images1 = VentImage./max(VentImage);
Images1 = (permute(Images1, [3, 1, 2]));
figure; Global.imslice(auto_segmented_mask)

auto_mask = double(auto_segmented_mask);
figure; Global.imslice(squeeze(auto_mask(1,:,:,:)))






