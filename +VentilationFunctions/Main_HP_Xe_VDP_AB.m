
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matlab script for processing VDP Analysis of Philips dcm data %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Writen by David Roach <David.Roach@cchmc.org>
% Edited & orginized by Abdullah Bdaiwi <Abdullah.Bdaiwi@cchmc.org> 1/16/2020
%
% Function called:  medFilter.m
%                   rfCorrection.m
%                   xform_nii.m
%                   DICOM_Load.m % to read all dicom files at ones : Abdullah 1/20/2020
%                   msgboxw.m 
%                   load_nii.m
%                   exportToPPTX.m  % to save result in powerpoint : Abdullah 1/25/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;% start with fresh variables, etc.
%% setting up analysis
%% 
%Original Analysis set up
data_type = 'multi'; %choose "single_dcm" if it's only one file
RFCorrection='yes';
Original_analysis='yes'; %will work only when Original_analysis='yes' & N4_bias_analysis='no' 

%N4 bias Analysis set up
N4_bias_analysis='no';  %will work only when N4_bias_analysis='yes' & Original_analysis='no';

%General set up (for both Original Analysis and N4 bias Analysis)
Median_Filter= 'yes';
savedata='yes';
Incomplete_threshold = 60; % Defects Thresholds
Complete_threshold = 15;
calculate_SNR='yes';       

%%  Subject Info.
%% 
fprintf( 'Patients information being entered, please wait...\n' );
Subject_Info = inputdlg({'Subject No.','Study Group', 'Age', 'Gender', 'Scan Date','Notes'},...
              'Subject_Info', [1 35; 1 35; 1 35; 1 35; 1 35; 3 35],{'IRC186H_240',' Asthma,',' 16 yr,',' M,', ' 20190430', ' Note:none'}); 
Subject_num=Subject_Info(1);
ScanDate=Subject_Info(5);
Notes=Subject_Info(6);

Subject_num=strjoin(Subject_num);
ScanDate=strjoin(ScanDate);
Notes=strjoin(Notes);

Subject_Info=strjoin(Subject_Info);
%% load data
%% 
if strcmp(Original_analysis, 'yes')==1 && strcmp(N4_bias_analysis, 'no')==1 
    if strcmp(data_type, 'single_dcm')==1
       msgboxw('Please, select only one file...',14);
       [filename, parentPath] = uigetfile('.dcm','Select File');
       cd(parentPath)
       info=dicominfo([parentPath filename]);
       A = dicomread(info);
       A = double(squeeze(A));
       MR=A;
    else
        msgboxw('Please, select all files...',14);
        [MR, parentPath, FileNames] = DICOM_Load; % now allow for import of extensionless file names and list of .dcm
        info=dicominfo([parentPath char(FileNames(1))]);
    end
% parentPath = parentPath;
cd(parentPath)
% Normalize the intensity of the original image to fall between [0,1].
MR2=MR-min(MR(:));
MR2=MR/max(MR(:));
%% Segment lung parenchyma
%% 
fprintf( 'Segmeting Lung Parenchyma...\n' );
    for p = 1:size(MR2,3)
        scaledImage(:,:,p) = imadjust(MR2(:,:,p));
        scaledImage2(:,:,p) = adapthisteq(scaledImage(:,:,p));
        imshow(scaledImage2(:,:,p),[]);   
        HPmask(:,:,p) = roipoly(scaledImage2(:,:,p)); 
%         set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
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
                                  HPmask(:,:,p) = roipoly(scaledImage2(:,:,p));                                   
                                  break
                          end
                      end                  
        HPmask2(:,:,p) = roipoly(scaledImage2(:,:,p));
%         set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        choice = questdlg('Poly ROI acceptable?', ...
                          'Save Poly ROI', ...
                          'Yes','No','Yes');
                        while ~isequal(choice,'Yes')
                            switch choice
                                case 'Yes'
                                    disp(' Outline Next Lung')
                                    break
                                case 'No'
                                    disp(' Redraw ROI')
                                    HPmask2(:,:,p) = roipoly(scaledImage2(:,:,p));
                                    break
                           end
                       end
        mask(:,:,p) = imadd(HPmask(:,:,p),HPmask2(:,:,p));
        maskarray(:,:,p) = medFilter(mask(:,:,p));
        imshow(maskarray(:,:,p));
    end 
%     load('E:\Lab\VPD CODE\VDP_TEST1_AB\TEST1_VEN_Data\maskarray.mat');
end 
%% N4_Bias
%% 
if strcmp(N4_bias_analysis, 'yes')==1  && strcmp( Original_analysis, 'no')==1
    msgboxw('Please, Select VentilationN4 Nifti file...',14);
    [parentFile1,parentPath1] = uigetfile('*.nii.gz', 'Select Ventilation Nifti file');
    
    %read in the dicom images
    A1 = load_nii([parentPath1 parentFile1]); % Original
    A=A1.img;
    A = double(squeeze(A));
    I90=imrotate(A,-90);
    Ifv=flipdim(I90,1); % Original
    MR=Ifv;
    msgboxw('Please, Select Ventilation Mask Nifti file...',14);
    [parentFile2,parentPath2] = uigetfile('*.nii.gz', 'Select Ventilation Mask Nifti file');
    A2 = load_nii([parentPath2 parentFile2]); % Original
    A3=A2.img;
    A3 = double(squeeze(A3));
    I290=imrotate(A3,-90);
    Ifv2=flipdim(I290,1); % Original
    maskarray=Ifv2;
    parentPath=parentPath1;
    cd(parentPath);
    % Normalize the intensity of the original image to fall between [0,1].
    MR2=MR-min(MR(:));
    MR2=MR/max(MR(:));
end 
%% Create a folder
%% 
if strcmp(Original_analysis, 'yes')==1
    c_parentPath=parentPath;
    mkdir([c_parentPath 'Orignal Analysis']);
    sub_parentPath=[c_parentPath,'\Orignal Analysis\'];
else
     if strcmp(N4_bias_analysis, 'yes')==1
        c_parentPath=parentPath;
        n4_foldername=['N4 Bias Analysis_',num2str(Incomplete_threshold)];
        mkdir([c_parentPath n4_foldername]);
        sub_parentPath=[c_parentPath,'\',n4_foldername,'\'];
     end
end 

%% SNR 
%% 
if strcmp(calculate_SNR, 'yes')==1 
%Created a structuring element needed to dilate the Xe mask array
%The structuring element size '28' can be changed as needed to remove
%partial volume pixels and the trachea 
SE=strel('square',28);
%Dilates the mask arrary according to the structuring element
maskarraydilated=imdilate(maskarray,SE);
%This takes the complement of the dilated mask array and creates a background noise array 
backgroundmask=imcomplement(maskarraydilated);
background1=MR.*(backgroundmask);
%Any 0s in the background noise array are labeled as NaNs
background1(background1==0)=NaN;
%rf correction of the Xe intensity within the Xe vent maskarray
NormMR=rfCorrection(MR,maskarray);
NormMR2=MR.*(maskarray>0);
%Any 0s in the Xe vent array are labeled as NaNs 
NormMR2(NormMR2==0)=NaN;
%calculates the column size needed to reshape the background noise array for accurate StdDev calculation 
x=size(background1,1)*size(background1,2);
%Calculates the SNR for each slice 
for r = 1:size(maskarray,3)
   SNRarray(:,:,r)=((((nanmean(nanmean(NormMR2(:,:,r)))))/(nanstd(reshape(background1(:,:,r),x,1))))*sqrt(2-(pi/2)));
end 

%List the overall average SNR of the image set in the Matlab workspace 
Overall_SNRarray=nanmean(SNRarray);
%Generates a text file with the SNR of each slice
C=num2cell(SNRarray);
if strcmp(Original_analysis, 'yes')==1
   mkdir([sub_parentPath 'SNRvent_Original']);
   outputPath8 = ([sub_parentPath 'SNRvent_Original']);
   cd(sub_parentPath);
   filePh = fopen('SNR_Original.txt','w');
   fprintf(filePh,'%s\n',C{:});
   fclose(filePh);
elseif strcmp(N4_bias_analysis, 'yes')==1
        mkdir([sub_parentPath 'SNRvent_N4']);
        outputPath8 = ([sub_parentPath 'SNRvent_N4']);
        cd(sub_parentPath);
        filePh = fopen('SNR_N4.txt','w');
        fprintf(filePh,'%s\n',C{:});
        fclose(filePh);
end

fprintf( 'Writing out SNR, please wait...\n' );
%Generates Xe ventilation images with the SNR listed in the top lefthand cover of each image slice
figure
for m = 1:size(MR2,3)
      outputFileName8 = sprintf('%d.png',m);
      fullFileName = fullfile(outputPath8, outputFileName8);
      imshow(MR2(:,:,m),[]);
      hText = text(5,10, num2str(SNRarray(:,:,m)),'Color',[1 1 1],'FontSize',10);
      hFrame = getframe(gcf); %// Get content of figure
      imwrite(hFrame.cdata,fullFileName,'png');
      pause(0);
end
if strcmp(savedata,'yes') == 1
    if strcmp(Original_analysis, 'yes')==1
        mkdir([parentPath 'SNR_Background_Original']);
        outputPath8 = ([parentPath 'SNR_Background_Original']);
    elseif strcmp(N4_bias_analysis, 'yes')==1
        mkdir([parentPath 'SNR_Background_N4']);
        outputPath8 = ([parentPath 'SNR_Background_N4']);
    end
end

%Generates background noise images with the SNR listed in the top lefthand cover of each image slice
figure
for m = 1:size(background1,3)
      outputFileName8 = sprintf('%d.png',m);
      fullFileName = fullfile(outputPath8, outputFileName8);
      imshow(background1(:,:,m),[]);
      hText = text(5,10, num2str(SNRarray(:,:,m)),'Color',[1 1 1],'FontSize',10);
      hFrame = getframe(gcf); %// Get content of figure
      imwrite(hFrame.cdata,fullFileName,'png');
      pause(1);
end
close all
end
%% Eroding edge of masks
% for p = 1:size(maskarray,3)
%     se = strel('square',3);
%     maskarray2(:,:,p) = imerode(maskarray(:,:,p),se);
% end
% figure; imslice(MR)
%% Write out ventilation images and masks as png files
%% 
fprintf( 'Writing out ventilation images, please wait...\n' );
cd(sub_parentPath)
mkdir([sub_parentPath 'scaledVent']);
outputPath = ([sub_parentPath 'scaledVent']);

%will a adjust the intensity values in grayscale image original to new
%values in processed such that 1% of data is saturated at low and high
%intensities of the original image and print out png images of the
%processed image to parentPath\masks folder
for m = 1:size(MR2,3)
      scaledImage(:,:,m) = imadjust(MR2(:,:,m));
      scaledImage2(:,:,m) = adapthisteq(scaledImage(:,:,m));
      outputFileName = sprintf('Img_%d.png',m);
      fullFileName = fullfile(outputPath, outputFileName);
      imwrite(scaledImage2(:,:,m),fullFileName,'png');
      imshow(scaledImage2(:,:,m),[]);
      title(sprintf('Processed Image # %d',m));
      pause(1);
end

fprintf( 'Writing out mask, please wait...\n' );
%create a folder labeled 'masks' in the parentPath folder to output ventilation images saved as png files 
mkdir([sub_parentPath 'masks']);
outputPath = ([sub_parentPath 'masks']);
%write out png files of created masks for validation of Matlab code vs Robby's R code 
for m = 1:size(maskarray,3)
      outputFileName = sprintf('%d.png',m);
      fullFileName = fullfile(outputPath, outputFileName);
      imwrite(maskarray(:,:,m),fullFileName,'png');
      imshow(maskarray(:,:,m),[]);
      title(sprintf('Processed Mask # %d',m));
      pause(1);
end
close; 
%% Applies rf correction and median filter as well as quantifying VDP
%% 
%Applies a median filter(medFilter) to mask, next it will calculate rf correction and apply it to the
%MR images and maskarray to produce a corrected and segmented NormMR array
if strcmp(RFCorrection, 'yes')==1 && strcmp( Original_analysis, 'yes')==1
    NormMR=rfCorrection(MR,maskarray);
    NormMR2=NormMR.*(maskarray>0);
else
    NormMR2=MR.*(maskarray>0);
end

if strcmp(Median_Filter, 'yes')==1
    %calculate the defect ventilation percentage with median filter
    Incomplete_defect = medFilter(double(NormMR2<(mean(NormMR2(NormMR2>0))*(Incomplete_threshold/100))*(maskarray>0)));
    Complete_defect = medFilter(double(NormMR2<(mean(NormMR2(NormMR2>0))*(Complete_threshold/100))*(maskarray>0)));
    Hyper_defect = medFilter(double(NormMR2>(mean(NormMR2(NormMR2>0))*2)*(maskarray>0)));
    defectArray = Incomplete_defect+Complete_defect+4*Hyper_defect;
else
    %calculate the defect ventilation percentage without median filter
    Incomplete_defect = NormMR2<(mean(NormMR2(NormMR2>0))*Incomplete_threshold)*(maskarray>0);
    Complete_defect = NormMR2<(mean(NormMR2(NormMR2>0))*Complete_threshold)*(maskarray>0);
    Hyper_defect = NormMR2>(mean(NormMR2(NormMR2>0))*2)*(maskarray>0);
    defectArray = Incomplete_defect+Complete_defect+4*Hyper_defect;
end
%bin defect signal for histogram
d_Incomplete = NormMR2.*(defectArray==1); %Incomplete
d_Complete = NormMR2.*(defectArray==2);%Complete
d_Hyper = NormMR2.*(defectArray==4);%Hyper
d1 = (NormMR2.*(defectArray==0));%Normal

%Calculating precent of defects
VDP=100*(sum(sum(sum(defectArray==1)))+(sum(sum(sum(defectArray==2)))))/sum(sum(sum(maskarray>0)));
Incomplete= (100*(sum(sum(sum(defectArray==1)))/sum(sum(sum(maskarray>0)))));
Complete = 100*(sum(sum(sum(defectArray==2)))/sum(sum(sum(maskarray>0))));
Hyper = 100*(sum(sum(sum(defectArray==4)))/sum(sum(sum(maskarray>0))));
Normal = 100*((sum(sum(sum(maskarray>0))))-(sum(sum(sum(defectArray==1))))-(sum(sum(sum(defectArray==2))))-(sum(sum(sum(defectArray==4)))))/(sum(sum(sum(maskarray>0))));

%Find mean and standard deviation of defects
md1 = mean(d1(d1>0));
sd1 = std(d1(d1>0));
md_Incomplete = mean(d_Incomplete(d_Incomplete>0));
sd_Incomplete = std(d_Incomplete(d_Incomplete>0));
md_Complete = mean(d_Complete(d_Complete>0));
sd_Complete = std(d_Complete(d_Complete>0));
md_Hyper = mean(d_Hyper(d_Hyper>0));
sd_Hyper = std(d_Hyper(d_Hyper>0));

%equalizing histogram bin counts for each defect
acnts = 50; 
a = hist(d1(d1>0),acnts);
bcnts = length(a)*(max(d_Incomplete(d_Incomplete>0))-min(d_Incomplete(d_Incomplete>0)))/(max(d1(d1>0))-min(d1(d1>0)));
ccnts = length(a)*(max(d_Complete(d_Complete>0))-min(d_Complete(d_Complete>0)))/(max(d1(d1>0))-min(d1(d1>0)));
dcnts = length(a)*(max(d_Hyper(d_Hyper>0))-min(d_Hyper(d_Hyper>0)))/(max(d1(d1>0))-min(d1(d1>0)));
%% Output histogran and montages of ventilation, mask, and VDP images 
%% Create histogram figure of ventilation defects and print/save figure
figure('position',[350 350 750 350]);
if strcmp(Original_analysis, 'yes')==1
    hist(d1(d1>0),acnts), xlim([0 10]);
    hold on
    hist(d_Incomplete(d_Incomplete>0),bcnts);
    hist(d_Complete(d_Complete>0),ccnts);
    hist(d_Hyper(d_Hyper>0),dcnts);
else
    hist(d1(d1>0),acnts);
    hold on
    hist(d_Incomplete(d_Incomplete>0),bcnts);
    hist(d_Complete(d_Complete>0));
    hist(d_Hyper(d_Hyper>0),dcnts);
end

h = findobj(gcf,'Type','patch'); 
set(h(1),'FaceColor','r','EdgeColor','k','facealpha',0.5);
set(h(2),'FaceColor','b','EdgeColor','k','facealpha',0.5);
set(h(3),'FaceColor','y','EdgeColor','k','facealpha',0.5);
set(h(4),'FaceColor','w','EdgeColor','k','facealpha',0.5);
legend1 = sprintf('Normal: %0.2f±%0.2f (%0.1f%%)',md1, sd1, Normal);
legend2 = sprintf('Incomplete: %0.2f±%0.2f (%0.1f%%)',md_Incomplete,sd_Incomplete,Incomplete);
legend3 = sprintf('Complete: %0.2f±%0.2f (%0.1f%%)',md_Complete,sd_Complete, Complete);
legend4 = sprintf('Hyper: %0.2f±%0.2f (%0.1f%%)',md_Hyper,sd_Hyper, Hyper);
legend({legend1 legend2 legend3 legend4});
% title1 = sprintf('Ventilation Defect Percentage %0.1f%%',VDP,'\n',Subject_Info);
title1=sprintf(Subject_Info);
title2 = sprintf('\nVentilation Defect Percentage %0.1f%%',VDP);
titlemain=[title1 ,title2,'%'];
titlemain=sprintf(titlemain);
title({titlemain},'Fontweight','bold','FontSize',12);
set(gcf,'PaperPosition',[0 0 7.5 3.5]);
print('Ventilation_Histogram','-dpng','-r300');
close all
%% Creating color maps for ventilation defects
%% 
d_Incomplete_rgb=(defectArray==1)*255;
d_Complete_rgb=(defectArray==2)*255;
d_Hyper_rgb=(defectArray==4)*255;
d1rgb=(defectArray==0)*255;

cm60=[0 0 0
      1 1 0];
array4d60 = zeros(size(d_Incomplete_rgb,1), size(d_Incomplete_rgb,2),3, size(d_Incomplete_rgb,3));
for slice = 1:size(d_Incomplete_rgb,3)
    array4d60(:,:,:,slice) =  ind2rgb(d_Incomplete_rgb(:,:,slice),cm60);
end
%imshow(array4d60(:,:,:,6),[]);

cm15=[0 0 0
      0 0 1];
array4d15 = zeros(size(d_Complete_rgb,1), size(d_Complete_rgb,2),3, size(d_Complete_rgb,3));
for slice = 1:size(d_Complete_rgb,3)
    array4d15(:,:,:,slice) =  ind2rgb(d_Complete_rgb(:,:,slice),cm15);
end
%imshow(array4d15(:,:,:,6),[]);

cm250=[0 0 0
      1 0 0];
array4d250 = zeros(size(d_Hyper_rgb,1), size(d_Hyper_rgb,2),3, size(d_Hyper_rgb,3));
for slice = 1:size(d_Hyper_rgb,3)
    array4d250(:,:,:,slice) =  ind2rgb(d_Hyper_rgb(:,:,slice),cm250);
end
%imshow(array4d250(:,:,:,6),[]);

cm1=[0 0 0
     1 1 1];
array4d1 = zeros(size(d1rgb,1), size(d1rgb,2),3, size(d1rgb,3));
for slice = 1:size(d1rgb,3)
    array4d1(:,:,:,slice) =  ind2rgb(d1rgb(:,:,slice),cm1);
end
%imshow(array4d1(:,:,:,6),[]);
%% Turn 3D array of images into 4D array for overlaying colored ventilation defects on ventilation images
array4dMRrgb = zeros(size(scaledImage2,1), size(scaledImage2,2), 3,size(scaledImage2,3));
array4dMRrgb(:,:,1,:)=scaledImage2;
array4dMRrgb(:,:,2,:)=scaledImage2;
array4dMRrgb(:,:,3,:)=scaledImage2;
figure, imshow(array4dMRrgb(:,:,:,3),[]);

%Turn 3D array of masks into 4D array to eliminate the background from the
%normal ventilation defect array mask
array4dMaskrgb = zeros(size(maskarray,1), size(maskarray,2), 3,size(maskarray,3));
array4dMaskrgb(:,:,1,:)=maskarray;
array4dMaskrgb(:,:,2,:)=maskarray;
array4dMaskrgb(:,:,3,:)=maskarray;

%Eliminate background from 'normal' ventilation defect array mask and add all defect array color
%masks into one array
array4dMRrgb2=array4d1.*array4dMRrgb;
defectArray_rgb2 = array4d60+array4d15+4*array4d250;
midslice=round((size(defectArray_rgb2,4)/2));
array4dMRrgb3=defectArray_rgb2+array4dMRrgb2;
array4dMRrgb4=uint8(255*array4dMRrgb3);
figure,imshow(array4dMRrgb4(:,:,:,midslice),[]);

fprintf( 'display ventilation image slice with an overlay of colored ventilation masks, please wait...\n' );
%Will create a folder for ventilation defect images, display ventilation image slice 
%with an overlay of colored ventilation masks and print out overlayed images 
mkdir([sub_parentPath 'ventDefect']);
outputPath2 = ([sub_parentPath 'ventDefect']);

for m = 1:size(array4dMRrgb4,4)
      outputFileName = sprintf('%d.png',m);
      fullFileName = fullfile(outputPath2, outputFileName);
      imwrite(array4dMRrgb4(:,:,:,m),fullFileName,'png');
      imshow(array4dMRrgb4(:,:,:,m),[]);
      title(sprintf('Processed Mask # %d',m));
      pause(1);
end
close all;
%% Add noise back into the segmented images
%% 
maskarray3=maskarray*255;
cmap=[0 0 0
      1 1 0];
array4dcoloredmask = zeros(size(maskarray3,1), size(maskarray3,2),3, size(maskarray3,3));
for slice = 1:size(maskarray3,3)
    array4dcoloredmask(:,:,:,slice) =  ind2rgb(maskarray3(:,:,slice),cmap);
end

cmap2=[1 1 1
       0 0 0];
array4dcoloredmask2 = zeros(size(maskarray3,1), size(maskarray3,2),3, size(maskarray3,3));
for slice = 1:size(maskarray3,3)
    array4dcoloredmask2(:,:,:,slice) =  ind2rgb(maskarray3(:,:,slice),cmap2);
end

array4dcoloredmask3=array4dMRrgb.*array4dcoloredmask;
array4dcoloredmask4=array4dMRrgb.*array4dcoloredmask2;
array4dcoloredmask5=array4dcoloredmask3+array4dcoloredmask4;
imshow(array4dcoloredmask5(:,:,:,6),[])

mkdir([sub_parentPath 'colorMasks'])
outputPath4 = ([sub_parentPath 'colorMasks']);

fprintf( 'write out png files of created masks for validation, please wait...\n' );
%write out png files of created masks for validation of Matlab code vs Robby's R code 
for m = 1:size(array4dcoloredmask5,4)
      outputFileName4 = sprintf('%d.png',m);
      fullFileName = fullfile(outputPath4, outputFileName4);
      imwrite(array4dcoloredmask5(:,:,:,m),fullFileName,'png');
      imshow(array4dcoloredmask5(:,:,:,m),[]);
      title(sprintf('Processed Mask # %d',m));
      pause(1);
end
close all;
fprintf( 'Creating and displaying montages of ventilation images, ventilation defects, and mask arrays, please wait...\n' );
%Creating and displaying montages of ventilation images, ventilation defects, and mask arrays
figure
subplot(4,1,1)
HPXeMontage = montage(reshape(scaledImage,[size(scaledImage,1), size(scaledImage,2), 1, size(scaledImage,3)]),'Size',[1 size(scaledImage,3)],'DisplayRange',[]);
title('Xenon Montage')
savedMontage=get(HPXeMontage,'CData');
savedMontage2 = imadjust(savedMontage);
outputFileName2= 'HPXeMontage.png';
fullFileName2 = fullfile(sub_parentPath,outputFileName2);
imwrite(savedMontage2,fullFileName2,'png');
subplot(4,1,2)
MaskMontage=montage(reshape(maskarray,[size(maskarray,1), size(maskarray,2), 1, size(maskarray,3)]),'Size',[1 size(maskarray,3)],'DisplayRange',[]);
title('Mask Array')
savedMaskMontage=get(MaskMontage,'CData');
savedMaskMontage2 = imadjust(savedMaskMontage);
outputFileName2= 'MaskMontage.png';
fullFileName2 = fullfile(sub_parentPath,outputFileName2);
imwrite(savedMaskMontage2,fullFileName2,'png');
subplot(4,1,3)
DefectArrayMontage=montage(reshape(array4dMRrgb4,[size(array4dMRrgb4,1), size(array4dMRrgb4,2),size(array4dMRrgb4,3),size(array4dMRrgb4,4)]),'Size',[1 size(array4dMRrgb4,4)],'DisplayRange',[]);
title('Defect Array')
savedDefectArrayMontage=get(DefectArrayMontage,'CData');
outputFileName2= 'DefectArrayMontage.png';
fullFileName2 = fullfile(sub_parentPath,outputFileName2);
imwrite(savedDefectArrayMontage,fullFileName2,'png');
subplot(4,1,4)
ColoredMaskMontage=montage(reshape(array4dcoloredmask5,[size(array4dcoloredmask5,1), size(array4dcoloredmask5,2),size(array4dcoloredmask5,3),size(array4dcoloredmask5,4)]),'Size',[1 size(array4dcoloredmask5,4)],'DisplayRange',[]);
title('Colored Mask')
savedColoredMaskMontage=get(ColoredMaskMontage,'CData');
outputFileName4= 'ColoredMaskMontage.png';
fullFileName4 = fullfile(sub_parentPath,outputFileName4);
imwrite(savedColoredMaskMontage,fullFileName4,'png');
close all;


%% save result in a powerpoint file // Abdullah 1/27/2020
%% 
% %how many slices to incluse in the final powerpoint file
% sl_include=zeros(1,size(maskarray,3));
% for sl=1:size(maskarray,3)
%     choose_slice=maskarray(:,:,sl);       
%     logica_test2=sum(choose_slice(:));
%     if logica_test2>1
%          sl_include(sl)=sl;
%     else
%          sl_include(sl)=0;
%     end 
% end
% sl_include( :, ~any(sl_include,1) ) = [];  %columns
% sl_1=sl_include(1);
% sl_end=sl_include(end);
% if you want the un-segemted slices not to show up on the powepoint, run
% the above code and add this line (,'Indices', sl_1:sl_end) to the end of
% each Montage. 
%%% get rid of the gray frame in each montage
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(reshape(scaledImage,[size(scaledImage,1), size(scaledImage,2), 1, size(scaledImage,3)]),...
    'Size',[1 size(scaledImage,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

MaskMontage = figure('Name','Mask');set(MaskMontage,'WindowState','minimized');
montage(reshape(maskarray,[size(maskarray,1), size(maskarray,2), 1, size(maskarray,3)]),...
    'Size',[1 size(maskarray,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

DefectArrayMontage = figure('Name','Defects');set(DefectArrayMontage,'WindowState','minimized');
montage(reshape(array4dMRrgb4,[size(array4dMRrgb4,1), size(array4dMRrgb4,2),size(array4dMRrgb4,3),...
    size(array4dMRrgb4,4)]),'Size',[1 size(array4dMRrgb4,4)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

% ColoredMaskMontage = figure('Name','ColoredMask');set(ColoredMaskMontage,'WindowState','minimized');
% montage(reshape(array4dcoloredmask5,[size(array4dcoloredmask5,1), size(array4dcoloredmask5,2),...
%   size(array4dcoloredmask5,3),size(array4dcoloredmask5,4)]),'Size',[1 size(array4dcoloredmask5,4)],'DisplayRange',[]);
% set(gca,'units','pixels'); % set the axes units to pixels
% x = get(gca,'position'); % get the position of the axes
% set(gcf,'units','pixels'); % set the figure units to pixels
% y = get(gcf,'position'); % get the figure position
% set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
% set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
ReportTitle = [Subject_num,'_',ScanDate,'_','_Final_VDP_Summary'];
if strcmp(Original_analysis, 'yes')==1
    cd(parentPath);
    %Start new presentation
    isOpen  = exportToPPTX();
    if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
    %Add slides
    exportToPPTX('addslide'); % Image/mask/VDP
    exportToPPTX('addtext',Subject_Info,'Position',[4 0.1 10 1],'Color','k','FontSize',30);
    exportToPPTX('addpicture',VentMontage,'Position',[0 1 16 1.5]);
    exportToPPTX('addpicture',MaskMontage,'Position',[0 2.3 16 1.5]);
    exportToPPTX('addpicture',DefectArrayMontage,'Position',[0 3.5 16 1.5]);
    % exportToPPTX('addpicture',ColoredMaskMontage,'Position',[0 3.8 16 1.5]);
    exportToPPTX('addtext','Vent Images','Position',[0 1 1.6 0.4],'Color','r','FontSize',18,'BackgroundColor','k');
    exportToPPTX('addtext','Mask','Position',[0 2.3 1 0.4],'Color','r','FontSize',18,'BackgroundColor','k');
    vdp_title = sprintf('Original VDP= %0.1f%%',VDP);
    exportToPPTX('addtext',vdp_title,'Position',[0 3.5 3 0.4],'Color','r','FontSize',18,'BackgroundColor','k');
    
    exportToPPTX('addslide'); %Histogram
    cd(sub_parentPath);
    Original_hist=imread('Ventilation_Histogram.png');    
    exportToPPTX('addpicture',Original_hist,'Position',[0 2.5 5.5 3]);
    org_hist_title=['Original %',num2str(Incomplete_threshold)];
    exportToPPTX('addtext',org_hist_title,'Position',[2 5.5 3 1]);
    cd(parentPath);
    exportToPPTX('save',fullfile(parentPath, ReportTitle));
    exportToPPTX('close');
    fprintf('PowerPoint file has been saved\n');
elseif strcmp(N4_bias_analysis, 'yes')==1
        cd(parentPath);
        ppt_file_name=[ReportTitle,'.pptx'];
        if Incomplete_threshold==60 && exist(ppt_file_name, 'file') == 2            
            %Start new presentation
            isOpen  = exportToPPTX();
            if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
                exportToPPTX('close');
            end
            exportToPPTX('open',ppt_file_name);
            exportToPPTX('switchslide',1);
            exportToPPTX('addpicture',DefectArrayMontage,'Position',[0 5 16 1.5]);
            % exportToPPTX('addpicture',ColoredMaskMontage,'Position',[0 3.8 16 1.5]);
            vdp_title = sprintf('N4 Defects_60 VDP= %0.1f%%',VDP);
            exportToPPTX('addtext',vdp_title,'Position',[0 5 3 0.4],'Color','r','FontSize',18 ,'BackgroundColor','k');
            
            exportToPPTX('switchslide',2); %Histogram
            cd(sub_parentPath);
            N4_hist=imread('Ventilation_Histogram.png');           
            exportToPPTX('addpicture',N4_hist,'Position',[5 2.5 5.5 3]);
            org_hist_title=['N4 %',num2str(Incomplete_threshold)];
            exportToPPTX('addtext',org_hist_title,'Position',[7 5.5 3 1]);
            cd(parentPath);
            exportToPPTX('save',fullfile(parentPath, ReportTitle));
            exportToPPTX('close');
            fprintf('PowerPoint file has been updated\n');
        elseif Incomplete_threshold==60
                     %Start new presentation
            isOpen  = exportToPPTX();
            if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
                exportToPPTX('close');
            end
            exportToPPTX('new','Dimensions',[16 9], ...
                'Title',ReportTitle, ...
                'Author','CPIR @ CCHMC');
            %Add slides
            exportToPPTX('addslide'); % Image/mask/VDP
            exportToPPTX('addtext',Subject_Info,'Position',[4 0.1 10 1],'Color','k','FontSize',30);
            exportToPPTX('addpicture',VentMontage,'Position',[0 1 16 1.5]);
            exportToPPTX('addpicture',MaskMontage,'Position',[0 2.3 16 1.5]);
            exportToPPTX('addpicture',DefectArrayMontage,'Position',[0 3.5 16 1.5]);
            exportToPPTX('addtext','Vent Images','Position',[0 1 1.6 0.4],'Color','r','FontSize',18,'BackgroundColor','k');
            exportToPPTX('addtext','Mask','Position',[0 2.2 1 0.4],'Color','r','FontSize',18,'BackgroundColor','k');
            vdp_title = sprintf('N4 Defects_60 VDP= %0.1f%%',VDP);
            exportToPPTX('addtext',vdp_title,'Position',[0 3.5 3 0.4],'Color','r','FontSize',18,'BackgroundColor','k');

            exportToPPTX('addslide'); %Histogram
            cd(sub_parentPath);
            N4_hist=imread('Ventilation_Histogram.png');           
            exportToPPTX('addpicture',N4_hist,'Position',[5 2.5 5.5 3]);
            org_hist_title=['N4 %',num2str(Incomplete_threshold)];
            exportToPPTX('addtext',org_hist_title,'Position',[7 5.5 3 1]);
            cd(parentPath);
            exportToPPTX('save',fullfile(parentPath, ReportTitle));
            exportToPPTX('close');
            fprintf('PowerPoint file has been saved\n');            
        elseif Incomplete_threshold==75 && exist(ppt_file_name, 'file') == 2
                %Start new presentation
                isOpen  = exportToPPTX();
                if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
                    exportToPPTX('close');
                end
                exportToPPTX('open',ppt_file_name);

                exportToPPTX('switchslide',1);
                exportToPPTX('addpicture',DefectArrayMontage,'Position',[0 6.5 16 1.5]);
                vdp_title = sprintf('N4 Defects_75 VDP= %0.1f%%',VDP);
                exportToPPTX('addtext',vdp_title,'Position',[0 6.5 3 0.4],'Color','r','FontSize',18 ,'BackgroundColor','k');

                exportToPPTX('switchslide',2); %Histogram
                cd(sub_parentPath);
                N4_hist=imread('Ventilation_Histogram.png');           
                exportToPPTX('addpicture',N4_hist,'Position',[10 2.5 5.5 3]);
                org_hist_title=['N4 %',num2str(Incomplete_threshold)];
                exportToPPTX('addtext',org_hist_title,'Position',[12 5.5 5 1]);
                cd(parentPath);
                exportToPPTX('save',fullfile(parentPath, ReportTitle));
                exportToPPTX('close');
                fprintf('PowerPoint file has been updated\n');
        elseif Incomplete_threshold==75
                         %Start new presentation
                isOpen  = exportToPPTX();
                if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
                    exportToPPTX('close');
                end
                exportToPPTX('new','Dimensions',[16 9], ...
                    'Title',ReportTitle, ...
                    'Author','CPIR @ CCHMC');
                %Add slides
                exportToPPTX('addslide'); % Image/mask/VDP
                exportToPPTX('addtext',Subject_Info,'Position',[4 0.1 10 1],'Color','k','FontSize',30);
                exportToPPTX('addpicture',VentMontage,'Position',[0 1 16 1.5]);
                exportToPPTX('addpicture',MaskMontage,'Position',[0 2.3 16 1.5]);
                exportToPPTX('addpicture',DefectArrayMontage,'Position',[0 3.5 16 1.5]);
                exportToPPTX('addtext','Vent Images','Position',[0 1 1.6 0.4],'Color','r','FontSize',18 ,'BackgroundColor','k');
                exportToPPTX('addtext','Mask','Position',[0 2.2 1 0.4],'Color','r','FontSize',18 ,'BackgroundColor','k');
                vdp_title = sprintf('N4 Defects_75 VDP= %0.1f%%',VDP);
                exportToPPTX('addtext',vdp_title,'Position',[0 3.5 3 0.4],'Color','r','FontSize',18,'BackgroundColor','k');

                exportToPPTX('addslide'); %Histogram
                cd(sub_parentPath);
                N4_hist=imread('Ventilation_Histogram.png');           
                exportToPPTX('addpicture',N4_hist,'Position',[10 2.5 5.5 3]);
                org_hist_title=['N4 %',num2str(Incomplete_threshold)];
                exportToPPTX('addtext',org_hist_title,'Position',[12 5.5 5 1]);
                cd(parentPath);
                exportToPPTX('save',fullfile(parentPath, ReportTitle));
                exportToPPTX('close');
                fprintf('PowerPoint file has been saved\n');               
         end
end
close all;

%% save mask 
if strcmp(Original_analysis, 'yes')==1
    [a,b,c] = size(maskarray); dicomwrite(reshape(maskarray,[a,b,1,c]),'maskarray.dcm', info, 'CreateMode', 'Copy','MultiframeSingleFile', true);
end 
%% save data 
if strcmp(savedata,'yes') == 1
    if strcmp(Original_analysis, 'yes')==1
        cd(sub_parentPath);
        save 'HP_Xe_gas_VDP_Orignal_analysis.mat';
    elseif strcmp(N4_bias_analysis, 'yes')==1
        cd(sub_parentPath);
        save 'HP_Xe_gas_VDP_N4_analysis.mat';
    end
end
%% All Analysis is done
fprintf( 'Analysis Is Done...\n' );
msgboxw('Analysis Is Done...',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end VDP Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%

