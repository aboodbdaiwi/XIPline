
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matlab script for processing Philips raw data into ADC maps %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10/10/2017 : Written by Zackary I. Cleveland (ZIC) and Abdullah S. Bdaiwi
% 10/11/2017 : ZIC, airway segmentation added and pathway bug fixed.
% 10/12/2017 : ZIC, edited to allow for fitting of selected b values
% 11/15/2017 : AB,  add 2D filters 
% 11/15/2017 : AB,  edited to add SNR and weighting vector calculation 
% 11/15/2017 : AB,  edited to add weighting to the Polynomial fit to calculate the ADC 
% 11/29/2017 : AB,  edited to add weighting to the Correlation coefficients fit to calculate the Rmap 
% 05/21/2018 : AB,  edited to add Bayesain fitting
% 02/04/2020 : AB,  add maxlikelihood fit
%
% Functions Called:     load_philips_extr1_2D.m
%                       ADC_thresh.m
%                       ADC_Airway_Segment.m
%                       fermi_filter_2D_AB.m
%                       Gaussian_filter_2D_AB.m
%                       hann_filter_2D_AB.m
%                       matbugs.m
%                       wpolyfit.m or polyfitweighted.m    
%                       weightedcorrs.m                    
%                       DICOM_Load.m   
%                       fmaxlikelihood_ADC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for any questions, please contact Abdullah Bdaiwi at Abdullah.Bdaiwi@cchmc.org

%% Things you need to do before you hit run
%% 
clc;clear;close all;% start with fresh variables, etc.
tic
subject_info='yes';
savedata = 'yes'; % if yes, save workspace to current directory
filterim = 'Fermi'; % choose 'Fermi','Gaussian', 'Han' or 'None'
segment = 'theshold';  % choose 'manual' to segment lung parenchyma manually otherwise will use threshold
airway_segment = 'manual';  % segment airways mannually (default) or "quicknEZ"  with 
                            % a quick and easy signal ratio involving
                            % awcof. Manual is the default   
Linear_fit='yes'; % if yes,  linear fit will be performed
W_Linear_fit = 'yes'; % if yes, weighted linear fit will be performed
ML_fit='yes'; % if yes,  maxlikelihood fit will be performed
Bayesain_fit = 'yes'; % if yes, Bayesain fit will be performed
%% User imput needed for recon and processing
%% 
awcof = 15; % cuttoff factor for removing airways using b0/bmax ratio
ch_range = [1 1]; % assume one channel for data import by default
Wf = 10; % GE Heathcare default Fermi Filter
SE = 1; % structuring element size for erode/dilate
nzcof = 15; % noise floor cuttoff factor for setting threshold
nbValues=5;
max_bValues=25;
bValues = linspace(0,max_bValues,nbValues);  % default to xenon
bValues=bValues';
% bval = length(bValues)-1; % index of 2nd bvalue selected for 2 bvalue fit
bs2fit = 1:length(bValues); % b values used in ADC fitting
%If you want to change it to fit only the 1st, 3rd, and 4th bvalues, simply edit to
% bs2fit = [1 2 3]; % b values used in ADC fitting

% R2_threshold=.1; %R-squared threshold
%% Subject Info.
%% 
if strcmp(subject_info, 'yes')==1
    Subject_Info = inputdlg({'Subject No.','Study Group', 'Age', 'Gender', 'Visit Date'},...
              'Subject_Info', [1 30; 1 30; 1 30; 1 30; 1 30],{'IRC186H-174','-Control-','20 yr-','M-', '04/30/19'}); 
        Subject_num=Subject_Info(1);
        Subject_num=strjoin(Subject_num);
        Subject_Info=strjoin(Subject_Info);
else
    Subject_Info='No Subject Info';
end


%% Read in Data
%% 
home_path = pwd;
data_type = 'multi_dcm';
[filename, path] = uigetfile('.data','Select File');
% [filename, path] = uigetfile('Select File');
[filepath,name,ext] = fileparts(filename);
if  strcmp(ext, '.data')   
    cd(path)
    filename(end-4:end)=[];
    [ImgK,NoiK,kx_oversample_factor] = load_philips_extr1_2D(filename,ch_range);
    diffK = ImgK(:,:,:,:);
    cd(home_path);

    data_size = size(diffK);
    ky = data_size(1);
    kx =data_size(2);
    slices = data_size(3);
    nbs = data_size(4);
    diffimgcmplx = zeros(ky,kx/kx_oversample_factor,slices, nbs);
else
    if strcmp(data_type, 'single_dcm')==1
         cd(path)
         info=dicominfo([path filename]);
         A = dicomread(info);
         A = double(squeeze(A));
         diffimg_dicom=A;
    else
         diffimg_dicom = DICOM_Load; % now allow for import of extensionless file names and list of .dcm
    end
end
%% Reconstruct Data
%% image filter
if strcmp(filterim,'Fermi') == 1
        filter= fermi_filter_2D_AB(ky,kx,10);
elseif strcmp(filterim,'Gaussian') == 1
        filter= Gaussian_filter_2D_AB(ky,kx,1.5);
elseif strcmp(filterim,'Han') == 1
        filter= hann_filter_2D_AB(ky,kx);
else
        filter = 1;
end

for sl =1:slices
    for b = 1:nbs
        data_slice = diffK(:,:,sl,b).*filter;
        recon_slice = fftshift(fft2(data_slice),2);
        recon_slice(:,1:kx/(2*kx_oversample_factor))=[]; % crop extra pixels from oversample factor
        recon_slice(:, (1+kx/kx_oversample_factor):end)=[]; 
        diffimgcmplx(:,:,sl,b) = recon_slice;
    end
end
diffimg = abs(diffimgcmplx); % make magnitude data
diffimg = rot90(rot90(diffimg));

% figure; imslice(Ndiffimg)
%% Segment lung parenchyma and airways
%% 
if strcmp(segment,'manual') == 1
        lung_mask =ADC_lung_parenchyma_Segment(airway_segment,awcof,diffimg);
else
        lung_mask = ADC_thresh(NoiK,diffimg,SE,nzcof);
end
%Segment the airway 
airway_mask = ADC_Airway_Segment(airway_segment,awcof,diffimg);
% find the final mask 
final_mask =lung_mask;
final_mask(airway_mask==1)=0;

% flip the binary mask to create the noise mask(1 to 0, 0 to 1)
noise_mask = ((~lung_mask)+(~airway_mask))-1;
noise_mask=noise_mask>0;
% if there's artifacts, segemnt them. 
% artifact_mask = Segment_artifact(airway_segment,awcof,diffimg); 
% noise_mask(artifact_mask==1)=0;
%% get rid of any slice that has no signal (un-segmented)
slice_checker=zeros(1,slices);
for check_slice=1:slices
    sl=final_mask(:,:,check_slice);
    if sum(sl(:))>1
        slice_checker(check_slice)=1;
    elseif sum(sl(:))==0
        slice_checker(check_slice)=0;
    end
end
find_slice_checker=find(slice_checker);
Ndiffimg= zeros(ky,kx/kx_oversample_factor,length(find_slice_checker), nbs);
Nfinal_mask= zeros(ky,kx/kx_oversample_factor,length(find_slice_checker));
Nnoise_mask= zeros(ky,kx/kx_oversample_factor,length(find_slice_checker));
for find_sl=1:length(find_slice_checker)
    choose_sl=find_slice_checker(find_sl);
    Ndiffimg(:,:,find_sl,:)=diffimg(:,:,choose_sl,:);
    Nfinal_mask(:,:,find_sl)=final_mask(:,:,choose_sl);
    Nnoise_mask(:,:,find_sl)=noise_mask(:,:,choose_sl);
end
% figure; imslice(Ndiffimg)
% figure; imslice(Nfinal_mask)
diffimg=Ndiffimg;
final_mask=Nfinal_mask;
noise_mask=Nnoise_mask;
slices=length(find_slice_checker);
%% Normalize the intensity of the original image to fall between [0,100].
Ndiffimg=diffimg-min(diffimg(:));
Ndiffimg=diffimg/max(diffimg(:));
Ndiffimg=Ndiffimg*100;

final_mask_size=size(final_mask);
%convert 4D to 2D array (lenghth pixel inside the final mask and number of b values
Data_vec=reshape(Ndiffimg.*final_mask,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),length(bValues)]);
Data_vec = Data_vec(any(Data_vec,2),:);
pixel_location=find(final_mask); %to find pixel location to convert the array back to 3D image
num_ones=sum(final_mask(:)); % for loopping 
%% Write out Diffusion images and masks as png files
%% 
% cd(path)
% mkdir([path 'Diffusion Images']);
% outputPath = ([path 'Diffusion Images']);
% cd(outputPath);
% fprintf( 'Write out Diffusion images...\n' );
% %The following will print out each individual ADC map for each slice as a png file with 600 dpi resolution 
% for xi = 1:size(Ndiffimg,3)
%     for yi= 1:size(Ndiffimg,4)
%       outputFileName = sprintf('Img-slice%d-b%d.png',xi,yi);
%       imagesc(Ndiffimg(:, :, xi,yi)); colormap(gray); colorbar; caxis([0 100]);
%       title(outputFileName);
%       print(outputFileName,'-r150','-dpng');
%       pause(1);
%       
%     end
% end
% %write out png files of created masks for Diffusion 
% fprintf( 'Write out Mask...\n' );
% mkdir([path 'Masks']);
% outputPath = ([path 'Masks']);
% cd(outputPath);
% 
% for m = 1:size(Ndiffimg,3)
%       outputFileName = sprintf('slice%d.png',m);
%       imagesc(final_mask(:, :, m));  colormap(gray);
%       title(outputFileName);
%       print(outputFileName,'-r150','-dpng');
%       pause(1);
% end
% close all;
%%  Write out Diffusion images and masks
%% 
cd(path)
mkdir([path 'Diffusion Images']);
outputPath = ([path 'Diffusion Images']);

for m = 1:size(Ndiffimg,3)
    for n=1:size(Ndiffimg,4)
        chosen_img=uint16(Ndiffimg(:,:,m,n));
        scaledImage = imadjust(chosen_img);
        scaledImage2 = adapthisteq(scaledImage);
        outputFileName = sprintf('Img_slice%d-b%d.png',m,n);
        fullFileName = fullfile(outputPath, outputFileName);
        imwrite(scaledImage2,fullFileName,'png');
        imshow(scaledImage2,[]);
        title(sprintf('Processed Image slice%d-b%d',m,n));
        pause(1);
    end
end

%create a folder labeled 'masks' in the parentPath folder to output Diffusion images saved as png files 
mkdir([path 'Masks']);
outputPath = ([path 'Masks']);
%write out png files of created masks for Diffusion 
for m = 1:size(Ndiffimg,3)
      outputFileName = sprintf('slice%d.png',m);
      fullFileName = fullfile(outputPath, outputFileName);
      imwrite(final_mask(:,:,m),fullFileName,'png');
      imshow(final_mask(:,:,m),[]);
      title(sprintf('Processed Mask slice %d',m));
      pause(1);
end
close
%% Calcualting the SNR and the Weighting Vector
%% 
%SNR=signal_n/noise_n
% imag_airway= zeros(size(Ndiffimg));
SNR_vec = zeros(1, length(bValues));
signal_n_avg = zeros(1, length(bValues));
std_noise_n = zeros(1, length(bValues));
mean_noise= zeros(1, length(bValues));
for n = 1:length(bValues)
        %Calculating the siganl vector
        signal_vec = Ndiffimg(:,:,:,n);
        signal_vec(final_mask==0)=[];
        signal_n_avg(n) = mean(signal_vec);
%         imag_airway(:,:,:,n) =Ndiffimg(:,:,:,n).*airway_mask;
        %Calculating the noise vector
        noise_vec=(Ndiffimg(:,:,:,n));      
        noise_vec(noise_mask==0)=[];       
        mean_noise(n) = mean(noise_vec);            %the mean of the noise               
        std_noise_n(n) = std(noise_vec);            %the standard deviation of the noise
        SNR_vec(n) = signal_n_avg(n) / std_noise_n(n); %signal to noise ratio
end
SNR_vec=SNR_vec';
weight_vec = zeros(1, length(SNR_vec));
for n = 1:length(SNR_vec)
    weight_vec(n) = SNR_vec(n) / SNR_vec(1); %the weighting vector
end
 weight_vec =  weight_vec';
cd(path);
SNR_table=table( bValues,SNR_vec,weight_vec);
f=figure('position',[350 350 300 250]);
t = uitable(f,'Data',SNR_table{:,:},'ColumnName',...
    SNR_table.Properties.VariableNames,'Position',[20 20 262 204]);
print('SNR Table','-dpng','-r300');
close all; 
%% Calculating ADC & R maps using linear fit
%% 
if strcmp(W_Linear_fit,'yes') == 1
    fprintf( 'performing W.Linear Fitting...\n' );
ADCwvec=zeros(1,num_ones);
Rw_vec=zeros(1,num_ones);
ADCmapw=zeros(final_mask_size);
Rmapw=zeros(final_mask_size);
for wlin_loop=1:num_ones
    x=bValues';
    y=log(Data_vec(wlin_loop,:));              
    w=weight_vec';            
    m=1;   % m is the polynomial degree
    %Calculating the ADC map
    p_w = polyfitweighted(x,y,m,w); %or we can use p = wpolyfit(x,y,m,w); it does the same thing  
    ADCwvec(wlin_loop) = -p_w(1);
    %Calculating the R map           
    R_Matrix1 = zeros (length(x),2);
    R_Matrix1(:, 1) =x;
    R_Matrix1(:, 2) =y;
    R_weighted1 = weightedcorrs(R_Matrix1,w); % Weighted Correlation Matrix
    Rw_vec(wlin_loop) = R_weighted1(2);      %for weithting R map  
    %convert vector to 2D image 
    choosen_pixel=pixel_location(wlin_loop);
    ADCmapw(choosen_pixel) = ADCwvec(wlin_loop);
    Rmapw(choosen_pixel) = Rw_vec(wlin_loop);

end 
    fprintf( 'done...\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove background noise and airways from ADC map for linear fit
ADCbgrw = ADCmapw.*final_mask;
ADCw_R2=(Rmapw.*final_mask).^2;

% limit ADC to physically allowable values
ADCbgrw(ADCbgrw<0)=0;
ADCbgrw(ADCbgrw>0.14)=0; %Xe  % default to xenon

meanADCw = mean(ADCbgrw(ADCbgrw>0));       %ADC 
stdADCw = std(ADCbgrw(ADCbgrw>0));         %standard deviation
CVADCw = stdADCw/meanADCw;                 %Coefficient of Variation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write out ADC maps and histogram for linear fitting 
mkdir([path 'W.Linear ADC Analysis']);
outputPath = ([path 'W.Linear ADC Analysis']);
cd(outputPath);
colormap jet;
ljet = colormap;
ljet(1,:)=0;
fprintf( 'Write out ADC maps and Histogram...\n' );
%The following will print out each individual ADC map for each slice as a png file with 600 dpi resolution 
for xi = 1:size(Ndiffimg,3)
      outputFileName = sprintf('ADC map-W.Linear-slice%d.png',xi);
      imagesc(ADCbgrw(:, :, xi)); colormap(ljet);
      xx = colorbar;
      xx.Label.String = ' (cm^2.S^-^1)';caxis([0 0.14]);
      title(outputFileName);
%       fullFileName = fullfile(outputPath, outputFileName);
      print(outputFileName,'-r150','-dpng');
      pause(1);
end
%create histogram and print it out 
str_meanADCw=sprintf('Mean = %.4f cm^2/s', meanADCw);
str_stdADCw=sprintf(', STD = %.3f', stdADCw);
str_CVADCw=sprintf(', CV = %.3f', CVADCw);
str_legendw=[str_meanADCw,str_stdADCw,str_CVADCw];
figure('position',[350 350 1000 500]);
h=histogram(ADCbgrw(ADCbgrw>0),100); 
title_str1=['W_Linear ADC Histogram-',Subject_Info];
title(title_str1);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.1]);
set(gca,'FontSize',12)
legend(str_legendw);
print('Weighted ADC Histogram','-dpng','-r300');
close all;
wLiner_Elapsed_time=toc/60
fprintf( ' W-Linear Analysis done...\n' );
end
%% Calculating ADC & R maps using linear fit
%% 
if strcmp(Linear_fit,'yes') == 1
    fprintf( 'performing Linear Fitting...\n' );
ADCvec=zeros(1,num_ones);
R_fit=zeros(1,num_ones);
ADCmap = zeros(final_mask_size);
Rmap = zeros(final_mask_size);
for lin_loop=1:num_ones
    p = polyfit(bValues', log(Data_vec(lin_loop,:)), 1);
    ADCvec(lin_loop) = -p(1);
    R1 = corrcoef(bValues',log(Data_vec(lin_loop,:))); 
    R_fit(lin_loop)=R1(2);
    %convert vector to 2D image 
    choosen_pixel=pixel_location(lin_loop);
    ADCmap(choosen_pixel) = ADCvec(lin_loop);
    Rmap(choosen_pixel) = R_fit(lin_loop); % get R values from fit 
end
    fprintf( 'done...\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove background noise and airways from ADC map for linear fit
ADCbgr = ADCmap.*final_mask;
ADC_R2=(Rmap.*final_mask).^2;

% limit ADC to physically allowable values
ADCbgr(ADCbgr<0)=0;
ADCbgr(ADCbgr>0.14)=0; %Xe  % default to xenon

meanADC = mean(ADCbgr(ADCbgr>0));       %ADC 
stdADC = std(ADCbgr(ADCbgr>0));         %standard deviation
CVADC = stdADC/meanADC;                 %Coefficient of Variation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write out ADC maps and histogram for linear fitting 
mkdir([path 'Linear ADC Analysis']);
outputPath = ([path 'Linear ADC Analysis']);
cd(outputPath);
colormap jet;
ljet = colormap;
ljet(1,:)=0;
fprintf( 'Write out ADC maps and Histogram...\n' );
%The following will print out each individual ADC map for each slice as a png file with 600 dpi resolution 
for xi = 1:size(Ndiffimg,3)
      outputFileName = sprintf('ADC map-Linear-slice%d.png',xi);
      imagesc(ADCbgr(:, :, xi)); colormap(ljet);
      xx = colorbar;
      xx.Label.String = ' (cm^2.S^-^1)';caxis([0 0.14]);
      title(outputFileName);
%       fullFileName = fullfile(outputPath, outputFileName);
      print(outputFileName,'-r150','-dpng');
      pause(1);
end
%create histogram and print it out 
str_meanADC=sprintf('Mean = %.4f cm^2/s', meanADC);
str_stdADC=sprintf(', STD = %.3f', stdADC);
str_CVADC=sprintf(', CV = %.3f', CVADC);
str_legend=[str_meanADC,str_stdADC,str_CVADC];
figure('position',[350 350 1000 500]);
h=histogram(ADCbgr(ADCbgr>0),100); 
title_str1=['Linear ADC Histogram-',Subject_Info];
title(title_str1);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.1]);
set(gca,'FontSize',12)
legend(str_legend);
print('Linear ADC Histogram','-dpng','-r300');
close all;
Liner_Elapsed_time=toc/60
fprintf( ' Linear Analysis done...\n' );
end
%% Calculating ADC & R maps using ML fit
%% 
if strcmp(ML_fit,'yes') == 1
    tic
    fprintf( 'performing ML Fitting...\n' );
ADCmlvec=zeros(1,num_ones);
Sovec=zeros(1,num_ones);
ADCmlmap = zeros(final_mask_size);
Somapml = zeros(final_mask_size);
for ml_loop=1:num_ones
    ydata=Data_vec(ml_loop,:);
    x0 = [ydata(1),0.05];
    fun2 = @(x)fmaxlikelihood_ADC(x,bValues',ydata,std_noise_n);
    bestx2 = fminsearch(fun2,x0);
    Sovec(ml_loop)=bestx2(1);
    ADCmlvec(ml_loop)=bestx2(2);
    %convert vector to 2D image 
    choosen_pixel=pixel_location(ml_loop);
    ADCmlmap(choosen_pixel) = ADCmlvec(ml_loop);
    Somapml(choosen_pixel) = Sovec(ml_loop); % get R values from fit 
end
    fprintf( 'done...\n' );
%Remove background noise and airways from ADC map for linear fit
ADCmlbgr = ADCmlmap.*final_mask;

% limit ADC to physically allowable values
ADCmlbgr(ADCmlbgr<0)=0;
ADCmlbgr(ADCmlbgr>0.14)=0; %Xe  % default to xenon

meanADCml = mean(ADCmlbgr(ADCmlbgr>0));       %ADC 
stdADCml = std(ADCmlbgr(ADCmlbgr>0));         %standard deviation
CVADCml = stdADCml/meanADCml;                 %Coefficient of Variation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write out ADC maps and histogram for linear fitting 
mkdir([path 'Maxlikelihood ADC Analysis']);
outputPath = ([path 'Maxlikelihood ADC Analysis']);
cd(outputPath);
colormap jet;
ljet = colormap;
ljet(1,:)=0;
fprintf( 'Write out ADC maps and Histogram...\n' );
%The following will print out each individual ADC map for each slice as a png file with 600 dpi resolution 
for xi = 1:size(Ndiffimg,3)
      outputFileName = sprintf('ADC map-Maxlikelihood-slice%d.png',xi);
      imagesc(ADCmlbgr(:, :, xi)); colormap(ljet);
      xx = colorbar;
      xx.Label.String = ' (cm^2.S^-^1)';caxis([0 0.14]);
      title(outputFileName);
%       fullFileName = fullfile(outputPath, outputFileName);
      print(outputFileName,'-r150','-dpng');
%       pause(1);
end
%create histogram and print it out 
str_meanADC=sprintf('Mean = %.4f cm^2/s', meanADCml);
str_stdADC=sprintf(', STD = %.3f', stdADCml);
str_CVADC=sprintf(', CV = %.3f', CVADCml);
str_legend=[str_meanADC,str_stdADC,str_CVADC];
figure('position',[350 350 1000 500]);
h=histogram(ADCmlbgr(ADCmlbgr>0),100); 
title_str1=['Maxlikelihood ADC Histogram-',Subject_Info];
title(title_str1);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.1]);
set(gca,'FontSize',12)
legend(str_legend);
print('Maxlikelihood ADC Histogram','-dpng','-r300');
close all;
ML_Elapsed_time=toc/60
fprintf( ' Maxlikelihood Analysis done...\n' );
end
%% Bayesain ADC
%% 
if strcmp(Bayesain_fit,'yes') == 1    
fprintf( 'performing Bayesain Fitting...\n' );
cd(path);
ADC_Model='ADC_Bayesian_Model_1D.txt';
ADC_Bayesian=zeros(final_mask_size);
So_Bayesian=zeros(final_mask_size);
sigma_n=zeros(slices,nbs);
Bayes_nu0 = zeros(final_mask_size(3),final_mask_size(1)*final_mask_size(2));
Bayes_alpha = zeros(final_mask_size(3),final_mask_size(1)*final_mask_size(2));
for slice_n=1:slices
        M = zeros(final_mask_size(1),final_mask_size(2),nbs);
        for n = 1:length(bValues)
            noise_vec=Ndiffimg(:,:,slice_n,n);
            noise_vec(noise_mask(:,:,slice_n)==0)=[];
            std_slice_n = std(noise_vec(:));
            sigma_n(slice_n,n)=std_slice_n;
            M(:,:,n)=((Ndiffimg(:,:,slice_n,n).*final_mask(:,:,slice_n)).^2)/(std_slice_n^2);
        end
        M_Data_vec=reshape(M,[final_mask_size(1)*final_mask_size(2),length(bValues)]);
        M = M_Data_vec(any(M_Data_vec,2),:);
        dataStruct= struct('M', M, 'b', bValues, 'sigma', sigma_n(slice_n,:));
        num_ones=final_mask(:,:,slice_n);
        num_ones=sum(num_ones(:)); % for loopping 
        Bayes_Model_3D=ADC_Bayesian_Model_1D(num_ones,bValues); % Write Bayes Model

        %initializing
        a=0; b_int=0.14;
        init_adc = a + (b_int-a).*rand(num_ones,1); %adc prior  

        nu0_prior=reshape(Ndiffimg(:,:,slice_n,1).*final_mask(:,:,slice_n),[final_mask_size(1)*final_mask_size(2),1]); % So %adc prior
        nu0_prior = nu0_prior(any(nu0_prior,2),:);  %So prior 
        initStructs =struct('nu0', nu0_prior','alpha', init_adc');

        % Call WindBuges
        nchains  = 1; % How Many Chains?
        nburnin  = 3000; % How Many Burn-in Samples?
        nsamples = 2000;  % How Many Recorded Samples?
        fprintf( 'Running WinBUGS for Bayesian Fit...\n' );
        [samples, stats] = matbugs(dataStruct, ...
            fullfile(path, ADC_Model), ...
            'init', initStructs, ...
            'nChains', nchains, ...
            'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', 1, 'DICstatus', 0, 'refreshrate',100,  ...
            'monitorParams', {'nu0','alpha'}, ...
            'Bugdir', 'C:\Users\BDA5IK\Downloads\winbugs143_unrestricted\winbugs14_full_patched\WinBUGS14');

Bayes_nu0(slice_n,1:num_ones) =stats.mean.nu0;
Bayes_alpha(slice_n,1:num_ones) =stats.mean.alpha;
end

pixel_location=find(final_mask); %to find pixel location to convert the array back to 3D image
num_ones=sum(final_mask(:)); % for loopping 

Bayes_nu0_2=Bayes_nu0';
Bayes_nu0_2=reshape(Bayes_nu0_2,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),1]);
Bayes_nu0_2 = Bayes_nu0_2(any(Bayes_nu0_2,2),:);

Bayes_alpha_2=Bayes_alpha';
Bayes_alpha_2=reshape(Bayes_alpha_2,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),1]);
Bayes_alpha_2 = Bayes_alpha_2(any(Bayes_alpha_2,2),:);

%convert vector to 2D image 
Bayes_nu0_map=zeros(final_mask_size);
Bayes_ADC_map=zeros(final_mask_size);
for bayes_loop=1:num_ones
    choosen_pixel=pixel_location(bayes_loop);
    Bayes_nu0_map(choosen_pixel) = Bayes_nu0_2(bayes_loop);
    Bayes_ADC_map(choosen_pixel) = Bayes_alpha_2(bayes_loop);
end 

%    figure; imslice(Bayes_ADC_map)
Bayes_adc=Bayes_ADC_map.*final_mask;
Bayes_adc(Bayes_adc<0)=0;
Bayes_adc(Bayes_adc>0.14)=0; %Xe  % default to xenon

Mean_ADC_Bayes=mean(Bayes_adc(Bayes_adc>0));
stdADC_Bayes = std(Bayes_adc(Bayes_adc>0));         %standard deviation
CVADC_Bayes = stdADC_Bayes/Mean_ADC_Bayes;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write out ADC maps and histogram for Bayesain fitting 
mkdir([path 'Bayesain ADC Analysis']);
outputPath = ([path 'Bayesain ADC Analysis']);
cd(outputPath);
colormap jet;
ljet = colormap;
ljet(1,:)=0;
fprintf( 'Write out ADC maps and Histogram...\n' );
%The following will print out each individual ADC map for each slice as a png file with 600 dpi resolution 
for xi = 1:size(Ndiffimg,3)
      outputFileName = sprintf('ADC map-Bayesain-slice%d.png',xi);
      imagesc(Bayes_adc(:, :, xi)); colormap(ljet);
      xx = colorbar;
      xx.Label.String = ' (cm^2.S^-^1)';caxis([0 0.14]);
      title(outputFileName);
%       fullFileName = fullfile(outputPath, outputFileName);
      print(outputFileName,'-r150','-dpng');
      pause(1);
end
%create histogram and print it out 
str_mean_Bayes=sprintf('Mean = %.4f cm^2/s', Mean_ADC_Bayes);
str_std_Bayes=sprintf(', STD = %.3f', stdADC_Bayes);
str_CVADC_Bayes=sprintf(', CV = %.3f', CVADC_Bayes);
str_legendw=[str_mean_Bayes,str_std_Bayes,str_CVADC_Bayes];
figure('position',[350 350 1000 500]);
h=histogram(Bayes_adc(Bayes_adc>0),100); 
title_str2=['Bayesain ADC Histogram-',Subject_Info];
title(title_str2);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.1]);
set(gca,'FontSize',12)
legend(str_legendw);
print('Bayesain ADC Histogram','-dpng','-r300');
close all;
Bayes_Elapsed_time=toc/60
fprintf( ' Bayesain Analysis done...\n' );
end

%% save workspace
%% 
if strcmp(savedata,'yes') == 1
   save_data=[path Subject_num,'_ADC_Analysis_workspace.mat'];
   save(save_data)
%    save('\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images\IRC186H-XenonGeneral\IRC186H Data and Analyses\IRC186H-240\20190417\RAW IRC186H-240\ADC RESULT IRC186H_240')
end
%% 





