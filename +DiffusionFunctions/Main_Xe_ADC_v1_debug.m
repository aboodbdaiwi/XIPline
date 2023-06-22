
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matlab script for processing Philips raw data into ADC maps %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10/10/2017 : Written by Zackary I. Cleveland (ZIC) and Abdullah S. Bdaiwi
% 10/11/2017 : ZIC, airway segmentation added and pathway bug fixed.
% 10/12/2017 : ZIC, edited to allow for fitting of selected b values 
% 11/15/2017 : AB,  edited to add SNR and weighting vector calculation 
% 11/15/2017 : AB,  edited to add weighting to the Polynomial fit to calculate the ADC 
% 11/29/2017 : AB,  edited to add weighting to the Correlation coefficients fit to calculate the Rmap 
% 05/21/2018 : AB,  edited to add Bayesain fitting
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for any questions, please contact Abdullah Bdaiwi at Abdullah.Bdaiwi@cchmc.org
clc;clear;close all;% start with fresh variables, etc.
%% Thing you need to do before you hit run
%% 
tic
subject_info='no';
savedata = 'no'; % if yes, save workspace to current directory
filterim = 'Fermi'; % choose 'Fermi','Gaussian', 'Han' or 'None'
segment = 'theshold';  % choose 'manual' to segment lung parenchyma manually otherwise will use threshold
airway_segment = 'manual';  % segment airways mannually (default) or "quicknEZ"  with 
                            % a quick and easy signal ratio involving
                            % awcof. Manual is the default   
Linear_fit='no'; % if Yes,  linear fit will be performed
W_Linear_fit = 'yes'; % if Yes, weighted linear fit will be performed
Bayesain_fit = 'yes'; % if Yes, Bayesain fit will be performed
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
% bval = length(bValues)-1; % index of 2nd bvalue selected for 2 bvalue fit
bs2fit = 1:length(bValues); % b values used in ADC fitting
%If you want to change it to fit only the 1st, 3rd, and 4th bvalues, simply edit to
% bs2fit = [1 2 3]; % b values used in ADC fitting

% R2_threshold=.1; %R-squared threshold
%% Subject Info.
%% 
if strcmp(subject_info, 'Yes')==1
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

% Normalize the intensity of the original image to fall between [0,100].
Ndiffimg=diffimg-min(diffimg(:));
Ndiffimg=diffimg/max(diffimg(:));
Ndiffimg=Ndiffimg*100;
figure; imslice(Ndiffimg)
%% Segment lung parenchyma and airways
%% 
if strcmp(segment,'manual') == 1
        lung_mask =ADC_lung_parenchyma_Segment(airway_segment,awcof,diffimg);
else
        lung_mask = ADC_thresh(NoiK,diffimg,SE,nzcof);
end
%Segment the airway 
airway_mask = ADC_Airway_Segment(airway_segment,awcof,diffimg);

final_mask =lung_mask;
final_mask(airway_mask==1)=0;

% flip the binary mask to create the noise mask(1 to 0, 0 to 1)
noise_mask = ((~lung_mask)+(~airway_mask))-1;
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
imag_airway= zeros(size(Ndiffimg));
SNR_vec = zeros(1, length(bValues));
signal_n_avg = zeros(1, length(bValues));
std_noise_n = zeros(1, length(bValues));
mean_noise= zeros(1, length(bValues));
for n = 1:length(bValues)
        %Calculating the siganl vector
        signal_vec = Ndiffimg(:,:,:,n);
        signal_vec(final_mask==0)=[];
        signal_n_avg(n) = mean(signal_vec);
        imag_airway(:,:,:,n) =Ndiffimg(:,:,:,n).*airway_mask;
        %Calculating the noise vector
        noise_vec=(Ndiffimg(:,:,:,n));      
        noise_vec(noise_mask==0)=[];       
        mean_noise(n) = mean(noise_vec);            %the mean of the noise               
        std_noise_n(n) = std(noise_vec);            %the standard deviation of the noise
        SNR_vec(n) = signal_n_avg(n) / std_noise_n(n); %signal to noise ratio
end

weight_vec = zeros(1, length(SNR_vec));
for n = 1:length(SNR_vec)
    weight_vec(n) = SNR_vec(n) / SNR_vec(1); %the weighting vector
end
 weight_vec =  weight_vec';
cd(path);
SNR_table=table( bValues',SNR_vec',weight_vec);
f=figure('position',[350 350 300 250]);
t = uitable(f,'Data',SNR_table{:,:},'ColumnName',...
    SNR_table.Properties.VariableNames,'Position',[20 20 262 204]);
print('SNR Table','-dpng','-r300');
close all; 
% disp(SNR_vec)
% disp(weight_vec)
%% Calculating ADC & R maps using linear fit
%% 
if strcmp(W_Linear_fit,'Yes') == 1
    fprintf( 'performing W.Linear Fitting...\n' );
ADCmapw = zeros(ky,kx/kx_oversample_factor ,slices);
Rmapw = zeros(ky,kx/kx_oversample_factor ,slices);
ADC_airways = zeros(ky,kx/kx_oversample_factor ,slices);
Rmap_airways = zeros(ky,kx/kx_oversample_factor ,slices);
    for sl=1:slices
        for kys=1:ky
            for kxs=1:kx/kx_oversample_factor     
                pixel_value = squeeze(Ndiffimg(kys,kxs,sl,:))'; 
                s_2=squeeze(imag_airway(kys,kxs,sl,:))';
                %%%%%%%%%%%%%%% airway  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                p_air = polyfit(bValues(bs2fit),log(s_2(bs2fit)),1); 
                ADC_airways(kys,kxs,sl) = -p_air(1);
                R1_air = corrcoef(bValues(bs2fit),log(s_2(bs2fit))); 
                Rmap_airways(kys,kxs,sl) = R1_air(2); % get R values from fit 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%  Weighted Linear %%%%%%%%%%%%%%%%%%
                x=bValues(bs2fit);
                y=log(pixel_value(bs2fit));              
                w=weight_vec(bs2fit)';            
                m=1;   % m is the polynomial degree
                %Calculating the ADC map
                p_w = polyfitweighted(x,y,m,w); %or we can use p = wpolyfit(x,y,m,w); it does the same thing  
                ADCmapw(kys,kxs,sl) = -p_w(1);
                %Calculating the R map           
                R_Matrix1 = zeros (length(x),2);
                R_Matrix1(:, 1) =x;
                R_Matrix1(:, 2) =y;
                R_Traditional1 = corrcoef(R_Matrix1);     % Traditional Correlation Matrix 
                R_weighted1 = weightedcorrs(R_Matrix1,w); % Weighted Correlation Matrix
                %Rmap(kys,kxs,sl) =  R_Traditional1(2);  %without weighting 
                Rmapw(kys,kxs,sl) = R_weighted1(2);      %for weithting R map          
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            end
        end
    end
    fprintf( 'done...\n' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove background noise and airways from ADC map for linear fit
ADCbgrw = ADCmapw.*final_mask;
ADC_airbgr=ADC_airways;
ADC_airbgr(airway_mask==0)=[];

Airway_R2_airway=(Rmap_airways.*airway_mask).^2;
Airway_R2_airway(airway_mask==0)=[];

ADCw_R2=(Rmapw.*final_mask).^2;

% limit ADC to physically allowable values
ADCbgrw(ADCbgrw<0)=0;
ADC_airbgr(ADC_airbgr<0)=0;

ADCbgrw(ADCbgrw>0.14)=0; %Xe  % default to xenon
ADC_airbgr(ADC_airbgr>0.14)=0;

meanADCw = mean(ADCbgrw(ADCbgrw>0));       %ADC 
stdADCw = std(ADCbgrw(ADCbgrw>0));         %standard deviation
CVADCw = stdADCw/meanADCw;                 %Coefficient of Variation 

meanADC_airway = mean(ADC_airbgr(ADC_airbgr>0));       %ADC 
stdADC_airway = std(ADC_airbgr(ADC_airbgr>0));         %standard deviation
CVADC_airway = stdADC_airway /meanADC_airway;    %Coefficient of Variation 
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
Liner_Elapsed_time=toc/60
fprintf( ' W-Linear Analysis done...\n' );

%% Calculating ADC & R maps using linear fit
%% 
if strcmp(Linear_fit,'Yes') == 1
    fprintf( 'performing Linear Fitting...\n' );
ADCmap = zeros(ky,kx/kx_oversample_factor ,slices);
Rmap = zeros(ky,kx/kx_oversample_factor ,slices);
ADC_airways = zeros(ky,kx/kx_oversample_factor ,slices);
Rmap_airways = zeros(ky,kx/kx_oversample_factor ,slices);
    for sl=1:slices
        for kys=1:ky
            for kxs=1:kx/kx_oversample_factor     
                pixel_value = squeeze(Ndiffimg(kys,kxs,sl,:))'; 
                s_2=squeeze(imag_airway(kys,kxs,sl,:))';
                %%%%%%%%%%%%%%% Without Weighting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                p = polyfit(bValues(bs2fit),log(pixel_value(bs2fit)),1); 
                ADCmap(kys,kxs,sl) = -p(1);
                R1 = corrcoef(bValues(bs2fit),log(pixel_value(bs2fit))); 
                Rmap(kys,kxs,sl) = R1(2); % get R values from fit 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%% airways  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
                p_air = polyfit(bValues(bs2fit),log(s_2(bs2fit)),1); 
                ADC_airways(kys,kxs,sl) = -p_air(1);
                R1_air = corrcoef(bValues(bs2fit),log(s_2(bs2fit))); 
                Rmap_airways(kys,kxs,sl) = R1_air(2); % get R values from fit 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
    fprintf( 'done...\n' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove background noise and airways from ADC map for linear fit
ADCbgr = ADCmap.*final_mask;
ADC_airbgr=ADC_airways;
ADC_airbgr(airway_mask==0)=[];

Airway_R2_airway=(Rmap_airways.*airway_mask).^2;
Airway_R2_airway(airway_mask==0)=[];

ADC_R2=(Rmap.*final_mask).^2;

% limit ADC to physically allowable values
ADCbgr(ADCbgr<0)=0;
ADC_airbgr(ADC_airbgr<0)=0;

ADCbgr(ADCbgr>0.14)=0; %Xe  % default to xenon
ADC_airbgr(ADC_airbgr>0.14)=0;

meanADC = mean(ADCbgr(ADCbgr>0));       %ADC 
stdADC = std(ADCbgr(ADCbgr>0));         %standard deviation
CVADC = stdADC/meanADC;                 %Coefficient of Variation 

meanADC_airway = mean(ADC_airbgr(ADC_airbgr>0));       %ADC 
stdADC_airway = std(ADC_airbgr(ADC_airbgr>0));         %standard deviation
CVADC_airway = stdADC_airway /meanADC_airway;    %Coefficient of Variation 
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
if strcmp(W_Linear_fit,'Yes') == 1    
title_str1=['W_Linear ADC Histogram-',Subject_Info];
else
title_str1=['Linear ADC Histogram-',Subject_Info];
end
title(title_str1);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.1]);
set(gca,'FontSize',12)
legend(str_legend);
print('Weighted ADC Histogram','-dpng','-r300');
close all;
Liner_Elapsed_time=toc/60
fprintf( ' Linear Analysis done...\n' );
% 
%% Bayesain ADC
%% 
fprintf( 'performing Linear Fitting...\n' );
if strcmp(Bayesain_fit,'Yes') == 1    
tic
% Write Bayes Model
cd(path)
diff_size=size(Ndiffimg);
Bayes_Model_3D=ADC_Bayesian_Model_3D(diff_size(1),diff_size(2),bValues);
ADC_Model='ADC_Bayes_Model_3D.txt';
ADC_Bayesian=zeros(ky,kx/kx_oversample_factor ,slices);
So_Bayesian=zeros(ky,kx/kx_oversample_factor ,slices);
sigma_n=zeros(slices,nbs);
for slice_n=1:slices   
        for n = 1:length(bValues)
            noise_vec=Ndiffimg(:,:,slice_n,n);
            noise_vec(noise_mask(:,:,slice_n)==0)=[];
            std_slice_n = std(noise_vec(:));  
            sigma_n(slice_n,n)=std_slice_n;
        end
        sigma=sigma_n(slice_n,:);
        M=Ndiffimg(:,:,slice_n,:);
        for n=1:nbs
            M(:,:,n)=(Ndiffimg(:,:,slice_n,n).^2)/(sigma(n)^2);
        end
        M=reshape(M,[ky,kx/kx_oversample_factor,length(bValues)]);
        dataStruct= struct('M', M, 'b', bValues, 'sigma', sigma);
        %initializing
        a=0;
        b_int=0.14;
        init_adc = a + (b_int-a).*rand(ky,kx/kx_oversample_factor); %adc prior
        
        nu0_prior=Ndiffimg(:,:,slice_n,1); % So %adc prior

        initStructs =struct('nu0', nu0_prior,'alpha', init_adc);
% figure; imslice(nu0_prior)
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
       %
        Bayes_nu0=stats.mean.nu0;
        So_Bayesian(:,:,slice_n)=Bayes_nu0;

        Bayes_alpha=stats.mean.alpha;
        ADC_Bayesian(:,:,slice_n)=Bayes_alpha;  
end
fprintf( 'Done...\n' );
%    figure; imslice(Bayes_adc)
Bayes_adc=ADC_Bayesian.*final_mask;
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
print('Weighted ADC Histogram','-dpng','-r300');
close all;
Bayes_Elapsed_time=toc/60
end
fprintf( ' Bayesain Analysis done...\n' );

%% save workspace
%% 
if strcmp(savedata,'yes') == 1
   save_data=[path Subject_num,'ADC_Analysis_workspace.mat'];
   save(save_data)
%    save('\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images\IRC186H-XenonGeneral\IRC186H Data and Analyses\IRC186H-240\20190417\RAW IRC186H-240\ADC RESULT IRC186H_240')
end
%% 





