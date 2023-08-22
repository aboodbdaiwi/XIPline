
function [ADCmap,ADCcoloredmap,SNR_table,meanADC,stdADC,ADC_hist] = ADC_Analysis (diffimg, lung_mask, airway_mask, bvalues, ADCFittingType,ADCAnalysisType, WinBUGSPath, outputpath,PatientAge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% ADC Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Inputs:
%       diffimg = 4D images (x,y,slices,bvalue)
%       lung_mask = 3D lung mask, (x,y,slices)
%       airway_mask = 3D airway mask (x,y,slices)
%       bvalues = string: 
%       ADCFittingType = string:  'Log Weighted Linear' | 'Log Linear' | 'Non-Linear' | 'Bayesian'
%       ADCAnalysisType = string: human | animals
%       WinBUGSPath = string: WinBUGS .exe file location
%       outputpath = string containing the path of the folder to store the output
%

%   Package: 
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%% create a final mask
final_mask =lung_mask;
final_mask(airway_mask==1)=0;

% flip the binary mask to create the noise mask(1 to 0, 0 to 1)
noise_mask = ((~lung_mask)+(~airway_mask))-1;
noise_mask=noise_mask>0;

%% get rid of any slice that has no signal (un-segmented)
slices = size(diffimg,3);
Nbvalues = length(bvalues);
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
Ndiffimg= zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker), Nbvalues);
Nfinal_mask= zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nnoise_mask= zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
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
% Ndiffimg=diffimg-min(diffimg(:));
Ndiffimg=diffimg/max(diffimg(:));
Ndiffimg=Ndiffimg*100;

final_mask_size=size(final_mask);
%convert 4D to 2D array (lenghth pixel inside the final mask and number of b values
Data_vec=reshape(Ndiffimg.*final_mask,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),length(bvalues)]);
Data_vec = Data_vec(any(Data_vec,2),:);
pixel_location=find(final_mask); %to find pixel location to convert the array back to 3D image
num_ones=sum(final_mask(:)); % for loopping 

%% Calcualting the SNR and the Weighting Vector
close all; 
SNR_vec = zeros(1, length(bvalues));
signal_n_avg = zeros(1, length(bvalues));
std_noise_n = zeros(1, length(bvalues));
mean_noise= zeros(1, length(bvalues));
for n = 1:length(bvalues)
        %Calculating the siganl vector
        signal_vec = Ndiffimg(:,:,:,n);
        signal_vec(final_mask==0)=[];
        signal_n_avg(n) = mean(signal_vec);
        
        %Calculating the noise vector
        noise_vec=(Ndiffimg(:,:,:,n));      
        noise_vec(noise_mask==0)=[];       
        mean_noise(n) = mean(noise_vec);            %the mean of the noise               
        std_noise_n(n) = std(noise_vec);            %the standard deviation of the noise
        SNR_vec(n) = round(signal_n_avg(n) / std_noise_n(n),2)*sqrt(2 - (pi/2)); %signal to noise ratio
end
SNR_vec=SNR_vec';
weight_vec = zeros(1, length(SNR_vec));
for n = 1:length(SNR_vec)
    weight_vec(n) = SNR_vec(n) / SNR_vec(1); %the weighting vector
end
 weight_vec =  weight_vec';
cd(outputpath);
SNR_table=table( bvalues',SNR_vec,weight_vec);
SNRFig = figure('Name','SNR','units','normalized','position',[350 350 300 250]);
set(SNRFig,'WindowState','minimized');
set(SNRFig,'color','white','Units','inches','Position',[0.25 0.25 3 2.5])

% SNRFig=figure('position',[350 350 300 250]);
t = uitable(SNRFig,'Data',SNR_table{:,:},'ColumnName',...
    SNR_table.Properties.VariableNames,'Position',[20 20 260 200]);
% print('SNR Table','-dpng','-r300');
saveas(gca,'SNR_Table.png')
close all; 
%% Calculating ADC 
%% 
ADCvec = zeros(1,num_ones);
Rvec = zeros(1,num_ones);
Sovec = zeros(1,num_ones);
ADCmap = zeros(final_mask_size);
Rmap = zeros(final_mask_size);
Somap = zeros(final_mask_size);


switch ADCFittingType
    case 'Log Linear'
        fprintf( 'performing Linear Fitting...\n' );
        for lin_loop=1:num_ones
            p = polyfit(bvalues, log(Data_vec(lin_loop,:)), 1);
            ADCvec(lin_loop) = -p(1);
            R1 = corrcoef(bvalues',log(Data_vec(lin_loop,:))); 
            Rvec(lin_loop)=R1(2);
            %convert vector to 2D image 
            choosen_pixel=pixel_location(lin_loop);
            ADCmap(choosen_pixel) = ADCvec(lin_loop);
            Rmap(choosen_pixel) = Rvec(lin_loop); % get R values from fit 
        end
        fprintf( 'done...\n' );
        %Remove background noise and airways from ADC map for linear fit
        ADCmap = ADCmap.*final_mask;
        Rmap=(Rmap.*final_mask).^2;
    case 'Log Weighted Linear'
        fprintf( 'performing W.Linear Fitting...\n' );
        for wlin_loop = 1:num_ones
            x = bvalues;
            y = log(Data_vec(wlin_loop,:));              
            w = weight_vec';            
            m = 1;   % m is the polynomial degree
            %Calculating the ADC map
            p_w = DiffusionFunctions.polyfitweighted(x,y,m,w); %or we can use p = wpolyfit(x,y,m,w); it does the same thing  
            ADCvec(wlin_loop) = -p_w(1);
            %Calculating the R map           
            R_Matrix1 = zeros (length(x),2);
            R_Matrix1(:, 1) = x;
            R_Matrix1(:, 2) = y;
            R_weighted1 = DiffusionFunctions.weightedcorrs(R_Matrix1,w); % Weighted Correlation Matrix
            Rvec(wlin_loop) = R_weighted1(2);      %for weithting R map  
            %convert vector to 2D image 
            choosen_pixel = pixel_location(wlin_loop);
            ADCmap(choosen_pixel) = ADCvec(wlin_loop);
            Rmap(choosen_pixel) = Rvec(wlin_loop);
        end 
        fprintf( 'done...\n' );
        %Remove background noise and airways from ADC map for linear fit
        ADCmap = ADCmap.*final_mask;
        Rmap = (Rmap.*final_mask).^2;

    case 'Non-Linear'
        fprintf( 'performing None-Linear Fitting...\n' );
        for nonl_loop=1:num_ones
            ydata=Data_vec(nonl_loop,:);
            x0 = [ydata(1),0.35];
            fun2 = @(x)sseval(x,bvalues,ydata);
            bestx2 = fminsearch(fun2,x0);
            Sovec(nonl_loop)=bestx2(1);
            ADCvec(nonl_loop)=bestx2(2);
            %convert vector to 2D image 
            choosen_pixel=pixel_location(nonl_loop);
            ADCmap(choosen_pixel) = ADCvec(nonl_loop);
            Somap(choosen_pixel) = Sovec(nonl_loop); % get R values from fit 
        end
        fprintf( 'Completed...\n' );      
        %Remove background noise and airways from ADC map for linear fit
        ADCmap = ADCmap.*final_mask;
        Rmap=(Rmap.*final_mask).^2;        
    case 'Bayesian'
        clc
        fprintf( 'performing Bayesain Fitting...\n' );
        cd(outputpath);
        if strcmp(ADCAnalysisType, 'human') == 1  
            ADC_Model='ADC_Model_Human_1D.txt';
        elseif strcmp(ADCAnalysisType, 'animals') == 1
            ADC_Model='ADC_Model_nimal_1D.txt';
        end        
        sigma_n=zeros(slices,Nbvalues);
        Bayes_nu0 = zeros(final_mask_size(3),final_mask_size(1)*final_mask_size(2));
        Bayes_alpha = zeros(final_mask_size(3),final_mask_size(1)*final_mask_size(2));
        for slice_n=1:slices
            M = zeros(final_mask_size(1),final_mask_size(2),Nbvalues);
            for n = 1:Nbvalues
                noise_vec=Ndiffimg(:,:,slice_n,n);
                noise_vec(noise_mask(:,:,slice_n)==0)=[];
                std_slice_n = std(noise_vec(:));
                sigma_n(slice_n,n)=std_slice_n;
                M(:,:,n)=((Ndiffimg(:,:,slice_n,n).*final_mask(:,:,slice_n)).^2)/(std_slice_n^2);
            end
            M_Data_vec=reshape(M,[final_mask_size(1)*final_mask_size(2),length(bvalues)]);
            M = M_Data_vec(any(M_Data_vec,2),:);
            dataStruct= struct('M', M, 'b', bvalues, 'sigma', sigma_n(slice_n,:));
            num_ones=final_mask(:,:,slice_n);
            num_ones=sum(num_ones(:)); % for loopping 
            if strcmp(ADCAnalysisType, 'human') == 1  
                Bayes_Model_1D = DiffusionFunctions.ADC_Model_Human_1D(num_ones,Nbvalues); % Write Bayes Model
                b_int=0.14;
            elseif strcmp(ADCAnalysisType, 'animals') == 1
                Bayes_Model_1D = DiffusionFunctions.ADC_Model_animal_1D(num_ones,Nbvalues); % Write Bayes Model  
                b_int=0.039;
            end
            %initializing
            a=0; 
            init_adc = a + (b_int-a).*rand(num_ones,1); %adc prior  

            nu0_prior=reshape(Ndiffimg(:,:,slice_n,1).*final_mask(:,:,slice_n),[final_mask_size(1)*final_mask_size(2),1]); % So %adc prior
            nu0_prior = nu0_prior(any(nu0_prior,2),:);  %So prior 
            initStructs =struct('nu0', nu0_prior','alpha', init_adc');
            % Call WindBuges
            nchains  = 1; % How Many Chains?
            nburnin  = 3000; % How Many Burn-in Samples?
            nsamples = 2000;  % How Many Recorded Samples?
            fprintf( 'Running WinBUGS for Bayesian Fit...\n' );
            [samples, stats] =  DiffusionFunctions.matbugs(dataStruct, ...
                                fullfile(outputpath, ADC_Model), ...
                                'init', initStructs, ...
                                'nChains', nchains, ...
                                'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
                                'thin', 1, 'DICstatus', 0, 'refreshrate',100,  ...
                                'monitorParams', {'nu0','alpha'}, ...
                                'Bugdir', WinBUGSPath);
        Bayes_nu0(slice_n,1:num_ones) =stats.mean.nu0;
        Bayes_alpha(slice_n,1:num_ones) =stats.mean.alpha;
        end

        pixel_location=find(final_mask); % find pixel location to convert the array back to 3D image
        num_ones=sum(final_mask(:)); % for loopping 
        Bayes_nu0_2=Bayes_nu0';
        Bayes_nu0_2=reshape(Bayes_nu0_2,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),1]);
        Bayes_nu0_2 = Bayes_nu0_2(any(Bayes_nu0_2,2),:);
        Bayes_alpha_2=Bayes_alpha';
        Bayes_alpha_2=reshape(Bayes_alpha_2,[final_mask_size(1)*final_mask_size(2)*final_mask_size(3),1]);
        Bayes_alpha_2 = Bayes_alpha_2(any(Bayes_alpha_2,2),:);
        %convert vector to 2D image 
        for bayes_loop=1:num_ones
            choosen_pixel=pixel_location(bayes_loop);
            Somap(choosen_pixel) = Bayes_nu0_2(bayes_loop);
            ADCmap(choosen_pixel) = Bayes_alpha_2(bayes_loop);
        end 

end

% limit ADC to physically allowable values
ADCmap(ADCmap<0)=0;
ADCmap(ADCmap>0.14)=0; %Xe  % default to xenon
meanADC = mean(ADCmap(ADCmap>0));       %ADC 
stdADC = std(ADCmap(ADCmap>0));         %standard deviation
CVADC = stdADC/meanADC;                 %Coefficient of Variation 
%create histogram and print it out 
str_meanADC=sprintf('Mean = %.4f cm^2/s', meanADC);
str_stdADC=sprintf(', STD = %.3f', stdADC);
str_CVADC=sprintf(', CV = %.3f', CVADC);
str_legend=[str_meanADC,str_stdADC,str_CVADC];
figure('position',[350 350 1000 500]);

h=histogram(ADCmap(ADCmap>0),100); 
histtitle_str1=['ADC Histogram-',ADCFittingType];
title(histtitle_str1);
xlabel('ADC cm^2/sec','FontSize',12,'Color','k');
ylabel('Pixel Count.','FontSize',12,'Color','k');
xlim([0 0.14]);
set(gca,'FontSize',12)
% legend(str_legend);
% print(histtitle_str1,'-dpng','-r300');

% plot healthy dist.
hold on
ADCLB_RefMean = 0.0002*PatientAge+0.029; % 
ADCLB_RefSD = 5e-5*PatientAge+0.0121; 
Mean_Edges = linspace(0,0.14,100);
ADC_vec = ADCmap;
ADC_vec(final_mask==0)=[];
[bin_count,~] = histcounts(ADC_vec,Mean_Edges);
y = 0:0.001:1;
f = exp(-(y-ADCLB_RefMean).^2./(2*ADCLB_RefSD^2))./(ADCLB_RefSD*sqrt(2*pi));
f = f./max(f(:));
f = f.*max(bin_count(:));
plot(y,f,'k-.','LineWidth',2);
str_MeanADC_Ref = sprintf('Ref: Mean = %.3f', ADCLB_RefMean);
str_STDADC_Ref = sprintf(' ± %.3f', ADCLB_RefSD);
str_legend2=[str_MeanADC_Ref,str_STDADC_Ref];
legend(str_legend,str_legend2);
saveas(gca,'ADC_Histogram.png')
close all; 

cmap = colormap (jet);
cmap(1,:)=0;
close all
%% %% write tiff and read back Binneddiff maps
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis 
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving ADCmap Tiff...')
%ADCmap Binned
Proton_Image = zeros(size(ADCmap));
for slice=1:size(ADCmap,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(ADCmap(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); clim([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(ADCmap,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[outputpath,'ADCmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'ADCmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('Saving ADCmap Tiff Completed.')
close all;
% read tiff
cd(outputpath)
tiff_info = imfinfo('ADCmap.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
ADCcoloredmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('ADCmap.tif', ii);
    ADCcoloredmap(:,:,:,ii) = temp_tiff;
end
% S = orthosliceViewer((BinnedVentmap)); %colormap(SixBinMap);


%% 

% print out ADC maps
ADC_map_Montage = figure('Name','ADC maps');set(ADC_map_Montage,'WindowState','minimized');
montage(reshape(ADCmap,[size(ADCmap,1), size(ADCmap,2), 1, size(ADCmap,3)]),...
    'Size',[1 size(ADCmap,3)],'DisplayRange',[]);
% set(gca,'units','pixels'); % set the axes units to pixels
% x = get(gca,'position'); % get the position of the axes
% set(gcf,'units','pixels'); % set the figure units to pixels
% y = get(gcf,'position'); % get the figure position
% set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
% set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap(cmap); clim([0 0.14]);
xx = colorbar;
xx.Label.String = ' (cm^2/s)';
title_str1=['ADC maps-',ADCFittingType];
title(title_str1);
set(gca,'FontSize',12)
% print(title_str1,'-dpng','-r300');
saveas(gca,'ADCmap.png');
close all;
fprintf( ' ADC Analysis completed...\n' );

%% Save Report in a powerpoint file
%Diffusion Images
Diffimg_Montage = figure('Name','Diffusion Images');set(Diffimg_Montage,'WindowState','minimized');
montage(reshape(Ndiffimg,[size(Ndiffimg,1), size(Ndiffimg,2), size(Ndiffimg,4)* size(Ndiffimg,3)]),...
    'Size',[size(Ndiffimg,4) size(Ndiffimg,3)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
Diffusion_MontagePosition=get(gcf,'position');

%ADC maps
ADC_map_Montage = figure('Name','ADC maps');set(ADC_map_Montage,'WindowState','minimized');
montage(reshape(ADCmap,[size(ADCmap,1), size(ADCmap,2), 1, size(ADCmap,3)]),...
    'Size',[1 size(ADCmap,3)],'DisplayRange',[]);
% set(gca,'units','pixels'); % set the axes units to pixels
% x = get(gca,'position'); % get the position of the axes
% set(gcf,'units','pixels'); % set the figure units to pixels
% y = get(gcf,'position'); % get the figure position
% set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
% set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap(cmap); clim([0 0.14]);
xx = colorbar;
xx.Label.String = ' (cm^2/s)';
title_str1=['ADC maps-',ADCFittingType];
title(title_str1);
set(gca,'FontSize',12)
ADC_map_MontagePosition=get(gcf,'position'); 

nslices = size(ADCmap,3);
if nslices > 8
    nslices = 8;
end
scale_fac = 2*nslices;
% save ppt 
ReportTitle='Diffusion_Analysis';
cd(outputpath);
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
ppt_file_name = 'Diffusion_Analysis.pptx';
if exist(ppt_file_name, 'file') == 2            
    Global.exportToPPTX('open',ppt_file_name);
    Global.exportToPPTX('switchslide',1);
else            
    Global.exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
end 
%Add slides
%%%%%%%%%%%%%%%%% Subject/scan date/Processing Date  %%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); %Title Slide
Global.exportToPPTX('addtext',sprintf('Diffusion Analysis'), ...
    'Position',[0 0 16 4.5], ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'FontSize',72, ...
    'FontWeight','bold');
Global.exportToPPTX('addtext',sprintf(['Scan Date: ','','\nProcessing Date: ',datestr(date,29)]), ...
    'Position',[0 4.5 16 4.5], ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','center', ...
    'FontSize',36);
%%%%%%%%%%%%%%%%%%%%%%%%%%% b0 image and  ADC map %%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % diffusion images
Global.exportToPPTX('addpicture',Diffimg_Montage,...
    'Position',[1 0.5 scale_fac scale_fac*(Diffusion_MontagePosition(4)/Diffusion_MontagePosition(3))]);

Global.exportToPPTX('addtext',sprintf('Diffusion Images'),'Position',[0 0 5 1],'Color','b','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Slices'),'Position',[6 0.1 5 1],'Color','r','FontSize',20);
Global.exportToPPTX('addtext',sprintf(['bvalues (SNR= ',num2str(SNR_vec'),')']),'Position',[0 3 1.1 8],'Color','r','FontSize',20);

Global.exportToPPTX('addslide'); % ADC maps
Global.exportToPPTX('addpicture',ADC_map_Montage,...
    'Position',[0 1 scale_fac scale_fac*(ADC_map_MontagePosition(4)/ADC_map_MontagePosition(3))]);
Global.exportToPPTX('addtext',sprintf(title_str1),'Position',[6 0 5 1],'Color','b','FontSize',25);

%Histogram 
ADC_hist=imread('ADC_Histogram.png');
Global.exportToPPTX('addpicture',ADC_hist,'Position',[0 3.5 10 10*(size(ADC_hist,1)/size(ADC_hist,2))]);

% save all
Global.exportToPPTX('save',fullfile(outputpath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved'); 
close all;
%% save maps in mat file
save_data=[outputpath,'ADC_Analysis','.mat'];
save(save_data);   

%% 
fprintf('ADC Analysis compeleted\n'); 
Run_ADC_Time=toc/60

end


% Please add updates here
% 
% 
% 
% 