function [Ventilation] = calculate_LB_VDP(MR, maskarray, maskarraytrachea, ventmean, ventstd, parentPath, N4_bias_analysis, Overall_SNR,Ventilation,Proton,MainInput)
%% MATLAB script to perform Linear Binning VDP calculation

% Authors: Joseph Plummer
%         Abdullah Bdaiwi
% Date: 05/09/2021.

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% 1) Define inputs of lung data:

% if nargin == 7
%     ventmean = ventmean;
%     ventstd = ventstd;
%     outputPath = parentPath;
%     N4_bias_analysis = N4_bias_analysis;
% elseif nargin == 4
%     ventmean = 0.5421; % Default (non-N4 corrected, based on data up to Jan 2021)
%     ventstd = 0.1509; % Default (non-N4 corrected, based on data up to Jan 2021)
% else
%     disp('User must supply inputs for image data, masked-trachea data (1s and 0s), ventilation distribution mean, and standard deviation, the file path name, and whether or not N4 bias correction was applied.')
%     return
% end

% Define output folder names based on whether the images are corrected or not:
% if N4_bias_analysis == "no"
%     foldername = ["Linear Binning Ventilation Analysis\"];
% elseif N4_bias_analysis == "yes"
%     foldername = ["N4 Bias Linear Binning Ventilation Analysis\"];
% else
%     disp('User must declare whether or they are analyzing N4 corrected images. Check settings/input_params(), and make sure N4_bias_analysis = yes or no.')
% end
cd(parentPath)
foldername = "VDP Analysis\";
mkdir(foldername)
LB_outputpath = char(foldername);
LB_outputpath = [parentPath, LB_outputpath];
cd(LB_outputpath)

%% 2) Perform N4 bias correction OR RF correction:
if N4_bias_analysis == "no"
    % Perform RF correction on the lungs using the trachea as the max signal:
    NormMR = VentilationFunctions.rfCorrection(MR,maskarraytrachea);
    NormMR2 = NormMR.*(maskarraytrachea > 0); % Normalize
    VentImage = MR;
elseif N4_bias_analysis == "yes"
    % Run the N4 bias correction on the images WITH the trachea:
    [N4, Bias] = VentilationFunctions.N4_bias_correction(MR, maskarraytrachea, parentPath);
    NormMR2 = N4.*(maskarraytrachea > 0); % Normalize
    VentImage = N4;
end

% Make a logical array of the trachea:
maskarray4 = logical(maskarraytrachea);

%% 3) Generate bins for the linear binned data:
% Divide the lungs up into percentiles based on the trachea signal:
VentMax2 = prctile((NormMR2(maskarray4(:))),99.0); % 99.5th percentile

% Re-scale the images:
ScaledVentImage2 = NormMR2/VentMax2;
ScaledVentImage2(ScaledVentImage2 > 1) = 1;
ScaledVentImagenan2 = ScaledVentImage2;
ScaledVentImagenan2(ScaledVentImagenan2 == 0) = nan;

% % Display the histogram of the results:
% figure, histogram(ScaledVentImagenan2, 75);
% disp('Scaling Ventilation Image Completed.')

% Method for calculating mean and std of an image:
% mean(ScaledVentImage2(ScaledVentImage2 > 0))        
% std(ScaledVentImage2(ScaledVentImage2 > 0))  

SixBinMap = [1 0 0  %Red %Used for Vent and RBC
             1 0.7143 0  %Orange
             0.4 0.7 0.4  %Green1
             0 1 0  %Green2
             0 0.57 0.71  %Blue1
             0 0 1  %Blue2
             ];

% Define thresholds based on standard deviations:         
VentThresh2 = [ventmean-2*ventstd ,ventmean-1*ventstd ,ventmean,ventmean+ventstd,ventmean+2*ventstd];

% Bin the scaled image data accordingly:
VentBinMap2 = zeros(size(ScaledVentImage2));
for ii = 1:size(ScaledVentImage2,3)
    for jj = 1:size(ScaledVentImage2,2)
        for kk = 1:size(ScaledVentImage2,1)
            if maskarray(kk,jj,ii) == 0
                VentBinMap2(kk,jj,ii) = 0; %0 - black
            elseif abs(ScaledVentImage2(kk,jj,ii)) <= VentThresh2(1)
                VentBinMap2(kk,jj,ii) = 1; %1 - red
            elseif abs(ScaledVentImage2(kk,jj,ii)) <= VentThresh2(2)
                VentBinMap2(kk,jj,ii) = 2; %2 - orange
            elseif abs(ScaledVentImage2(kk,jj,ii)) <= VentThresh2(3)
                VentBinMap2(kk,jj,ii) = 3; %3 - green 1
            elseif abs(ScaledVentImage2(kk,jj,ii)) <= VentThresh2(4)
                VentBinMap2(kk,jj,ii) = 4; %4 - green 2
            elseif abs(ScaledVentImage2(kk,jj,ii)) <= VentThresh2(5)
                VentBinMap2(kk,jj,ii) = 5; %5 - blue 1
            elseif abs(ScaledVentImage2(kk,jj,ii)) > VentThresh2(5)
                VentBinMap2(kk,jj,ii) = 6; %6 - blue 2
            end
        end
    end
end

% Apply the lung & trachea masks to the binned data:
VentBinMask2 = logical((VentBinMap2-1).*maskarray4);

disp('Quantifying Images...')
% find mean and std of ADC in each bin 
Sum_all_pixels=sum(maskarray(:)==1);

Low1VentBinMap=zeros(size(VentBinMap2));
Low1VentBinMap(VentBinMap2==1)=1;
VentLow1Percent=(sum(Low1VentBinMap(:))/Sum_all_pixels)*100;

Low2VentBinMap=zeros(size(VentBinMap2));
Low2VentBinMap(VentBinMap2==2)=1;
VentLow2Percent=(sum(Low2VentBinMap(:))/Sum_all_pixels)*100;

Normal1VentBinMap=zeros(size(VentBinMap2));
Normal1VentBinMap(VentBinMap2==3)=1;
VentNormal1Percent=(sum(Normal1VentBinMap(:))/Sum_all_pixels)*100;

Normal2VentBinMap=zeros(size(VentBinMap2));
Normal2VentBinMap(VentBinMap2==4)=1;
VentNormal2Percent=(sum(Normal2VentBinMap(:))/Sum_all_pixels)*100;

High1VentBinMap=zeros(size(VentBinMap2));
High1VentBinMap(VentBinMap2==5)=1;
VentHigh1Percent=(sum(High1VentBinMap(:))/Sum_all_pixels)*100;

High2VentBinMap=zeros(size(VentBinMap2));
High2VentBinMap(VentBinMap2==6)=1;
VentHigh2Percent=(sum(High2VentBinMap(:))/Sum_all_pixels)*100;


%% 4) Quantify the ventilation via linear binning method:
LB_mean = mean(ScaledVentImage2(maskarray4(:)));
LB_std = std(ScaledVentImage2(maskarray4(:)));
LB_VDP = sum(VentBinMap2(:) == 1)/sum(maskarray4(:) == 1) * 100;
VentLowPercent = sum(VentBinMap2(:) == 2)/sum(maskarray4(:) == 1) * 100;
VentHighPercent = sum(VentBinMap2(:) > 4)/sum(maskarray4(:) == 1) * 100;

%% 5) Visualize the ventilation linear binning maps:
% Match each lung bin to their color:
VentDefect = (VentBinMap2 == 1) * 255;
LowVent = (VentBinMap2 == 2) * 255;
NormalVent1 = (VentBinMap2 == 3) * 255;
NormalVent2 = (VentBinMap2 == 4) * 255;
HighVent1 = (VentBinMap2 == 5) * 255;
HighVent2 = (VentBinMap2 == 6) * 255;

% Generate color arrays for each signal bin:
cmVentDefect=[0 0 0
              1 0 0];
array4dVentDefect = zeros(size(VentDefect,1), size(VentDefect,2),3, size(VentDefect,3));
for slice = 1:size(VentDefect,3)
    array4dVentDefect(:,:,:,slice) =  ind2rgb(VentDefect(:,:,slice),cmVentDefect);
end
%imshow(array4dVentDefect(:,:,:,6),[]);

cmLowVent=[0 0 0
          1 0.7143 0];
array4dLowVent = zeros(size(LowVent,1), size(LowVent,2),3, size(LowVent,3));
for slice = 1:size(LowVent,3)
    array4dLowVent(:,:,:,slice) =  ind2rgb(LowVent(:,:,slice),cmLowVent);
end
%imshow(array4dLowVent(:,:,:,6),[]);

cmNormalVent1=[0 0 0
              0.4 0.7 0.4];
array4dNormalVent1 = zeros(size(NormalVent1,1), size(NormalVent1,2),3, size(NormalVent1,3));
for slice = 1:size(NormalVent1,3)
    array4dNormalVent1(:,:,:,slice) =  ind2rgb(NormalVent1(:,:,slice),cmNormalVent1);
end
%imshow(array4dNormalVent1(:,:,:,6),[]);


cmNormalVent2=[0 0 0
               0 1 0];
array4dNormalVent2 = zeros(size(NormalVent2,1), size(NormalVent2,2),3, size(NormalVent2,3));
for slice = 1:size(NormalVent2,3)
    array4dNormalVent2(:,:,:,slice) =  ind2rgb(NormalVent2(:,:,slice),cmNormalVent2);
end
%imshow(array4dNormalVent2(:,:,:,6),[]);

cmHighVent1=[0 0 0
             0 0.57 0.71];
array4dHighVent1 = zeros(size(HighVent1,1), size(HighVent1,2),3, size(HighVent1,3));
for slice = 1:size(HighVent1,3)
    array4dHighVent1(:,:,:,slice) =  ind2rgb(HighVent1(:,:,slice),cmHighVent1);
end
%imshow(array4dHighVent1(:,:,:,6),[]);

cmHighVent2=[0 0 0
             0 0 1];
array4dHighVent2 = zeros(size(HighVent2,1), size(HighVent2,2),3, size(HighVent2,3));
for slice = 1:size(HighVent2,3)
    array4dHighVent2(:,:,:,slice) =  ind2rgb(HighVent2(:,:,slice),cmHighVent2);
end
%imshow(array4dHighVent2(:,:,:,6),[]);

% Add all parts together to make the full colormap:
colorVentmap = array4dVentDefect + array4dLowVent + array4dNormalVent1 + ...
    array4dNormalVent2 + array4dHighVent1 + array4dHighVent2;

% Make a montage:
LinearBinningVDPMontage = montage(reshape(colorVentmap,[size(colorVentmap,1),...
    size(colorVentmap,2),size(colorVentmap,3),size(colorVentmap,4)]),...
    'Size',[1 size(colorVentmap,4)],'DisplayRange',[]);
title('Linear Binning VDP')
savedLinearBinningVDPMontage = get(LinearBinningVDPMontage,'CData');

% Write out the montage using the input details above:
FileName = [foldername + "Linear_Binning_Ventilation_Defect_Montage_with_Trachea.png"];
fullFileName = fullfile(parentPath,FileName);
imwrite(savedLinearBinningVDPMontage,fullFileName,'png');
%% Output mask images with proton overlays:
cd(LB_outputpath);
cus_colormap = zeros(100,3);
% colorgradient = 0:0.01:0.99;
cus_colormap(:,1) = 0;
cus_colormap(:,2) = 1;
cus_colormap(:,3) = 1;

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')
%Mask images
for slice=1:size(maskarray,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton.ProtonRegistered(:,:,slice))),squeeze(Ventilation.LungMask(:,:,slice)),[1,100],[0,0.99*max(Proton.ProtonRegistered(:))],cus_colormap,0.6,gca);
    colormap(gca,cus_colormap)   
%     X = print('-RGBImage',['-r',num2str(size(VentBinMap2,2)/2)]);%2 inches
%      hImage = get( gca, 'Children' ); 
%      X = hImage.CData;   
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[LB_outputpath,'MaskRegistered.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[LB_outputpath,'MaskRegistered.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('MaskRegistered Tiff Completed.')
close all;
% read tiff
cd(LB_outputpath)
tiff_info = imfinfo('MaskRegistered.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
MaskRegistered = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('MaskRegistered.tif', ii);
    MaskRegistered(:,:,:,ii) = temp_tiff;
end
MaskRegistered = permute(MaskRegistered,[1 2 4 3]);

S = orthosliceViewer((MaskRegistered)); %colormap(SixBinMap);


%% 6) Calculate Histograms:
disp('Calculating Histograms...')
Vent_vec = ScaledVentImage2;
Vent_vec(maskarray==0)=[];
Vent_vec(Vent_vec == 0) = NaN;

healthy_mean_Osc=ventmean;
healthy_std_Osc=ventstd;
norm_lim_high = healthy_mean_Osc + healthy_std_Osc;
norm_lim_low = healthy_mean_Osc - healthy_std_Osc;
Mean_Edges = linspace(0,1,100);
[bin_count,edges] = histcounts(Vent_vec,Mean_Edges);
Norm_Factor_1 = sum(bin_count);

Healthy_Mean_Dist = normpdf(Mean_Edges,healthy_mean_Osc,healthy_std_Osc);
NormFactor = sum(Healthy_Mean_Dist(:));
Healthy_Mean_Dist = Healthy_Mean_Dist./NormFactor;
Healthy_Mean_Dist = Healthy_Mean_Dist*Norm_Factor_1;

% Ventilation Histogram
% VentHistFig = figure('Name','Ventilation Histogram');
% VentEdges = [-0.5, linspace(0,1,50) 1.5];
% set(VentHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
All_Bins=VentThresh2;
Bin1 = All_Bins(1);
Bin2 = All_Bins(2);
Bin3 = All_Bins(3);
Bin4 = All_Bins(4);
Bin5 = All_Bins(5);


pos = [1 1 10 6];
HistogramFig_Health = figure('Name','Siemens 2101 Histogram');
set(HistogramFig_Health,'color','white','Units','inches','Position',pos)
histogram(Vent_vec,Mean_Edges,'FaceColor','w','FaceAlpha',1);
hold on
histogram(Vent_vec(Vent_vec<=Bin1),Mean_Edges,'FaceColor',[1 0 0],'FaceAlpha',1)
histogram(Vent_vec(Vent_vec>Bin1 & Vent_vec<=Bin2),Mean_Edges,'FaceColor',[1 0.7143 0],'FaceAlpha',1)
histogram(Vent_vec(Vent_vec>Bin2 & Vent_vec<=Bin3),Mean_Edges,'FaceColor',[0.4 0.7 0.4],'FaceAlpha',1)
histogram(Vent_vec(Vent_vec>Bin3 & Vent_vec<=Bin4),Mean_Edges,'FaceColor',[0 1 0],'FaceAlpha',1)
histogram(Vent_vec(Vent_vec>Bin4 & Vent_vec<=Bin5),Mean_Edges,'FaceColor',[0 0.57 0.71],'FaceAlpha',1)
histogram(Vent_vec(Vent_vec>Bin5 ),Mean_Edges,'FaceColor',[0 0 1],'FaceAlpha',1)
hold on

% histogram((ScaledVentImage2(maskarray4(:))),VentEdges,'Normalization','probability','EdgeColor',[0 0.65 0],'FaceColor',[1 1 1],'FaceAlpha',0.75);

plot([norm_lim_low norm_lim_low],[0 5500],'k-.',[norm_lim_high,norm_lim_high],[0 5500],'k-.','LineWidth',2)
if max(bin_count)>= max(Healthy_Mean_Dist)
    ylim([0 max(bin_count)+100])
elseif max(bin_count)< max(Healthy_Mean_Dist)
    ylim([0  max(Healthy_Mean_Dist)+100])
end
% plot healthy dist.
y = 0:0.001:1;
f = exp(-(y-ventmean).^2./(2*ventstd^2))./(ventstd*sqrt(2*pi));
f = f./max(f(:));
f = f.*max(bin_count(:));
plot(y,f,'k-.','LineWidth',2);

xlim([0 1])
set(gca,'linewidth',1.5)
set(gca,'FontSize',12)

% Bins:
ylims = ylim;
xlims = xlim;
pat = patch([xlims(1) xlims(1) VentThresh2(1) VentThresh2(1)],[0 ylims(2) ylims(2) 0],SixBinMap(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(VentThresh2,2)
    pat = patch([VentThresh2(bin-1) VentThresh2(bin-1) VentThresh2(bin) VentThresh2(bin)],[0 ylims(2) ylims(2) 0],SixBinMap(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([VentThresh2(5) VentThresh2(5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinMap(6,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
% axis([-0.01 1.01 -inf inf])
legend1 = sprintf('VDP:%0.1f%%  ',LB_VDP);
legend2 = sprintf('Vent Low Percent:%0.1f%%  ' ,VentLowPercent);
legend3 = sprintf('Vent High Percent:%0.1f%%  ' ,VentHighPercent);
legend4 = sprintf(' Mean:%0.4f%    '  ,LB_mean);
legend5 = sprintf(' StdDev:%0.4f%    '  ,LB_std);
legend([legend1  legend2  legend3  legend4   legend5]);
title1 = sprintf('Ventilation Defect Percentage %0.1f%%',LB_VDP);
title({title1},'Fontweight','bold','FontSize',12)
set(gca,'FontSize',14)
% print([foldername + 'Ventilation_Histogram'],'-dpng','-r300');
cd(LB_outputpath)
saveas(gca,'LB_Ventilation_Histogram.png');
hold off
%% %% write tiff and read back BinnedVent maps
%Tiffs (Binned Images)
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')
%Vent Binned
for slice=1:size(VentBinMap2,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton.ProtonRegistered(:,:,slice))),squeeze(VentBinMap2(:,:,slice)),[1,6],[0,0.99*max(Proton.ProtonRegistered(:))],SixBinMap,1,gca);
    colormap(gca,SixBinMap)   
%     X = print('-RGBImage',['-r',num2str(size(VentBinMap2,2)/2)]);%2 inches
%      hImage = get( gca, 'Children' ); 
%      X = hImage.CData;   
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[LB_outputpath,'BinnedVent.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[LB_outputpath,'BinnedVent.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('Saving Vent Tiff Completed.')
close all;
% read tiff
cd(LB_outputpath)
tiff_info = imfinfo('BinnedVent.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
BinnedVentmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedVent.tif', ii);
    BinnedVentmap(:,:,:,ii) = temp_tiff;
end
% S = orthosliceViewer((BinnedVentmap)); %colormap(SixBinMap);
%% save result in a powerpoint file // Abdullah 1/27/2020
%% 
cd(parentPath)
close all;
VentImage = VentImage./max(VentImage(:));
for m = 1:size(VentImage,3)
    VentImage2(:,:,m) = imadjust(VentImage(:,:,m));
    VentImage2(:,:,m) = adapthisteq(VentImage2(:,:,m));
end
% %how many slices to incluse in the final powerpoint file
sl_include=zeros(1,size(maskarray,3));
for sl=1:size(maskarray,3)
    choose_slice=maskarray(:,:,sl);       
    logica_test2=sum(choose_slice(:));
    if logica_test2>1
         sl_include(sl)=sl;
    else
         sl_include(sl)=0;
    end 
end
sl_include( :, ~any(sl_include,1) ) = [];  %columns
sl_1=sl_include(1);
sl_end=sl_include(end);
NumSliceView = length(sl_include);
if NumSliceView > 16
    NumSliceView = 16;
end
% if you want the un-segemted slices not to show up on the powepoint, run
% the above code and add this line (,'Indices', sl_1:sl_end) to the end of
% each Montage. 
%%% get rid of the gray frame in each montage
VentscaledImage = VentImage2(:,:,sl_1:sl_end);
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(reshape(VentscaledImage,[size(VentscaledImage,1), size(VentscaledImage,2), 1, size(VentscaledImage,3)]),...
    'Size',[1 size(VentscaledImage,3)],'DisplayRange',[0 max(VentscaledImage(:))]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
VentMontagePosition=get(gcf,'position'); 

ProtonMontage = figure('Name','Proron Image');set(ProtonMontage,'WindowState','minimized');
if strcmp(MainInput.NoProtonImage,'yes') == 1 
    ProtonImage = Proton.ProtonRegisteredColored(:,:,sl_1:sl_end);
    montage(reshape(ProtonImage,[size(ProtonImage,1), size(ProtonImage,2), 1, size(ProtonImage,3)]),...
    'Size',[1 size(ProtonImage,3)],'DisplayRange',[0 1]);
else 
    ProtonImage = permute(Proton.ProtonRegisteredColored, [1 2 4 3]);
    ProtonImage = ProtonImage(:,:,:,sl_1:sl_end);
    montage(reshape(ProtonImage,[size(ProtonImage,1), size(ProtonImage,2),size(ProtonImage,3),...
    size(ProtonImage,4)]),'Size',[1 size(ProtonImage,4)],'DisplayRange',[]);
end
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

MaskRegistered = permute(MaskRegistered,[1 2 4 3]);
maskarrayimg = MaskRegistered(:,:,:,sl_1:sl_end);
MaskMontage = figure('Name','Mask');set(MaskMontage,'WindowState','minimized');
montage(reshape(maskarrayimg,[size(maskarrayimg,1), size(maskarrayimg,2), size(maskarrayimg,3), size(maskarrayimg,4)]),...
    'Size',[1 size(maskarrayimg,4)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
% MaskMontagePosition = get(gcf,'position'); 

% BinnedVentmap = permute(BinnedVentmap,[1 2 4 3]);
Ventcolormap = BinnedVentmap(:,:,:,sl_1:sl_end);
DefectArrayMontage = figure('Name','Defects');set(DefectArrayMontage,'WindowState','minimized');
montage(reshape(Ventcolormap,[size(Ventcolormap,1), size(Ventcolormap,2),size(Ventcolormap,3),...
    size(Ventcolormap,4)]),'Size',[1 size(Ventcolormap,4)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
% DefectArrayMontagePosition=get(gcf,'position'); 

% save ppt 
ReportTitle = 'Ventilation_Analysis';
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
ppt_file_name = 'Ventilation_Analysis.pptx';
if isfile(ppt_file_name)           
    disp('file existed');
    Global.exportToPPTX('open',ppt_file_name);
    Global.exportToPPTX('switchslide',1);
else            
    Global.exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
end 
%Add slides
Global.exportToPPTX('addslide'); % Image/mask/VDP
Global.exportToPPTX('addtext',sprintf(char(foldername)),'Position',[5 0 7 1],'Color','b','FontSize',25);
Global.exportToPPTX('addpicture',VentMontage,'Position',...
    [0 0.4 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);
Global.exportToPPTX('addpicture',ProtonMontage,'Position',...
    [0 1.6 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);
Global.exportToPPTX('addpicture',MaskMontage,'Position',...
    [0 2.8 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);
Global.exportToPPTX('addpicture',DefectArrayMontage,'Position',...
    [0 4 NumSliceView NumSliceView*(VentMontagePosition(4)/VentMontagePosition(3))]);
LB_Venthist = imread([parentPath , char(foldername) ,'LB_Ventilation_Histogram.png']);    
Global.exportToPPTX('addpicture',LB_Venthist,'Position',[0 5 6.8 6.8*(size(LB_Venthist,1)/size(LB_Venthist,2))]);

Global.exportToPPTX('addtext',['Vent Images (SNR=',num2str(Overall_SNR),')'],'Position',[0 0 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Registered Proton Images','Position',[0 1.2 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
Global.exportToPPTX('addtext','Mask','Position',[0 2.4 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
vdp_title = sprintf('Linear Binning VDP= %0.1f%%',LB_VDP);
Global.exportToPPTX('addtext',vdp_title,'Position',[0 3.6 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');

% add table/bins summary
BinTable   = { ...
    '',{'LVR1','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'LVR2','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Normal1','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Normal2','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'HVR1','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'HVR2','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Bins','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    'Percent (%)',num2str(VentLow1Percent,'%1.3f'),num2str(VentLow2Percent,'%1.3f'),num2str(VentNormal1Percent,'%1.3f'),num2str(VentNormal2Percent,'%1.3f'),num2str(VentHigh1Percent,'%1.3f'),num2str(VentHigh2Percent,'%1.3f')};
Global.exportToPPTX('addtable',BinTable,'Position',[7.5 6 7.75 2],'Vert','middle','Horiz','center','FontSize',12);


Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');               

close all;
%% save excel file
% save mean and std in a table and export as excel file in data location 
headers={'summary';'Bin1';'Bin2'; 'Bin3'; 'Bin4'; 'Bin5'; 'Bin6'}';
% MeanSTDNames = {'Mean(cm2/s)';'Std';'Percent(%)'};
Vent_BinTable = table({'Percent(%)'},...
                     [VentLow1Percent],...
                     [VentLow2Percent],...
                     [VentNormal1Percent],...
                     [VentNormal2Percent],...
                     [VentHigh1Percent],...
                     [VentHigh2Percent]);
Vent_BinTable.Properties.VariableNames = headers;
excel_file_name = 'LB_VDP_Summary.xlsx';
writetable(Vent_BinTable,[LB_outputpath,excel_file_name],'Sheet',1)


%% 7) Save workspace data:
Ventilation.LB_VDP = LB_VDP;
Ventilation.LB_mean = LB_mean;
Ventilation.LB_std = LB_std;
Ventilation.VentscaledImage = VentscaledImage;
Ventilation.Ventcolormap = Ventcolormap;
Ventilation.LB_Venthist = LB_Venthist;
Ventilation.VentLow1Percent = round(VentLow1Percent,2);
Ventilation.VentLow2Percent = round(VentLow2Percent,2);
Ventilation.VentNormal1Percent = round(VentNormal1Percent,2);
Ventilation.VentNormal2Percent = round(VentNormal2Percent,2);
Ventilation.VentHigh1Percent = round(VentHigh1Percent,2);
Ventilation.VentHigh2Percent = round(VentHigh2Percent,2);

save_title = [LB_outputpath, 'LinearBinningAnalysis.mat'];
save(save_title);
end % end of function
