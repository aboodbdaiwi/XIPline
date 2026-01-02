function [Diffusion] = LinearBinningADC(Diffusion, MainInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% Linear Binning ADC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Inputs:
%
%   Outputs:
%       success - 1 if ran successfully, 0 if not

%   Package: https://github.com/aboodbdaiwi/XIPline
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 
%
% This code is based on the linear Binning Method by Ziyi Wang
% (https://pubmed.ncbi.nlm.nih.gov/30024058/)
% PMID: 30024058 PMCID: PMC6318005 DOI: 10.1002/mrm.27377
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Images = Diffusion.diffimg;
ADCmap = Diffusion.ADCmap;
maskarray = Diffusion.lung_mask;
age = Diffusion.PatientAge;
DataLocation = Diffusion.outputpath;
HealthyMean = Diffusion.ADCLB_RefMean;
HealthySD = Diffusion.ADCLB_RefSD;

success = 1;
disp('Linear Binning ADC Processing Beginning...')
tic;
%check folder is valid
if DataLocation == 0 %exit if no folder selected
    disp('Error. Must select a folder')
    success = 0;
    return;
end
if isfolder(char(DataLocation)) == 0  %exit if input is not a folder
   disp('Error. Input must be a folder')
   success = 0;
   return;
end
%% Normalize data
Images = (Images./max(Images(:))).*100;

%% Determine healthy distribution based on the subject's age 
% HealthyADCmean=0.0002*age+0.029;
% HealthyADCstd=5e-5*age+0.0121;

try
    HealthyADCmean = eval(char(HealthyMean));
catch
    HealthyADCmean = HealthyMean;
end
HealthyADCmean = HealthyADCmean(1);
try
    HealthyADCstd = eval(HealthySD);
catch
    HealthyADCstd = HealthySD;
end
HealthyADCstd = HealthyADCstd(1);
ADCThresh = [HealthyADCmean-2*HealthyADCstd ,HealthyADCmean-1*HealthyADCstd ,HealthyADCmean, HealthyADCmean+ HealthyADCstd, HealthyADCmean+2*HealthyADCstd];

DiffBinMap2 = zeros(size(ADCmap));
for ii = 1:size(ADCmap,3)
    for jj = 1:size(ADCmap,2)
        for kk = 1:size(ADCmap,1)
            if maskarray(kk,jj,ii) == 0
                DiffBinMap2(kk,jj,ii) = 0; %0 - black
            elseif (ADCmap(kk,jj,ii)) <= ADCThresh(1)
                DiffBinMap2(kk,jj,ii) = 1; %1 - red
            elseif (ADCmap(kk,jj,ii)) <= ADCThresh(2)
                DiffBinMap2(kk,jj,ii) = 2; %2 - orange
            elseif (ADCmap(kk,jj,ii)) <= ADCThresh(3)
                DiffBinMap2(kk,jj,ii) = 3; %3 - green 1
            elseif (ADCmap(kk,jj,ii)) <= ADCThresh(4)
                DiffBinMap2(kk,jj,ii) = 4; %4 - green 2
            elseif (ADCmap(kk,jj,ii)) <= ADCThresh(5)
                DiffBinMap2(kk,jj,ii) = 5; %5 - blue 1
            elseif (ADCmap(kk,jj,ii)) > ADCThresh(5)
                DiffBinMap2(kk,jj,ii) = 6; %6 - blue 2
            end
        end
    end
end

DiffBinMask2 = logical((DiffBinMap2-1).*maskarray);

disp('Quantifying Images...')
%Diffusion Quantification
ADC_vec=ADCmap;
ADC_vec(maskarray==0)=[];
DiffMean = mean(ADC_vec);
DiffStd = std(ADC_vec);

% find mean and std of ADC in each bin 
Sum_all_pixels=sum(maskarray(:)==1);

Low1DiffBinMap=zeros(size(DiffBinMap2));
Low1DiffBinMap(DiffBinMap2==1)=1;
DiffLow1Percent=(sum(Low1DiffBinMap(:))/Sum_all_pixels)*100;

Low2DiffBinMap=zeros(size(DiffBinMap2));
Low2DiffBinMap(DiffBinMap2==2)=1;
DiffLow2Percent=(sum(Low2DiffBinMap(:))/Sum_all_pixels)*100;

Normal1DiffBinMap=zeros(size(DiffBinMap2));
Normal1DiffBinMap(DiffBinMap2==3)=1;
DiffNormal1Percent=(sum(Normal1DiffBinMap(:))/Sum_all_pixels)*100;

Normal2DiffBinMap=zeros(size(DiffBinMap2));
Normal2DiffBinMap(DiffBinMap2==4)=1;
DiffNormal2Percent=(sum(Normal2DiffBinMap(:))/Sum_all_pixels)*100;

High1DiffBinMap=zeros(size(DiffBinMap2));
High1DiffBinMap(DiffBinMap2==5)=1;
DiffHigh1Percent=(sum(High1DiffBinMap(:))/Sum_all_pixels)*100;

High2DiffBinMap=zeros(size(DiffBinMap2));
High2DiffBinMap(DiffBinMap2==6)=1;
DiffHigh2Percent=(sum(High2DiffBinMap(:))/Sum_all_pixels)*100;

DiffLowPercent =DiffLow1Percent+DiffLow2Percent;
DiffNormalPercent = DiffNormal1Percent+DiffNormal2Percent;
DiffHighPercent = DiffHigh1Percent+DiffHigh2Percent;

ADC_Bin1=ADCmap;
ADC_Bin1(Low1DiffBinMap==0)=[];
Mean_ADC_Bin1=mean(ADC_Bin1);
Std_ADC_Bin1=std(ADC_Bin1);

ADC_Bin2=ADCmap;
ADC_Bin2(Low2DiffBinMap==0)=[];
Mean_ADC_Bin2=mean(ADC_Bin2);
Std_ADC_Bin2=std(ADC_Bin2);

ADC_Bin3=ADCmap;
ADC_Bin3(Normal1DiffBinMap==0)=[];
Mean_ADC_Bin3=mean(ADC_Bin3);
Std_ADC_Bin3=std(ADC_Bin3);

ADC_Bin4=ADCmap;
ADC_Bin4(Normal2DiffBinMap==0)=[];
Mean_ADC_Bin4=mean(ADC_Bin4);
Std_ADC_Bin4=std(ADC_Bin4);

ADC_Bin5=ADCmap;
ADC_Bin5(High1DiffBinMap==0)=[];
Mean_ADC_Bin5=mean(ADC_Bin5);
Std_ADC_Bin5=std(ADC_Bin5);

ADC_Bin6=ADCmap;
ADC_Bin6(High2DiffBinMap==0)=[];
Mean_ADC_Bin6=mean(ADC_Bin6);
Std_ADC_Bin6=std(ADC_Bin6);

% save mean and std in a table and export as excel file in data location 
headers={'summary';'Bin1';'Bin2'; 'Bin3'; 'Bin4'; 'Bin5'; 'Bin6'}';
% MeanSTDNames = {'Mean(cm2/s)';'Std';'Percent(%)'};
ADC_BinTable = table({'ADC(cm2/s)';'std'; 'Percent(%)'},...
                     [Mean_ADC_Bin1;Std_ADC_Bin1;DiffLow1Percent],...
                     [Mean_ADC_Bin2;Std_ADC_Bin2;DiffLow2Percent],...
                     [Mean_ADC_Bin3;Std_ADC_Bin3;DiffNormal1Percent],...
                     [Mean_ADC_Bin4;Std_ADC_Bin4;DiffNormal2Percent],...
                     [Mean_ADC_Bin5;Std_ADC_Bin5;DiffHigh1Percent],...
                     [Mean_ADC_Bin6;Std_ADC_Bin6;DiffHigh2Percent]);
ADC_BinTable.Properties.VariableNames = headers;
excel_file_name = ['\','LB_ADC_Summary.xlsx'];
writetable(ADC_BinTable,[DataLocation,excel_file_name],'Sheet',1)

% saveas(gca,'MeanSTDPercent_ADCBins.png');
% close all;
%% create Histogram
healthy_mean_Osc=HealthyADCmean;
healthy_std_Osc=HealthyADCstd;
norm_lim_high = healthy_mean_Osc + healthy_std_Osc;
norm_lim_low = healthy_mean_Osc - healthy_std_Osc;
Mean_Edges = linspace(0,0.14,100);
[bin_count,edges] = histcounts(ADC_vec,Mean_Edges);
Norm_Factor_1 = sum(bin_count);

Healthy_Mean_Dist = normpdf(Mean_Edges,healthy_mean_Osc(1),healthy_std_Osc(1));
NormFactor = sum(Healthy_Mean_Dist(:));
Healthy_Mean_Dist = Healthy_Mean_Dist./NormFactor;
Healthy_Mean_Dist = Healthy_Mean_Dist*Norm_Factor_1;

% Bin1 = healthy_mean_Osc - 2*healthy_std_Osc;
% Bin2 = healthy_mean_Osc - healthy_std_Osc;
% Bin3 = healthy_mean_Osc;
% Bin4 = healthy_mean_Osc + healthy_std_Osc;
% Bin5 = healthy_mean_Osc + 2*healthy_std_Osc;
% Bin6 = healthy_mean_Osc + 3*healthy_std_Osc;
% Bin7 = healthy_mean_Osc + 4*healthy_std_Osc;
% Bin8 = healthy_mean_Osc + 5*healthy_std_Osc;
% All_Bins=[Bin1,Bin2,Bin3,Bin4,Bin5,Bin6,Bin7,Bin8];
All_Bins=ADCThresh;
Bin1 = All_Bins(1);
Bin2 = All_Bins(2);
Bin3 = All_Bins(3);
Bin4 = All_Bins(4);
Bin5 = All_Bins(5);

pos = [1 1 15 6];
HistogramFig_Health = figure('Name','LB Histogram');
set(HistogramFig_Health,'color','white','Units','inches','Position',pos)
histogram(ADC_vec,Mean_Edges,'FaceColor','w','FaceAlpha',1);
hold on

histogram(ADC_vec(ADC_vec<=Bin1),Mean_Edges,'FaceColor',[0 0 1],'FaceAlpha',1)
histogram(ADC_vec(ADC_vec>Bin1 & ADC_vec<=Bin2),Mean_Edges,'FaceColor',[0 0.57 0.71],'FaceAlpha',1)
histogram(ADC_vec(ADC_vec>Bin2 & ADC_vec<=Bin3),Mean_Edges,'FaceColor',[0.4 0.7 0.4],'FaceAlpha',1)
histogram(ADC_vec(ADC_vec>Bin3 & ADC_vec<=Bin4),Mean_Edges,'FaceColor',[0 1 0],'FaceAlpha',1)
histogram(ADC_vec(ADC_vec>Bin4 & ADC_vec<=Bin5),Mean_Edges,'FaceColor',[1 0.7143 0],'FaceAlpha',1)
histogram(ADC_vec(ADC_vec>Bin5 ),Mean_Edges,'FaceColor',[1 0 0],'FaceAlpha',1)

% histogram(ADC_vec(ADC_vec>Bin4 & ADC_vec<Bin5),Mean_Edges,'FaceColor',[184/255 226/255 145/255],'FaceAlpha',1)
% histogram(ADC_vec(ADC_vec>Bin5 & ADC_vec<Bin6),Mean_Edges,'FaceColor',[243/255 205/255 213/255],'FaceAlpha',1)
% histogram(ADC_vec(ADC_vec>Bin6 & ADC_vec<Bin7),Mean_Edges,'FaceColor',[225/255 129/255 162/255],'FaceAlpha',1)
% histogram(ADC_vec(ADC_vec>Bin7),Mean_Edges,'FaceColor',[197/255 27/255 125/255],'FaceAlpha',1)
% plot(Mean_Edges,Healthy_Mean_Dist,'k--','LineWidth',2);

plot([norm_lim_low norm_lim_low],[0 5500],'k-.',[norm_lim_high,norm_lim_high],[0 5500],'k-.','LineWidth',2)
if max(bin_count)>= max(Healthy_Mean_Dist)
    ylim([0 max(bin_count)+50])
elseif max(bin_count)< max(Healthy_Mean_Dist)
    ylim([0  max(Healthy_Mean_Dist)+50])
end

% plot healthy dist.
% y = 0:0.001:1;
% f = exp(-(y-HealthyADCmean).^2./(2*HealthyADCstd^2))./(HealthyADCstd*sqrt(2*pi));
% f = f./max(f(:));
% f = f.*max(bin_count(:));

f = DiffusionFunctions.ADC_RefGamma(age);
f = f ./ max(f(:));   % normalize to peak = 1
f = f .* max(bin_count(:));      % scale to peak = max bin count
x = linspace(0,0.14,400);
plot(x, f, 'k-.', 'LineWidth', 1.5);

xlim([0 0.14])
set(gca,'linewidth',1.5)
set(gca,'FontSize',12)

ylims = ylim;
xlims = xlim;
SixBinMap = [0 0 1  %Blue1
             0 0.57 0.71  %Blue2
             0.4 0.7 0.4  %Green1
             0 1 0  %Green2
             1 0.7143 0   %Orange  
             1 0 0   %Red 
             ];
         
% SixBinMap = [1 0 0  %Red
%              1 0.7143 0  %Orange
%              0.4 0.7 0.4  %Green1
%              0 1 0  %Green2
%              184/255 226/255 145/255  %Green3
%              243/255 205/255 213/255  %Light Pink
%              225/255 129/255 162/255  %Med Pink
%              197/255 27/255 125/255  %Dark Pink
%              1 1 1 % White (For regions of RBC Defect)
%              ];
         
pat = patch([xlims(1) xlims(1) All_Bins(1) All_Bins(1)],[0 ylims(2) ylims(2) 0],SixBinMap(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(All_Bins,2)
    pat = patch([All_Bins(bin-1) All_Bins(bin-1) All_Bins(bin) All_Bins(bin)],[0 ylims(2) ylims(2) 0],SixBinMap(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([All_Bins(5) All_Bins(5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinMap(6,:),'HandleVisibility','off');
% pat = patch([All_Bins(8) All_Bins(8) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinMap(8,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
xlim([0 0.14])
% axis([Mean_Edges(2) Mean_Edges(end-1) -inf inf])
legend1 = sprintf('LDP:%0.1f%%  ' ,DiffLowPercent);
legend2 = sprintf('NDP:%0.1f%%  ' ,DiffNormalPercent);
legend3 = sprintf('HDP:%0.1f%%  ' ,DiffHighPercent);
legend4 = sprintf('Mean: %f±%f ',round(DiffMean,4),round(DiffStd,4));

Mean_Healthy_predicted=sprintf('Healthy: %f±%f ',HealthyADCmean,HealthyADCstd);
lgd = legend([legend1  legend2  legend3  legend4],Mean_Healthy_predicted);
lgd.FontSize = 16;   % <-- bigger legend
title('Linear Binning ADC Histogram','Fontweight','bold','FontSize',15)
xlabel('ADC (cm^2/s)')
set(gca,'FontSize',20, 'FontWeight','bold')
hold off
saveas(gca,'LB_ADC_Histogram.png');
close all;

%% Creating color maps for Diffusion 

LowDiff1=(DiffBinMap2==1)*255;
LowDiff2=(DiffBinMap2==2)*255;
NormalDiff1=(DiffBinMap2==3)*255;
NormalDiff2=(DiffBinMap2==4)*255;
HighDiff1=(DiffBinMap2==5)*255;
HighDiff2=(DiffBinMap2==6)*255;

cmDiffDefect=[0 0 0
              0 0 1];
array4dLowDiff1 = zeros(size(LowDiff1,1), size(LowDiff1,2),3, size(LowDiff1,3));
for slice = 1:size(LowDiff1,3)
    array4dLowDiff1(:,:,:,slice) =  ind2rgb(LowDiff1(:,:,slice),cmDiffDefect);
end
cmLowDiff2=[0 0 0
          0 0.57 0.71];
array4dLowDiff2 = zeros(size(LowDiff2,1), size(LowDiff2,2),3, size(LowDiff2,3));
for slice = 1:size(LowDiff2,3)
    array4dLowDiff2(:,:,:,slice) =  ind2rgb(LowDiff2(:,:,slice),cmLowDiff2);
end
cmNormalDiff1=[0 0 0
              0.4 0.7 0.4];
array4dNormalDiff1 = zeros(size(NormalDiff1,1), size(NormalDiff1,2),3, size(NormalDiff1,3));
for slice = 1:size(NormalDiff1,3)
    array4dNormalDiff1(:,:,:,slice) =  ind2rgb(NormalDiff1(:,:,slice),cmNormalDiff1);
end
cmNormalDiff2=[0 0 0
               0 1 0];
array4dNormalDiff2 = zeros(size(NormalDiff2,1), size(NormalDiff2,2),3, size(NormalDiff2,3));
for slice = 1:size(NormalDiff2,3)
    array4dNormalDiff2(:,:,:,slice) =  ind2rgb(NormalDiff2(:,:,slice),cmNormalDiff2);
end
cmHighDiff1=[0 0 0
             1 0.7143 0]; 
array4dHighDiff1 = zeros(size(HighDiff1,1), size(HighDiff1,2),3, size(HighDiff1,3));
for slice = 1:size(HighDiff1,3)
    array4dHighDiff1(:,:,:,slice) =  ind2rgb(HighDiff1(:,:,slice),cmHighDiff1);
end
cmHighDiff2=[0 0 0
             1 0 0];
array4dHighDiff2 = zeros(size(HighDiff2,1), size(HighDiff2,2),3, size(HighDiff2,3));
for slice = 1:size(HighDiff2,3)
    array4dHighDiff2(:,:,:,slice) =  ind2rgb(HighDiff2(:,:,slice),cmHighDiff2);
end
BinnedADCmap=array4dLowDiff1+array4dLowDiff2+array4dNormalDiff1+array4dNormalDiff2+array4dHighDiff1+array4dHighDiff2;
LinearBinningADCMontage = montage(reshape(BinnedADCmap,[size(BinnedADCmap,1), size(BinnedADCmap,2),size(BinnedADCmap,3),size(BinnedADCmap,4)]),'Size',[1 size(BinnedADCmap,4)],'DisplayRange',[0 6]);
title('Linear Binning ADC','Fontweight','bold','FontSize',15)
colormap_Bins = [0 0 1  %Blue1
             0 0.57 0.71  %Blue2
             0.4 0.7 0.4  %Green1
             0 1 0  %Green2
             1 0.7143 0   %Orange  
             1 0 0   %Red 
             ];
colormap(colormap_Bins); colorbar('XTickLabel',{'0','LADC','LADC','NADC','NADC','HADC','HADC'},'XTick', 0:1:6);
saveas(gca,'LinearBinningADCMontage.png');
close all;

%% %% write tiff and read back Binneddiff maps
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis 
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving ADCcolormap Tiff...')
%ADCmap Binned
Proton_Image = zeros(size(DiffBinMap2));
for slice=1:size(DiffBinMap2,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(DiffBinMap2(:,:,slice)),[1,6],[0,0.99*max(Proton_Image(:))],colormap_Bins,1,gca);
    colormap(gca,colormap_Bins); %caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(DiffBinMap2,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedADCmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedADCmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end
disp('Saving BinnedADCmap Tiff Completed.')
close all;
% read tiff
cd(DataLocation)
tiff_info = imfinfo('BinnedADCmap.tif'); % return tiff structure, one element per image
% tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
BinnedADCmaptif = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedADCmap.tif', ii);
    BinnedADCmaptif(:,:,:,ii) = temp_tiff;
end
% BinnedADCmaptif = permute(BinnedADCmaptif,[1 2 4 3]);
% S = orthosliceViewer((test)); %colormap(SixBinMap);

%% Save Report in a powerpoint file
%Diffusion Images
Diffusion_Montage = figure('Name','Diffusion Images');set(Diffusion_Montage,'WindowState','minimized');
montage(reshape(Images(:,:,:,1),[size(Images,1), size(Images,2), 1, size(Images,3)]),...
    'Size',[1 size(Images,3)],'DisplayRange',[0 50]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
Diffusion_MontagePosition=get(gcf,'position');

%LBADC maps
LBADC_Montage = figure('Name','LBADC Images');set(LBADC_Montage,'WindowState','minimized');
montage(reshape(BinnedADCmap,[size(BinnedADCmap,1), size(BinnedADCmap,2),size(BinnedADCmap,3),...
    size(BinnedADCmap,4)]),'Size',[1 size(BinnedADCmap,4)],...
    'DisplayRange',[0 6]);set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap_Bins = [0 0 1  %Blue1
             0 0.57 0.71  %Blue2
             0.4 0.7 0.4  %Green1
             0 1 0  %Green2
             1 0.7143 0   %Orange  
             1 0 0   %Red 
             ];
colormap(colormap_Bins); colorbar('XTickLabel',{'0','LADC','LADC','NADC','NADC','HADC','HADC'},'XTick', 0:1:6);
LBADC_MontagePosition=get(gcf,'position');

nslices = size(BinnedADCmap,3);
if nslices > 8
    nslices = 8;
end
scale_fac = 2*nslices;

ReportTitle='Diffusion_Analysis';
cd(DataLocation);
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
% exportToPPTX('addslide'); %Title Slide
% exportToPPTX('addtext',sprintf('Subject'), ...
%     'Position',[0 0 16 4.5], ...
%     'VerticalAlignment','bottom', ...
%     'HorizontalAlignment','center', ...
%     'FontSize',72, ...
%     'FontWeight','bold');
% exportToPPTX('addtext',sprintf(['Scan Date: ','','\nProcessing Date: ',datestr(date,29)]), ...
%     'Position',[0 4.5 16 4.5], ...
%     'VerticalAlignment','top', ...
%     'HorizontalAlignment','center', ...
%     'FontSize',36);
%%%%%%%%%%%%%%%%%%%%%%%%%%% b0 image and LB ADC map %%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); 
Global.exportToPPTX('addpicture',Diffusion_Montage,...
    'Position',[0.5 0.5 scale_fac scale_fac*(Diffusion_MontagePosition(4)/Diffusion_MontagePosition(3))]);
Global.exportToPPTX('addtext',sprintf('ADC-Linear Binning Analysis'),'Position',[6 0 5 1],'Color','b','FontSize',25);
Global.exportToPPTX('addtext',sprintf('b0 Image'),'Position',[0.5 0.15 5 3],'Color','r','FontSize',20);

Global.exportToPPTX('addpicture',LBADC_Montage,...
    'Position',[0.5 2.5 scale_fac scale_fac*(LBADC_MontagePosition(4)/LBADC_MontagePosition(3))]);
Global.exportToPPTX('addtext',sprintf('ADC Color Map'),'Position',[0.5 2.2 5 3],'Color','r','FontSize',20);

%Histogram 
LBADC_hist=imread('LB_ADC_Histogram.png');    
Global.exportToPPTX('addpicture',LBADC_hist,'Position',[0 4.5 7.5 7.5*(size(LBADC_hist,1)/size(LBADC_hist,2))]);

%ADC Bins Summary
[~,pm,~,zm] = ztest(DiffMean,HealthyADCmean,HealthyADCstd);
Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(DiffMean,'%1.4f'),'±',num2str(DiffStd,'%1.4f'),' (Predicted Healthy = ',num2str(HealthyADCmean,'%1.4f'),'±',num2str(HealthyADCstd,'%1.4f'),'); Z = ',num2str(zm,'%1.4f'),', P = ',num2str(pm,'%1.4f')]),'Position',[6.5 5.5 9.5 2],'HorizontalAlignment','center');
BinTable   = { ...
    '',{'LDR1','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'LDR2','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Normal1','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Normal2','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'HDR1','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'HDR2','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Bins','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    'Percent (%)',num2str(DiffLow1Percent,'%1.3f'),num2str(DiffLow2Percent,'%1.3f'),num2str(DiffNormal1Percent,'%1.3f'),num2str(DiffNormal2Percent,'%1.3f'),num2str(DiffHigh1Percent,'%1.3f'),num2str(DiffHigh2Percent,'%1.3f');...
    'ADC (cm^2/s)',[num2str(Mean_ADC_Bin1,'%1.3f'),'±',num2str(Std_ADC_Bin1,'%1.3f')],...
                     [num2str(Mean_ADC_Bin2,'%1.3f'),'±',num2str(Std_ADC_Bin2,'%1.3f')],...
                     [num2str(Mean_ADC_Bin3,'%1.3f'),'±',num2str(Std_ADC_Bin3,'%1.3f')],...
                     [num2str(Mean_ADC_Bin4,'%1.3f'),'±',num2str(Std_ADC_Bin4,'%1.3f')],...
                     [num2str(Mean_ADC_Bin5,'%1.3f'),'±',num2str(Std_ADC_Bin5,'%1.3f')],...
                     [num2str(Mean_ADC_Bin6,'%1.3f'),'±',num2str(Std_ADC_Bin6,'%1.3f')],;...
    };
Global.exportToPPTX('addtable',BinTable,'Position',[7.5 6 7.75 2],'Vert','middle','Horiz','center','FontSize',12);

% save all
cd(DataLocation);
Global.exportToPPTX('save',fullfile(DataLocation, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n'); 
% close all;
BinnedADCmap =  permute(BinnedADCmap,[1 2 4 3]); 

% store result    
Diffusion.LB_colorDiffmap = BinnedADCmap;
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

% Convert to numeric
DL1 = str2double(Diffusion.DiffLow1Percent);
DL2 = str2double(Diffusion.DiffLow2Percent);

DN1 = str2double(Diffusion.DiffNormal1Percent);
DN2 = str2double(Diffusion.DiffNormal2Percent);

DH1 = str2double(Diffusion.DiffHigh1Percent);
DH2 = str2double(Diffusion.DiffHigh2Percent);

% Compute diffusion ranges
Diffusion.LDR = DL1 + DL2;       % Low range diffusion
Diffusion.NDR = DN1 + DN2;       % Normal diffusion
Diffusion.HDR = DH1 + DH2;       % High diffusion

%% save .m file
save_data=[DataLocation,'\','ADC_LinearBinningAnalysis','.mat'];
save(save_data); 
%% 

Processing_Time=toc/60

end

