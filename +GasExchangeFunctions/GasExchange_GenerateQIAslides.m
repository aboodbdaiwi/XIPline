function [] = GasExchange_GenerateQIAslides(DataLocation)

%
% workspace (string or character vector) -- is the full file path to the gas
%   exhange analysis workspace .m file
%

workspace = [DataLocation '\GasExchangeWorkspace.mat']; 

load(workspace, ...
'NumPlotSlices',...
'ProtonImageRegistered',...
'Slices_Co', ...
'ProtonMax',...
'VentBinMap', ...
'SixBinMap', ...
'Slices_Ax',...
'HealthyData', ...
'RBCTransferMean',...
'BarrierUptakeMean',...
'RBCBarrierMean', ...
'Subject', ...
'Notes', ...
'ScanDate', ...
'AppendedDissolvedNMRFit',...
'AppendedDissolvedFit',...
'DissolvedKSpace', ...
'XePulses', ...
'CorrectedDissKSpace_SS', ...
'GasKSpace',...
'HPhaseDynamics', ...
'motion_ind', ...
'HPulses',...
'HDynFit',...
'HDynRMSE',...
'SS_ind',...
'ProtonImage',...
'H_RecMatrix',...
'VentImage',...
'ProtonMaskRegistered', ...
'snr_kspace_end', ...
'snr_kspace_beginning', ...
'snr_spectrum_end'...
);

Registration = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registration,'WindowState','minimized');
set(Registration,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices
    nexttile
    imshowpair(abs(ProtonImage(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix+Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    if slice == round(NumPlotSlices/2)
        title('Before Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(fliplr(rot90(squeeze(abs(ProtonImage(H_RecMatrix+Slices_Ax(slice),H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(abs(ProtonImageRegistered(:,:,Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    if slice == round(NumPlotSlices/2)
        title('After Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
end

%% Signal Dynamics + Dissolved NMR Spectrum

%View H and Xe Dynamic and post dissolvedNMR spectrum data
DissolvedNMR_SigDynamics_fhandle = figure('Name','Signal Dynamics and Dissolved NMR');set(DissolvedNMR_SigDynamics_fhandle,'WindowState','minimized');
set(DissolvedNMR_SigDynamics_fhandle,'color','white','Units','inches','Position',[0.25 0.25 13.33 7.5])

%PostDissolvedNMR spectrum
subplot(1,3,1);
set(gca,'color','white','Units','normalized','Position',[0.036 0.1 0.35 0.8])%set(gca,'color','white','Units','inches','Position',[0.75 1 6 5])
hold on
%data
plot(AppendedDissolvedNMRFit.f, abs(AppendedDissolvedNMRFit.spectralDomainSignal),'b')%mag
plot(AppendedDissolvedNMRFit.f, real(AppendedDissolvedNMRFit.spectralDomainSignal),'r')%real
plot(AppendedDissolvedNMRFit.f, imag(AppendedDissolvedNMRFit.spectralDomainSignal),'Color',[0,0.6,0.2]')%imag
%fits
plot(AppendedDissolvedNMRFit.f,abs(sum(AppendedDissolvedFit,2)),'Color',[0 0 1 0.33],'LineWidth',3)%mag
plot(AppendedDissolvedNMRFit.f,real(sum(AppendedDissolvedFit,2)),'Color',[1 0 0 0.33],'LineWidth',3)%real
plot(AppendedDissolvedNMRFit.f,imag(sum(AppendedDissolvedFit,2)),'Color',[0,0.6,0.2,0.33]','LineWidth',3)%imag
%components
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,1)),'FaceColor','r','FaceAlpha',0.33,'LineStyle','none')%rbc
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,2)),'FaceColor','b','FaceAlpha',0.33,'LineStyle','none')%barrier
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,3)),'FaceColor',[0,0.6,0.2]','FaceAlpha',0.33,'LineStyle','none')%gas
%settings
lgnd = legend({'Spectrum - Magnitude','Spectrum - Real','Spectrum - Imaginary','Fit - Magnitude','Fit - Real','Fit - Imaginary','RBC - Real','Membrane - Real','Gas - Real'},'Location','best','NumColumns',3);
set(lgnd,'color','none', ...
    'FontSize',9);
xlim([-8000 2000])
set(gca, 'XDir','reverse')
title(['Dissolved Phase Spectra and Fit: RBC/Membrane = ',num2str(AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2))])
xlabel('Frequency (Hz)')
ylabel('NMR Signal (a.u.)')
%SNR text box
spec_SNR_str = sprintf('Spectrum SNR = %0.2f',snr_spectrum_end);
annotation('textbox',[0.17 0.35 0.3 0.3],'String',spec_SNR_str,'FitBoxToText','on','FontSize',12,'FontWeight','bold')
hold off

% Xe Dynamic data setup  
XePulses = (1:size(CorrectedDissKSpace_SS,2))'+SS_ind;%start after steady state


%Xe - SS and Detrending

subplot(1,3,2);
hold on
set(gca,'FontSize',11,'Units','normalized','Position',[0.44 0.15 0.25 0.70])%set(gca,'FontSize',12,'Units','inches','Position',[6.75 1 4 5])
%Dissolved
plot(abs(DissolvedKSpace(1,:)),'Color',[0.5 0 0],'LineWidth',1);
plot(XePulses,abs(CorrectedDissKSpace_SS(1,:)),'Color',[1.0 0 0],'LineWidth',1);
%Gas
plot(abs(GasKSpace(1,:)),'Color',[0 0 1],'LineWidth',1);
%SS
h1 = patch([0 SS_ind SS_ind 0], [max(ylim) max(ylim) 0 0], [0.66 0.66 0.66],'EdgeColor','none');%Remove due to steady state
xl = xlim; yl = ylim;
uistack(h1, 'bottom')
xlim(xl);ylim(yl);
%Add SNR text boxes around points x = XePulses(1) and y = CorrectedDissKSpace_SS(1,1)
x1 = XePulses(1);
y1 = abs(CorrectedDissKSpace_SS(1,1));
ax = gca;
ax_pos = ax.Position;
x_lim = ax.XLim;
y_lim = ax.YLim;
norm_x1 = (x1-x_lim(1))/(x_lim(2)-x_lim(1))*ax_pos(3)+ax_pos(1);
norm_y1 = (y1-y_lim(1))/(y_lim(2)-y_lim(1))*ax_pos(4)+ax_pos(2);
kspacebegin_SNR_str = sprintf('FID SNR = %0.2f',snr_kspace_beginning);
annotation('textarrow', [norm_x1 + 0.1, norm_x1], [norm_y1 + 0.1, norm_y1], 'String', kspacebegin_SNR_str);
%Add SNR text boxes around points x = XePulses(end) and y = CorrectedDissKSpace_SS(1,end)
x2 = XePulses(end);
y2 = abs(CorrectedDissKSpace_SS(1,end));
norm_x2 = (x2-x_lim(1))/(x_lim(2)-x_lim(1))*ax_pos(3)+ax_pos(1);
norm_y2 = (y2-y_lim(1))/(y_lim(2)-y_lim(1))*ax_pos(4)+ax_pos(2);
kspaceend_SNR_str = sprintf('FID SNR = %0.2f',snr_kspace_end);
annotation('textarrow', [norm_x2 - 0.04, norm_x2], [norm_y2 + 0.2, norm_y2], 'String', kspaceend_SNR_str);
%settings
title('Xe Signal Intensity Dynamics')
box on
xlabel('Projection')
ylabel('Signal Intensity')
legend({'Steady State Discard','Dissolved','Corrected Dissolved','Gas'},'Location','best','NumColumns',2,'FontSize',9)
hold off
%H - Motion
subplot(1,3,3);
hold on
set(gca,'FontSize',11,'Units','normalized','Position',[0.74 0.15 0.25 0.70])%set(gca,'FontSize',12,'Units','inches','Position',[11 1 4 5])
h1 = plot(HPhaseDynamics,'k','LineWidth',1);
xl = xlim; yl = ylim;
patch([motion_ind length(HPulses) length(HPulses) motion_ind], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.5 0.5 0.5],'EdgeColor','none')%Remove due to motion
uistack(h1, 'top');
plot(HDynFit,'Color',[0 1 0 0.5],'LineWidth',3);
plot(HDynFit+3*HDynRMSE,'Color',[0 1 0 0.5])
plot(HDynFit-3*HDynRMSE,'Color',[0 1 0 0.5])
xlim(xl);ylim(yl);
title('H Signal Phase Dynamics')
box on
xlabel('Projection')
ylabel('Signal Phase (deg)')
legend({'Unused - Motion','Data','Fit'},'Location','best','FontSize',10)
hold off


%% binning map figures

SumVentFig = figure('Name','Ventitlation Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumVentFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0 0 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','none');
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(VentBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Ventilation Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(VentBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end

set(SumVentFig,'WindowState','minimized');

%%
load(workspace,...
'RBCBinMap'...
)

SumRBCFig = figure('Name','RBC Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0 0 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','none');
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCFig,'WindowState','minimized');
%%
load(workspace,...
    'RBCBarrierBinMap', ...
    'SixBinRBCBarMap'...
    )

SumRBCBarFig = figure('Name','RBC:Membrane Ratio Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCBarFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0 0 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','none');
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBarrierBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
    colormap(gca,SixBinRBCBarMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC:Membrane Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
    colormap(gca,SixBinRBCBarMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCBarFig,'WindowState','minimized');

%%
load(workspace,...
    'BarrierBinMap',...
    'EightBinMap'...
    )

SumBarrFig = figure('Name','Membrane Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumBarrFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0 0 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','none');
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(BarrierBinMap(:,:,Slices_Co(slice),1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Membrane Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(BarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumBarrFig,'WindowState','minimized');
%%

isOpen  = GasExchange_exportToPPTX.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    GasExchange_exportToPPTX.exportToPPTX('close');
end

% Create a new presentation
GasExchange_exportToPPTX.exportToPPTX('new','Dimensions',[16 9]);

% Add a slide for the spectra/dynamics
rose = [1 0.773 0.796];
GasExchange_exportToPPTX.exportToPPTX('addslide','BackgroundColor',rose);

% Add figures 
GasExchange_exportToPPTX.exportToPPTX('addpicture',DissolvedNMR_SigDynamics_fhandle,'Scale','maxfixed');

% Add title for the signal dynamics slide
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Gas Exchange | %s%s, age,sex,disease | scan date: %s | page 1/2',Notes,Subject,ScanDate),'Position',[2 0 13 3],'FontSize',22);




% Add a slide for the montages
GasExchange_exportToPPTX.exportToPPTX('addslide','BackgroundColor',rose);

% Title for binning map montage slide
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Gas Exchange | %s%s, age,sex,disease | scan date: %s | page 2/2',Notes,Subject,ScanDate),'Position',[2 0 13 3],'FontSize',22);

% add montages
sfact = 1.7;
shift_left = 0.25;
shift_up = 0.5;
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumVentFig,'Position',[-shift_left shift_up 14/sfact 4/sfact]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumRBCFig,'Position',[14/sfact-shift_left shift_up 14/sfact 4/sfact]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumBarrFig,'Position',[-shift_left shift_up+4/sfact 14/sfact 4/sfact]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumRBCBarFig,'Position',[14/sfact-shift_left shift_up+4/sfact 14/sfact 4/sfact]);

% add table

meanSummary = {...
    '',{'Membrane','FontWeight','bold'},{'RBC','FontWeight','bold'},{'RBC:Membrane','FontWeight','bold'};...
    {'Mean','FontWeight','bold'},num2str(BarrierUptakeMean,'%1.5f'),num2str(RBCTransferMean,'%1.3e'),num2str(RBCBarrierMean,'%1.3f');...
    {'Healthy Mean','FontWeight','bold'},num2str(HealthyData.HealthyMeans.MeanBarrier,'%1.5f'),num2str(HealthyData.HealthyMeans.MeanRBC,'%1.3e'),num2str(HealthyData.HealthyMeans.MeanRBCBar,'%1.3f');...
    };
GasExchange_exportToPPTX.exportToPPTX('addtable',meanSummary,'Position',[0.1 7.25 5.5 1.5],'Vert','middle','Horiz','center','FontSize',16)

% add registration
sf = 2;
GasExchange_exportToPPTX.exportToPPTX('addpicture',Registration,'Position',[6 5 Registration.Position(3)/sf Registration.Position(4)/sf]);

close all;
% Save the presentation
GasExchange_exportToPPTX.exportToPPTX('save',[DataLocation '\' 'ShortGasExchangeSummary_' Notes Subject '_' ScanDate '.pptx']);


