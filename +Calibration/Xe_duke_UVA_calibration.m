function [GasExResults, CalResults] = Xe_duke_UVA_calibration(MainInput)

% 5/10/17 - updated to report key values to scanner tech more clearly
% 5/11/17 - Add red warnings for frequency, TE90, SNR, and Vref
% 6/29/17 - Added the option to fit barrier to voigt lineshape instead
% 6/30/17 - fix bug in TE90 calculation associated with negative deltaF
% 7/27/17 - fix additional TE90 bug if inital delta phase is negative
% 10/10/18 - update to use combined dynamic/calibration data
% 01/08/21 - updated to accept either Duke or UVA twix files
% 01/14/21 - updated to accept Duke, Duke/UVA, or UVA twix files
% clc; clear all; close all;


%% Define number of calibration, number to skip, number of gas, target flip
FlipTarget = 20; % target flip angle from calibration
nDis = 500;
nSkip = 100; % number of dissolved fids to skip to ensure steady-state
nAvg = 50; % Dissolved fids to average. Set to 50 to match Elly's 1-sec static. Set to >600 to use all available
nCal = 20; % number of flip angle calibration frames following, normally 20Voigt = 1; % Barrier lorentzian for 0, Voigt if 1. Requires Elly's NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1
Voigt = 1;
disp_fids = 0; % set to 1 to display first and last fids for noise checks.

%% define tolerances for wanrings
freq_tol = 200; % frequency adjustment tolerance
deltaPhase1_tol = 90; % if calibration phase>90 use minTE
TrueRefScale_tol = 1.33; 
SNR_tol = 25; % allowable

%% first read in the file and create the twix object
% if exist('DataFile','var') == 0
%   [DataFile, DataPath] = uigetfile('*.dat', 'Select a .dat file');  % starts in current dir 
%   DataFile = fullfile(DataPath, DataFile);
% end
% DataFiles = dir([DataPath,'\*.dat']);
file_name = MainInput.XeFileName;
file_folder = MainInput.XeDataLocation;

% file_with_path = strcat(file_folder,'\',file_name);  % join path and filename to open

twix_obj = Calibration.mapVBVD(MainInput.XeFullPath); 
filename = file_name(15:end-4); % return filename back to something short
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems
file_loc = regexp(MainInput.XeDataLocation,filesep,'split');
file_loc = file_loc{end-1};

%%
% extract favorite variables from the header
Seq_name = twix_obj.hdr.Config.SequenceDescription;
weight = twix_obj.hdr.Dicom.flUsedPatientWeight;  %
ref_amp = twix_obj.hdr.Spice.TransmitterReferenceAmplitude; % Reference amplitude entered
te = twix_obj.hdr.Phoenix.alTE{1}; 
tr = twix_obj.hdr.Config.TR;
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
nFids = twix_obj.hdr.Config.NRepMeas;
nPts = twix_obj.hdr.Config.VectorSize;  % present in cal and FID but not imaging sequences.

% get scan date
scanDate = twix_obj.hdr.Phoenix.tReferenceImage0; 
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];

if isfield(twix_obj.hdr.Phoenix, 'sWiPMemBlock')
    % Duke twix file
    if isfield(twix_obj.hdr.Phoenix.sWiPMemBlock,'adFree')
        VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{4};
        rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;  % seems to be in all sequences
        rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;  % in cali this is the 3rd flip angle (calibration, usually)
    %rf_amp3=1; %dummy override for when using hard pulse
    elseif isfield(twix_obj.hdr.Phoenix.sWiPMemBlock,'alFree')
        % Duke twix file using UVA sequence
        VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1};
        rf_amp1 = [];
        rf_amp2 = [];
    else
        disp('WARNING: twix file type not supported, cannot determine reference voltage')
    end 
elseif isfield(twix_obj.hdr.Phoenix, 'sWipMemBlock')
    % UVA twix file
    VRef = twix_obj.hdr.Phoenix.sWipMemBlock.alFree{1};
    rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{2}.flAmplitude;  % dissolved phase
    rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;  % calibration
else
    disp('WARNING: twix file type not supported, cannot determine reference voltage')
end 

% print out key variable values
fprintf('\n\n');
fprintf('\tFile Location = %s\n',file_loc);
fprintf('\tName = %s\n',Seq_name');
fprintf('\tScan Date = %s\n',scanDateStr');
fprintf('\tWeight = %0.1f kg \n',weight);
fprintf('\tTE = %0.0f us\n',te);
fprintf('\tTR = %0.0f us\n',tr);
fprintf('\tFrequency = %0.0f Hz\n',freq);
fprintf('\tReference Voltage = %0.1f V\n',VRef);
fprintf('\tPulse Voltage (Dissolved) = %0.1f V\n',rf_amp1');
fprintf('\tPulse Voltage (Gas) = %0.1f V\n',rf_amp2);
fprintf('\tDwell Time = %0.0f ns\n',dwell_time);
fprintf('\tnPts = %0.0f\n',nPts);
fprintf('\tnFrames = %0.0f\n\n',nFids);

dwell_time=dwell_time*1E-9; % convert to seconds for calculations

% parse out the various fids - dissolved, gas, and calibration
theFID = squeeze(double(twix_obj.image()));
nDis1 = nDis - nSkip; % if skipping first dissolved FIDs
disData = theFID(:,1:nDis); % all dissolved data
gasData = theFID(:,end-nCal+1);
if nAvg > nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = theFID(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = theFID(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end %if
disData_avg = mean(disData,2);  % average of all dissolved
disData1_avg = mean(disData1,2);  % average of dissolved after skipping
calData = theFID(:,end-nCal+1:end); % the data left for flip angle calculations.
t = double((0:(length(disData)-1))*dwell_time');
noise_est = max(abs(disData(end-10:end,end)));
%% Fit gas Spectrum

fprintf('Analysis of Gas FID\n');
gasfitObj = Calibration.NMR_TimeFit(gasData,t,1e-4,-84,30,0,0,10000);
gasfitObj.fitTimeDomainSignal();
%PJN - Give Figure a name and a position
% figure('Name','Gas Phase Analysis');
% ax = gasfitObj.plotTimeAndSpectralFit;
gasfitObj.describe();
% set(gcf,'Units','Normalized','Position',[0 0.04 0.5 0.85]);

% figure; hold on
% plot(gasfitObj.fMat,real(gasfitObj.individualSpectrums(:,1)),'r','LineWidth',2)
% plot(gasfitObj.fMat,real(gasfitObj.individualSpectrums(:,2)),'b','LineWidth',2)
% plot(gasfitObj.fMat,real(gasfitObj.individualSpectrums(:,3)),'Color',[0 0.5 0],'LineWidth',2)

%% Fit dissolved Spectrum
fprintf('\nAnalysis of Averaged Dissolved FID, skipping %0.0f\n',nSkip);
%disfitObj = NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[200 200 40],[0 0 0],0,10000);
if Voigt == 1
    fprintf('Using Voigt Lineshape to fit Barrier\n');
    disfitObj = Calibration.NMR_TimeFit_v(disData1_avg,t,[1 1 1],[0 -700  -7400],[250 200 30],[0 200 0],[0 0 0],0,length(t)); % first widths lorenzian, 2nd are gauss
else
    fprintf('Using Lorentzian Lineshape to fit Barrier\n');
    disfitObj = Calibration.NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[240 240 40],[0 0 0],0,length(t));
end % if

disfitObj = disfitObj.fitTimeDomainSignal();
%PJN - Give a name and position to this plot
% figure('Name','Dissolved Phase Analysis')
% disfitObj.plotTimeAndSpectralFit;
% set(gcf,'Units','Normalized','Position',[0.5 0.04 0.5 0.85]); 

disp('            Area     Freq (Hz)   Linewidths(Hz)   Phase(degrees)');
peakName = {'    RBC:','Barrier:','    Gas:'};
for iComp = 1:length(disfitObj.area)
    disp([peakName{iComp} '  '...
        sprintf('%8.3e',disfitObj.area(iComp)) ' ' ...
        sprintf('%+8.2f',disfitObj.freq(iComp))  '  ' ...
        sprintf('%8.2f',abs(disfitObj.fwhm(iComp))) '  ' ...
        sprintf('%8.2f',abs(disfitObj.fwhmG(iComp))) '  ' ...
        sprintf('%+9.2f',disfitObj.phase(iComp))]);
end


%% calculate target transmit receive frequency from gas fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n',freq_target); % report new target frequency
if abs(freq_target-freq)>freq_tol
    fprintf(2,'WARNING! Frequency adjust exceeds tolerances; Check system\n');
end 

%% Calculate and report TE90 NEW
deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-barrier freq difference
deltaTe90= (90-deltaPhase)/(360*deltaF); % how far off are we?
te90 =te + deltaTe90*1e6; % in usec
fprintf('TE90 = %3.2f ms\n',te90/1000); % report TE90 in ms
if abs(deltaPhase)> 90
    fprintf(2,'WARNING! Phi_cal = %3.0f%c!; Use min TE!\n',deltaPhase,char(176));
end %if

%%
% Find Amplitudes faster by using max amplitudes from each FID in calibration
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1) + noise_est;   % cos theta decay
% fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1)+coefs(3);   % cos theta decay
guess(1)=max(flipCalAmps);
guess(2)=20*pi/180;       % just guess 10 degrees
% guess(3)= 0;      

xdata=1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);
ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
param_err=fitparams-ci(:,1)';
flip_angle=abs(fitparams(2)*180/pi);
flip_err=param_err(2)*180/pi;
FlipScaleFactor = flip_angle/FlipTarget;

%PJN - give figure a name
% figure('Name','Flip Angle Calibration');  % plot up the calibration data
% %PJN - 
% plot(xdata,ydata,'bo','MarkerFaceColor','b');
% hold on;
% plot(xdata,fitfunct(fitparams,xdata),'-r');
% legend('Acquired','Fit');
% xlabel('Frame Number');
% ylabel('Magnitude');
% %PJN - Show Ref Voltage instead of Pulse Voltage
% %a=sprintf('Flip Cal: V_p = %0.1f yields flip = %0.1f±%0.1f%c', rf_amp2,flip_angle,flip_err,char(176));
% a=sprintf('Flip Cal: V_p = %0.1f yields flip = %0.1f±%0.1f%c', VRef,flip_angle,flip_err,char(176));
% title(a); % note degree symbol is char(176)

%% provide reference amplitude and warnings
VRefScaleFactor=FlipTarget/flip_angle; % How much should Ref Voltage Be Scaled?
fprintf(['For ', num2str(VRef), 'V calibration, True_Ref = %3.0f V\n'],VRef*VRefScaleFactor); % report TE90 in ms
fprintf(['(Estimated Weight Reference Voltage = %3.0f V)\n'],558 + 1.4*weight); % report TE90 in ms
if VRefScaleFactor>TrueRefScale_tol
    fprintf(2,'WARNING! Excessive voltage calibration scale factor; Check system\n');
end %if

%% Calculate and report dissolved phase spectral SNR
timeFit=disfitObj.calcTimeDomainSignal(disfitObj.t);  % calculate fitted data in time domain
timeRes=timeFit - disfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_dis=disfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('\nSNR Check for dissolved and gas peaks...\n'); % 
fprintf('SNR for RBC peak = %3.1f \n',SNR_dis(1)); 
fprintf('SNR for Barrier peak = %3.1f \n',SNR_dis(2)); 
fprintf('SNR for Gas peak = %3.1f \n',SNR_dis(3)); 

%% Calculate and report final gas phase FID spectral SNR
timeFit=gasfitObj.calcTimeDomainSignal(gasfitObj.t);  % calculate fitted data in time domain
timeRes=timeFit - gasfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_gas=gasfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('SNR for final gas peak = %3.1f \n',SNR_gas); %
if SNR_gas<SNR_tol
    fprintf(2,'WARNING! Gas FID SNR below minimums; Check Coil Plug\n');
end %if

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3)/sum(disfitObj.area(1:2));
fprintf('\nGasDisRatio = %3.3f \n',GasDisRatio); 

%% Quantify RBC:Barrier ratio
RbcBarRatio = disfitObj.area(1)/disfitObj.area(2);
fprintf('RbcBarRatio = %3.3f\n',RbcBarRatio); % 

%% display representative fids if requested
% if disp_fids == 1 % display representative fids only if desired
%     figure();  % plot up representative first FID for time-domain noise check
%     plot(real(theFID(:,1)),'-b');
%     hold on;
%     plot(imag(theFID(:,1)),'-r');
%     legend('Real','Imag');
%     xlabel('Points');
%     ylabel('Real Signal');
%     a=sprintf('First Dissolved Frame');
%     title(a);
%     
%     figure();  % plot up representative last FID for time-domain noise check
%     plot(real(theFID(:,nDis)),'-b');
%     hold on;
%     plot(imag(theFID(:,nDis)),'-r');
%     legend('Real','Imag');
%     xlabel('Points');
%     ylabel('Real Signal');
%     a=sprintf('Last Dissolved Frame');
%     title(a);
% end % if

%% create old fashioned matrix of fit parameters of old
% if Voigt == 1
%     fitparams_array=[disfitObj.area' (disfitObj.freq'-round(gasfitObj.freq)) disfitObj.fwhm' disfitObj.fwhmG' disfitObj.phase'];  % modified for Voigt
% else
%     fitparams_array=[disfitObj.area' (disfitObj.freq'-round(gasfitObj.freq)) disfitObj.fwhm' disfitObj.phase'];  % standard
% end %if

%% PJN - Creat Pop-up box of parameters that need to be set in real time
% CreateStruct.Interpreter = 'tex';
% CreateStruct.WindowStyle = 'non-modal';
% 
% textstr = {['\bf Gas Center Frequency: \rm' num2str(freq_target,8) ' Hz'];' '; [' \bf TE_{90}: \rm' num2str(te90/1000,2) ' ms']; ' '; ['\bf Reference Voltage: \rm' num2str(VRef*VRefScaleFactor,3) 'V']};
% outbox = msgbox(textstr,'Imaging Parameters',CreateStruct);
% 
% th = findall(outbox, 'Type', 'Text');                   %get handle to text within msgbox
% th.FontSize = 18;                                   %Change the font size
% 
% deltaWidth = sum(th.Extent([1,3]))-outbox.Position(3) + th.Extent(1);
% deltaHeight = sum(th.Extent([2,4]))-outbox.Position(4) + 10;
% outbox.Position([3,4]) = outbox.Position([3,4]) + [deltaWidth, deltaHeight];
% 
% outbox.Resize = 'on';

%% Store Results for Additional Analysis - Abood 
%Results to Pass to Gas Exchange Recon
GasExResults.RbcBarRatio = RbcBarRatio;
GasExResults.GasDisRatio = GasDisRatio;
GasExResults.DisFit = disfitObj;

%Results to Pass to Calibration information
CalResults.flip_angle = flip_angle;
CalResults.flip_err = flip_err;
CalResults.FlipScaleFactor = FlipScaleFactor;
CalResults.Pulses = xdata;
CalResults.GasDecay = ydata;
CalResults.DecayFit = fitparams;
CalResults.GasFit = gasfitObj;
CalResults.te90 = te90/1000; %to ms
CalResults.freq_target = freq_target;
CalResults.Reference_Voltage = VRef*VRefScaleFactor;
CalResults.dwell_time = dwell_time;  
CalResults.noise_est = noise_est;  
CalResults.freq = freq; 
CalResults.FlipTarget = FlipTarget;
CalResults.te = te;  
CalResults.nDis = nDis;
CalResults.nCal = nCal;
CalResults.VRefScaleFactor = VRefScaleFactor;
CalResults.VRef = VRef;

end




