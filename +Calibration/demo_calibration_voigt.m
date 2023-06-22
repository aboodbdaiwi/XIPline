% 5/10/17 - updated to report key values to scanner tech more clearly
% 5/11/17 - Add red warnings for frequency, TE90, SNR, and Vref
% 6/29/17 - Added the option to fit barrier to voigt lineshape instead
% 6/30/17 - fix bug in TE90 calculation associated with negative deltaF
% 7/27/17 - fix additional TE90 bug if inital delta phase is negative
% 10/10/18 - update to use combined dynamic/calibration data
clc; clear all; close all;


%% Define number of calibration, number to skip, number of gas, target flip
FlipTarget = 20; % target flip angle from calibration
nSkip = 100; % number of dissolved fids to skip to ensure steady-state
nAvg = 50; % Dissolved fids to average. Set to 50 to match Elly's 1-sec static. Set to >600 to use all available
nCal = 20; % number of flip angle calibration frames following, normally 20
nGas = 1; % set to 1 for original 3T calibration, set to 0 for new dynamic/calibration combination scan
Voigt = 1; % Barrier lorentzian for 0, Voigt if 1. Requires Elly's NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1
disp_fids = 0; % set to 1 to display first and last fids for noise checks.

%% define tolerances for wanrings
freq_std = 34091550; % roughly expected frequency
freq_tol = 200; % frequency adjustment tolerance
deltaPhase1_tol = 90; % if calibration phase>90 use minTE
TrueRefScale_tol = 1.17; % corresponds to 700V ref voltage
SNR_tol = 25; % allowable
RbcBarAdjust = 80; % adjust RBC:barrier to scale from 73ms to 15 ms in imaging - legacy

%% first read in the file and create the twix object
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open
twix_obj = mapVBVD(file_with_path); 
filename = file(15:end-4); % return filename back to something short
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems

% extract favorite variables from the header
Seq_name = twix_obj.hdr.Config.SequenceDescription;
weight = twix_obj.hdr.Dicom.flUsedPatientWeight;  %
ref_amp = twix_obj.hdr.Spice.TransmitterReferenceAmplitude; % Reference amplitude entered
te = twix_obj.hdr.Meas.alTE(1); 
tr = twix_obj.hdr.Config.TR;
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
nFids = twix_obj.hdr.Config.NRepMeas;
nPts = twix_obj.hdr.Config.VectorSize;  % present in cal and FID but not imaging sequences.
rf_amp1 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude;  % seems to be in all sequences
rf_amp2 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,2}.flAmplitude;  % in cali this is the 2nd flip angle (gas-phase usually)
rf_amp3 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,3}.flAmplitude;  % in cali this is the 3rd flip angle (calibration, usually)
%rf_amp3=1; %dummy override for when using hard pulse

% print out key variable values
fprintf('\n\n');
fprintf('\tName = %s\n',Seq_name');
fprintf('\tWeight = %0.1f kg \n',weight);
fprintf('\tTE=%0.0f us\n',te);
fprintf('\tTR=%0.0f us\n',tr);
fprintf('\tFrequency=%0.0f Hz\n',freq);
fprintf('\tvoltage=%0.1f V\n',rf_amp1);
fprintf('\tdwell_time=%0.0f ns\n',dwell_time);
fprintf('\tnPts=%0.0f\n',nPts);
fprintf('\tnFrames=%0.0f\n\n',nFids);

dwell_time=dwell_time*1E-9; % convert to seconds for calculations

% parse out the various fids - dissolved, gas, and calibration
theFID = squeeze(double(twix_obj.image()));
nDis=nFids - nCal - nGas; 
nDis1 = nDis - nSkip; % if skipping first dissolved FIDs
disData = theFID(:,1:nDis); % all dissolved data
gasData = theFID(:,(nDis+1));
if nAvg >nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = theFID(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = theFID(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end %if
disData_avg = mean(disData,2);  % average of all dissolved
disData1_avg = mean(disData1,2);  % average of dissolved after skipping
calData = theFID(:,(nDis+2):(nDis+nCal+1)); % the data left for flip angle calculations.

%% Fit gas Spectrum
%   t = double((0:nPts-1)*dwell_time');
if nGas >0
    fprintf('Analysis of Gas FID\n');
    t = double((0:(length(gasData)-1))*dwell_time');
    gasfitObj = NMR_TimeFit(gasData,t,1e-4,-84,30,0,0,10000);
    gasfitObj.fitTimeDomainSignal();
    figure();
    gasfitObj.plotTimeAndSpectralFit;
    gasfitObj.describe();
end %if

%% Fit dissolved Spectrum
fprintf('\nAnalysis of Averaged Dissolved FID, skipping %0.0f\n',nSkip);
%disfitObj = NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[200 200 40],[0 0 0],0,10000);
if Voigt == 1
    fprintf('Using Voigt Lineshape to fit Barrier\n');
    disfitObj = NMR_TimeFit_v(disData1_avg,t,[1 1 1],[0 -700  -7400],[250 200 30],[0 200 0],[0 0 0],0,length(t)); % first widths lorenzian, 2nd are gauss
else
    fprintf('Using Lorentzian Lineshape to fit Barrier\n');
    disfitObj = NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[240 240 40],[0 0 0],0,length(t));
end % if
disfitObj = disfitObj.fitTimeDomainSignal();
figure()
disfitObj.plotTimeAndSpectralFit;
disfitObj.describe();


%% calculate target transmit receive frequency from gas fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n',freq_target); % report new target frequency
if abs(freq_target-freq_std)>freq_tol
    fprintf(2,'WARNING! Frequency adjust exceeds tolerances; Check system\n');
end %if

%% Calculate and report TE90 NEW
deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-barrier freq difference
deltaTe90= (90-deltaPhase)/(360*deltaF); % how far off are we?
te90 =te + deltaTe90*1e6; % in usec
fprintf('TE90=%3.2f ms\n',te90/1000); % report TE90 in ms
if abs(deltaPhase)> 90
    % DGM 3/3/20: Adding specific minimum TE90 value of 0.45ms to warning message
    fprintf(2,'WARNING! Phi_cal =%3.0f%c!; Use min TE of 0.45ms!\n',deltaPhase,char(176));
end %if

%% Calculate and report TE90 OLD
% minte = 0;
% deltaF = disfitObj.freq(2)-disfitObj.freq(1); % RBC-barrier freq difference
% deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
% 
% time180 = abs(1/(2*deltaF));
% deltaTe90= (90-deltaPhase)/(360*deltaF);
% while(deltaTe90>(minte+time180))
%     % This te is too high, so subtract 180 deg of phase
%     deltaTe90 = deltaTe90 - time180;
% end
% while(deltaTe90<minte)
%     % This TE is too low, so add 180 deg of phase
%     deltaTe90 = deltaTe90 + time180;
% end
% te90 =te + deltaTe90*1e6; % in usec
% fprintf('TE90=%3.2f ms\n',te90/1000); % report TE90 in ms
% % calculate initial phase difference in calibration as double check
% deltaPhase1 = deltaPhase;
% if abs(deltaPhase)>180
%     deltaPhase1=360-deltaPhase;
% end %if
% if abs(deltaPhase1)>deltaPhase1_tol
%     fprintf(2,'WARNING! Phi_cal =%3.0f%c!; Use min TE!\n',abs(deltaPhase1),char(176));
% end %if

%% Find amplitudes (slow! full fitting only for super noisy data)
% for iSpec = 1:nCal    
%     nmrFit = NMR_TimeFit(calData(:,iSpec), 1:1024, 1,1,1,0, 0, 1024);
%     nmrFit = nmrFit.fitTimeDomainSignal();
%     flipCalAmps(iSpec) = nmrFit.area;
% end

% Find Amplitudes faster by using max amplitudes from each FID in calibration
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1);   % cos theta decay
guess(1)=max(flipCalAmps);
guess(2)=20*pi/180;       % just guess 10 degrees

xdata=1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);
ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
param_err=fitparams-ci(:,1)';
flip_angle=abs(fitparams(2)*180/pi);
flip_err=param_err(2)*180/pi;
%fprintf('flip angle=%3.1f \n',flip_angle); % report TE90 in ms

figure();  % plot up the calibration data
plot(xdata,ydata,'-b');
hold on;
plot(xdata,fitfunct(fitparams,xdata),'-r');
legend('Acquired','Fit');
xlabel('Frame Number');
ylabel('Magnitude');
a=sprintf('Flip Cal: V_p=%0.1f yields flip=%0.1f±%0.1f%c', rf_amp3,flip_angle,flip_err,char(176));
title(a); % note degree symbol is char(176)

%% provide reference amplitude and warnings
VRefScaleFactor=FlipTarget/flip_angle; % How much should Ref Voltage Be Scaled?
fprintf('TrueRef scale factor = %3.3f \n',VRefScaleFactor); %
fprintf('For 600V calibration, True_Ref = %3.0f V\n',600*VRefScaleFactor); % report TE90 in ms
if VRefScaleFactor>TrueRefScale_tol
    fprintf(2,'WARNING! Excessive calibration scale factor; Check system\n');
end %if

%% Calculate and report dissolved phase spectral SNR
timeFit=disfitObj.calcTimeDomainSignal(disfitObj.t);  % calculate fitted data in time domain
timeRes=timeFit - disfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_dis=disfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('\nSNR Check for dissolved and gas peaks...\n'); % 
fprintf('SNR for dissolved peaks =%3.1f \n',SNR_dis); % 

%% Calculate and report final gas phase FID spectral SNR
timeFit=gasfitObj.calcTimeDomainSignal(gasfitObj.t);  % calculate fitted data in time domain
timeRes=timeFit - gasfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_gas=gasfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('SNR for final gas peak =%3.1f \n',SNR_gas); %
if SNR_gas<SNR_tol
    fprintf(2,'WARNING! Gas FID SNR below minimums; Check Coil Plug\n');
end %if

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3)/sum(disfitObj.area(1:2));
fprintf('\nGasDisRatio=%3.3f \n',GasDisRatio); 

%% Quantify RBC:Barrier ratio
RbcBarRatio = disfitObj.area(1)/disfitObj.area(2);
fprintf('RbcBarRatio=%3.3f\n',RbcBarRatio); % 
%fprintf('RbcBar_%2.0f = %3.3f\n',RbcBarAdjust,RbcBarRatio*RbcBarAdjust/100)

%% display representative fids if requested
if disp_fids == 1 % display representative fids only if desired
    figure();  % plot up representative first FID for time-domain noise check
    plot(real(theFID(:,1)),'-b');
    hold on;
    plot(imag(theFID(:,1)),'-r');
    legend('Real','Imag');
    xlabel('Points');
    ylabel('Real Signal');
    a=sprintf('First Dissolved Frame');
    title(a);
    
    figure();  % plot up representative last FID for time-domain noise check
    plot(real(theFID(:,nDis)),'-b');
    hold on;
    plot(imag(theFID(:,nDis)),'-r');
    legend('Real','Imag');
    xlabel('Points');
    ylabel('Real Signal');
    a=sprintf('Last Dissolved Frame');
    title(a);
end % if

%% create old fashioned matrix of fit parameters of old
if Voigt == 1
    fitparams_array=[disfitObj.area' (disfitObj.freq'-round(gasfitObj.freq)) disfitObj.fwhmL' disfitObj.fwhmG' disfitObj.phase'];  % modified for Voigt
else
    fitparams_array=[disfitObj.area' (disfitObj.freq'-round(gasfitObj.freq)) disfitObj.fwhm' disfitObj.phase'];  % standard
end %if