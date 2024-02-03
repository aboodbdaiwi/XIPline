function [GasExResults, CalResults] = XeCTC_Calibration_MRD(MainInput)
%XeCTC_ibration - Analysis of Xe CTC calibration scan
%   Will process the resulting XeCTC Calibration data from the CPIR complete
%    patch 
%
%   Syntax:  results = XeCTC_Calibration(DataFile)
%
%   Inputs:
%      DataFile - (optional) path to data file to process
%
%   Outputs:
%      GasExResults - structure containing the results from the analysis
%           pertaining to GasEx
%      CalResults - structure containing the results for the calibrations
%
%   Package: https://https://github.com/cchmc-cpir/XeCTC_Calibration
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/

%   Please add updates at the end. Ex: 3/10/24 - ASB: updated flip angle caluclation for spiral diffusion
% clc; 
% clearvars -except DataFile;
% close all;

%% Fixed Parameters
%CCHMC Parameters 
freq = 2463;
dwell_time = 4.88*1E-6;% micro-s to s
te = 450; %us

%Basic Parameters
FlipTarget = 20; % target flip angle from calibration
nDis = 500;
nSkip = 100; % number of dissolved fids to skip to ensure steady-state
nAvg = 50; % Dissolved fids to average. Set to 50 to match Elly's 1-sec static. Set to >600 to use all available
nCal = 20; % number of flip angle calibration frames following, normally 20
Voigt = 1; % Barrier lorentzian for 0, Voigt if 1. Requires Elly's NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1

%% Define Tolerances for Warnings
freq_tol = 200; % frequency adjustment tolerance
deltaPhase1_tol = 90; % if calibration phase>90 use minTE
FAScale_tol = 1.33; 
SNR_tol = 25; % allowable

%% Import Data File
cd(MainInput.XeDataLocation)
if exist(MainInput.XeFileName, 'file')
    dsetCal = LoadData.ismrmrd.Dataset(MainInput.XeFileName, 'dataset');
else
    error(['File ' MainInput.XeFileName ' does not exist.  Please generate it.'])
end
hdrCal = LoadData.ismrmrd.xml.deserialize(dsetCal.readxml);
DCal = dsetCal.readAcquisition();

cellSize = cell2mat(DCal.data(1,1));
theFID = zeros(size(cellSize,1),size(DCal.data,2));
for i = 1:size(DCal.data,2)
    theFID(:,i) = cell2mat(DCal.data(1,i));
end

%% Parse Out the Various FIDs
if size(theFID,2) < 500
    nDis = size(theFID,2);
end
disData = theFID(:,1:nDis); % all dissolved data
gasData = theFID(:,end-nCal+1);
if nAvg > nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = theFID(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = theFID(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end
disData1_avg = mean(disData1,2);  % average of dissolved after skipping
calData = theFID(:,end-nCal+1:end); % the data left for flip angle calculations.
t = double((0:(length(disData)-1))*dwell_time');
noise_est = max(abs(disData(end-10:end,end)));
% figure; plot(abs(disData(:,end)))
% figure; plot(abs(gasData(:,end)))
% figure; imslice(abs(theFID))
%% calculate frequency 

% gasdata_fa = theFID(:,end-nCal+1:end);
% gasdata_fa_avg = mean(gasdata_fa,2);
% % figure; plot(abs(gasdata_fa_avg));
% gasdata_fa_avg_fft = fftshift(ifft((gasdata_fa_avg)));
% ts = 4.88*1E-6;% dwell time  - micro-s to s
% fs = 1/(ts); % ferquency rate
% N = 2048; % number of samples 
% f=(-N/2:N/2-1)/N*fs;
% figure; plot(f,abs(gasdata_fa_avg_fft))
% % figure; plot(gasfitObj.f,abs(gasdata_fa_avg_fft))
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% startPoints = [0 0 10000 0.1];
% x = f;
% y = abs(gasdata_fa_avg_fft)';
% crop = 900:1100;
% f1 = fit(x(1,crop)',y(1,crop)',gaussEqn,'Start', startPoints);
% plot(f1,x,y)

%% Fit gas Spectrum
fprintf('Analysis of Gas FID\n');
gasfitObj = Calibration.NMR_TimeFit(gasData,t,1e-4,-84,30,0,0,10000);
gasfitObj.fitTimeDomainSignal();
% figure('Name','Gas Phase Analysis');
% gasfitObj.plotTimeAndSpectralFit;
gasfitObj.describe();
% set(gcf,'Units','Normalized','Position',[0 0.04 0.5 0.85]);

%% Fit Dissolved Spectrum
fprintf('\nAnalysis of Averaged Dissolved FID, skipping %0.0f\n',nSkip);
if Voigt == 1
    fprintf('Using Voigt Lineshape to fit Barrier\n');
    disfitObj = Calibration.NMR_TimeFit_v(disData1_avg,t,[1 1 1],[0 -720  -7700],[250 200 30],[0 200 0],[0 0 0],0,length(t)); % first widths lorenzian, 2nd are gauss
else
    fprintf('Using Lorentzian Lineshape to fit Barrier\n');
    disfitObj = Calibration.NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -720  -7700],[240 240 40],[0 0 0],0,length(t));
end
disfitObj = disfitObj.fitTimeDomainSignal();
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

%% Calculate Target TR Frequency from Gas Fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n',freq_target); % report new target frequency
if abs(freq_target-freq)>freq_tol
    fprintf(2,'WARNING! Frequency adjust exceeds tolerances; Check system\n');
end 

%% Calculate and Report TE90
deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-barrier freq difference
deltaTe90 = (90-deltaPhase)/(360*deltaF); % how far off are we?
te90 = te + deltaTe90*1e6; % in usec
fprintf('TE90 = %3.2f ms\n',te90/1000); % report TE90 in ms
if abs(deltaPhase)> deltaPhase1_tol
    fprintf(2,'WARNING! Phi_cal = %3.0f%c!; Use min TE!\n',deltaPhase,char(176));
end

%% Calculate True Flip Angle
flipCalAmps = max(abs(calData)); % Add code to avoid saturation?

% calculate flip angle
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1) + noise_est;   % cos theta decay
% fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1)+coefs(3);   % cos theta decay
guess(1)=max(flipCalAmps);
guess(2)=20*pi/180;
% guess(3)=0;

xdata = 1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,~,residual,~,~,~,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);
ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
param_err = fitparams-ci(:,1)';
flip_angle = abs(fitparams(2)*180/pi);
flip_err = param_err(2)*180/pi;
FlipScaleFactor = flip_angle/FlipTarget;

% figure('Name','Flip Angle Calibration');  % plot up the calibration data
% plot(xdata,ydata,'bo','MarkerFaceColor','b');
% hold on;
% plot(xdata,fitfunct(fitparams,xdata),'-r');
% legend('Acquired','Fit');
% xlabel('Acquisition Number');
% ylabel('Magnitude');
% a=sprintf('Flip Cal: %0.1f%c yields flip = %0.1f±%0.1f%c', FlipTarget,char(176),flip_angle,flip_err,char(176));
% title(a); %note degree symbol is char(176)

%% Provide Flip Angle Factor and Warnings
fprintf('For flip angle calibration of %0.1f%c, the flip angle factor is %0.3f \n', FlipTarget,char(176),FlipScaleFactor);
if FlipScaleFactor>FAScale_tol
    fprintf(2,'WARNING! Excessive flip angle calibration scale factor; Check system\n');
end

%% Calculate and Report Dissolved Phase Spectral SNR
timeFit = disfitObj.calcTimeDomainSignal(disfitObj.t);  % calculate fitted data in time domain
timeRes = timeFit - disfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_dis = disfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('\nSNR Check for dissolved and gas peaks...\n'); % 
fprintf('SNR for RBC peak = %3.1f \n',SNR_dis(1)); 
fprintf('SNR for Barrier peak = %3.1f \n',SNR_dis(2)); 
fprintf('SNR for Gas peak = %3.1f \n',SNR_dis(3)); 

%% Calculate and Report Final Gas Phase FID spectral SNR
timeFit=gasfitObj.calcTimeDomainSignal(gasfitObj.t);  % calculate fitted data in time domain
timeRes=timeFit - gasfitObj.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_gas=gasfitObj.area/std25; % will be a matrix for multiple peak fit
fprintf('SNR for final gas peak = %3.1f \n',SNR_gas); %
if SNR_gas<SNR_tol
    fprintf(2,'WARNING! Gas FID SNR below minimums; Check Coil Plug\n');
end

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3)/sum(disfitObj.area(1:2));
fprintf('\nGasDisRatio = %3.3f \n',GasDisRatio); 

%% Quantify RBC:Barrier ratio
RbcBarRatio = disfitObj.area(1)/disfitObj.area(2);
fprintf('RbcBarRatio = %3.3f\n',RbcBarRatio);

%% Store Results for Additional Analysis
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
CalResults.Reference_Voltage = 0;
CalResults.dwell_time = dwell_time;  
CalResults.noise_est = noise_est;
CalResults.freq = freq; 
CalResults.FlipTarget = FlipTarget;
CalResults.te = te;  
CalResults.nDis = nDis;
CalResults.nCal = nCal;
CalResults.VRefScaleFactor = 0;
CalResults.VRef = 0;

end