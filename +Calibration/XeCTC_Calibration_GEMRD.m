function [GasExResults, CalResults] = XeCTC_Calibration_GEMRD(MainInput)
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

%% Define number of calibration, spectral fit and display toggles,number to skip, number of gas, target flip

% user control of  outputs
disp_fids = 1; % set to 1 to display first (noise) frame and first FID
save_csv =1; % export derived variables as csv file

% user control of skipping and averaging
seconds2skip = 2; %number of seconds to skip in breath hold
seconds2avg = 1;
nDis = 500; % Assume consortium standard
FlipTarget = 10; % target flip angle from calibration

% frequency guesses in ppm for dissolved phase fits
rbc_freq_ppm = 217.2;
mem_freq_ppm = 197.7;
gas_freq_ppm = 0;

%% define tolerances for warnings
freq_tol = 100; % frequency adjustment tolerance
deltaPhase1_tol = 90; % if calibration phase>90 use minTE
TrueRefScale_tol = 1.33;
SNR_tol = 50; % Ideal is 100 for RBC SNR, but allow 50

%% first read in the file and create the twix object
% [file, path] = uigetfile('*.*', 'Select file'); % starts in current dir
% file_with_path = strcat(path, file); % join path and filename to open
% [~,filename,ext] = fileparts(file_with_path); % get file extension
cali_struct = LoadData.ismrmrd.GE.Functions.readRawCali(MainInput.CalFileName);
cali_struct.data = cali_struct.data(:,1:515);
filename = MainInput.CalFileName;
ext = MainInput.CalDataext;
if strcmp(ext,'.dat') 
   filename = file(15:end-4); % return filename back to something short
else
   filename = [filename,'_cali'];
end   
cal_vars_export = [filename,'.csv']; % generate the filename for data comparison
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems
file_loc = regexp(path, filesep, 'split');
file_loc = file_loc{end-1};

%% Prepare variables
% extract favorite variables from the data struct
weight = cali_struct.weight;
te = cali_struct.te;
tr = cali_struct.tr;
dwell_time = cali_struct.dwell_time;
freq = cali_struct.freq;
xeFreqMHz = cali_struct.xeFreqMHz;
theFID = cali_struct.data;
nFids = size(theFID, 2);
nCal = nFids-nDis; % assume remaining FIDS past dissolved are cal
nPts = size(theFID, 1);
VRef = cali_struct.vref;
scanDateStr = cali_struct.scan_date;
rf_excitation_ppm = cali_struct.rf_excitation_ppm;
disp(['YESSIR, rf_excitation_ppm = ', num2str(rf_excitation_ppm)]);

% % print out key variable values
% fprintf('\n\n');
% fprintf('\tFile Location = %s\n', file_loc);
% fprintf('\tName = %s\n', cali_struct.seq_name);
% fprintf('\tScan Date = %s\n', scanDateStr');
% fprintf('\tWeight = %0.1f kg \n', weight);
% fprintf('\tTE = %0.0f us\n', te);
% fprintf('\tTR = %0.0f us\n', tr(1)); % some sequences now have an array of TRs for bonus spectra
% fprintf('\tFrequency = %0.0f Hz\n', freq);
% fprintf('\tRF Excitation Offset = %0.0f ppm\n', rf_excitation_ppm);
% fprintf('\tReference Voltage = %0.1f V\n', VRef);
% fprintf('\tDwell Time = %0.0f ns\n', dwell_time*1e9);
% fprintf('\tnPts = %0.0f\n', nPts);
% fprintf('\tnFrames = %0.0f\n', nFids);
% fprintf('\tnDissolved Frames = %0.0f\n', nDis);
% fprintf('\tnGas Frames = %0.0f\n\n', nCal);

%% parse out the various fids - dissolved, gas, and calibration
calData = theFID(:, end-nCal+1:end); % the data left for flip angle calculations
gasData = theFID(:, end-nCal+1); % the first gas frame for frequency calculations

tr_s = tr(1) * 1e-6; %tr in seconds
t_tr = tr_s * (1:nFids);
nSkip = round(seconds2skip/tr_s);
nAvg = round(seconds2avg/tr_s);
t = dwell_time * (0:(nPts - 1))';
disData = theFID(:, nSkip:nSkip+nAvg); % just the dissolved data of interest
noise_est = max(abs(disData(end-10:end,end)));
disData_avg = mean(disData,2); % average the dissolved data

%% Fit gas Spectrum
fprintf('Analysis of Gas FID\n');
% gasfitObj = NMR_TimeFit(gasData, t, 1e-4, -84, 30, 0, 0, 10000);
gasfitObj = Calibration.NMR_TimeFit(gasData, t, 1e-10, 1000, 30, 0, 0, 10000);
gasfitObj.fitTimeDomainSignal();
% h1 = figure('Name', 'Gas Phase Analysis');
% gasfitObj.plotTimeAndSpectralFit;
% xlim([0, 0.01]);
gasfitObj.describe();

%% calculate target transmit receive frequency from gas fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n', freq_target); % report new target frequency
if abs(freq_target-freq) > freq_tol
    fprintf(2, 'WARNING! Frequency adjust exceeds tolerances; Check system\n');
end

%% Perform Flip Angle Calibration on Gas Phase Signals
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs, xdata)coefs(1) * cos(coefs(2)).^(xdata - 1) + noise_est; % cos theta decay
guess(1) = max(flipCalAmps);
guess(2) = 20 * pi / 180; % just guess 10 degrees
% guess(3) = noise_est; % with absolute values there will be a baseline offset per Matt Wilmering 4/12/21

xdata = 1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit', 'Display', 'off');
[fitparams, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(fitfunct, guess, xdata, ydata, [], [], fitoptions);
ci = nlparci(fitparams, residual, jacobian); % returns 95% conf intervals on fitparams by default
param_err = fitparams - ci(:, 1)';
flip_angle = abs(fitparams(2)*180/pi);
flip_err = param_err(2) * 180 / pi;
FlipScaleFactor = flip_angle/FlipTarget;

% h2=figure('Name', 'Flip Angle Calibration'); % plot up the calibration data
% plot(xdata, ydata, 'bo', 'MarkerFaceColor', 'b');
% hold on;
% plot(xdata, fitfunct(fitparams, xdata), '-r');
% legend('Acquired', 'Fit');
% xlabel('Frame Number');
% ylabel('Magnitude');
% a = sprintf('Flip Cal: V_{ref} = %0.1f yields flip = %0.1f�%0.1f%c', VRef, flip_angle, flip_err, char(176));
% title(a); % note degree symbol is char(176)
% pause(1); % had to insert pause to avoid occasional plotting conflicts - don't know why

%% Provide New Reference Amplitude and Warnings

if isnan(VRef)
   hint_refV = '\bf Reference Voltage: \rmNot found';
   fprintf(2,'WARNING! Reference Voltage not found.\n');
else   
   fprintf(['For ', num2str(VRef), 'V calibration, True_Ref = %3.0f V\n'], VRef*FlipScaleFactor); 
   fprintf('(Estimated Weight Reference Voltage = %3.0f V)\n', 395.3+1.1431*weight); % Ari Bechtel estimate, new fit
   hint_refV = ['\bf Reference Voltage: \rm', num2str(VRef*FlipScaleFactor, 3), 'V'];
end   
if FlipScaleFactor > TrueRefScale_tol
    fprintf(2, 'WARNING! Excessive voltage calibration scale factor; Check system\n');
end 

%% Fit dissolved Spectrum - Make this last so user sees it
fprintf('\nAnalysis of Averaged Dissolved FID, skipping %0.0f\n', nSkip);
fprintf('Using Voigt Lineshape to fit Membrane\n');

% Adjust frequency guesses based on rf_excitation then to Hz
rbc_freq_adj = rbc_freq_ppm - rf_excitation_ppm;
mem_freq_adj = mem_freq_ppm - rf_excitation_ppm;
gas_freq_adj = gas_freq_ppm - rf_excitation_ppm;

freq_guess = [rbc_freq_adj, mem_freq_adj, gas_freq_adj] * xeFreqMHz; % in Hz

%set all other initial parameter guesses
area_guess =  [1, 1, 1]; % no benefit seen in tailoring these guesses
fwhmL_guess = [8.8, 5.0, 2] * xeFreqMHz;
fwhmG_guess = [0, 6.1, 0] * xeFreqMHz;
phase_guess = [0, 0, 0]; % no benefit seen in tailoring these guesses

% first widths lorenzian, 2nd are gaussian
disfitObj = Calibration.NMR_TimeFit_v(disData_avg, t, area_guess, freq_guess, fwhmL_guess, fwhmG_guess, phase_guess, [],[]);
disfitObj = disfitObj.fitTimeDomainSignal();
% h3 = figure('Name', 'Dissolved Phase Analysis'); 
% disfitObj.plotTimeAndSpectralFit;
% % update axis limits to see gas signal
% ax = findall(h3, 'type', 'axes'); % find all the axes in figure (ax1 defined in NMR_TimeFit)
% set(ax(1),'xlim',[0 .01]); % narrow gas plot limits for meaningful time
% set(ax(5),'xlim',[-8000 2000])  % expand plot limits to show gas peak

disp('            Area     Freq (Hz)   Linewidths(Hz)   Phase(degrees)');
peakName = {'     RBC:', 'Membrane:', '     Gas:'};
for iComp = 1:length(disfitObj.area)
    disp([peakName{iComp}, '  ', ...
        sprintf('%8.3e', disfitObj.area(iComp)), ' ', ...
        sprintf('%+8.2f', disfitObj.freq(iComp)), '  ', ...
        sprintf('%8.2f', abs(disfitObj.fwhm(iComp))), '  ', ...
        sprintf('%8.2f', abs(disfitObj.fwhmG(iComp))), '  ', ...
        sprintf('%+9.2f', disfitObj.phase(iComp))]);
end

%% Calculate and report TE90 and Warnings
deltaPhase = disfitObj.phase(2) - disfitObj.phase(1); % RBC-Membrane phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase), 180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-membrane freq difference
deltaTe90 = (90 - deltaPhase) / (360 * deltaF); % how far off are we?
te90 = (te + deltaTe90 * 1e6)/1000; % in usec

% Hint the actual TE90 for usage
if te90 < 0.445 
   hint_te90 = '\color[rgb]{1,0,0} WARNING! Low TE90; Use 0.45ms.';
elseif  te90 >= 0.505
   hint_te90 = '\color[rgb]{1,0,0} WARNING! High TE90; Use 0.50ms.';
else
   hint_te90 = ['\color[rgb]{0,0,0} Use ',num2str(te90, 2), 'ms.'];
end    
fprintf('\nTE90 = %3.2f ms\n', te90); % report TE90 in ms
fprintf(2,[hint_te90(20:end),' \n']);

%% Calculate and report dissolved phase spectral SNR
pause(0.1);
disData_tail = disData(nPts/2+1:end,:); % just take the second half of FIDs
mean_disData_tail = mean(disData_tail,2);
diff = mean_disData_tail - disData_tail; % array of difference fids
diff = diff - mean (diff); % subtract off any DC bias
noiseFrame = diff(:,1); % take first frame as representative noise frame
noiseDis = std(real(diff));
SNR_frames_d = disfitObj.area'./noiseDis; % SNR of each frame
SNRsnf_d = mean(SNR_frames_d,2)*sqrt(nAvg); % simulated noise frame SNR
noise_mean = mean(noiseDis,2); % save mean noise frame for gas calculation

fprintf('\nSNR Check for dissolved and gas peaks...\n'); 
fprintf('SNR of dissolved spectra (Dai method):\n');
fprintf('SNR for RBC peak = %3.1f \n', SNRsnf_d(1));
fprintf('SNR for Membrane peak = %3.1f \n', SNRsnf_d(2));
fprintf('SNR for Gas peak = %3.1f \n', SNRsnf_d(3));

% Calculate and report SNR of analyzed gas FID using the Dai method
SNR_gas_frame = gasfitObj.area'./noise_mean; % use noise calculated from dissolved
fprintf('SNR for dedicated gas peak = %3.1f\n',SNR_gas_frame); 

%% Quantify RBC:Membrane ratio
RbcMemRatio = disfitObj.area(1) / disfitObj.area(2);
fprintf('\nRbcMemRatio = %3.3f\n', RbcMemRatio); %

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3) / sum(disfitObj.area(1:2));
fprintf('GasDisRatio = %3.3f \n', GasDisRatio);

%% Store Results for Additional Analysis
%Results to Pass to Gas Exchange Recon
GasExResults.RbcBarRatio = RbcMemRatio;
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