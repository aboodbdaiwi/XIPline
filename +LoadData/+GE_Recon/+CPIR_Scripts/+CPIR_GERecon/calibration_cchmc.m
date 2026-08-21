function [GasExResults, CalResults] = calibration_cchmc(MainInput)

% calibration_cchmc
% Process GE calibration data for XIPline.
%
% Inputs:
%   MainInput.CalFullPath - path to GE calibration P-file
%   MainInput.Cal_h       - optional GE header structure
%   MainInput.Cal_sv      - optional save flag / filename
%
% Outputs:
%   GasExResults
%   CalResults

if nargin < 1
    help(mfilename);
    return;
end

d = MainInput.CalFullPath;
fname = MainInput.CalFileName;
if isfield(MainInput,'Cal_h')
    h = MainInput.Cal_h;
else
    h = [];
end

if isfield(MainInput,'Cal_sv')
    sv = MainInput.Cal_sv;
else
    sv = false;
end
%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = LoadData.GE_Recon.fidall.read_p(d,true);
    d=conj(d);
end

%% filenames
if isempty(h), error('isempty(h)'); end

if ischar(sv)
    % if filename given in sv field
    fname = sv;
    if ~isempty(regexpi(fname,'\.7$', 'once')),  fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
    sv = true;
else
    % remove suffix
    if sv
        tmp = regexpi(fname,'\.');
        if ~isempty(tmp)
            fname = fname(1:(tmp(end)-1));
        else
            fname = fname;
        end
    end
end

%% reading frequencies and flip angles
% use h.rdb_hdr.user14 (rf pulse shape) to get cv number to get correct array
freq_array_name = sprintf('vap_freq%d.fdl',h.rdb_hdr.user14);
freq_array_path = fullfile('~/fidall/parameters',freq_array_name);

% Get the folder containing this calibration_cchmc.m function
thisFunctionPath = mfilename('fullpath');
thisFunctionDir = fileparts(thisFunctionPath);

% calibration_cchmc.m is located in:
% CPIR_Scripts/+CPIR_GERecon/
%
% Go up one folder to CPIR_Scripts, then enter:
% CPIR_Scripts/cal/
CPIRScriptsDir = fileparts(thisFunctionDir);
calDir = fullfile(CPIRScriptsDir, 'cal');

% Select the appropriate frequency array based on institute
if strcmp(MainInput.Institute, 'CCHMC')
    freq_array_path = fullfile( ...
        calDir, ...
        'cal_cpir_nexc520_ndiss500_dfreq202_freq.fdl');
elseif strcmp(MainInput.Institute, 'XeCTC')
    freq_array_path = fullfile( ...
        calDir, ...
        'cal_xectc_nexc520_ndiss500_dfreq218_freq.fdl');
else
    error('Unsupported institute: %s', MainInput.Institute);
end

% Make sure the frequency array file actually exists
if ~isfile(freq_array_path)
    error('Calibration frequency array file not found:\n%s', ...
        freq_array_path);
end

freq_array = LoadData.GE_Recon.fidall.read_fdl(freq_array_path);

%% Separate spectra
dps = freq_array ~= 0;
kdata_dissolved=d(dps==1,:);

skipped_frames_diss = 50;
dps(1:skipped_frames_diss) = 0;
data_dissolved = d(dps == 1,:);
data_avg = mean(data_dissolved);

gps = freq_array == 0;
data_gas = d(gps == 1,:);

%% Do the spectral fitting
rbc_idx = 1;
barrier_idx = 2;
gas_idx = 3;
dwell = 1/h.rdb_hdr.spectral_width;
ts = dwell*(1:size(d,2));

%initial guesses
area_guess = double([0.27*abs(data_avg(1)) abs(data_avg(1)) 0.1*abs(data_avg(1))]); 
area_lowerBounds = double([0 0 0]); % Always positive
area_upperBounds = double(realmax/10*[1 1 1]); % Make it finite, but huge
freq_guess = double([500 -200 -7200]);
freq_lowerBounds = double([0 -2000 -9000]);
freq_upperBounds = double([2000 0 -5000]);
fwhm_guesses = double([200 200 25]);
fwhm_lowerBounds = double([0 0 0]);
fwhm_upperBounds = double(inf*[1 1 1]);
% fwhmG_guess = double([0 500 0]);
% fwhmG_lowerBounds = double([0 250 0]);
% fwhmG_upperBounds = double([0 inf 0]);
phase_guesses = double([-75 -140 -70]);
phase_lowerBounds = double(-inf*[1 1 1]);
phase_upperBounds = double(inf*[1 1 1]);

% dissolvedFit = NMR_TimeFit_v(double(data_avg),double(ts),area_guess,freq_guess,fwhm_guesses,fwhmG_guess,phase_guesses,[],[]);
dissolvedFit = Calibration.NMR_TimeFit(double(data_avg),double(ts),area_guess,freq_guess,fwhm_guesses,phase_guesses,[],[]);
dissolvedFit.fitTimeDomainSignal();%initial fit
dissolvedFit.setBounds(area_lowerBounds,area_upperBounds,...
    freq_lowerBounds,freq_upperBounds,...
    fwhm_lowerBounds,fwhm_upperBounds,...
    phase_lowerBounds,phase_upperBounds);%add bounds
dissolvedFit.fitTimeDomainSignal();%refit

%% RBC:Barrier Ratio
meanRbc2barrier = dissolvedFit.area(rbc_idx)/dissolvedFit.area(barrier_idx);

%% Compute TE90
te = h.rdb_hdr.te*1e-6;

deltaF = dissolvedFit.freq(barrier_idx)-dissolvedFit.freq(rbc_idx);
deltaPhase = dissolvedFit.phase(barrier_idx)-dissolvedFit.phase(rbc_idx);
time180 = abs(1/(2*deltaF));
te90 = (90-deltaPhase*(1+360/deltaPhase*(deltaPhase < -90)))/(360*deltaF);

while(te90<0)
    % This TE is too low, so add 180 deg of phase
    te90 = te90 + time180;
end
while(te90>time180)
     % This te is too high, so subtract 180 deg of phase
    te90 = te90 - time180;
end

% make te90 a multiple of 4us + 2us
te90 = ceil(te90*1e6);
if mod(te90,4)
    te90 = te90 + 4 - mod(te90,4);
end
te90 = (te90 + 2)*1e-6;

te90 = te90 + te;

%% Estimate flip angle
% Compute the gas signal magnitude for each calibration spectrum
Ng = size(data_gas,1);
area_g = zeros(1,Ng);
area_g(:) = abs(data_gas(:,4));

% Target flip angle prescribed on the scanner
targetFA = h.image.mr_flip;

% -------------------------------------------------------------------------
% Fit the GE gas-decay data using the SAME model expected by XIPline:
%
%     S(n) = A * cos(alpha)^(n-1) + noise
%
% This makes CalResults.DecayFit directly compatible with the existing
% XIPline plotting code.
% -------------------------------------------------------------------------

xdata = 1:Ng;
ydata = reshape(area_g,1,[]);

% Estimate noise floor from the end of the decay.
% Keep it non-negative and small relative to the measured signal.
nTail = max(3,round(0.15*Ng));
noise_est = median(ydata(end-nTail+1:end));

% The tail still contains real signal, so avoid treating the entire tail
% magnitude as noise if the decay has not reached the noise floor.
noise_est = min(noise_est,0.05*max(ydata));

% Fit parameters:
%   fitparams(1) = initial signal amplitude
%   fitparams(2) = flip angle in radians
fitfunct = @(coefs,x) ...
    coefs(1).*cos(coefs(2)).^(x-1) + noise_est;

guess = [max(ydata)-noise_est, targetFA*pi/180];

% Physically reasonable bounds
lb = [0, 1*pi/180];
ub = [Inf, 45*pi/180];

fitoptions = optimoptions('lsqcurvefit','Display','off');

[fitparams,~,residual,~,~,~,jacobian] = lsqcurvefit( ...
    fitfunct, ...
    guess, ...
    xdata, ...
    ydata, ...
    lb, ...
    ub, ...
    fitoptions);

% Actual flip angle
actualFA = abs(fitparams(2)*180/pi);

% Estimate 95% confidence interval / uncertainty
try
    ci = nlparci(fitparams,residual,'jacobian',jacobian);
    param_err = fitparams-ci(:,1)';
    flip_err = abs(param_err(2)*180/pi);
catch
    flip_err = NaN;
end

% Flip-angle scale factor
FlipScaleFactor = actualFA/targetFA;

% Dixon TG - compute the TG that provides the target FA
targetTG = ceil( ...
    h.rdb_hdr.ps_mps_tg + ...
    200*log10(targetFA./actualFA));

%% Compute gas frequency / AX

center_freq = h.rdb_hdr.ps_mps_freq/10;

% -------------------------------------------------------------------------
% First retain the original GE method for determining the robust frequency
% correction: fit the first five gas spectra individually and use their
% median frequency.
% -------------------------------------------------------------------------

nGasFreqFits = min(5,Ng);

freq_g = zeros(1,nGasFreqFits);

area_orig = max(abs(data_gas(:)));
freq_orig = 0;
fwhm_orig = 25;
phase_orig = 0;

for s = 1:nGasFreqFits

    geGasFit = ...
        LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.NMR_Fit( ...
        double(data_gas(s,:)), ...
        double(ts), ...
        double(area_orig), ...
        double(freq_orig), ...
        double(fwhm_orig), ...
        double(phase_orig), ...
        [],[]);

    [~,freq_g(s),~,~] = ...
        LoadData.GE_Recon.fidall.calcTimeDomainSignalFit(geGasFit);

end

% Robust frequency offset from the GE calibration
gas_freq_offset = median(freq_g);

% Absolute scanner frequency
targetAX = round(center_freq + gas_freq_offset);

% -------------------------------------------------------------------------
% Generate a GasFit object compatible with the existing XIPline plotting
% code.
%
% Use the gas spectrum closest to the median frequency so that the displayed
% spectrum is representative and is not dependent on whichever spectrum
% happened to be processed last.
% -------------------------------------------------------------------------

[~,representativeGasIndex] = ...
    min(abs(freq_g-gas_freq_offset));

gasDataForFit = data_gas(representativeGasIndex,:);

% Initial estimates
gas_area_guess = max(abs(gasDataForFit));
gas_freq_guess = gas_freq_offset;
gas_fwhm_guess = 25;
gas_phase_guess = angle(gasDataForFit(1))*180/pi;

% Use the same NMR_TimeFit object type expected elsewhere in XIPline
gasfitObj = Calibration.NMR_TimeFit( ...
    double(gasDataForFit), ...
    double(ts), ...
    double(gas_area_guess), ...
    double(gas_freq_guess), ...
    double(gas_fwhm_guess), ...
    double(gas_phase_guess), ...
    0, ...
    length(ts));

gasfitObj.fitTimeDomainSignal();

% Use this object for XIPline
gasFit = gasfitObj;

%% Store Results for XIPline

% Results to Pass to Gas Exchange Recon
RbcBarRatio = meanRbc2barrier;

GasDisRatio = dissolvedFit.area(gas_idx) / ...
    sum(dissolvedFit.area([rbc_idx barrier_idx]));

GasExResults.RbcBarRatio = RbcBarRatio;
GasExResults.GasDisRatio = GasDisRatio;
GasExResults.DisFit = dissolvedFit;

% Results to Pass to Calibration information
CalResults.flip_angle = actualFA;
CalResults.flip_err = flip_err;
CalResults.FlipScaleFactor = FlipScaleFactor;

CalResults.Pulses = xdata;
CalResults.GasDecay = ydata;

% IMPORTANT:
% These are now [amplitude, flip-angle-radians], exactly what the
% existing XIPline fitfunct expects.
CalResults.DecayFit = fitparams;

% Proper NMR_TimeFit gas object for XIPline plotting
CalResults.GasFit = gasfitObj;

CalResults.te90 = te90*1e3;            % ms
CalResults.freq_target = targetAX;

CalResults.Reference_Voltage = 0;
CalResults.dwell_time = dwell;
CalResults.freq = center_freq;
CalResults.FlipTarget = targetFA;
CalResults.te = te*1e6;                % us
CalResults.nDis = size(kdata_dissolved,1);
CalResults.nCal = Ng;
CalResults.VRefScaleFactor = 0;
CalResults.VRef = 0; 
CalResults.targetTG = targetTG;
% Same fixed additive term used by the fitted decay model
CalResults.noise_est = noise_est;

%% saving result as mat file
if ~isempty(fname)
    save([fname '.mat'],'-mat','-v7.3',...
        'targetTG','targetAX','te90','meanRbc2barrier','data_dissolved','data_gas','data_avg','dissolvedFit','area_g','fname');
end

%% saving data for external processing
kdata_dissolved = kdata_dissolved.';
tr = h.image.tr*1e-6;
spec_tsp = ts(2);
date = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];

[fp,fn] = fileparts(fname);

save([fp '/Spect_' fn '.mat'],'-mat','-v7.3',...
        'kdata_dissolved','tr','center_freq','spec_tsp','date');

%% plotting
CalFig = figure;
set(CalFig,'DefaultAxesFontSize',14);

fa_text = sprintf( ...
    'Target FA = %.1f\nActual FA = %.1f\nTarget TG = %d', ...
    targetFA,actualFA,targetTG);

fit_line = fitfunct(fitparams,xdata);

ax1 = subplot(2,1,1,'Parent',CalFig);
plot(ax1,xdata,ydata,'bo','MarkerFaceColor','b');
hold(ax1,'on');
plot(ax1,xdata,fit_line,'-r','LineWidth',2);
hold(ax1,'off');
axis(ax1,[1, Ng, 0.8*min(area_g), 1.2*max(area_g)]);
text(ax1,round(Ng*.65),max(area_g),fa_text);

F = dissolvedFit.f;
dspec = fftshift(fft(data_avg));  % measured average spectrum

fitSpectrum = calcSpectralDomainSignal(dissolvedFit,F);
spec_scale = max(abs(dspec)) ./ max(abs(fitSpectrum));

ds_text = sprintf( ...
    'RBC:Barrier = %1.2f\nGas Frequency = %dHz\nTE90=%1.3fms\n', ...
    meanRbc2barrier,targetAX,te90*1e3);

% Plot measured and fitted spectra
ax2 = subplot(2,1,2,'Parent',CalFig);
plot(ax2,F,abs(dspec),'r', ...
         F,abs(fitSpectrum .* spec_scale),'b');

title(ax2,'Magnitude');
legend(ax2,'Measured','Estimated','Location','northeast');
xlim(ax2,[F(1),F(end)]);
set(ax2,'XDir','reverse');

text(ax2,F(end-1),max(abs(dspec))*.5,ds_text);

% Save figure to file
figstr = sprintf('P%05d Exam%d Series%d', ...
    h.image.rawrunnum, ...
    h.exam.ex_no, ...
    h.series.se_no);

set(CalFig,'Name',figstr);

if ~isempty(fname)

    exportgraphics(CalFig,[fname '.png'],'Resolution',600);

    savefig(CalFig,[fname '.fig']);

    inp.InstanceNumber = h.rdb_hdr.user4;
    LoadData.GE_Recon.fidall.write_scdicom([fname '.dcm'],CalFig,h,inp);

end
close all;
end      % main function recon_calibration_cchmc.m


