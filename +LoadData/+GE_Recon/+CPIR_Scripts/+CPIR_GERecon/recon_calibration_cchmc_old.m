function [meanRbc2barrier,te90,targetAX,targetTG] = recon_calibration_cchmc(d,h,sv)
%recon_calibration_cchmc  Display spectra as stacked plot
%  recon_calibration_cchmc(d,h,lb,how,sv)
%   d  raw P-file data
%      optionally filename to P-file or image archive
%   h  header structure
%  sv  Print pix and save data to file   [logical]    (false)
%      alternatively: filename for printing
%
% 12/2023 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% input parameters
if ~exist('sv','var'),  sv = []; end
if isempty(sv),         sv = false; end

%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
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
freq_array_path = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Scripts\cal\cal_cpir_nexc520_ndiss500_dfreq202_freq.fdl';
freq_array = read_fdl(freq_array_path);

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
% Compute the gas peak area for each specta
Ng = size(data_gas,1);
area_g = zeros(1,Ng);
area_g(1,:) = abs(data_gas(:,4));

flip_fit = fit((1:Ng)',area_g','exp1');
targetFA = h.image.mr_flip;                     % This is the target we wanted to hit
actualFA = acosd(exp(flip_fit.b));                          % This is the estimate

% Dixon TG - compute the TG that provides the targetFA
targetTG = ceil(h.rdb_hdr.ps_mps_tg + 200*log10(targetFA./actualFA));

%% Compute AX
freq_g = zeros(1,5);
center_freq = h.rdb_hdr.ps_mps_freq/10;
area_orig = max(abs(data_gas(:)));
freq_orig = 0;
fwhm_orig = 1;
phase_orig = 0;
for s = 1:5
    gasFit = Calibration.NMR_Fit(double(data_gas(s,:)), double(ts), double(area_orig),...
        double(freq_orig),double(fwhm_orig),double(phase_orig),[],[]);
    [~,freq_g(s),~,~] = calcTimeDomainSignalFit(gasFit);
end
targetAX = round(center_freq + median(freq_g));

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
figure; set(gcf,'DefaultAxesFontSize',14);
fa_text = sprintf('Target FA = %.1f\nActual FA = %.1f\nTarget TG = %d',targetFA,actualFA,targetTG);
fit_line = flip_fit(1:Ng);
subplot(211),plot(1:Ng,area_g,'b',1:Ng,fit_line,'r'); axis([1, Ng, 0.8*min(area_g), 1.2*max(area_g)]);
text(round(Ng*.65), max(area_g), fa_text);

F = dissolvedFit.f;
dspec = fftshift(fft(data_avg));             % measured average spectrum
spec_scale = max(abs(dspec))./max(abs(calcSpectralDomainSignal(dissolvedFit,F)));   % scaling factor for display
ds_text = sprintf('RBC:Barrier = %1.2f\nGas Frequency = %dHz\nTE90=%1.3fms\n',meanRbc2barrier,targetAX,te90*1e3);

% Plot measured and fitted spectra
subplot(212),plot(F,abs(dspec),'r',F,abs(calcSpectralDomainSignal(dissolvedFit,F).*spec_scale),'b'); title('Magnitude');legend('Measured','Estimated','Location','northeast');
xlim([F(1), F(end)])
set ( gca, 'xdir', 'reverse' )
text(F(end-1),max(abs(dspec))*.5,ds_text);

% Save figure to file
figstr = sprintf('P%05d Exam%d Series%d',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
set(gcf,'name',figstr);
if ~isempty(fname)
    print(fname,'-dpng','-r600');
    saveas(gcf,[fname '.fig'],'fig');
    if true
        inp.InstanceNumber = h.rdb_hdr.user4;
        write_scdicom([fname '.dcm'],gcf,h,inp);
    end
end


end      % main function recon_calibration_cchmc.m


