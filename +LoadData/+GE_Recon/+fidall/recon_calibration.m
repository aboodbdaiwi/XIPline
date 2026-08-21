function [meanRbc2barrier,te90,targetAX,targetTG] = recon_calibration(d,h,wfn,delay,fname)
%RECON_CALIBRATION Reconstruct 129Xe calibration
% [bb,bbabs] = recon_calibration(d,h,wfn,mtx,delay,lb,fname,comb_time)
%                                                               (default)
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%       wfn  Location of waveform .mat file
%     delay  Gradient-acquisition delay                   [pts] (0)
%     fname  Print <fname>.png and save reco as <fname>.mat     ([]) 
%            also export dicom if template.dcm is found
%
%        bb  Reconstructed data       (mtx,mtx,mtx,#timesteps,#coils)
%     bbabs  RMS coil-combined data   (mtx,mtx,mtx,#timesteps)
%
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
if ~exist('delay','var'), delay = []; end
if isempty(delay),      delay = 0; end
if length(delay)~=1, error('length(delay)(=%g)~=1',length(delay)); end
if ~exist('fname','var'), fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    end
end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end

%% reading frequencies and flip angles
freq_array = read_fdl([wfn(1:end-4) '_freq.fdl']);
flip_array = read_fdl([wfn(1:end-4) '_flip.fdl']);

%% check fields required for reconstruction
fina = {'ts','n_ind'};
for l=1:length(fina)
    if ~isfield(wf,fina{l}), error('wf.%s not existing',fina{l}); end
end

%% Separate spectra

skipped_frames_diss = 50;                        % Number of spectra to ignore for downstream saturation
skipped_frames_gas = 20;
%bw = h.rdb_hdr.user0;                       % full acq bandwidth [Hz]
%[nexc,~,nphases,nechoes,nslices,ncoils] = size(d); % data size

dps = freq_array ~= 0;
kdata_dissolved = d(dps==1,2:end);
dps(1:skipped_frames_diss) = 0;
gps = freq_array == 0;

skipped_frames_gas = min(skipped_frames_gas,ceil(sum(gps)/3));
skipped_frames_gas = 0;
data_dissolved = d(dps == 1,2:end);
data_gas = d(gps(1:size(d,1)) == 1,2:end);
data_gas(1:skipped_frames_gas,:) = [];
end_frame = min(size(data_gas,1),80);
data_gas = data_gas(1:end_frame,:);


if mean(angle(mean(data_dissolved(1:2:end,:)./data_dissolved(2:2:end,:)))) < -3
    data_dissolved(2:2:end,:) = data_dissolved(2:2:end,:).*exp(1i*pi);
    data_gas(2:2:end,:) = data_gas(2:2:end,:).*exp(1i*pi);
elseif mean(angle(mean(data_dissolved(1:2:end,:)./data_dissolved(2:2:end,:)))) > 3
    data_dissolved(2:2:end,:) = data_dissolved(2:2:end,:).*exp(-1i*pi);
    data_gas(2:2:end,:) = data_gas(2:2:end,:).*exp(-1i*pi);
end


data_avg = mean(data_dissolved);

%% Do the spectral fitting
rbc_idx = 3;
barrier_idx = 2;
gas_idx = 1;

area_orig = [.1 1 1];
freq_orig = [max(abs(freq_array)) 800 -100];
fwhm_orig = [1 250 250];
phase_orig = pi*randn(1,3);

dissolvedFit = NMR_Fit(data_avg, wf.ts(2:end), area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
[area_f, freq_f, fwhm_f, phase_f] = calcTimeDomainSignalFit(dissolvedFit);
dissolvedFit.area = area_f;
dissolvedFit.freq = freq_f;
dissolvedFit.fwhm = fwhm_f;
dissolvedFit.phase = phase_f;

if dissolvedFit.freq(rbc_idx) > 0
    rbc_idx = 2;
    barrier_idx = 3;
end

%% RBC:Barrier Ratio
meanRbc2barrier = dissolvedFit.area(rbc_idx)/dissolvedFit.area(barrier_idx);

%% Compute TE90
te = h.rdb_hdr.te*1e-6;

deltaF = freq_f(barrier_idx)-freq_f(rbc_idx);
deltaPhase = phase_f(barrier_idx)-phase_f(rbc_idx);
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
area_orig = 1000;
freq_orig = 0;
fwhm_orig = 1;
phase_orig = 0;

Ng = size(data_gas,1);
area_g = zeros(1,Ng);
freq_g = zeros(1,5);

for s = 1:5
    gasFit = NMR_Fit(data_gas(s,:), wf.ts(2:end), area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
    [~,freq_g(s),~,~] = calcTimeDomainSignalFit(gasFit);
end

%% Compute AX
center_freq = h.rdb_hdr.ps_mps_freq/10;
targetAX = round(center_freq - median(freq_g));
area_g(1,:) = abs(data_gas(:,4));

flip_fit = fit((1:Ng)',area_g','exp1');
targetFA = flip_array(end);                     % This is the target we wanted to hit
actualFA = acosd(exp(flip_fit.b));                          % This is the estimate


% Dixon TG - compute the TG that provides the targetFA
targetTG = ceil(h.rdb_hdr.ps_mps_tg + 200*log10(targetFA./actualFA));

%% saving result as mat file
if ~isempty(fname)
    save([fname '.mat'],'-mat','-v7.3',...
        'targetTG','targetAX','te90','meanRbc2barrier','data_dissolved','data_gas','data_avg','area_f','freq_f','fwhm_f','phase_f','area_g','fname');
end

%% saving data for external processing
kdata_dissolved = kdata_dissolved.';
tr = h.image.tr*1e-6;
spec_tsp = wf.ts(2);
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
subplot(212),plot(F,abs(dspec),'r',F,abs(calcSpectralDomainSignal(dissolvedFit,F).*spec_scale),'b'); title('Magnitude');legend('Measured','Estimated','Location','NorthWest');
text(F(10),max(abs(dspec))*.5,ds_text);

% Save figure to file
figstr = sprintf('P%05d Exam%d Series%d',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
set(gcf,'name',figstr);
if ~isempty(fname)
    print(fname,'-dpng','-r600');
end


end      % main function recon_calibration.m


