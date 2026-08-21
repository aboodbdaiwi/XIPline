function [fitres,resnames,params,pk] = fit_oxsa(spec,h,lb,fname,...
    plt,pk_fname,peaks,offset,time_delay,figoff,reverse_freq)
%FIT_OXSA  Fit reconstructed spectra using OXSA AMARES
%[fitres,resnames,params,pk] = fit_oxsa(spec,h,lb,fname,plt,...
%                    pk_fname,peaks,offset,time_delay,figoff,reverse_freq)
%        spec   Reconstructed spectra [ns,nx,ny,nz] or [nexc,ns]
%           h   GE P-file header structure
%          lb   Applied Gaussian linebroadening -> deapodisation [Hz]  ([])
%       fname   Save to file (if fname given; or true->fname from fnin)('')
%         plt   Plotting: plt(1): 0=off,1=final,2=all,3=random15%
%               plt(2) = save amplitude as dicom
%               plt(3) 1=imagesc_row, 2=imagesc_ind3d             ([3 0 1])
%    pk_fname   OXSA prior knowledge file name                         ('')
%       peaks   Peaks to include in fitting      []->all               ([])
%      offset   Frequency offset                 []->from max   [ppm]  ([])
%  time_delay   Acquisition delay (linear phase) []->from hdr   [s]     (0)
%      figoff   Figure ID offset                                      (300)
%reverse_freq   Invert frequency axis (fid=conj(fid))               (false)
%               for recon_mrsi/csi prior to 5/2024
%
%      fitres   Fit result maps (nx,ny,nz,npeaks,nres)
%    resnames = {'amplitude','chemShift','linewidth','crbamp'}
%      params   OXSA parameter structure
%          pk   Prior knowledge structure
%
%  2/2025 Rolf Schulte (input from Michael Vaeggemose and Alixander Khan)
if nargin<1, help(mfilename); return; end

if ~exist('lb','var'),       lb = []; end
if ~exist('fname','var'),    fname = ''; end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = 3; end
if length(plt)<2,            plt(2) = 0; end
if length(plt)<3,            plt(3) = 1; end
if ~exist('pk_fname','var'), pk_fname = ''; end
if ~exist('offset','var'),   offset = []; end
if ~exist('time_delay','var'), time_delay = 0; end
if ~exist('figoff','var'),   figoff = []; end
if isempty(figoff),          figoff = 300; end
if ~exist('peaks','var'),    peaks = []; end
if ~exist('reverse_freq','var'), reverse_freq = []; end
if isempty(reverse_freq),    reverse_freq = false; end


%% loading, if filename given
if isempty(spec)
    [spec,dname] = uigetfile('*.mat','Select reconstructed MRSI mat file');
    if isequal(spec,0), error('No file selected'); end
    spec = [dname spec];
end
if ~isnumeric(spec)
    if isempty(regexpi(spec,'\.mat$','once')), spec = [spec '.mat']; end
    if exist(spec,'file')
        fnin = spec;
        xx = load(spec);
        if isfield(xx,'spec')
            spec = xx.spec;
        else
            error('no field spec');
        end
        if isfield(xx,'h')
            h = xx.h;
        else
            error('no field h');
        end
        if isfield(xx,'lb') && isempty(lb)
            lb = xx.lb;
        end
    else
        error('wrong input: spec');
    end
end


%% determine whether data is MRS or MRSI data
is_mrsi = true;    % by default data treated as MRSI
if size(spec,1)==h.rdb_hdr.user4
    is_mrsi = false;
    spec = permute(spec,[2 1 3]);
end



%% create output filename
if islogical(fname)
    if fname
        if exist('fnin','var')
            fname = fnin;
        else
            warning('no input filename (as spec) fnin; skipping saving');
            fname = '';
        end
    else
        fname = '';
    end
end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
    if isempty(regexpi(fname,'_oxsa$','once')), fname = [fname '_oxsa']; end
    fprintf('Saving results to ''%s''\n',fname);
end


%% invert frequency axis for recon_mrsi/csi prior to 5/2024
if reverse_freq
    fprintf('Inverting frequency axis\n');
    spec = ifftshift(fft(conj(ifft(fftshift(spec,1),[],1)),[],1),1);
end


%% averaging, if acquired with multiple time steps
if size(spec,5)>1
    fprintf('Averaging along dim 5\n');
    spec = mean(spec,5);
end
if length(size(spec))>4, error('length(size(spec))(=%g)>4',length(size(spec))); end
[ns,nx,ny,nz] = size(spec);


%% extract oxsa parameter structure
params = sub_header2oxsa(h,ns,time_delay);
nucleus = params.nucleus;


%% calculate mask for MRSI data
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
if is_mrsi
    bb = mrsi2pseuimage(spec,0);

    figure(figoff+5);
    if plt(3)==1
        imagesc_row(bb);
    else
        imagesc_ind3d(bb);
    end
    set(figoff+5,'name',[figstr ' - pseudo image']);
    if ~isempty(fname), print('-dpng','-r300',[fname '_pseudoimage.png']); end

    mask = calc_mask(bb,2,[],[],[],[],[10 90]);

    figure(figoff+6);
    if plt(3)==1
        imagesc_row(mask);
    else
        imagesc_ind3d(mask);
    end
    set(figoff+6,'name',[figstr ' - mask']);
    if ~isempty(fname), print('-dpng','-r300',[fname '_mask.png']); end
    pause(2);
else
    mask = true(1,nx,ny,nz);
end
nmask = sum(mask(:));
specvec = spec(:,mask);


%% add oxsa to path
if ~isdeployed
    % if ~contains(lower(path),'oxsa')    % not existing in V8.1
    if isempty(regexpi(path,'oxsa', 'once'))
        fprintf('Attention: adding oxsa to matlab path\n');
        addpath_oxsa;
    end
end


%% options
switch nucleus
    case 2,  xzm = [-10 15];
    case 13, xzm = [-20 20];
    case 31, xzm = [-25 15];
    otherwise, xzm = [min(params.ppmAxis) max(params.ppmAxis)];
end


%% prior-knowledge files
if isempty(pk_fname)
    % call directly --> PK* included in compiler
    switch nucleus
        case 2,   pk = PK_3T_2H_5Peaks;
        case 13,  pk = PK_3T_1C13_pyr;
        case 129, pk = PK_3T_129Xe_3Peaks;
        otherwise
            error('No pk file for nucleus=''%s''',num2str(nucleus))
    end
else
    try
        pk = feval(pk_fname);  %#ok<FVAL>
    catch ME
        warning('''pk = feval(pk_fname);'' failed; trying ''load(pk_fname,''pk'')''');
        fprintf('Error message: %s\n',ME.message);
        load(pk_fname,'pk')
    end
end
% add ppm_off
if abs(params.ppm_off)>0.01
    for ll=1:length(pk.bounds)
        pk.bounds(1,ll).chemShift = pk.bounds(1,ll).chemShift + params.ppm_off;
        pk.initialValues(1,ll).chemShift = pk.initialValues(1,ll).chemShift + ...
            params.ppm_off;
    end
end
% select peaks
if ~isempty(peaks)
    pk.bounds = pk.bounds(peaks);
    pk.initialValues = pk.initialValues(peaks);
    pk.priorKnowledge = pk.priorKnowledge(peaks);
end
[params.met_names,params.met_freqs,params.peak_names] = sub_pk2meta(pk);


%% determine offset from peak maximum
% assumes that reference peak has highest signal
if isempty(offset)
    s1d = sqrt(sum(specvec.*conj(specvec),2));
    [~,ii] = max(s1d);
    offset = -params.ppmAxis(ii);
end


%% shift spectra
if abs(offset)>1d-4
    fid = ifft(fftshift(specvec,1),[],1);
    phafu = exp(1i*2*pi*params.timeAxis.'*offset*params.imagingFrequency);
    fid = bsxfun(@times,fid,phafu);
    specvec = ifftshift(fft(fid,[],1),1);
end


%% plot rms spec
sub_plot_rms_spec(specvec,params,xzm,figoff+7,[figstr ' - RMS spectrum']);
if ~isempty(fname), print('-dpng','-r300',[fname '_rms_spec.png']); end


%% set experimental parameters
switch plt(1)
    case {0,1}     % off,final
     showPlot = false(1,nmask);
    case 2         % all
        showPlot = true(1,nmask);
    case 3         % random-15% 
        if nmask>100
            showPlot = (randn(1,nmask)>1);
        else
            showPlot = true(1,nmask);
        end
    otherwise, error('plt(1)(=%g) ~=0,1,2,3',plt(1));
end
npeaks = length(params.peak_names);
resnames = {'amplitude','chemShift','linewidth','crbamp'};
nres   = length(resnames);
resvec = zeros(nmask,npeaks,nres);


%% transform back spectra to FID
fid = double(ifft(fftshift(specvec,1),[],1));
if ~isempty(lb)
    if abs(lb(1))>0.01
        fprintf('De-apodising with lb=%g[Hz]\n',lb(1));
        fil = lb_fun(params.samples,params.bandwidth,-lb(1)).';
        fid = bsxfun(@times,fid,fil);
    end
end


%% actual amares fitting
for ll=1:nmask
    fstr = sprintf('ll=%d/%d',ll,nmask);
    fprintf('%s\n',fstr);
    
    [fitResults, fitStatus, ~, CRBResults] = ...
        AMARES.amaresFit(fid(:,ll),params,pk,0);

    if showPlot(1,ll)
        [spectrum,fittedSpectrum,separateSpectra] = ...
            sub_amares2spec(fitStatus);
        spec_plt = [spectrum,fittedSpectrum,separateSpectra];
        sub_plot_spec_result(spec_plt,params,xzm,figoff+8,figstr,fstr)
    end

    for lr=1:nres
        if lr<4
            resvec(ll,:,lr) = fitResults.(resnames{lr});
        else
            resvec(ll,:,lr) = CRBResults.amplitude;
        end
    end
end


%% reshape vector into 3D spatial
fitres = zeros(nx,ny,nz,npeaks,nres);
fitres(repmat(mask,[1 1 1 npeaks nres])) = resvec;



%% plotting
if ((plt(1)>0) && (any([nx,ny,nz]>1)))
    for lr=1:nres
        figure(figoff+lr); clf;
        cm = [];
        if lr==4
            tmp = fitres(:,:,:,:,lr); 
            cm = [0 0.1*mean(tmp(tmp>0))];
        end
        if plt(3)==1
            imagesc_row(fitres(:,:,:,:,lr),cm);
        else
            imagesc_ind3d(fitres(:,:,:,:,lr),[],cm);
        end
        colorbar
        str = [figstr ' - ' resnames{lr} ': '];
        for ll=1:npeaks, str = sprintf('%s-%s',str,params.peak_names{ll}); end
        set(figoff+lr,'name',str);
        drawnow;
        if ~isempty(fname), print('-dpng','-r300',[fname '_' resnames{lr} '.png']); end
    end
end


%% print info
resstr = sprintf('OXSA Amares fit: \nnucleus=%g\n',nucleus);
resstr = sprintf('%sns=%g; nx=%d; ny=%g; nz=%g; nmask=%d\n',resstr,ns,nx,ny,nz,nmask);
resstr = sprintf('%soffset=%g[ppm]; lb=%g; time_delay=%g\n',...
    resstr,offset,lb,time_delay);
str = '          ';
for lp=1:npeaks, str = sprintf('%s\t%s',str,params.peak_names{lp}); end
resstr = sprintf('%s%s\n',resstr,str);
for lr=1:nres
    str = '';
    for lp=1:npeaks
        str = sprintf('%s\t%.3g+-%.3g',str,...
            mean(resvec(:,lp,lr)),std(resvec(:,lp,lr)));
    end
    resstr = sprintf('%s%s %s\n',resstr,resnames{lr},str);
end


fprintf('%s',resstr);
if ~isempty(fname)
    fileid = fopen([fname '.txt'],'wt');
    fprintf(fileid,resstr);
    fclose(fileid);
end


%% saving
if ~isempty(fname)
    save([fname '.mat'],'fitres','resnames','params','pk','resstr');
end


%% dicom export
if ~isempty(fname) && plt(2)>0
    inp.MRAcquisitionType = '3';
    inp.SeriesNumber = h.series.se_no + figoff;
    write_dicom(fname,fitres(:,:,:,:,1),h,inp);
end

end      % fit_oxsa.m


%% sub-functions
%% plot (spatial) rms spectrum with initial values of peaks positions
function sub_plot_rms_spec(spec,params,xzm,figid,figstr)

% calculate spatial RMS spectrum
s1d = sqrt(sum(sum(sum(spec.*conj(spec),4),3),2));
s1d = s1d/max(s1d);                    % normalise 
lgnd = {'RMS spectrum',params.met_names{1:end}};
ppm = params.ppmAxis - params.ppm_off;
fmet = repmat(params.met_freqs,[2 1]) - params.ppm_off;

% actual plotting
figure(figid); clf
set(figid,'name',figstr);
plot(ppm,s1d,'k',...
    fmet,[0 1],':','LineWidth',1.5)
axis([xzm min(s1d) max(s1d)]);
set(gca,'xdir','reverse');
xlabel('freq [ppm]');
% set(gca,'ytick','');
grid on
legend(lgnd);
drawnow

end      % sub_plot_rms_spec


%% plot results
function sub_plot_spec_result(spec_plt,params,xzm,figid,figstr,tstr)

spec_plt = real(spec_plt);
spec_res = spec_plt(:,1) - spec_plt(:,2);
axs = [xzm min(spec_plt(:)) max(spec_plt(:))];
ppm = params.ppmAxis - params.ppm_off;
lgnd = {'orig','fit','residual',params.met_names{1:end}};

% actual plotting
figure(figid); clf
set(figid,'name',figstr);
plot(ppm,spec_plt(:,1),'k',ppm,spec_plt(:,2),'g:',...
    ppm,spec_res,'r:',ppm,spec_plt(:,3:end),'',...
    'LineWidth',1.5)
legend(lgnd);
title(tstr)
axis(axs)
set(gca,'xdir','reverse');
grid on
drawnow

end      % sub_plot_spec


%% convert GE p-file header to oxsa amares parameters structure
function params = sub_header2oxsa(h,ns,time_delay)

% get parameters
nucleus = h.image.specnuc;
bw = h.rdb_hdr.spectral_width;
switch nucleus
    case {1,2}, ppm_off = -4.7;
    otherwise,  ppm_off = 0;
end
f0 = h.rdb_hdr.ps_mps_freq/10;
hz = (-ns/2:ns/2-1)/ns*bw;
ppm = hz(:)/(f0*1d-6);
t = (0:ns-1)/bw;

B0 = 2*pi*f0/gyrogamma(nucleus);
B0 = round(B0*10)/10;

% echo time -> linear phase
if ~isempty(time_delay)
    beginTime = time_delay;
else
    % if contains(lower(h.image.psd_iname),'fidall')
    if ~isempty(regexpi(h.image.psd_iname,'fidall', 'once'))
        beginTime = h.rdb_hdr.te*1d-6;
    end
    % if contains(lower(h.image.psd_iname),'fidcsi')
    if ~isempty(regexpi(h.image.psd_iname,'fidcsi', 'once'))
        beginTime = 0.5d-3;
    end
end

% set OXSA structure
params.samples = ns;                             % #samples
params.bandwidth = bw;                           % spectral BW      [Hz]
params.dwellTime = 1/bw;                         % dwellTime        [s]
params.ppmAxis = ppm;                            % ppmAxis          [ppm]
params.timeAxis = t;                             % timeAxis         [s]
params.imagingFrequency = f0*1d-6;               % centre frequency [MHz]
params.ppm_off = ppm_off;                        % f0 offset        [ppm]
params.B0 = B0;                                  % field strength   [T]
params.beginTime = beginTime;                    % linear phase     [s]
params.nucleus = nucleus;                        % nucleus

end      % sub_header2oxsa


%% convert pk naming to usable names+freqs
function [met_names,met_freqs,peak_names] = sub_pk2meta(pk)


nmet = length(pk.initialValues);       % #metabolites
met_freqs = zeros(1,nmet);             % freqs of metabolites in [ppm]

lp = 0;
for lm=1:nmet
    met_freqs(1,lm) = pk.initialValues(lm).chemShift;
    peakName = pk.initialValues(lm).peakName;
    if iscell(peakName)
        for l2=1:length(peakName)
            lp = lp+1;
            peak_names{lp} = peakName{l2};
        end
        met_names{lm} = strrep(peakName{1}(1:end-1),'_',' ');
    else
        lp = lp+1;
        peak_names{lp} = peakName;
        met_names{lm} = peakName;
    end
end

end      % sub_pk2met


%% calculate resulting spectra
function [spectrum,fittedSpectrum,separateSpectra] = sub_amares2spec(varargin)

% Process options
options = processVarargin(varargin{:});

% Required "options"
exptParams = options.exptParams;
inputFid = options.inputFid;
xFit = options.xFit;
constraintsCellArray = options.constraintsCellArray;

% Defaults
if ~isfield(options,'apodization')
    options.apodization = 30;
end

% Apply prior knowledge to yield full solution
fitResults = AMARES.applyModelConstraints(xFit,constraintsCellArray);

% Calculate the zero and first order phase correction
if ~isfield(options,'firstOrder') || options.firstOrder % On by default
    if ~isfield(exptParams,'freqAxis')
        exptParams.freqAxis = exptParams.ppmAxis * exptParams.imagingFrequency;
    end
    % Set zero-order phase relative to the reference peak
    actualRefPeak = AMARES.getActualRefPeakDx(options.pkWithLinLsq);
    if ~isempty(actualRefPeak)
        zeroOrderPhaseRad = fitResults.phase(actualRefPeak) * pi / 180;
    else
        zeroOrderPhaseRad = 0;
    end
    fprintf('zeroOrderPhaseRad = %.3f\n',zeroOrderPhaseRad)
    hz_axis = (-exptParams.samples/2:exptParams.samples/2-1).' / ...
        exptParams.dwellTime/exptParams.samples;
    firstOrderCorrection = exp(-1i*(zeroOrderPhaseRad + 2*pi*hz_axis*exptParams.beginTime));
else
    firstOrderCorrection = 1;
end


% Calculate the fixed spectra to plot
timeAxis = exptParams.dwellTime*(0:exptParams.samples-1).';
spectrum = specApodize(timeAxis, specFft(inputFid,1).* ...
    firstOrderCorrection,options.apodization);
[modelFid,~,modelFids] = AMARES.makeModelFidAndJacobianReIm(xFit,...
    constraintsCellArray,exptParams.beginTime,exptParams.dwellTime,...
    exptParams.imagingFrequency,exptParams.samples, 'complexOutput', true);
fittedSpectrum =  specApodize(timeAxis, specFft(modelFid).* ...
    firstOrderCorrection,options.apodization);

% Individual peaks
if true    
    nMetab = length(options.pkWithLinLsq.bounds);    
    peakDx = 0;
    for lMetab = 1:nMetab
        if iscell(options.pkWithLinLsq.bounds(lMetab).peakName)
            indivFID = modelFids(:,peakDx+1);
            for l=2:length(options.pkWithLinLsq.bounds(lMetab))
                indivFID = indivFID + modelFids(:,peakDx+l);
            end
            peakDx = peakDx+length(options.pkWithLinLsq.bounds(lMetab).peakName);
        else
            peakDx = peakDx+1;
            indivFID = modelFids(:,peakDx);
        end
        indivSpectrum = specApodize(timeAxis,specFft(indivFID).* ...
            firstOrderCorrection,options.apodization);
        separateSpectra(:,lMetab) = indivSpectrum; % MV 2023-01-03
    end 
end


end      % sub_amares2spec
