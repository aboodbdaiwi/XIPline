function [snr,lw,sz,dur,mdcf,psf_vol] = analyse_trajectory(wfn,r2s,lb,...
    zf_fac,verb,fwhi,fov)
%ANALYSE_TRAJECTORY  Determine real resolution and SNR of given trajectory
%[snr,lw,sz,dur,mdcf,psf_vol] = analyse_trajectory(wfn,r2s,lb,zf_fac,...
%                                 verb,fwhi,fov)
%    wfn   k-space waveform name (*.mat)
%    r2s   Decay rate due to T2*                           [Hz]        (15)
%     lb   Apodisation applied in reconstruction           [Hz,spat,Hz]([])
% zf_fac   Zero-filling factor                                          (2)
%   verb   Verbose: 0=off, 1=verbose, 2=verbose+plotting                (2)
%   fwhi   Full-width at fwhi-max: 0.64=Raleigh criterion in optics  (0.64)
%    fov   Field of view (default from wf.fov)             [mm]
%
%    snr   Relative to optimally sampled SNR               [%]
%     lw   Linewidth factor from nominal to real size
%     sz   Nominal/real pixel/voxel sizes (nom/real,dim)   [mm]
%    dur   Trajectory duration                             [s]
%   mdcf   mean(dcf)
%psf_vol   PSF volume (normalised) (use zf_fac>=2)
%
%  3/2023 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('r2s','var'),    r2s = []; end
% if isempty(r2s),           r2s = 15; end
if ~exist('lb','var'),     lb = []; end
if isempty(lb),            lb = 0; end
if length(lb)<2,           lb = [lb 0]; end
if length(lb)<3,           lb = [lb 0]; end
if ~exist('zf_fac','var'), zf_fac = []; end
if isempty(zf_fac),        zf_fac = 2; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end
if ~exist('fwhi','var'),   fwhi = []; end
if isempty(fwhi),          fwhi = 0.64; end
if fwhi>0.999, error('fwhi(=%g)>1',fwhi); end
if fwhi<0.001, error('fwhi(=%g)<0',fwhi); end
if ~exist('fov','var'),    fov = []; end
zf = 1024;         % zf for 1D cross-section plot
nnoiseall = 5d4;
do_vol = true;     % calculate volume-based linewidth factor


%% reading waveform file
wf = load_waveform(wfn);
dim = size(wf.k,3);
mtx = wf.mtx;
nexc = wf.nviews;
npts = max(sum(wf.ind,2));
nk = size(wf.k,2);
t = wf.t;
ksp = wf.k;
dcf = wf.dcf;


%% parameters
if isempty(fov)
    fov = wf.fov(1)*1d3;    % FOV [mm]
end
% dim = size(ksp,3);
switch dim
    case 2
        if length(mtx)==1, mtx = [mtx mtx]; end
        if length(fov)==1, mtx = [fov fov]; end
    case 3
        if length(mtx)==1
            mtx = [mtx mtx mtx]; 
        elseif length(mtx)==2
            mtx = [mtx 1];
        end
        if length(fov)==1
            fov = [fov fov fov]; 
        elseif length(fov)==2
            fov = [fov 1];
        end
    otherwise, error('dim(=%d)~=2 or 3',dim); 
end
mtx = mtx(:).';
if nexc*npts~=nk
    if isfield(wf,'scheme')
        fprintf('MRSI acquisition\n');
        nexc = nk;
        npts = 1;
        dcf = ones(size(dcf));
        fprintf('Setting nexc=nk=%d; npts=1; dcf=1s\n',nk);
    else
        warning('nexc(=%d)*npts(=%d)~=nk(=%d)',nexc,npts,nk);
        npts = round(nk/nexc);
        fprintf('Setting npts=nk/nexc=%d\n',npts);
    end
end
if any(diff(mtx)~=0), warning('not isotropic'); end
nnoise = ceil(nnoiseall/prod(mtx));
fprintf('dim = %g; nexc = %d; npts = %d; nk = %d\n',dim,nexc,npts,nk);


%% calculate point-spread-function
if verb>0, fprintf('Calculating point spread function\n'); end
if ~isempty(r2s)
    if verb>0 
        fprintf('Assumed R2* = %g [Hz]; T2* = %g [ms]\n',...
            r2s,1d3/(pi*r2s));
    end
    dt = round((t(1,2)-t(1,1))*1d6);
    bw = 1d6/dt;
    if isinf(bw), warning('isinf(bw)'); end
    if bw>5d5, warning('bw(=%g[Hz])>500[kHz]',bw); end
    if bw<4999, warning('bw(=%g[Hz])<5[kHz]',bw); end
    fprintf('Acq. BW = %g [Hz]\n',bw);
    
    if any(diff(sum(wf.ind,2))~=0)
        dd = [];
        for l1=1:size(wf.ind,1)
            n2 = sum(wf.ind(l1,:));
            decay = lb_fun(n2,bw,[0 r2s],-t(1,1));
            dd = [dd decay];
        end
    else
        decay = lb_fun(npts,bw,[0 r2s],-t(1,1));
        dd = repmat(decay,[nexc 1]).';
        dd = dd(:).';
    end
else
    dd = ones(1,nk);
end
if nk~=size(dd,2)
    warning('nk(=%d)~=size(dd,2)(=%d)',nk,size(dd,2));
    % nk = size(dd,2);
end
dn = randn(nnoise,nk)+1i*randn(nnoise,nk);
dn = dn/sqrt(2);                       % scale to std=1
dd = [dd ; dn];
if any(abs(lb)>0)
    if verb>0
        fprintf('Apodisation in reconstruction = %g %g %g [Hz,spat,Hz]\n',...
            lb(1),lb(2),lb(3));
    end
    
    % dd = ak_apodise(dd,[],[],t,lbt,false);
    dd = ak_apodise(dd,ksp,lb(2),t,[lb(1) lb(3)],false);
end
% dd = [dd ; dn];

if zf_fac>0
    ksp = ksp/zf_fac; 
    mtx_reco = round(mtx*zf_fac);
else
    mtx_reco = mtx;
end

switch dim
    case 2
        %bb = gridding(dd,ksp,dcf,mtx_reco); % windows mex not shared
        ksp_3d = cat(3,ksp, zeros(1,size(ksp,2),1));
        bb = gridding3dsingle(ksp_3d,dd,dcf,[mtx_reco, 1]);
        b1{1} = bb(:,ceil(mtx_reco(2)/2+1),1);
        b1{2} = bb(ceil(mtx_reco(1)/2+1),:,1).';
        tmp1 = bb(:,:,1+(1:nnoise));
        % tmp2 = bb(:,:,1+nnoise+(1:nnoise));
    case 3
        if size(dcf,2)~=size(ksp,2) % account for dcf only calculated for 1 
            factor = round(size(ksp,2)/size(dcf,2));
            dcf = repmat(dcf,1,factor,1);
        end
        bb = gridding3dsingle(ksp,dd,dcf,mtx_reco);
        b1{1} = bb(:,ceil(mtx_reco(2)/2+1),ceil(mtx_reco(3)/2+1),1);
        b1{2} = bb(ceil(mtx_reco(1)/2+1),:,ceil(mtx_reco(3)/2+1),1).';
        b1{3} = squeeze(bb(ceil(mtx_reco(1)/2+1),ceil(mtx_reco(2)/2+1),:,1));
        tmp1 = bb(:,:,:,1+(1:nnoise));        
end
bn = tmp1(:).';
for l=1:dim
    if ~isempty(zf)
        d1 = ifftshift(truma(fftshift(ifft(fftshift(b1{l}))),false,[zf 1]));
        b1{l} = abs(ifftshift(fft(d1)));
    else
        b1{l} = abs(b1{l});
    end
end


%% determine fwhm
lw = zeros(1,dim);
sz = zeros(2,dim);
for l=1:dim
    [lw(1,l),out{l}] = fwhm(b1{l}.',mtx(l),[],false,false,fwhi);
    sz(1,l) = fov(l)/mtx(1,l);            % nominal pixel/voxel size
    sz(2,l) = fov(l)/mtx(1,l)*lw(1,l);    % actual size
end
mdcf = mean(dcf(:));


%% determine PSF factor in 2D/3D
if do_vol
    if zf_fac<2
        warning('zf_fac(=%g)<2: psf_vol might be imprecise',zf_fac);
    end
    if dim==2
        btmp = bb(:,:,1);
    else
        btmp = bb(:,:,:,1);
    end
    btmp = btmp/max(abs(btmp(:)));
    mask = btmp>0.005; btmp = btmp.*mask;
    psf_vol = abs(sum(btmp(:)))/(zf_fac^dim);
    % psf_vol = abs(sum(btmp(:)));
    % fprintf('-- %g %g\n',psf_vol,sum(abs(btmp(:)))/(zf_fac^dim));
    clear btmp
end


%% plotting 
if verb>1
    figure(51); clf; 
    switch dim
        case 2, imagesc(abs(bb(:,:,1)));
        case 3, imagesc_ind3d(abs(bb(:,:,:,1)));
    end
    
    figure(52); clf; set(gcf,'name','log');
    switch dim
        case 2, imagesc(log(abs(bb(:,:,1))));
        case 3, imagesc_ind3d(log(abs(bb(:,:,:,1))));
    end
    colorbar;
    imax = log(max(abs(b1{1})));
    caxis([-8 0]+imax);

    figure(53); clf; clr = {'b','r','g'};
    for l=1:dim
        x_ax = (-zf/2:zf/2-1)/zf*mtx(l);
        maxval = b1{l}(out{l}.ind_max,1);
        plot(-x_ax,b1{l},clr{l},...
             out{l}.freqs,maxval,[clr{l} '*'],...
             out{l}.fp,maxval*fwhi,[clr{l} '*'],...
             out{l}.fm,maxval*fwhi,[clr{l} '*']);
        hold on
    end
    hold off
    grid on; axis tight;
end


%% calculate real resolution
if verb>0
    lwstr = sprintf('%.3g %.3g',lw(1,1),lw(1,2));
    for l=1:2
        resistr{l} = sprintf('%.3g %.3g',sz(l,1),sz(l,2));
    end
    nosistr = sprintf('%.3g %.3g',fov/mtx(1,1),fov/mtx(1,2));
    switch dim
        case 2, vstr = 'pixel';            
        case 3
            vstr = 'voxel';
            lwstr = sprintf('%s %.3g',lwstr,lw(1,3));
            for l=1:2
                resistr{l} = sprintf('%s %.3g',resistr{l},sz(l,3));
            end
    end
    fprintf('FOV = %g [mm]\n',fov);
    fprintf('Nominal matrix size = %s\n',num2str(mtx));
    fprintf('Linewidth factor   = %s (full width at %g%% max)\n',...
        lwstr,fwhi*100);
    if do_vol, fprintf('PSF volume         = %.3g\n',psf_vol); end
    fprintf('Nominal %s size = %s [mm^%d]\n',vstr,resistr{1},dim);
    fprintf('Real %s size    = %s [mm^%d]\n',vstr,resistr{2},dim);
    fprintf('mean(dcf)          = %.3g\n',mdcf);
end


%% calculate theoretical SNR of FID
if verb>2
    if abs(lb(2))>0, warning('not considering spatial apodisation'); end
    calc_snr_fid(wf.t(end),r2s,[lb(3) lb(1)],wf.t(1));
end


%% calculate real SNR
signal = max(abs(b1{1}));
noise = std(bn,'',2);
noise = noise*sqrt(nk);                % scale to std=1
snr = signal/noise(1)*100;             % apodised signal-to-noise ratio
if verb>0, fprintf('SNR = %.3g/%.3g = %.3g [%%]\n',signal,noise(1),snr); end
dur = max(t(:));


end      % analyse_trajectory.m

