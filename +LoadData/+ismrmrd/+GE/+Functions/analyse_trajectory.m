function [snr,lw,sz,dur] = analyse_trajectory(wfn,r2s,lbt,zf_fac,verb,fwhi)
%ANALYSE_TRAJECTORY  Determine real resolution and SNR of given trajectory
%[snr,sz,dur] = analyse_trajectory(wfn,r2s,lbt,zf_fac,verb,fwhi)
%    wfn   k-space waveform name (*.mat)
%    r2s   Decay rate due to T2*                           [Hz]      (15)
%    lbt   Apodisation applied in reconstruction           [Hz]      ([])
% zf_fac   Zero-filling factor                                       (2)
%   verb   Verbose: 0=off, 1=verbose, 2=verbose+plotting             (2)
%   fwhi   Full-width at fwhi-max: 0.64=Raleigh criterion in optics  (0.64)
%    snr   Relative to optimally sampled SNR               [%]
%     lw   Linewidth factor from nominal to real size
%     sz   Nominal/real pixel/voxel sizes (nom/real,dim)   [mm]
%    dur   Trajectory duration                             [s]
%
% 12/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('r2s','var'),    r2s = []; end
if isempty(r2s),           r2s = 15; end
if ~exist('lbt','var'),    lbt = []; end
if ~exist('zf_fac','var'), zf_fac = []; end
if isempty(zf_fac),        zf_fac = 2; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end
if ~exist('fwhi','var'),   fwhi = []; end
if isempty(fwhi),          fwhi = 0.64; end
if fwhi>0.999, error('fwhi(=%g)>1',fwhi); end
if fwhi<0.001, error('fwhi(=%g)<0',fwhi); end
zf = 1024;         % zf for 1D cross-section plot
nnoiseall = 5d4;


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
if isfield(wf,'mtx')
    mtx = wf.mtx; 
else
    mtx = wf.npix;
end
if isfield(wf,'k')
    dim = size(wf.k,3);
else
    if isfield(wf,'ks')
        dim = size(wf.ks,3);
    else
        error('field wf.k and wf.ks not found');
    end
end

switch dim
    case 1
        if ~isreal(wf.k)
            ksp(1,:,1) = real(wf.k);
            ksp(1,:,2) = imag(wf.k);
        else
            error('(size(ksp,3)==1) && isreal(ksp)');
        end
        dcf = wf.dcf;
        t = wf.t.'; t = t(:).';
        [nexc,npts] = size(wf.t);
        dim = 2;
    case 2
        ksp = wf.k;
        dcf = wf.dcf;
        t = wf.t.'; t = t(:).';
        [nexc,npts] = size(wf.t);
    case 3
        [ksp,dcf,t,nexc,npts] = sub_wf2traj(wf);
    otherwise, error('dim(=%d)~=2 or 3',dim);
end
fov = wf.fov(1)*1d3;    % FOV [mm]
nk = size(ksp,2);


%% parameters
% dim = size(ksp,3);
switch dim
    case 2, if length(mtx)==1, mtx = [mtx mtx]; end
    case 3, if length(mtx)==1, mtx = [mtx mtx mtx]; end
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
    % bw = 1/(wf.t(1,2)-wf.t(1,1));
    bw = 1/(t(1,2)-t(1,1));
    if isinf(bw), warning('isinf(bw)'); end
    if bw>5d5, warning('bw(=%g[Hz])>500[kHz]',bw); end
    if bw<4999, warning('bw(=%g[Hz])<5[kHz]',bw); end
    fprintf('Acq. BW = %g [Hz]\n',bw);
    
    decay = lb_fun(npts,bw,[0 r2s],-t(1,1));
    dd = repmat(decay,[nexc 1]).';
    dd = dd(:).';
else
    dd = ones(1,nk);
end
dn = randn(nnoise,nk)+1i*randn(nnoise,nk);
dn = dn/sqrt(2);                       % scale to std=1
dd = [dd ; dn];
if ~isempty(lbt)
    if verb>0
        fprintf('Apodisation in reconstruction = %g [Hz]\n',lbt);
    end
    
    dd = ak_apodise(dd,[],[],t,lbt,false);
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
        bb = gridding(dd,ksp,dcf,mtx_reco);
        b1{1} = bb(:,ceil(mtx_reco(2)/2+1),1);
        b1{2} = bb(ceil(mtx_reco(1)/2+1),:,1).';
        tmp1 = bb(:,:,1+(1:nnoise));
        % tmp2 = bb(:,:,1+nnoise+(1:nnoise));
    case 3
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
    sz(1,l) = fov/mtx(1,l);            % nominal pixel/voxel size
    sz(2,l) = fov/mtx(1,l)*lw(1,l);    % actual size
end


%% plotting 
if verb>1
    figure(51); clf;
    switch dim
        case 2, imagesc(abs(bb(:,:,1)));
        case 3, imagesc_ind3d(abs(bb(:,:,:,1)));
    end
    
    figure(52); clf;
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
             out{l}.fp,maxval*0.64,[clr{l} '*'],...
             out{l}.fm,maxval*0.64,[clr{l} '*']);
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
    fprintf('Linewidth factor   = %s (full width at %g%% max)\n',lwstr,fwhi*100);
    fprintf('Nominal %s size = %s [mm^%d]\n',vstr,resistr{1},dim);
    fprintf('Real %s size    = %s [mm^%d]\n',vstr,resistr{2},dim);
end


%% calculate theoretical SNR of FID
if verb>2
    calc_snr_fid(wf.t(end),r2s,[0 lbt],wf.t(1));
end


%% calculate real SNR
signal = max(abs(b1{1}));
noise = std(bn,'',2);
noise = noise*sqrt(nk);                % scale to std=1
snr = signal/noise(1)*100;             % apodised signal-to-noise ratio
if verb>0, fprintf('SNR = %.3g/%.3g = %.3g [%%]\n',signal,noise(1),snr); end
dur = max(t(:));


end      % main function analise_trajectory.m


%% sub-functions
% generate 3D trajectory from reduced mat file
function [k,dcf,t,nexc,npts] = sub_wf2traj(wf)
do_single = true;
if isfield(wf,'nviews')
    % reduced mat size: only single view stored; regenerate info
    if ~isfield(wf,'k')
        if ~isfield(wf,'ks') && ~isfield(wf,'phi') && ~isfield(wf,'theta')
            error('wf.ks,phi,theta not existing');
        else
            ks = squeeze(wf.ks);
            nks = size(ks,1);
            if do_single, ks = single(ks); end
            k = zeros(nks,wf.nviews,3,class(ks));
            
            sinphi = sind(wf.phi);
            cosphi = cosd(wf.phi);
            sintheta = sind(wf.theta);
            costheta = cosd(wf.theta);
            
            k(:,:,1) =  ks(:,1)*(costheta.*cosphi) + ...
                ks(:,2)*(sinphi) - ...
                ks(:,3)*(sintheta.*cosphi);
            k(:,:,2) =  -ks(:,1)*(costheta.*sinphi) + ...
                ks(:,2)*(cosphi) + ...
                ks(:,3)*(sintheta.*sinphi);
            k(:,:,3) = ks(:,1)*(sintheta) + ...
                ks(:,3)*(costheta);
            % k = permute(k,[2 1 3]);
            k = reshape(k,[1 nks*wf.nviews 3]);
        end
    else
        k = wf.k;
    end
    if ~isfield(wf,'dcf')
        if ~isfield(wf,'dcfs')
            wfn_dcf = [wfn(1:end-4) '_dcf.mat'];
            if exist(wfn_dcf,'file')
                if verb>0
                    fprintf('Loading dcf from file ''%s''\n',wfn_dcf);
                end
                load(wfn_dcf);
            else
                if verb>0
                    fprintf('Calculating dcf via calc_sdc3\n');
                    fprintf('\tand storing to ''%s''\n',wfn_dcf);
                end
                dcf = calc_sdc3(k,wf.mtx,5);
                save(wfn_dcf,'dcf');
            end
        else
            dcfs = wf.dcfs(:).';
            if do_single, dcfs = single(dcfs); end
            dcf = repmat(dcfs,[1 wf.nviews]);
        end
    else
        dcf = wf.dcf;
        if isa(dcf,'uint16')
            if do_single
                dcf = single(dcf)/(2^16-1)*wf.dcf_range;
            else
                dcf = double(dcf)/(2^16-1)*wf.dcf_range;
            end
        end
    end
    if ~isfield(wf,'t')
        if ~isfield(wf,'ts')
            error('wf.ts not existing'); 
        else
            ts = wf.ts(:).';
            if do_single, ts = single(ts); end
            t = repmat(ts,[1 wf.nviews]);
        end
        npts = size(wf.ts,2);
    else
        t = wf.t;
        t = t.'; t = t(:).';
        npts = size(wf.t,2);
    end
    nexc = wf.nviews;
else
    if ~isfield(wf,'k'),   error('wf.k not existing'); end
    if ~isfield(wf,'dcf'), error('wf.dcf not existing'); end
    if ~isfield(wf,'t'),   error('wf.t not existing'); end
    k = wf.k;
    dcf = wf.dcf;
    t = wf.t;
    if size(t,1)==1, t = repmat(t,[size(k,1) 1]); end
    t = t.'; t = t(:).';
    [nexc,npts] = size(wf.t);
end

end      % sub_wf2traj


