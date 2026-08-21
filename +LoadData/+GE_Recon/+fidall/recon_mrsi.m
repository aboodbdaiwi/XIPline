function [spec,hz,s1d,ppm,par,h,bb] = recon_mrsi(d,h,wfn,zf,lb,fname,plt,...
    coco,ftmethod,ws,do_ecc,d_ref,freq,pseuim,nviews0)
%RECON_MRSI  Reconstruct MRSI data acquired with arbitrary spatial k-space
%[spec,hz,s1d,ppm,par,h,bb] = recon_mrsi(d,h,wfn,zf,lb,fname,plt,...
%                                  coco,slowft,ws,do_ecc,d_ref,freq,pseuim)
%       d   Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc] 
%           (ns=sampling pts; nc=coils)
%           Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%       h   Header structure
%     wfn   Location of waveform .mat file
%      zf   Zero fill to: [#s,#x,#y,#z]; if<0->multiplier
%      lb   Spectral line broadening [Hz]
%   fname   Print <fname>.png and save reco as <fname>.mat             ([])
%           also export max-image and plots as dicom
%     plt   Plotting: plt(1) 0=off, 1=normal, 2+=plot_mrsi              (2)
%           plt(2): save _plot_csi_sl1_ts1.png:              (3D->0, 2D->1)
%    coco   Coil combination [do percent_noise4decor]               ([1 0])
%ftmethod   Spatial Fourier transformation method                       (0)
%           0=gridding, 1=slowFT, 2=pinv
%      ws   Water suppression: 0=off; 1=HSVD; 2=subtraction             (0)
%  do_ecc   Phase deconvolution eddy current correction             (false)
%           requires H2O reference scan: 
%           if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%           if ~isempty(d_ref) -> default=true; 
%   d_ref   Reference scan (w/o water suppression)
%    freq   Reconstruct spectral domain with matrix inversion          ([])
%  pseuim   How to generate pseudo image                            ([1 4])
%           mrsi2pseuimage.m -> [how2ind npts]                  (2H->[2 4])
% nviews0   #initial view to zero d(1:nviews,:,:,:,:,:)=0              ([])
%
%    spec   Spectrum        [#s,#x,#y,#z,#t]
%      hz   Frequency axis  [Hz]
%     s1d   Pseudo spectra  [#s,3]
%     ppm   Frequency axis  [ppm]
%     par   Parameter structure for plot_mrsi
%
% 11/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;              % record execution time


%% misc input
if ~exist('wfn','var'),    wfn = ''; end
if ~exist('zf','var'),     zf = []; end
if ~exist('lb','var'),     lb = []; end
if isempty(lb),            lb = 0; end
if ~exist('fname','var'),  fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
    end
end
if ~exist('plt','var'),    plt = []; end
if isempty(plt),           plt = 2; end
if ~exist('coco','var'),   coco = []; end
if isempty(coco),          coco = 1; end
if length(coco)<2,         coco = [coco 0]; end
if ~exist('ftmethod','var'), ftmethod = []; end
if isempty(ftmethod),      ftmethod = 0; end
if ~exist('ws','var'),     ws = []; end
if isempty(ws),            ws = 0;  end
if length(ws)~=1,          warning('length(ws)~=1'); end
if ~exist('do_ecc','var'), do_ecc = []; end
if ~exist('d_ref','var'),  d_ref = []; end
if ~exist('freq','var'),   freq = []; end
if ~isempty(freq), if ~isnumeric(freq), warning('~isnumeric(freq)'); end; end
if ~exist('pseuim','var'), pseuim = []; end
if ~exist('nviews0','var'), nviews0 = []; end
do_ori = true;          % correct image orientation


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
end
d = single(d);
d = conj(d);                 % changing dir of axis (keep MRS convention)


%% null nviews0 initial data frames (for bSSFP)
if ~isempty(nviews0)
    fprintf('Nulling %g initial data frames\n',nviews0);
    d(1:nviews0,:,:,:,:,:) = 0;
end


%% reading waveform file
if iscell(wfn), wfn = wfn{h.image.user3}; end
if isempty(wfn), wfn = get_wfname(h); end
% wf = load_waveform(wfn,true);
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
    if  isempty(regexpi(wfn,'\.mat$', 'once')), wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end


%% set default for do_ecc
if isempty(do_ecc)
    do_ecc = false;
    if h.rdb_hdr.user23==1
        if h.image.numecho==2
            do_ecc = true;
        else
            warning('h.rdb_hdr.user23==1; nt(=%g)~=2',h.image.numecho);
        end
    end
    if ~isempty(d_ref), do_ecc = true; end 
end


%% if separate reference -> add as second time step
if ~isempty(d_ref)
    if any(size(d)~=size(d_ref)), error('size(d)~=size(d_ref)'); end
    if size(d,4)~=1,   warning('size(d,4)(=%g) ~= 1',size(d,4)); end
    tmp = complex(zeros(size(d,1),size(d,2),size(d,3),2,size(d,5),...
        size(d,6),class(d)));
    tmp(:,:,:,1,:,:) = d;
    tmp(:,:,:,2,:,:) = d_ref;
    d = tmp;
    clear tmp d_ref
    nt = 2;
end


%% misc parameters
ws_bounds = [-100 50];                      % bounds for HSVD water suppression
si = size(d); si = [si ones(1,6-length(si))];
mtx = wf.mtx(:).'; mtx = [mtx ones(1,3-length(mtx))];
nx = mtx(1);
ny = mtx(2);
nz = mtx(3);
ns = h.rdb_hdr.user1;                       % #sampled points (FID w/o ind)
nslices = h.rdb_hdr.nslices;
if ~exist('nt','var'), nt = h.image.numecho; end  % number of time steps
nc = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1; % number of Rx coil elements
nref = h.rdb_hdr.user19;                    % #reference frames (w/o WS)
if si(2)==1   % no spectral dimension; only spatial reconstruction
    fprintf('size(d,2)==1: only spatial reconstruction\n');
    nspec = 1;
    ind = true;
else
    nspec = wf.npts;                        % number of spectral points
    ind = wf.ind;
end
if isempty(pseuim)
    if h.image.specnuc==2, pseuim = 2;
    else, pseuim = 1; end
end
if length(pseuim)<2,       pseuim = [pseuim 4]; end



%% misc checks
if si(2)~=ns, warning('si(2)(=%g)~=ns(=%g)',si(2),ns); end
if si(3)~=1, warning('si(3)(=%g) ~= 1',si(3)); end
if nt~=si(4)
    warning('nt(=%g) ~= size(d,4)(=%g); setting nt=%g',nt,si(4),si(4));
    nt = si(4);
end
if si(5)~=nslices, warning('si(5)(=%g) ~= nslices(=%g)',si(5),nslices); end
if nc~=si(6)
    warning('nc(=%g) ~= size(d,6)(=%g); setting nc=%g',nc,si(6),si(6));
    nc = si(6);
end
if ((nz>1) && (nslices>1))
    warning('nz(=%g)>1 && nslices(=%g)>1; setting nslices=1',nz,nslices);
    nslices = 1;
end
dim = wf.dim;
switch dim
    case {1,2}, nn = [nspec mtx(1) mtx(2) nslices nt nc];
    case 3,     nn = [nspec mtx(1) mtx(2) mtx(3) nt nc];
    otherwise,  error('dim(=%g)~=1,2 or 3',dim);
end
if length(plt)<2
    if dim==3, plt(2) = 0;
    else,      plt(2) = 1; end
end
bw = h.rdb_hdr.spectral_width;              % sampling bandwidth [Hz]
t  = (0:nspec-1).'/bw;                      % time [s]
itmp = find(diff(wf.ind));
te = h.rdb_hdr.te*1d-6+(itmp(1)+1)/bw;
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
fov = h.rdb_hdr.fov;
kmax = max(abs(wf.k));
kmax = kmax(:).';
zf_pseuim = zf(5:end);
zf = zf(1:min(4,length(zf)));
[zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices,kmax);
fprintf('zf = %s; do_zf_spat = %d\n',num2str(zf),do_zf_spat);
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);


%% unchop data
if ischop(h), d(2:2:end,:) = -d(2:2:end,:); end


%% interleaved spectra from centre of k-space
if isfield(wf,'ind1')
    ind1 = wf.ind1;
else
    if isfield(wf,'ncok')
        ncok = wf.ncok;
        if isempty(ncok), ncok = 0; end
        if ncok(1)>1
            ind1 = true(wf.nexc,1);
            ind1(1:ncok(1):end,1) = false;
        else
            fprintf('ncok(1)==%g -> ind1 = [];\n',ncok(1));
            ind1 = [];
        end
    else
        ind1 = []; 
    end
end
if ~isempty(ind1)
    d = sub_reco_cok(d,h,ind1,wf.nexc,zf,lb,fname,plt,figstr);
    si(1) = size(d,1);
end


%% reference scan for eddy current compensation
if ((do_ecc) && (nref>0))
    if nt~=1, warning('nt(=%g)~=1',nt); end
    fprintf('Using reference scan for Eddy Current Compensation\n');
    tmp = complex(zeros(si(1)-nref,si(2),si(3),2,si(5),si(6),class(d)));
    tmp(:,:,:,1,:,:) = d(1:(si(1)-nref),:,:,1,:,:);
    tmp(1,:,:,2,:,:) = d((si(1)-nref+1),:,:,1,:,:);
    d = tmp;
    nt = 2;
    % nn(5) = 2;
end


%% k-space
if any(kmax>0.51) && do_zf_spat
    % density-weighted MRSI has kmax>0.5 for nominal=real resolution
    % to avoid information loss with zero-filling, k-space is rescaled
    % and reconstructed to a larger matrix
    mtx_reco = ceil(1.1*kmax.*mtx(1,1:dim))*2;
    for l=1:dim
        if mtx_reco(1,l)>zf(1,1+l), mtx_reco(1,l) = zf(1,1+l); end
    end
    kfac = mtx_reco./mtx(1,1:dim);
    kk = bsxfun(@rdivide,permute(wf.k,[2 3 1]),kfac);
    mtx_reco = [mtx_reco ones(1,3-dim)];
else
    kk = permute(wf.k,[2 3 1]);
    mtx_reco = mtx;
end
if size(d,1)>size(kk,1)
    rep = ceil(size(d,1)/size(kk,1));
    fprintf('size(d,1)(=%g)>size(kk,1)(=%g) -> replicating kk %g times\n',...
        size(d,1),size(kk,1),rep);
    kk = repmat(kk,[rep 1]);
end
if size(d,1)<size(kk,1)
    fprintf('size(d,1)(=%g)<size(kk,1)(=%g) -> truncating kk\n',...
        size(d,1),size(kk,1));
    kk = kk(1:size(d,1),:,:,:,:,:);
end
niso = [1,1,1];               % non-isotropic factor for off-iso shift
if isfield(wf,'fov')
    if length(wf.fov)==3
        if (any(diff(wf.fov)) || any(diff(wf.mtx)))
            niso = wf.fov(1)./wf.fov;
        end
    end
else
    warning('field wf.fov not found');
end
fprintf('niso = %s\n',num2str(niso));
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx_reco.*[xloc yloc zloc]./fov.*niso;
    % shft = mtx.*[xloc yloc zloc]./fov.*niso;
else
    shft = [];
end


%% misc info + warnings
fprintf('psdname=%s/%s\n',h.image.psdname,h.image.psd_iname);
fprintf('nspec=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
    nspec,nx,ny,nz,nslices,nt,nc);
fprintf('nucleus=%g\n',h.image.specnuc);

if wf.nexc~=si(1), warning('wf.nexc(=%g)~=si(1)(=%g)',wf.nexc,si(1)); end
if (do_ecc && (nt~=2))
    warning('do_ecc && (nt(=%g)~=2',nt);
    fprintf('reference required; setting ''do_ecc = false''\n');
    do_ecc = false;
end
if ((ws==2) && (nt~=2))
    warning('ws==2 && (nt(=%g)~=2; reference required',nt);
    fprintf('reference required; setting ''ws = 0''\n');
    ws = 0;
end


%% selecting indexed data
fprintf('Selecting indexed data\n');
d = d(:,ind,:,:,:,:);


%% spatial shift
% shft = []; warning('shft=[]');
if ~isempty(shft)
    fprintf('Shifting data by (%g,%g,%g) pixel\n',shft);
    shft = reshape(shft(1,1:dim),[1 dim]);
    phafu = exp(-sum(bsxfun(@times,1i*2*pi*kk,shft),2));
    d = bsxfun(@times,d,phafu);
    clear phafu
else
    fprintf('No off-iso shifting data\n');
end


%% apodisation + density compensation: one is inverse of other -> skipping
% d = bsxfun(@times,d,wf.dcf.');
% d = bsxfun(@times,d,wf.apo.');
if ~isempty(wf.apotune)
    % fine tune apodisation function apotune = apo.*dcf;
    d = bsxfun(@times,d,wf.apotune.');
end


%% spatial reconstruction
if dim~=3 && ~ftmethod
    warning('gridding only implemented for 3D; using slowFT');
    ftmethod = 1;
end
switch ftmethod
    case 0
        % size(d)  = [#exc, #spec, 1, #time, #slices, #coils]
        % size(dd) = [#spec, #x, #y, #z, #time, #coils]
        fprintf('Spatial reconstruction: gridding\n');
        si = size(d);
        dd = gridding3dsingle(kk,reshape(d,[si(1),prod(si(2:end))]).',...
            ones(1,si(1),'single'),mtx_reco);
        dd = permute(reshape(dd,[mtx_reco(1) mtx_reco(2) mtx_reco(3) nspec nt nc]),...
            [4 1 2 3 5 6]);

    case 1
        fprintf('Spatial reconstruction: slowFT\n');
        % size(d)  = [#exc,       #spec, 1, #time, #slices, #coils]
        % size(dd) = [(#x,#y,#z), #spec, 1, #time, #slices, #coils]
        % size(dd) = [#spec, #x, #y, #z, #time, #coils]
        dd = complex(zeros(prod(mtx_reco),nspec,1,nt,nslices,nc,class(d)));
        % reconstructing without density compensation function
        % -> apodisation with acquisition density -> SNR optimal
        xx = sub_grid(mtx_reco,dim);
        if nslices>1, mtx_reco(3) = nslices; end
        if isa(d,'single')
            xx = single(xx);
            kk = single(kk);
        end
        kk = 1i*2*pi*kk;    % before minus; wrong z-shift; spatial flip
        for lk=1:si(1)
            dd = dd + bsxfun(@times,d(lk,:,:,:,:,:),...
                exp(sum(bsxfun(@times,kk(lk,:),xx),2)));
            if (mod(lk,si(1)/20)<1), fprintf('.'); end
        end
        dd = dd/size(d,1);           % consistent scaling to recon_csi.m
        fprintf('\n');
        dd = reshape(permute(dd,[2 1 3 4 5 6]),...
            [nspec mtx_reco(1) mtx_reco(2) mtx_reco(3) nt nc]);

    case 2
        fprintf('Spatial reconstruction: pinvRecon\n');
        energy_threshold = 0.95;
        memE = size(kk,1)*prod(mtx_reco)*4*2*1d-6;
        fprintf('size(E) = [%d %d] = %g [MB]\n',size(kk,1),prod(mtx_reco),memE);
        si = size(d);
        % wfn_pinvE = [wfn(1:end-4) '_pinvE.mat'];
        wfn_pinvE = [wfn(1:end-4) '_pinvE.h5'];
        if (memE>100) && exist(wfn_pinvE,'file')
            fprintf('Loading pinvE from ''%s''\n',wfn_pinvE);
            % load(wfn_pinvE,'F');
            F = readh5complex(wfn_pinvE);
        else
            xx = sub_grid(mtx_reco,dim).';
            if nslices>1, mtx_reco(3) = nslices; end
            if isa(d,'single')
                xx = single(xx);
                kk = single(kk);
            end
            kk = -1i*2*pi*kk;
            switch dim
                case 2, E = exp(kk(:,1)*xx(1,:) + kk(:,2)*xx(2,:));
                case 3, E = exp(kk(:,1)*xx(1,:) + kk(:,2)*xx(2,:) + kk(:,3)*xx(3,:));
                otherwise, error('dim(=%g)~=2 or 3',dim);
            end
            [U,S,V] = svd(E,'econ');  % courtesy to Kylie Yeung
            clear E
            svdE = diag(S);
            cumulative_energy = cumsum(svdE) / sum(svdE);
            % figure;plot(cumulative_energy);title('cumulative energy')
            cutoff = find(cumulative_energy >= energy_threshold, 1);
            invS = 1./svdE;
            invS(cutoff:end) = 0;
            invS = diag(invS);
            F = V*invS*U';
            clear U S V svdE invS
            if (memE>100)
                if (memE>1000), warning('memE(=%g)>1000 [MB]',memE); end
                fprintf('Saving pinvE to ''%s''\n',wfn_pinvE);
                % save(wfn_pinvE,'F','-v7.3');
                saveh5complex(wfn_pinvE,F);
            end
        end
        d = bsxfun(@rdivide,d,wf.dcf.');    % dw-mrsi requires apodisation
        dd = F*reshape(d,[si(1) prod(si(2:end))]);
        dd = permute(reshape(dd,[mtx_reco(1) mtx_reco(2) mtx_reco(3) nspec nt nc]),...
            [4 1 2 3 5 6]);
        clear F

    otherwise
        error('ftmethod(=%g) not found',ftmethod);
end
clear d


%% (Gaussian) line broadening in spectral domain
if abs(lb(1))>1d-5
    fprintf('Gaussian apodisation\n');
    lb = [lb(:) ; zeros(4,1)];
    if any(lb(2:end,1)>0), warning('Only spectral apodisation possible'); end
    lb_opts = {'','h','h','h'};
    lb_bw = [bw nn(2:4)];
    lbf = lb_fun(nn(1),lb_bw(1),[0 lb(1)],lb_opts{1});
    % dd = bsxfun(@times,dd,reshape(lbf,[ones(1,l-1),nn(l),1]));
    dd = bsxfun(@times,dd,lbf.');
end


%% spatial zero filling (Fourier interpolation)
if do_zf_spat && ...
        ((size(dd,2)>zf(2)) || (size(dd,3)>zf(3)) || (size(dd,4)>zf(4)))
    fprintf('Zero-filling (spatial)\n');
    % spatial zf; manually, because spatial data symmetric
    zf_tmp = zf; zf_tmp(1) = nn(1);
    dd = truma(dd,2,zf_tmp);
end


%% frequency axes
hz  = (-zf(1)/2:zf(1)/2-1)/zf(1)*bw;
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;
if h.image.specnuc<3, ppm = ppm + 4.68; end   % relative to water


%% coil combination
if coco(1)>0 && nc>1
    fprintf('Coil combination\n');
    ind_svds = [];
    dd = csi_coil_combine(dd,[],ind_svds,[],coco(2));
end


%% eddy current correction - phase deconvolution
if do_ecc
    fprintf('Eddy current correction by phase deconvolution\n');
    dd = csi_ecc(dd);
end


%% water subtraction via HSVD
if ws==1
    fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
        ws_bounds(1),ws_bounds(2));
    [~,~,~,dd(:,:,:,:,1,:)] = hsvd(dd(:,:,:,:,1,:),25,ws_bounds,bw,1);
end


%% spectral reconstruction
% size(spec) = [#s,#x,#y,#z,#t]
if isempty(freq)
    % fft recon
    if zf(1)>1
        fprintf('Spectral reconstruction\n');
        spec = ifftshift(fft(dd,zf(1),1),1);
    else
        fprintf('Skipping spectral reconstruction\n');
        spec = dd;
    end
else
    % chemical shift modelling
    fprintf('Spectral reconstruction by chemical shift inversion\n');
    A = exp(-1i*2*pi*t*freq(:).');
    pinv_A = pinv(A);
    zf(1) = length(freq);
    spec = complex(zeros(zf),class(dd));
    for l2=1:zf(2)
        for l3=1:zf(3)
            for l4=1:zf(4)
                for l5=1:zf(5)
                    for l6=1:zf(6)
                        spec(:,l2,l3,l4,l5,l6) = pinv_A*dd(:,l2,l3,l4,l5,l6);
                    end
                end
            end
        end
    end
    hz = freq;
end
clear dd


%% Water subtraction
if ws==2
    fprintf('Water subtraction\n');
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws,true);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori
    spec = csi_orientation(spec,h.data_acq_tab.rotate,...
        h.data_acq_tab.transpose,true);
end


%% phase spectra
if pseuim(1)==3
    ii = ceil(size(spec,1)/2+0.5);
    fprintf('Phasing spec (ii=%d); for pseuim(1)==3\n',ii);
    spec = csi_phasing(spec,ii);
end


%% appr. 1D spectrum
if ((nargout>2) || (plt(1)>0)) || ~isempty(fname)
    % size(spec)=[#s,#x,#y,#z,#t,#c]
    s1d = zeros(zf(1),4);
    sss = sum(spec,5);
    s1d(:,3) = sqrt(sum(sum(sum(sum(abs(sss.*conj(sss)),6),4),3),2));
    s1d(:,4) = abs(sum(sum(sum(sss,2),3),4));
    sss = sqrt(sum(abs(sss.*conj(sss)),6));
    s1d(:,1) = sum(sum(sum(sss,2),3),4);
    s1d(:,2) = max(max(max(sss,[],2),[],3),[],4);
    s1d = bsxfun(@rdivide,s1d,max(s1d,[],1));
end
par.fax = hz; par.fig_str = figstr; par.lb = lb; % par structure for plot_mrsi
par.te = te;  % echo time [s]
par.niso = niso;
par.synthesizer_frequency = h.rdb_hdr.ps_mps_freq/10;
par.sample_frequency = h.rdb_hdr.spectral_width;


%% saving result as mat file
bb = mrsi2pseuimage(spec,zf_pseuim,pseuim(2),[],pseuim(1));
if ~isempty(fname)
    fprintf('Saving results to fname=''%s''\n',fname);
    spec = single(spec); % convert: double (64bit) to single (32bit) precision
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','plt',...
        'coco','freq','ftmethod','do_ori','nviews0',...
        'hz','s1d','nn','t','dim','bb',...
        'si','mtx','nslices','nspec','nt','nc','par');
end


%% plot check
if length(size(spec))>5
    warning('(length(size(spec))>5); skipping plotting');
    plt(1) = 0;
    fname = [];
end


%% plotting pseudo image
if (plt(1)>0) || ~isempty(fname)
    if plt(1)>0
        fprintf('Plotting pseudo image (visible)\n');
        fid = figure(20); clf
    else
        fprintf('Plotting pseudo image (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    %bs = sort(bb(:),'descend');
    %mm = 0.97*mean(bs(1:ceil(0.02*length(bs))));
    %imagesc_row(bb,[0 mm],'',size(bb,4)==1,true);
    imagesc_row(bb,[],'gl3',size(bb,4)==1,true);
    set(fid,'name',[figstr '- pseudo image']);
    colormap gray
    if ~isempty(fname)
        print(fid,[fname '_pseuimage'],'-dpng','-r300');
        if isdeployed
            % export to dicom database
            inp.SeriesDescription = sprintf('MRSI pseudo:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
                nspec,nx,ny,nz+nslices-1,nt,nc);
            if dim==3, inp.MRAcquisitionType = '3D'; end
            if any(abs(diff(niso))>1d-5)
                spacing = round(h.rdb_hdr.fov/nz/niso(3)*10)/10;
                inp.SpacingBetweenSlices = spacing;
            end
            write_dicom(fname,bb,h,inp,[],[],[],3);
        end
    end
    if plt(1)<=0, close(fid); end
end


%% pseudo spectrum
if ((plt(1)>0) || ~isempty(fname)) && ~isempty(s1d)
    if plt(1)>0
        fprintf('Plotting pseudo spectrum (visible)\n');
        fid = figure(21); clf
    else
        fprintf('Plotting pseudo spectrum (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    plot(hz,s1d(:,:,1).');
    set(gca,'XDir','reverse');
    legend('sum-abs','max','rms','sum-complex');
    set(fid,'name',[figstr '- spatially-combined spectra']);
    axis([min(hz) max(hz) -0.01 1.02]);
    grid on
    xlabel('freq [Hz]');
    if ~isempty(fname)
         print(fid,[fname '_s1d'],'-dpng','-r300');
        if isdeployed
             % export to dicom database
            inpscd.SeriesNumber = h.series.se_no+100;
            inpscd.SeriesDescription = ...
                sprintf('MRSI Se%d',h.series.se_no);
            inpscd.InstanceNumber = 0;
            write_scdicom([fname '_s1d.dcm'],fid,h,inpscd);
        end
    end
    if plt(1)<=0, close(fid); end
end


%% plotting csi grid (invisible + dicom)
if ~isempty(fname) && (plt(2)>0)
    if isdeployed, hh = h;
    else, hh = []; end
    plot_csi(spec,[],[],[figstr '- MRSI'],[fname '_plot_csi'],11,hh);
end


%% plotting 3-plane of pseudo images
if ((dim>2) && ((plt(1)>0) || ~isempty(fname)))
    % bb3 = mrsi2pseuimage(spec,[128 128 128],pseuim(2),[],pseuim(1));
    if isempty(zf_pseuim), zf_pseuim = [128 128 128]; end
    bb3 = mrsi2pseuimage(spec,zf_pseuim,pseuim(2),[],pseuim(1));
    if plt(1)>0
        fprintf('Plotting pseudo 3-plane (visible)\n');
        fid = figure(22); clf
    else
        fprintf('Plotting pseudo 3-plane (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    % bs = sort(bb3(:),'descend');
    % mm = 0.96*mean(bs(1:20d3));
    % imagesc_ind3d(bb3,[],[0 mm],[],false,true);
    rot = [];
    if (h.image.plane==2), rot = [0 1 1]; end
    imagesc_ind3d(bb3,[],[],[],false,true,'gl3',rot);
    set(fid,'name',[figstr '- pseudo 3-plane']);
    colormap gray
    if ~isempty(fname)
        print([fname '_pseu3plane'],'-dpng','-r300');
    end
    if plt(1)<=0, close(fid); end
end


%% plotting interactive MRSI viewer
if plt(1)>1
    fprintf('Starting plot_mrsi\n');
    plot_mrsi(spec,par);
end


%% finishing up
fprintf('recon_mrsi: runtime = %.2f s\n',toc(timerVal));


end      % recon_mrsi.m


%% sub-functions
% Cartesian grid
function xx = sub_grid(mtx,dim)

for l=1:dim, x{l} = (-mtx(l)/2:mtx(l)/2-1); end
% for l=1:dim, x{l} = (-mtx(l)/2:mtx(l)/2-1)+1; end

xx = zeros(1,prod(mtx),dim);
switch dim
    case 1, xx = x{l};
    case 2
        [x1,x2] = ndgrid(x{1},x{2});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
    case 3
        [x1,x2,x3] = ndgrid(x{1},x{2},x{3});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
        xx(1,:,3) = x3(:).';
    otherwise, error('dim(=%g)~=1,2 or 3');
end

xx = permute(xx,[2 3 1]);

end      % sub_grid


%% zero filling
function [zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices,kmax)
do_zf_spat = false;
if ~isempty(zf)
    zf = [zf(:).', nn((length(zf)+1):6)];
    if any(zf<0)             % negative numbers -> factor
        ii = (zf<0) & (nn>1);
        zf(ii) = -zf(ii).*nn(ii);
        zf(nn==1) = 1;
    end
    if any(abs(zf)<1d-10)    % 0 -> no zero-filling
        if any(abs(zf(2:4))<1d-10)
            zf_tmp = ceil(1.1*kmax.*nn(1,2:4))*2;
            ii = (abs(zf)<1d-10);
            ii(1,1) = false; ii(1,5:end) = false;
            zf(ii) = zf_tmp(1,ii(1,2:4));
        end
        ii = abs(zf)<1d-10;
        zf(ii) = nn(ii);
    end
    if any(zf<nn)            % no cropping
        warning('zf<nn');
        fprintf('zf = %s\n',num2str(zf));
        fprintf('nn = %s\n',num2str(nn));
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = %s\n',num2str(zf));
    end
    if nslices>1
        if zf(4)>nslices
            warning('zf(4)(=%g) > nslices(=%g); setting zf(4)=%g',...
                zf(4),nslices,nslices);
            zf(4) = nslices;
        end
    end
    if ((nslices==1) && (nz==1) && (zf(4)>1))
        warning('zf(4)(=%g)>1 and 2D single slice (nz=%g && nslices=%g)',...
            zf(4),nz,nslices);
        fprintf('  setting zf(4) to 1\n');
        zf(4) = 1;
    end
    zf = round(zf);
    
    if any(zf(2:4)~=nn(2:4)), do_zf_spat = true; end
else
    zf = nn;
end

end      % sub_zf


%% interleaved centre-of-k-space spectra
function d = sub_reco_cok(d,h,ind1,nexc,zf,lb,fname,plt,figstr)
fprintf('Interleaved spectra from centre of k-space\n');
% ind1 = true(nexc,1);
% ind1(1:ncok:end,1) = false;
nrep = size(d,1)/nexc;
if nrep>1
    warning('nrep=%g: replicating ind1',nrep);
    ind1 = repmat(ind1,[ceil(nrep) 1]);
end
if size(ind1,1)>size(d,1)
    warning('size(ind1,1)(=%d)>size(d,1)(=%d): truncating ind1',...
        size(ind1,1),size(d,1));
    ind1 = ind1(1:size(d,1),1);
end
dcok = d(~ind1,:,:,:,:,:);
d = d(ind1,:,:,:,:,:);
[spec,par,~,hz] = fid2spec(dcok,h,zf(1),lb,[],false);
if ~isempty(fname)
    save([fname '_cok.mat'],'-mat','spec','par','hz');
end

if (plt(1)>0) || ~isempty(fname)
    if plt(1)>0
        fid = figure(23); clf
    else
        fid = figure('Visible','off');
    end
    abspec = sqrt(mean(spec.*conj(spec),6));
    x = repmat(hz,[size(spec,1) 1]).';
    z = repmat((1:size(spec,1)).',[1 size(spec,2)]).';
    plot3(z,x,abspec);
    axis tight; grid on
    view(-85,23);
    ylabel('freq [Hz]');
    set(fid,'name',[figstr '- interleaved spectra']);
    if ~isempty(fname)
         print(fid,[fname '_cok'],'-dpng','-r300');
        if isdeployed
             % export to dicom database
            inpscd.SeriesNumber = h.series.se_no+100;
            inpscd.SeriesDescription = ...
                sprintf('MRSI Se%d',h.series.se_no);
            inpscd.InstanceNumber = 1;
            write_scdicom([fname '_cok.dcm'],fid,h,inpscd);
        end
    end
    if plt(1)<=0, close(fid); end
end

end      % sub_reco_cok
