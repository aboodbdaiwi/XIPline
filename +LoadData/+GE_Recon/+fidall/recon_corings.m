function [spec,hz,s1d,ppm,par] = recon_corings(d,h,wfn,zf,lb,fname,plt,...
    coco,ws,do_ecc,d_ref,pseuim)
%RECON_CORINGS  Reconstruct MRSI data acquired with concentric rings
%[spec,hz,s1d,ppm,par]=recon_corings()
%        d  Raw p-file data 
%           [(n_ring*n_intlv*n_zpe),ns,1,n_time,nslices,ncoils] 
%           Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%        h  Header structure
%      wfn  Location of waveform .mat file
%       zf  Zero fill to:   [#spec,#x,#y,#z]
%       lb  Line broadening [spec,x,y,z]         [Hz,mtx]               (0)
%    fname  Print <fname>.png and save reco as <fname>.mat             ([]) 
%           also export max-image and plots as dicom
%      plt  Plotting: 0=off, 1=normal, 2+=plot_mrsi                     (2)
%     coco  Coil combination [do percent_noise4decor]               ([1 0])
%       ws  Water suppression: 0=off; 1=HSVD; 2=subtraction             (0)
%   do_ecc  Phase deconvolution eddy current correction             (false)
%           requires H2O reference scan: 
%           if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%           if ~isempty(d_ref) -> default=true; 
%    d_ref  Reference scan (w/o water suppression)
%   pseuim  How to generate pseudo image                            ([1 4])
%           mrsi2pseuimage.m -> [how2ind npts]
%
%     spec  Spectrum        [#s,#x,#y,#z,#t]
%       hz  Frequency axis  [Hz]
%      s1d  Pseudo spectra  [#s,3]
%      ppm  Frequency axis  [ppm]
%      par  Parameter structure for plot_mrsi
%
%  7/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;               % record execution time


%% misc input
if ~exist('wfn','var'),    wfn = ''; end
if ~exist('zf','var'),     zf = []; end
if ~exist('lb','var'),     lb = []; end
if ~exist('fname','var'),  fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$', 'once')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$', 'once')), fname = fname(1:end-4); end
    end
end
if ~exist('plt','var'),    plt = []; end
if isempty(plt),           plt = 2; end
if length(plt)~=1,         warning('length(plt)~=1'); end
if ~exist('coco','var'),   coco = []; end
if isempty(coco),          coco = 1; end
if length(coco)<2,         coco = [coco 0]; end
if ~exist('ws','var'),     ws = []; end
if isempty(ws),            ws = 0;  end
if length(ws)~=1,          warning('length(ws)~=1'); end
if ~exist('do_ecc','var'), do_ecc = []; end
if ~exist('d_ref','var'),  d_ref = []; end
if ~exist('pseuim','var'), pseuim = []; end
if isempty(pseuim),        pseuim = 1;  end
if length(pseuim)<2,       pseuim = [pseuim 4]; end
do_ori = true;          % correct image orientation
dbg = false;            % print additional debugging info


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
end
d = single(d);
d = conj(d);                 % changing dir of axis (keep MRS convention)


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.k$', 'once'))
        if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
        wf = create_corings_wf(h,wfn);
    else
        if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
        if isempty(regexpi(wfn,'\.mat$', 'once')),  wfn = [wfn '.mat']; end
        if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
        wf = load(wfn);          % load waveform
    end
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
nslices = h.rdb_hdr.nslices;
if ~exist('nt','var'), nt = h.image.numecho; end      % #time steps
nc = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;           % #Rx coil elements
nref = h.rdb_hdr.user19;               % #reference frames (w/o WS)
ns = wf.npts;                          % #sampling points
ind = wf.ind;                          % index to select acquired data
nspec = wf.n_spec;                     % #spectral points
nk_circ = wf.nk_circ;                  % #sample points per 1 ring
bw_spec = wf.bw_spec;                  % spectral bandwidth [Hz]
n_ring = wf.n_ring;                    % #rings
n_intlv = wf.n_intlv;                  % #spectral interleaves
n_zpe = wf.n_zpe;                      % #z-phase encodings
if abs(rem(nk_circ,n_intlv))>1d-10
    warning('nk_circ(=%g) not multiple of n_intlv(=%g)',nk_circ,n_intlv);
end
if abs(rem(nspec,n_intlv))>1d-10
    warning('n_spec(=%g) not multiple of n_intlv(=%g)',nspec,n_intlv);
end

if size(d,1)<(n_ring*n_zpe*n_intlv)
    error('size(d,1)(=%d)<n_ring(=%d)*n_zpe(=%d)*n_intlv(=%d)',...
        size(d,1),n_ring,n_zpe,n_intlv);    
end
if size(d,1)>(n_ring*n_zpe*n_intlv)
    warning('size(d,1)(=%d)>n_ring(=%d)*n_zpe(=%d)*n_intlv(=%d): truncating d',...
        size(d,1),n_ring,n_zpe,n_intlv);
    ii = 1:(n_ring*n_zpe*n_intlv);
    d = d(ii,:,:,:,:,:);
end
if si(2)~=ns, warning('si(2)(=%g)~=ns(=%g)',si(2),ns); end
if si(3)~=1, warning('si(3)(=%g) ~= 1',si(3)); end
if nt~=si(4)
    warning('nt(=%g) ~= size(d,4)(=%g); setting nt=%g',nt,si(4),si(4));
    nt = si(4);
end
if si(5)~=nslices
    warning('si(5)(=%g) ~= nslices(=%g); setting n_slices=%d',...
        si(5),nslices,si(5));
    nslices = si(5);
end
if nc~=si(6)
    warning('nc(=%g) ~= size(d,6)(=%g); setting nc=%g',nc,si(6),si(6));
    nc = si(6);
end
if ((nz>1) && (nslices>1))
    warning('nz(=%g)>1 && nslices(=%g)>1; truncating slices',nz,nslices);
    nslices = 1;
    d = d(:,:,:,:,1,:);
end
nz = max([nz nslices]);                % use nz for 2D multi-slice and 3D
dim = wf.dim;                          % spatially 2D or 3D
nn = [nspec nx ny nz nt nc];           % final data size (6D)

bw = h.rdb_hdr.spectral_width;              % sampling bandwidth [Hz]
itmp = find(diff(wf.ind));
te = h.rdb_hdr.te*1d-6+(itmp(1)+1)/bw;
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
fov = h.rdb_hdr.fov;
zf_pseuim = zf(5:end);
zf = zf(1:min(4,length(zf)));


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
    shft = mtx.*[xloc yloc zloc]./fov.*niso;
    % shft = mtx.*[xloc yloc zloc]./fov;
else
    shft = [];
end


%% misc info + warnings
fprintf('psdname=%s/%s\n',h.image.psdname,h.image.psd_iname);
fprintf('nucleus=%g\n',h.image.specnuc);
fprintf('ns=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
    ns,nx,ny,nz,nslices,nt,nc);
fprintf('nspec=%g; nk_circ=%g; n_ring=%g; n_intlv=%g; n_zpe=%g\n',...
    nspec,nk_circ,n_ring,n_intlv,n_zpe);
fprintf('bw_spec=%g[Hz]; bw_spat=%g[kHz]\n',...
    bw_spec,2*h.rdb_hdr.bw);

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


%% unchop data
if ischop(h), d(2:2:end,:) = -d(2:2:end,:); end


%% selecting indexed data
fprintf('Selecting indexed data\n');
d = d(:,ind(1,:),:,:,:,:);


%% spatial shift
if ~isempty(shft)
    if ((dim>2) && (abs(shft(1,3))>0.01))
        zshft = true;
    else
        zshft = false;
    end

    kk = single(wf.k);
    if dim>2
        kk = repmat(kk,[1 n_zpe 1]);
        if zshft
            kk(1,:,3) = repelem(single(wf.kz),1,round(size(kk,2)/n_zpe));
        end
    end
    kk = 1i*2*pi*permute(reshape(kk,[size(d,2) size(d,1) (2+zshft)]),[2 1 3]);
    fprintf('Shifting data by (%g,%g,%g) pixel\n',shft);
    shft2 = reshape(shft(1,1:(2+zshft)),[1 1 (2+zshft)]);
    phafu = exp(-sum(bsxfun(@times,kk,shft2),3));
    d = bsxfun(@times,d,phafu);
    clear phafu kk
end


%% reshaping of data
% size(d)  = [(#rings*#zpe*#intlv), #acq, 1, #time, #slices, #coils]
% size(dd) = [#spec,  #k-space one ring, 1, #z-phase/slices, #time, #coils]
fprintf('Reshaping data\n');

dd = complex(zeros(nspec,n_ring*nk_circ,1,max([n_zpe nslices]),...
    nt,nc,class(d)));

if dbg, fprintf('l_spec  i2(1) i2(end) ii1(1) ii2(1) ii2(end)\n'); end
for lc=1:nc
    for lt=1:nt
        for l_slice=1:nslices
            for lz=1:n_zpe
                for l_spec=1:nspec
                    if n_zpe>1
                        l4 = lz;
                    else
                        l4 = l_slice;
                    end
                    if n_intlv==1
                        ii1 = (1:n_ring) + (lz-1)*n_ring;
                        ii2 = (1:nk_circ) + (l_spec-1)*nk_circ;
                        dd(l_spec,:,1,l4,lt,lc) = ...
                            reshape(d(ii1,ii2,1,lt,l_slice,lc).',...
                            [1 n_ring*nk_circ]);
                    else
                        for li=1:n_intlv
                            i2 = bsxfun(@plus,(1:(nk_circ/n_intlv)).',...
                                (0:(n_ring-1))*nk_circ+(li-1)*nk_circ/n_intlv);
                            i2 = i2(:).';
                            ss1 = mod(li-l_spec,n_intlv)+1;
                            ii1 = ss1:n_intlv:(n_intlv*n_ring);
                            ii2 = (1:(nk_circ/n_intlv)) + ...
                                (l_spec-1)*nk_circ/n_intlv;
                            if dbg
                                fprintf('%6d %6d %6d %6d %6d %6d\n',...
                                    l_spec,i2(1),i2(end),ii1(1),ii2(1),ii2(end));
                            end
                            tmp = reshape(d(ii1,ii2,1,lt,l_slice,lc).',...
                                [1 n_ring*nk_circ/n_intlv]);
                            dd(l_spec,i2,1,l4,lt,lc) = tmp;
                        end
                    end
                end
            end
        end
    end
end
clear d


%% correct for time evolution during one circle
fprintf('Correct for time evolution\n');
% phafu = repmat(exp(2*pi*1i*wf.t1*bw_spec),[1 n_ring*n_intlv]);
% dd = fft(bsxfun(@times,ifft(dd,[],2),phafu),[],2);
% dd = fft(bsxfun(@times,ifft(dd,[],1),phafu),[],1);
% dd = fft(bsxfun(@times,ifft(dd,256,1),phafu),[],1); dd = dd(1:10,:);
% dd = ifft(bsxfun(@times,fft(dd,256,1),phafu),[],1); dd = dd(1:10,:);
% dd = bsxfun(@times,dd,phafu);
warning('not yet implemented');


%%  spatial reconstruction along z
if n_zpe>1
    ntmp = size(dd); ntmp(4) = mtx(3);
    ddz = complex(zeros(ntmp,class(dd)));
    fprintf('Spatial reconstruction along z\n');
    kz = 1i*2*pi*wf.kz;
    xx = reshape((-mtx(3)/2):(mtx(3)/2-1),[1 1 1 mtx(3)]);
    for lkz=1:length(kz)
        ddz = ddz + bsxfun(@times,dd(:,:,:,lkz,:,:),...
            exp(kz(1,lkz)*xx));
    end
    if false
        [~,ii] = sort(wf.kz);
        ddz = ifftshift(ifft(fftshift(dd(:,:,:,ii,:,:),4),[],4),4);
    end
    dd = ddz;
    clear ddz
end


%%  spatial reconstruction along xy
% size(dd) = [#spec,  #k-space one ring, 1, #z, #time, #coils]
% size(spec) = [#spec,  #x, #y, #z, #time, #coils]
spec = complex(zeros(nn,'single'));
fprintf('Spatial reconstruction along xy\n');
for lz=1:max([nz nslices])
    for lt=1:nt
        for lc=1:nc
            spec(:,:,:,lz,lt,lc) = single(permute(gridding(dd(:,:,1,lz,lt,lc),...
                wf.k1,wf.dcf1,mtx(1)),[3 1 2]));
        end
    end
end
clear dd


%% (Gaussian) line broadening in spectral domain
if ~isempty(lb)
    fprintf('Gaussian apodisation\n');
    lb = [lb(:) ; zeros(4,1)];
    if any(lb(2:end,1)>0), warning('Only spectral apodisation possible'); end
    if abs(lb(1))>1d-10
        lbf = single(lb_fun(nn(1),bw_spec,[0 lb(1)]));
        spec = bsxfun(@times,spec,lbf.');
    end
end


%% zero filling (Fourier interpolation)
if ~isempty(zf)
    fprintf('Zero-filling\n');
    zf = [zf(:).', nn((length(zf)+1):6)];
    if any(zf<nn)
        warning('zf<nn'); 
        fprintf('zf = '); disp(zf)
        fprintf('nn = '); disp(nn)
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = '); disp(zf)
    end
    if ((nslices==1) && (nz==1) && (zf(4)>1))
        warning('zf(4)(=%g)>1 and 2D single slice (nz=%g && nslices=%g)',...
            zf(4),nz,nslices);
        fprintf('  setting zf(4) to 1\n');
        zf(4) = 1;
    end
    
    % spatial zf; manually, because spatial data symmetric
    zf_tmp = zf; zf_tmp(1) = nn(1);
    spec = truma(spec,true,zf_tmp);
else
    zf = nn;
end


%% frequency axes
hz  = (-zf(1)/2:zf(1)/2-1)/zf(1)*bw_spec;
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;
if h.image.specnuc<3, ppm = ppm + 4.68; end   % relative to water


%% coil combination
if coco(1)>0 && nc>1
    fprintf('Coil combination\n');
    spec = csi_coil_combine(spec,[],[],[],coco(2));
end


%% eddy current correction - phase deconvolution
if do_ecc
    fprintf('Eddy current correction by phase deconvolution\n');
    warning('must be adapted');
    spec = csi_ecc(spec);
    % dd = csi_ecc(dd,2);
end


%% HSVD water removal
if ws==1
    fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
        ws_bounds(1),ws_bounds(2));
    % spec = ifft(fftshift(spec,1),[],1);
    [~,~,~,spec] = hsvd(spec,25,ws_bounds,bw_spec,1);
    % spec = ifftshift(fft(spec,[],1),1);
end


%% spectral reconstruction
fprintf('Spectral reconstruction\n');
spec = ifftshift(fft(spec,zf(1),1),1);


%% Water subtraction
if ws==2
    fprintf('Water subtraction\n');
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori
    spec = csi_orientation(spec,h.data_acq_tab.rotate,...
        h.data_acq_tab.transpose,true);
end


%% appr. 1D spectrum
if ((nargout>2) || (plt>0)) || ~isempty(fname)
    % size(spec)=[#s,#x,#y,#z,#t,#c]
    s1d = zeros(zf(1),3);
    sss = sum(spec,5);
    s1d(:,3) = sqrt(sum(sum(sum(sum(abs(sss.*conj(sss)),6),4),3),2));
    sss = sqrt(sum(abs(sss.*conj(sss)),6));
    s1d(:,1) = sum(sum(sum(sss,2),3),4);
    s1d(:,2) = max(max(max(sss,[],2),[],3),[],4);
    s1d = bsxfun(@rdivide,s1d,max(s1d,[],1));
    clear sss
end
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
par.fax = hz; par.fig_str = figstr; par.lb = lb; % par structure for plot_mrsi
par.te = te;  % echo time [s]
par.niso = niso;


%% saving result as mat file
bb = mrsi2pseuimage(spec,zf_pseuim,pseuim(2),[],pseuim(1));
if ~isempty(fname)
    fprintf('Saving results to fname=''%s''\n',fname);
    spec = single(spec);
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','plt',...
        'coco','do_ori',...
        'hz','s1d','nn','dim',...
        'si','nx','ny','nz','nslices','ns','nt','nc');
end


%% plot check
if length(size(spec))>5
    warning('(length(size(spec))>5); skipping plotting');
    plt = 0;
    fname = [];
end


%% plotting pseudo image
if (plt>0) || ~isempty(fname)
    if plt>0
        fprintf('Plotting pseudo image (visible)\n');
        fid = figure(20); clf
    else
        fprintf('Plotting pseudo image (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    imagesc_row(bb,[],'gl3',size(bb,4)==1,true);
    set(fid,'name',[figstr '- pseudo image']);
    colormap gray
    if ~isempty(fname)
        if isdeployed
            print(fid,[fname '_pseuimage'],'-dpng','-r300','-painters');
            % export to dicom database
            inp.SeriesDescription = sprintf('MRSI pseudo:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
                nspec,nx,ny,nz+nslices-1,nt,nc);
            if dim==3, inp.MRAcquisitionType = '3D'; end
            if any(abs(diff(niso))>1d-5)
                spacing = round(h.rdb_hdr.fov/nz/niso(3)*10)/10;
                inp.SpacingBetweenSlices = spacing;
            end
            write_dicom(fname,bb,h,inp,[],[],[],3);
        else
            print(fid,[fname '_pseuimage'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% pseudo spectrum
if ((plt>0) || ~isempty(fname)) && ~isempty(s1d)
    if plt>0
        fprintf('Plotting pseudo spectrum (visible)\n');
        fid = figure(21); clf
    else
        fprintf('Plotting pseudo spectrum (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    plot(hz,s1d(:,:,1).');
    set(gca,'XDir','reverse');
    legend('sum','max','rms');
    set(fid,'name',[figstr '- pseudo 1D spectra']);
    axis([min(hz) max(hz) -0.01 1.02]);
    grid on
    xlabel('freq [Hz]');
    if ~isempty(fname)
        if isdeployed
            print(fid,[fname '_s1d'],'-dpng','-r300','-painters');
             % export to dicom database
            inpscd.SeriesNumber = h.series.se_no+100;
            inpscd.SeriesDescription = ...
                sprintf('MRSI Se%d',h.series.se_no);
            inpscd.InstanceNumber = 0;
            write_scdicom([fname '_s1d.dcm'],fid,h,inpscd);
        else
            print(fid,[fname '_s1d'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% plotting csi grid (invisible + dicom)
if ~isempty(fname)
    if isdeployed, hh = h;
    else, hh = []; end
    plot_csi(spec,[],[],[figstr '- MRSI'],[fname '_plot_csi'],11,hh);
end


%% plotting 3-plane of pseudo images
if ((dim>2) && ((plt>0) || ~isempty(fname)))
    % bb3 = mrsi2pseuimage(spec,[128 128 128],pseuim(2),[],pseuim(1));
    if isempty(zf_pseuim), zf_pseuim = [128 128 128]; end
    bb3 = mrsi2pseuimage(spec,zf_pseuim,pseuim(2),[],pseuim(1));
    if plt>0
        fprintf('Plotting pseudo 3-plane (visible)\n');
        fid = figure(22); clf
    else
        fprintf('Plotting pseudo 3-plane (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    % bs = sort(bb3(:),'descend');
    % mm = 0.96*mean(bs(1:20d3));
    % imagesc_ind3d(bb3,[],[0 mm],[],false,true);
    imagesc_ind3d(bb3,[],[],[],false,true,'gl3');
    set(fid,'name',[figstr '- pseudo 3-planee']);
    colormap gray
    if ~isempty(fname)
        if isdeployed
            print([fname '_pseu3plane'],'-dpng','-r300','-painters');
        else
            print([fname '_pseu3plane'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% plotting interactive MRSI viewer
if plt>1
    fprintf('Starting plot_mrsi\n');
    plot_mrsi(spec,par);
end


%% finishing up
fprintf('recon_corings: runtime = %.2f s\n',toc(timerVal));


end      % recon_corings.m

