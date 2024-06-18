function [spec,hz,s1d,ppm,par,h] = recon_mrsi(d,h,wfn,zf,lb,do_ecc,ws,fname,verb,...
    do_coco,abs_fid,freq,lcm_basis,d_ref,do_ori)
%RECON_MRSI  Reconstruct MRSI data acquired with arbitrary spatial k-space
%[spec,hz,s1d,ppm,par]=recon_mrsi(d,h,wfn,zf,lb,do_ecc,ws,fname,verb,...
%    do_coco,abs_fid,freq,lcm_basis,d_ref,do_ori)
%                                                               (default)
%        d  Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc] (ns=sampling pts;nc=coils)
%           Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%        h  Header structure
%      wfn  Location of waveform .mat file
%       zf  Zero fill to [#s,#x,#y,#z]; if<0->multiplier
%       lb  Spectral line broadening                    [Hz]
%   do_ecc  Phase deconvolution eddy current correction         (false)
%           requires H2O reference scan: 
%           if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%           if ~isempty(d_ref) -> default=true; 
%       ws  Water suppression: 0=off; 1=HSVD; 2=subtraction     (0)
%    fname  Print <fname>.png and save reco as <fname>.mat      ([])
%           also export max-image and plots as dicom
%     verb  Verbose mode + plotting:                            (2)
%           0=off; 1=verb; 2=verb+plot; 3=verb+invisble plot
%  do_coco  Do coil combination                                 (true)
%  abs_fid  Spectral Fourier transform of abs(fid)              (false)
%     freq  Reconstruct spectral domain with matrix inversion   ([])
%lcm_basis  Filename of LCModel basis -> save files             ([])
%    d_ref  Reference scan (w/o water suppression)
%   do_ori  Image orientation: x-y-z to AP-RL-SI                (true)
%
%     spec  Spectrum        [#s,#x,#y,#z,#t]
%       hz  Frequency axis  [Hz]
%      s1d  abs(spec) summed spatially+coil -> 1D spectral
%      ppm  Frequency axis  [ppm]
%      par  Parameter structure for plot_mrsi
%
% 12/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;               % record execution time


%% misc input
export_dcm = true;           % export max-image + plots as dicom
if ~exist('zf','var'),     zf = []; end
if ~exist('lb','var'),     lb = []; end
if ~exist('do_ecc','var'), do_ecc = []; end
if ~exist('ws','var'),     ws = []; end
if isempty(ws),            ws = 0;  end
if length(ws)~=1,          warning('length(ws)~=1'); end
if ~exist('fname','var'),  fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    end
end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end
if length(verb)~=1,        warning('length(verb)~=1'); end
if ~exist('do_coco','var'),do_coco = []; end
if isempty(do_coco),       do_coco = true; end
if length(do_coco)~=1,     warning('length(co_coco)~=1'); end
if ~exist('abs_fid','var'),abs_fid = []; end
if isempty(abs_fid),       abs_fid = false; end
if length(abs_fid)~=1,     warning('length(abs_fid)~=1'); end
if ~exist('freq','var'),   freq = []; end
if ~isempty(freq), if ~isnumeric(freq), warning('~isnumeric(freq)'); end; end
if ~exist('lcm_basis','var'), lcm_basis = []; end
if ~isempty(lcm_basis) 
    if ~islogical(lcm_basis) || ~ischar(lcm_basis)
        error('lcm_basis not logical or char');
    end
end
if ~exist('d_ref','var'),  d_ref = []; end
if ~exist('cor_intlv','var'), cor_intlv = []; end
if isempty(cor_intlv),     cor_intlv = true; end
if length(cor_intlv)~=1,   warning('length(cor_intlv)~=1'); end
if ~exist('do_ori','var'), do_ori = []; end
if isempty(do_ori),        do_ori = true; end
if length(do_ori)~=1,      warning('length(do_ori)~=1'); end


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
    tmp = complex(zeros(size(d,1),size(d,2),size(d,3),2,size(d,5),size(d,6)));
    tmp(:,:,:,1,:,:) = d;
    tmp(:,:,:,2,:,:) = d_ref;
    d = tmp;
    clear tmp d_ref
    nt = 2;
end


%% miscelaneous parameters
ws_bounds = [-100 50];                      % bounds for HSVD water suppression
si = size(d); si = [si ones(1,6-length(si))];
mtx = wf.mtx(:).'; mtx = [mtx ones(1,3-length(mtx))];
nx = mtx(1);
ny = mtx(2);
nz = mtx(3);
ns = h.rdb_hdr.user1;                       % #sampled points (FID)
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
bw = h.rdb_hdr.spectral_width;              % sampling bandwidth [Hz]
t  = (0:nspec-1).'/bw;                         % time [s]
itmp = find(diff(wf.ind));
te = h.rdb_hdr.te*1d-6+(itmp(1)+1)/bw;

if ((do_ecc) && (nref>0))
    if nt~=1, warning('nt(=%g)~=1',nt); end
    if verb>0, fprintf('Using reference scan for ECC\n'); end
    tmp = complex(zeros(si(1)-nref,si(2),si(3),2,si(5),si(6)));
    tmp(:,:,:,1,:,:) = d(1:(si(1)-nref),:,:,1,:,:);
    tmp(1,:,:,2,:,:) = d((si(1)-nref+1),:,:,1,:,:);
    d = tmp;
    nt = 2;
    % nn(5) = 2;
end
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
fov = h.rdb_hdr.fov; 
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx.*[xloc yloc zloc]./fov;
    warning('not sure about shft sign; 1H different from 2H!?');
    fprintf('check freq_inv; off-iso for 1H-2H, mrsi-spiral\n');
    fprintf('also check with different ak_iso in fidall\n');
else
    shft = [];
end


%% k-space
[zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices);
kmax = max(abs(wf.k));
kmax = kmax(:).';
if any(kmax>0.51) && do_zf_spat
    % density-weighted MRSI has kmax>0.5 for nominal=real resolution
    % to avoid information loss with zero-filling, k-space is rescaled
    % and reconstructed to a larger matrix
    mtx_reco = ceil(1.1*kmax.*mtx(1,1:dim))*2;
    for l=1:dim
        if mtx_reco(1,l)>zf(1,1+l), mtx_reco(1,l) = zf(1,1+l); end
    end
    kmax = mtx_reco./mtx(1,1:dim)/2;
    kk = bsxfun(@times,permute(wf.k,[2 3 1]),-1i*pi./kmax);
    mtx_reco = [mtx_reco ones(1,3-dim)];
else
    kk = permute(-1i*2*pi*wf.k,[2 3 1]);
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


%% misc info + warnings
if verb>0
    fprintf('psd_iname=%s\n',h.image.psd_iname);
    fprintf('ns=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
        nspec,nx,ny,nz,nslices,nt,nc);
    fprintf('nucleus=%g\n',h.rdb_hdr.user2);
end
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
if verb>0, fprintf('Selecting indexed data\n'); end
d = d(:,ind,:,:,:,:);


%% spatial shift
% shft = []; warning('shft=[]');
if ~isempty(shft)
    if verb>0, fprintf('Shifting data by (%g,%g,%g) pixel\n',shft); end
    shft = reshape(shft(1,1:dim),[1 dim]);
    phafu = exp(-sum(bsxfun(@times,kk,shft),2));
    d = bsxfun(@times,d,phafu);
end


%% apodisation + density compensation: one is inverse of other -> skipping
% d = bsxfun(@times,d,wf.dcf.');
% d = bsxfun(@times,d,wf.apo.');
if ~isempty(wf.apotune)
    % fine tune apodisation function apotune = apo.*dcf;
    d = bsxfun(@times,d,wf.apotune.');
end


%%  spatial reconstruction
% size(d)  = [#exc,       #spec, 1, #time, #slices, #coils]
% size(dd) = [(#x,#y,#z), #spec, 1, #time, #slices, #coils]
% size(dd) = [#spec, #x, #y, #z, #time, #coils]
if verb>0, fprintf('Reshaping data\n'); end
dd = complex(zeros(prod(mtx_reco),nspec,1,nt,nslices,nc));
% reconstructing without density compensation function
% -> apodisation with acquisition density -> SNR optimal
xx = sub_grid(mtx_reco,dim);


if verb>0, fprintf('Spatial reconstruction\n'); end
for lk=1:si(1)
    dd = dd + bsxfun(@times,d(lk,:,:,:,:,:),...
        exp(sum(bsxfun(@times,kk(lk,:),xx),2)));
    if verb>0 && (mod(lk,si(1)/20)<1)
        fprintf('.');
    end
end
dd = dd/size(d,1);           % consistent scaling to recon_csi.m

if verb>0, fprintf('\n'); end
dd = reshape(permute(dd,[2 1 3 4 5 6]),...
    [nspec mtx_reco(1) mtx_reco(2) mtx_reco(3) nt nc]);


%% (Gaussian) line broadening in spectral domain
if ~isempty(lb)
    if verb>0, fprintf('Gaussian apodisation\n'); end
    lb = [lb(:) ; zeros(4,1)];
    if any(lb(2:end,1)>0), warning('Only spectral apodisation possible'); end
    if lb(1)~=0
        lbf = lb_fun(nspec,bw,[0 lb(1)]);
        dd = bsxfun(@times,dd,lbf.');
    end
end


%% spatial zero filling (Fourier interpolation)
if do_zf_spat
    if verb>0, fprintf('Zero-filling (spatial)\n'); end
    % spatial zf; manually, because spatial data symmetric
    zf_tmp = zf; zf_tmp(1) = nn(1);
    dd = truma(dd,true,zf_tmp);
end
hz  = -(-zf(1)/2:zf(1)/2-1)/zf(1)*bw;        % frequency axis
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;


%% coil combination
if do_coco && nc>1
    if verb>0, fprintf('Coil combination\n'); end
    dd = csi_coil_combine(dd);
end


%% eddy current correction - phase deconvolution
if do_ecc
    if verb>0
        fprintf('Eddy current correction by phase deconvolution\n'); 
    end
    dd = csi_ecc(dd);
end


%% water subtraction via HSVD
if ws==1
    if verb>0
        fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
            ws_bounds(1),ws_bounds(2));
    end
    [~,~,~,dd(:,:,:,:,1,:)] = hsvd(dd(:,:,:,:,1,:),25,ws_bounds,bw,1);
end


%% spectral reconstruction
% size(spec) = [#s,#x,#y,#z,#t]
if isempty(freq)
    % fft recon
    if zf(1)>1
        if verb>0, fprintf('Spectral reconstruction\n'); end
        if abs_fid
            if verb>0, fprintf('\ttaking abs(dd) only\n'); end
            dd = abs(dd);
        end
        spec = ifftshift(fft(dd,zf(1),1),1);
    else
        if verb>0, fprintf('Skipping spectral reconstruction\n'); end
        spec = dd;
    end
else
    % chemical shift modelling
    if verb>0
        fprintf('Spectral reconstruction by chemical shift inversion\n'); 
    end
    A = exp(-1i*2*pi*t*freq(:).');
    pinv_A = pinv(A);
    zf(1) = length(freq);
    spec = complex(zeros(zf));
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


%% Water subtraction
if ws==2
    if verb>0, fprintf('Water subtraction\n'); end
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori
    spec = csi_orientation(spec,h);
end


%% appr. 1D spectrum
if ((nargout>2) || (verb>1)) || ~isempty(fname)
    s1d = zeros(zf(1),2,nt);
    for lt=1:nt
        if abs_fid
            sss = sum(real(spec(:,:,:,:,lt,:)),6);
        else
            sss = abs(sqrt(sum(spec(:,:,:,:,lt,:).*conj(spec(:,:,:,:,lt,:)),6)));
        end
        s1d(:,1,lt) = squeeze(sum(sum(sum(sum(sss,2),3),4),5));
        s1d(:,2,lt) = squeeze(max(max(max(max(sss,[],2),[],3),[],4),[],5));
    end
    s1d = s1d./repmat(max(s1d,[],1),[zf(1) 1]);
end
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
par.fax = hz; par.fig_str = figstr; par.lb = lb; % par structure for plot_mrsi
par.te = te;  % echo time [s]

%% saving result as mat file
if ~isempty(fname)
    if verb>0, fprintf('Saving results (fname=''%s'')\n',fname); end
    tmp = whos('spec');
    if verb>0, fprintf('Memory size of spec = %g [MB] (SI=Base10)\n',tmp.bytes*1d-6); end
    if tmp.bytes>100d6   % reduce size for data bigger than 100MB
        if verb>0
            fprintf('\t>100MB: converting spec (64bit) into single (32bit))\n'); 
        end
        spec = single(spec); % convert: double (64bit) to single (32bit) precision
    end
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','verb',...
        'do_coco','abs_fid','freq','lcm_basis','cor_intlv','do_ori',...
        'hz','s1d','nn','t','dim',...
        'si','mtx','nslices','nspec','nt','nc','par');
end


%% plotting
if ((verb>1) && (nspec>1))
    if (length(size(spec))>5) || ((verb==3)&&(isempty(fname)))
        warning('(length(size(spec))>5) || ((verb==3)&&(isempty(fname))); skipping plotting');
    else
        fprintf('Plotting\n');
        if abs_fid, pc = 0; else, pc = []; end
        if verb==3
            fid = figure('Visible','off');
        else
            figure(20); clf
        end
        imagesc_row(shiftdim(max(abs(spec),[],1),1),[],'row','',true);
        set(gcf,'name',[figstr '- image of max']);
        colormap gray
        if ~isempty(fname)
            print([fname '_maximage'],'-dpng','-r600','-painters');
        end
        if verb==3, close(fid); end
        
        if ~isempty(s1d)
            if verb==3
                fid = figure('Visible','off');
            else
                fid = figure(21); clf
            end
            plot(hz,s1d(:,:,1).');
            set(gca,'XDir','reverse');
            legend('rms','max');
            set(gcf,'name',[figstr '- pseudo 1D spectra']);
            grid on
            if ~isempty(fname)
                print([fname '_s1d'],'-dpng','-r600','-painters');
                if export_dcm>0
                    write_scdicom([fname '_s1d.dcm'],fid,h);
                end
            end
            if verb==3, close(fid); end
        end
        
        if ~isempty(fname)
            fname_plt_csi = [fname '_plot_csi'];
        else
            fname_plt_csi = [];
        end
        if export_dcm>0
            plot_csi(spec,pc,[],[figstr '- MRSI'],fname_plt_csi,11,h);
        else
            plot_csi(spec,pc,[],[figstr '- MRSI'],fname_plt_csi,11);
        end
        
        if verb~=3, plot_mrsi(spec,par); end
    end
end


%% exporting images of max as dicom into scanner database
if ~isempty(fname) && export_dcm
    if verb>0, fprintf('Storing dicom images\n'); end
    inp.SeriesDescription = sprintf('MRSI max:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
        ns,nx,ny,nz+nslices-1,nt,nc);
    if dim==3, inp.MRAcquisitionType = '3D'; end
    write_dicom(fname,shiftdim(max(abs(spec),[],1),1),h,inp);
end


%% write out to LCModel
if ~isempty(lcm_basis)
    if verb>0, fprintf('Writing LCModel files\n'); end
    if islogical(lcm_basis)
        if ~lcm_basis, return; end
        te = h.image.te*1d-3;
        switch te
            case 35,  lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press35/gamma_basis_TE35_B0_30_TE1_18_BW5000_N4096.basis';
            case 144, lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press144/gamma_basis_TE144_B0_30_TE1_20_BW5000_N4096.basis';
            otherwise
                lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press144/gamma_basis_TE144_B0_30_TE1_20_BW5000_N4096.basis';
                warning('TE(=%g)~=144[ms] (taking TE=144ms as basis)',h.image.te*1d-3); 
        end
    end
    if ~ischar(lcm_basis), error('~ischar(lcm_basis)'); end
    if verb>0, fprintf('Writing out raw data and control file for LCModel\n'); end
    if ~isempty(lb), warning('Data line-broadened'); end
    if (zf~=nspec),     warning('Data zero-filled');    end
    
    for l2=1:zf(2)
        for l3=1:zf(3)
            for l4=1:zf(4)
                for l5=1:zf(5)
                    lcm_fname = sprintf('lcm_matlab/P%.5d/P%.5d_l2_%g_l3_%g_l4_%g_l5_%g',...
                        h.image.rawrunnum,h.image.rawrunnum,l2,l3,l4,l5);
                    if verb>0, fprintf('%s\n',lcm_fname); end
                    write_lcmodel(dd(:,l2,l3,l4,l5),h,lcm_basis,[],lcm_fname);
                end
            end
        end
    end
end

if verb>0
    fprintf('recon_mrsi.m finished: ');
    fprintf('total reconstruction time = %.2f s\n',toc(timerVal));
end

end      % main function recon_mrsi.m


%% sub-functions
% Cartesian grid
function xx = sub_grid(mtx,dim)
for l=1:dim, x{l} = (-mtx(l)/2:mtx(l)/2-1); end

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
function [zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices)
do_zf_spat = false;
if ~isempty(zf)
    zf = [zf(:).', nn((length(zf)+1):6)];
    % if ~isempty(freq) && (zf(1)~=nspec)
    %    warning('zf(1)(=%g) ~= ns(=%g)',zf(1),nspec);
    % end
    if any(zf<0)             % negative numbers -> factor
        ii = zf<0;
        zf(ii) = -zf(ii).*nn(ii);
    end
    if any(abs(zf)<1d-10)    % 0 -> no zero-filling
        ii = abs(zf)<1d-10;
        zf(ii) = nn(ii);
    end
    if any(zf<nn)            % no cropping
        warning('zf<nn');
        fprintf('zf = '); disp(zf)
        fprintf('nn = '); disp(nn)
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = '); disp(zf)
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
    zf = nn(1);
end
end      %sub_zf
