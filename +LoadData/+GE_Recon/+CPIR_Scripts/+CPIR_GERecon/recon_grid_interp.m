function [bb,bbabs] = recon_grid_interp(d,h,wfn,mtx_reco,interp_factor,delay,lb,fname,...
    mask_fac,coco,plt,grdwrp)
%RECON_GRID Reconstruct 2D/3D non-Cartesian data
%[bb,bbabs] = recon_grid(d,h,wfn,mtx_reco,delay,lb,fname,mask_fac,coco,...
%                        plt,grdwrp)
%         d   Raw data (or P-File/ScanArchive fname)
%         h   Header from p-file (or empty)
%       wfn   Location of waveform .mat file (or wf structure)
%       mtx   Reconstruction matrix size (Fourier interp.)             ([])
%interp_factor Factor to interpolate by                                 (1)
%     delay   Gradient-acquisition delay                   [us]         (0)
%             1d or 2d=[WHOLE,ZOOM]
%        lb   Apodisation (Gauss time,spatial,exp time)[Hz,mtx,Hz]([0 0 0])
%     fname   Print <fname>.png and save reco as <fname>.mat           ([])
%             also export dicom if isdeployed
%  mask_fac   Threshold factor for masking                              (0)
%      coco   Coil combination via mri_coil_combine.m          ([1 0 1 16])
%             [method nper save_nococo mtx_coco]
%       plt   (1) plotting: 0=off,1=on  
%             (2) export as dicom          isdeployed->([0 1]) else ([1 0])
%    grdwrp   Apply gradwarp distortion correction (or gradcoil)     (true)
%
%        bb   Reconstructed data   (mtx(1),mtx(2),mtx(3),#timesteps,#coils)
%     bbabs   Magnitude data (gradwarp+coil combine)
%
%  2/2025  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
RFS_MAX_NSCANS = 16382;                % maximum number of scans
if ~exist('wfn','var'),   wfn = ''; end
if ~exist('mtx_reco','var'), mtx_reco = []; end
if isempty(mtx_reco),       mtx_reco = 0; end
if length(mtx_reco)<2,      mtx_reco = repmat(mtx_reco,[1 2]); end
if length(mtx_reco)<3,      mtx_reco = [mtx_reco,mtx_reco(end)]; end
mtx_reco = mtx_reco(:).';
if length(mtx_reco)~=3, error('length(mtx_reco)(=%d)~=3',length(mtx_reco)); end
if ~exist('interp_factor','var'), interp_factor = []; end
if isempty(interp_factor),       interp_factor = 1; end
if ~exist('delay','var'),   delay = []; end
if isempty(delay),          delay = 0; end
if length(delay)>1
    warning('length(delay)(=%d)>1; taking first delay',length(delay));
    delay = delay(1);
end
if ~exist('lb','var'),      lb = []; end
if isempty(lb),             lb = 0; end
if length(lb)<2,            lb = [lb 0]; end
if length(lb)<3,            lb = [lb 0]; end
if ~exist('fname','var'),   fname = []; end
if ~exist('mask_fac','var'), mask_fac = []; end
if isempty(mask_fac),       mask_fac = 0; end
if ~exist('coco','var'),    coco = []; end
if isempty(coco),           coco = 1; end
if length(coco)<2,          coco(2) = 0; end
if length(coco)<3,          coco(3) = 1; end
if length(coco)<4,          coco(4) = 16; end
if ~exist('plt','var'),     plt = []; end
if isempty(plt)
    if isdeployed,          plt = false;
    else,                   plt = true; end
end
if length(plt)<2
    if isdeployed,          plt(2) = true;
    else,                   plt(2) = false; end
end
if ~exist('grdwrp','var'),  grdwrp = []; end
if isempty(grdwrp),         grdwrp = true; end


%% other parameters
timerVal = tic;                        % record execution time
fidoff = 100;                          % figure ID offset
f0 = 0;                                % on-resonance reconstruction
fprintf('mtx=[%g,%g,%g], delay=%g, lb=[%g,%g,%g], plt=%g %g\n',...
    mtx_reco,delay,lb,plt(1),plt(2));


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
end
d = single(d);


%% create/check fname
if islogical(fname)
    if fname
        if isfield(h,'fname')
            fname = h.fname;
        else
            warning('no field h.fname');
            fname = [];
        end
    else
        fname = [];
    end
end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
    if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
end


%% reading waveform file
if iscell(wfn),  wfn = wfn{h.image.user3}; end
if isempty(wfn), wfn = get_wfname(h); end
if isstruct(wfn)
    wf = wfn;
else
    if isempty(regexpi(h.image.psdname,'spiral', 'once'))
        % load waveform structure from fidall.e
        if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
        if isempty(regexpi(wfn,'\.mat$', 'once')),  wfn = [wfn '.mat']; end
        if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
        wf = load_waveform(wfn,true,true);       % load waveform
    else
        if false
            % fidspiral.grad gradient waveform file
            wf = grad2waveform(wfn,h,true);
            fprintf('reducing delay by %g [us]\n',1d3/(4*h.rdb_hdr.bw));
            delay = delay - 1d3/(4*h.rdb_hdr.bw);
        else
            wf = calc_spiral_waveform(h);
        end
        delay = delay + 4;   % default in sonofrecon.py now -4; adding back
    end
end


%% extract fields required for reconstruction
if ~isfield(wf,'mtx'), error('wf.mtx not existing'); end
if ~isfield(wf,'k'),   error('wf.k not existing'); end
if ~isfield(wf,'dcf'), error('wf.dcf not existing'); end
if ~isfield(wf,'ind'), error('wf.ind not existing'); end
if ~isfield(wf,'t'),   error('wf.t not existing'); end
k = wf.k;
if abs(mtx_reco(1))<1
    k = wf.k/interp_factor; %X-fold interp
end
dcf = wf.dcf;
ind = wf.ind;
t = wf.t;
ns = sum(ind(end,:),2);
if size(t,1)==1
    if size(t,2)<=size(ind,2)
        t = repmat(t,[size(ind,1) 1]);
    end
end
t = t.'; t = t(:).';

if isempty(dcf),    error('dcf is empty'); end
if any(isnan(dcf)), error('dcf contains NaN'); end
if any(isinf(dcf)), error('dcf contains inf'); end

mtx_acq = wf.mtx(:).';       % nominal matrix resolution of acquisition
if size(mtx_acq,2)<2, mtx_acq = [mtx_acq mtx_acq]; end
if size(mtx_acq,2)<3, mtx_acq = [mtx_acq mtx_acq(1,end)]; end
niso = [1,1,1];              % non-isotropic factor for off-iso shift
if isfield(wf,'fov')
    if length(wf.fov)==3
        if (any(diff(wf.fov)) || any(diff(wf.mtx)))
            niso = wf.fov(1)./wf.fov;
        end
    end
end
clear wf
nslices = size(d,5);
dim = size(k,3);
if ((dim<2) || (dim>3)), error('dim(=%d) not 2D or 3D',dim); end

if (mtx_reco(1)<0)
    warning('mtx_reco(1)(=%g)<0; ignoring',mtx_reco(1)); 
    mtx_reco(1) = 0;
end
if abs(mtx_reco(1))<1
    mtx_reco = interp_factor*mtx_acq(:).'; %X-fold interp
else
    % scale k-space -> Fourier interpolation
    k = bsxfun(@times,k,shiftdim(mtx_acq(1,1:dim)./mtx_reco(1,1:dim),-1));
end
if dim==2
    mtx_reco(1,3) = nslices;
else
    if nslices>1, error('#slices(=%d)>1',nslices); end
end
fov = h.rdb_hdr.fov*1d-3;    % actual field-of-view [m]
k = single(k);
dcf = single(dcf);
t = single(t);


%% reshape data for vap_continue_loop
if ((sum(ind(:))>RFS_MAX_NSCANS) && (size(d,4)>1))
    fprintf('Reshaping data for vap_continue_loop\n');
    [n1,n2,n3,n4,n5,n6] = size(d);
    tmp = complex(zeros(n1*n4,n2,n3,1,n5,n6,class(d)));
    for l4=1:n4
        ii = (1:n1)+(l4-1)*n1;
        tmp(ii,:,:,1,:,:) = d(:,:,:,l4,:,:);
    end
    d = tmp;
    clear tmp
end


%% misc variables to checkings
bw = h.rdb_hdr.user0;
if ((abs(bw)<500) || ~isempty(regexpi(h.image.psdname,'spiral', 'once')))
    bw = h.rdb_hdr.bw*2d3; 
    fprintf('Setting bw=h.rdb_hdr.bw*2d3=%g[Hz]\n',bw);
end         
[nexc,n2,n3,n4,n5,ncoils] = size(d);
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
rot90fac = h.data_acq_tab.rotate;
trnsps = h.data_acq_tab.transpose;
if ~((trnsps(1)==0)||(trnsps(1)==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3: unknown image orientation',trnsps(1));
end
if n3~=1, warning('size(d,3)(=%g)~=1',n3); end


%% pre-processing of data
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx_reco.*[xloc yloc zloc]*1d-3/fov.*niso;
    if dim==2, shft = shft(1,1:2); end
else
    % shft = [];
    shft = zeros(1,dim);
end
for ld=1:dim
    % if ~isodd(mtx_reco(ld)), shft(ld) = shft(ld) + 0.5; end
    shft(ld) = shft(ld) + 0.5;
end

cart_down_fac = [];
dd = raw2grid(d,ischop(h),k,shft,cart_down_fac,[],ind,delay*1d-6,bw,false);
clear d ind
% size of reshaped data dd:
%  dim1=#coils; dim2=#indexed kspace data; dim3=#excitations; dim4=#slices
if any(abs(lb)>0), dd = ak_apodise(dd,k,lb(2),t,[lb(1) lb(3)],false); end
if abs(f0)>0.01
    fprintf('Off-resonance f0(=%g[Hz]) reconstruction\n',f0);
    phafu = exp(1i*2*pi*t*f0); 
    dd = bsxfun(@times,dd,phafu);
    clear phafu
end
clear t


%% actual reconstruction
ntimesteps = size(dd,3);
fprintf('Reconstructing/storing data as single (loop)\n');
if dim==3
    fprintf('3D gridding\n');
    if length(size(dd))>3, error('length(size(dd))(=%d)>3',length(size(dd))); end
    dd = reshape(permute(dd,[1 3 2]),[ncoils*ntimesteps,size(k,2)]);
    bb = gridding3dsingle(k,dd,dcf,mtx_reco);
    clear dd
    bb = permute(reshape(bb,[mtx_reco ncoils ntimesteps]),[1 2 3 5 4]);
else
    bb = complex(zeros([mtx_reco,ntimesteps,ncoils],'single'));
    for lt=1:ntimesteps
        fprintf('2D gridding %d/%d\n',lt,ntimesteps);
        for l3=1:nslices
            bb(:,:,l3,lt,:) = permute(...
                gridding(dd(:,:,lt,l3),k,dcf,mtx_reco(1,1:dim)),...
                [1 2 5 4 3]);
         end
    end
end
clear dd k dcf
bb = fliplr(bb);


%% separable 3D
if dim==2
    if isodd(bitget(h.rdb_hdr.data_collect_type,7))
        fprintf('Cartesian FFT along 3rd dimension\n');
        % bb = bb(:,:,end:-1:1,:,:,:);
        bb = fft(ifftshift(bb,3),[],3);
        bb = bb(:,:,end:-1:1,:,:,:);
        dim = 3;
    end
end


%% masking
if mask_fac>0
    mask = calc_mask(sqrt(mean(mean(bb.*conj(bb),5),4)),2,...
        0.1*mask_fac,[],[],true);
    if ncoils==1, bb = bsxfun(@times,bb,mask); end
else
    mask = [];
end


%% coil combination
if ncoils>1
    if (coco(3)>0) && ~isempty(fname)
        fprintf('Saving results to ''%s_nococo.mat''\n',fname);
        bb = single(bb);   % single (32bit) precision
        save([fname '_nococo.mat'],'-mat','-v7.3',...
            'bb','h','wfn','mtx_reco','delay','lb','fname',...
            'f0','cart_down_fac','dim','grdwrp');
    end
    if coco(1)>-1
        mtx_coco = [];
        if coco(4)>0, mtx_coco = coco(4); end
        bb = mri_coil_combine(bb,coco(2),'first',coco(1),mtx_coco,mask,true,dim);
    else
        warning('Skipping coil combination');
    end
end
clear mask
if any(isnan(bb(:)))
    warning('isnan(bb(:)): setting to 0');
    bb(isnan(bb)) = 0;
end
if any(isinf(bb(:)))
    warning('isinf(bb(:)): setting to 0');
    bb(isinf(bb)) = 0;
end


%% gradwarp distortion correction
if ischar(grdwrp)
    gradcoil = grdwrp;
    grdwrp = true;
else
    gradcoil = '';
end
if grdwrp
    if (dim==3)
        zfov = h.rdb_hdr.fov/niso(3);
    else
        zfov = [];
    end
    try
        bb = grad_warp(bb,h,gradcoil,(dim==3),'',zfov);
    catch ME
        warning(sprintf('grad_warp failed with Error\n\t%s',ME.message));
        grdwrp = -1;
    end
end


%% image orientation
bb = mri_orientation(bb,rot90fac,trnsps);


%% calculate real images
if size(bb,5)>1
    bbabs = abs(sqrt(sum(bb.*conj(bb),5)));
else
    if ~isreal(bb), bbabs = abs(bb);
    else,           bbabs = bb; end
end


%% saving result as mat file
if ~isempty(fname)
    fprintf('Saving results to ''%s.mat''\n',fname);
    bb = single(bb);   % single (32bit) precision
    bbabs = single(bbabs);
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','bbabs','h','wfn','mtx_reco','delay','lb','fname',...
        'f0','cart_down_fac','dim','grdwrp','niso');
end


%% plotting: 3 orientations
if (plt(1)>0) || ~isempty(fname)
    fprintf('Plotting results\n');
    if (plt(1)>0)
        fid = figure(fidoff);
    else
        fid = figure('visible','off');
    end
    if dim==3      % plot 3D
        imagesc_ind3d(bbabs,'','','',false,true,[]);
    else
        rshp = false;
        if (ntimesteps==1) && (nslices>8), rshp = true; end
        if (ntimesteps>8) && (nslices==1), rshp = true; end
        imagesc_row(bbabs,'',[],rshp,true);
    end
    colormap(gray);
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(fid,'name',figstr);
    drawnow
    if ~isempty(fname)
        print([fname '.png'],'-dpng','-r300');
    end
    if ~(plt(1)>0), close(fid); end
end


%% exporting images as dicom into scanner database
if ~isempty(fname) && (plt(2)>0)
    fprintf('Exporting images to dicom\n');
    inp.SeriesDescription  = sprintf('MNSRP;%dD%g;TS%g;%s',...
        dim,mtx_reco(1),ntimesteps,h.series.se_desc);
    inp.MRAcquisitionType = sprintf('%dD',dim);
    if any(abs(diff(niso))>1d-5)
        spacing = round(h.rdb_hdr.fov/mtx_reco(3)/niso(3)*10)/10;
        inp.SpacingBetweenSlices = spacing;
    end
    % matlab-based dicomwrite.m without template
    % dicom header generated from p-file header
    write_dicom_cchmc(fname,bbabs,h,inp,[],[],[],3);
end


%% print execution time
ttv = toc(timerVal);
fprintf('recon_grid: runtime = %g [m] %02g [s]\n',...
    floor(ttv/60),round(mod(ttv,60)));

if nargout<1, clear bb; end

end      % recon_grid.m
