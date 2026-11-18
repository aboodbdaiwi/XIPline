function [bb,bbabs,gx] = recon_grid3d_gx_XIPline(d,h,wfn,mtx_reco,interp_factor,delay,lb,fname,...
    mask_fac,coco,plt,grdwrp)
% Native GE 3D gridding with XIPline gas-exchange outputs.
% The imaging reconstruction follows the original recon_grid3d_gx.m.
% XIPline additions: extract two bonus FIDs, expose raw interleaves and
% trajectory, use the CPIR spatial shift, return metadata in gx.
%[bb,bbabs] = recon_grid3d_gx(d,h,wfn,mtx_reco,delay,lb,fname,mask_fac,coco,...
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
gx = struct();


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
if length(lb)<2,            lb = [lb lb]; end
if length(lb)<3,            lb = [lb lb]; end
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
    [d,h] = LoadData.GE_Recon.fidall.read_p(d,true);
end
d = single(conj(d));


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
        wf = LoadData.GE_Recon.fidall.load_waveform(wfn,true,true);       % load waveform
    else
        if false
            % fidspiral.grad gradient waveform file
            wf = LoadData.GE_Recon.fidall.grad2waveform(wfn,h,true);
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
% XIPline: save metadata needed by the downstream gas-exchange pipeline.
gx.mtx_acq = mtx_acq;
gx.niso = niso;
gx.wfn = wfn;
gx.header = h;
gx.WaveformFlips_deg = [];
gx.WaveformFreqs_Hz = [];
gx.ImagingFreqJump_Hz = NaN;
if isfield(wf,'flips'), gx.WaveformFlips_deg = double(wf.flips(:).'); end
if isfield(wf,'freqs'), gx.WaveformFreqs_Hz = double(wf.freqs(:).'); end
if numel(gx.WaveformFreqs_Hz)==2
    gx.ImagingFreqJump_Hz = abs(diff(gx.WaveformFreqs_Hz));
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


%% XIPline: extract the two appended spectroscopy FIDs (if present)
% Original imaging views come first, then gas FID, then dissolved FID.
% The last two whole DCF profiles must be zero. Keep unpadded FIDs for
% fitting and separate 512-point versions for display (no new information).
gx.GasBonusFID = [];
gx.DissolvedBonusFID = [];
gx.GasBonusFIDAcquired = [];
gx.DissolvedBonusFIDAcquired = [];
gx.GasBonusFIDPerCoil = [];
gx.DissolvedBonusFIDPerCoil = [];
gx.BonusTime_s = [];
gx.BonusTimePadded_s = [];
gx.BonusDwell_s = NaN;
gx.BonusAcquisitionIndices = [];
gx.BonusCoilIndex = NaN;
gx.BonusFIDStatus = 'not detected';

nViews = size(ind,1);
nRead = nnz(ind(1,:));
canCheck = size(d,1)==nViews && nViews>=4 && nRead>=32 && mod(nViews,2)==0 && ...
    size(d,3)==1 && size(d,4)==1 && size(d,5)==1 && ...
    all(sum(ind,2)==nRead) && size(d,2)>=size(ind,2) && ...
    numel(dcf)==nRead*nViews && size(k,2)==nRead*nViews && ...
    numel(t)==nRead*nViews;

if canCheck
    dcfByView = reshape(dcf,nRead,nViews);
    isBonus = all(dcfByView(:,end-1:end)==0,'all') && ...
        any(dcfByView(:,end-2)~=0) && ...
        numel(gx.WaveformFlips_deg)==2 && ...
        gx.WaveformFlips_deg(1)<gx.WaveformFlips_deg(2) && ...
        numel(gx.WaveformFreqs_Hz)==2 && ...
        gx.WaveformFreqs_Hz(1)~=gx.WaveformFreqs_Hz(2) && ...
        isequal(ind(end,:),ind(end-1,:));

    if isBonus
        adc = find(ind(end,:));
        nCoils = size(d,6);
        gasFID = reshape(d(end-1,adc,1,1,1,:),nRead,nCoils);
        disFID = reshape(d(end,adc,1,1,1,:),nRead,nCoils);
        [~,bestCoil] = max(sum(abs(disFID).^2,1));

        gx.GasBonusFIDPerCoil = gasFID;
        gx.DissolvedBonusFIDPerCoil = disFID;
        gx.GasBonusFIDAcquired = gasFID(:,bestCoil);
        gx.DissolvedBonusFIDAcquired = disFID(:,bestCoil);

        % Compatibility: display/export arrays zero-filled to >=512 points.
        Npad = max(512,nRead);
        gx.GasBonusFID = [gx.GasBonusFIDAcquired; complex(zeros(Npad-nRead,1,'single'))];
        gx.DissolvedBonusFID = [gx.DissolvedBonusFIDAcquired; complex(zeros(Npad-nRead,1,'single'))];
        gx.BonusCoilIndex = bestCoil;
        gx.BonusAcquisitionIndices = [nViews-1,nViews];
        tv = reshape(t,nRead,nViews);
        gx.BonusTime_s = double(tv(:,end));
        gx.BonusDwell_s = median(diff(gx.BonusTime_s));
        gx.BonusTimePadded_s = gx.BonusTime_s(1) + ...
            (0:Npad-1)' * gx.BonusDwell_s;
        gx.BonusFIDStatus = 'extracted: gas and dissolved bonus FIDs';
        fprintf('Bonus FIDs: views %d, %d; %d acquired samples; %d coil(s)\n', ...
            nViews-1,nViews,nRead,nCoils);

        % Remove bonus views from ALL imaging arrays before interleave split.
        d = d(1:end-2,:,:,:,:,:);
        ind = ind(1:end-2,:);
        kByView = reshape(k,nRead,nViews,size(k,3));
        k = reshape(kByView(:,1:end-2,:),1,nRead*(nViews-2),size(k,3));
        dcf = reshape(dcfByView(:,1:end-2),1,[]);
        t = reshape(tv(:,1:end-2),1,[]);
    end
end
if isempty(gx.BonusAcquisitionIndices)
    warning('GE bonus FIDs not identified. No automatic bonus subtraction.');
end

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

%% Split interleaves
d1 = d(1:2:end,:); % gas/dissolved interleave 1
d2 = d(2:2:end,:); % gas/dissolved interleave 2
gx.RawInterleave1 = d1;
gx.RawInterleave2 = d2;
d=cat(3,d1,d2);
d=reshape(d,size(d,1),size(d,2),1,size(d,3));
[nproj, nread, ~, ~] = size(d);

ind = ind(1:2:end,:);
nacq=sum(sum(ind,1)>0);

dcf_tmp = reshape(dcf,[],nproj);
dcfKeep = dcf_tmp(1:nacq,:);
dcf = dcfKeep(:)';

k_tmp = reshape(k,[],nproj,3);
kKeep = k_tmp(1:nacq,:,:);
k = reshape(kKeep,1,nacq*nproj,3);

t_tmp = reshape(t,[],nproj);
tKeep = t_tmp(1:nacq,:);
t = tKeep(:)';

% XIPline: preserve raw imaging k-space and matched trajectory.
gx.XeTraj = permute(kKeep,[3 1 2]);        % 3 x readout x profile
gx.Interleave1KSpace = transpose(d1(:,1:nacq));
gx.Interleave2KSpace = transpose(d2(:,1:nacq));
gx.nacq = nacq;
gx.nread = nread;
gx.nproj = nproj;
gx.dcf = dcfKeep;
gx.t = tKeep;


%% XIPline: same GE raw2grid step, with CPIR centering shift
shft = 0.5 * double(mod(mtx_reco(1:dim),2)==0);
shft(2) = shft(2) + 20;
fprintf('XIPline spatial shift: [%s] voxels\n',num2str(shft));

cart_down_fac = [];
dd = LoadData.GE_Recon.fidall.raw2grid(d,LoadData.GE_Recon.fidall.ischop(h),k,shft,cart_down_fac,[],ind,delay*1d-6,bw,false);
clear d ind
% size of reshaped data dd:
%  dim1=#coils; dim2=#indexed kspace data; dim3=#excitations; dim4=#slices
if any(abs(lb)>0), dd = LoadData.GE_Recon.fidall.ak_apodise(dd,k,lb(2),t,[lb(1) lb(3)],false); end
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
    bb = LoadData.GE_Recon.fidall.gridding3dsingle(k,dd,dcf,mtx_reco);
    clear dd
    bb = permute(reshape(bb,[mtx_reco ncoils ntimesteps]),[1 2 3 5 4]);
else
    bb = complex(zeros([mtx_reco,ntimesteps,ncoils],'single'));
    for lt=1:ntimesteps
        fprintf('2D gridding %d/%d\n',lt,ntimesteps);
        for l3=1:nslices
            bb(:,:,l3,lt,:) = permute(...
                LoadData.GE_Recon.fidall.gridding(dd(:,:,lt,l3),k,dcf,mtx_reco(1,1:dim)),...
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
        bb = LoadData.GE_Recon.fidall.grad_warp(bb,h,gradcoil,(dim==3),'',zfov);
    catch ME
        warning(sprintf('grad_warp failed with Error\n\t%s',ME.message));
        grdwrp = -1;
    end
end


%% image orientation
bb = LoadData.GE_Recon.fidall.mri_orientation(bb,rot90fac,trnsps);


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
        LoadData.GE_Recon.fidall.imagesc_ind3d(bbabs,'','','',false,true,'row');
    else
        rshp = false;
        if (ntimesteps==1) && (nslices>8), rshp = true; end
        if (ntimesteps>8) && (nslices==1), rshp = true; end
        LoadData.GE_Recon.fidall.imagesc_row(bbabs,'',[],rshp,true);
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


%% XIPline metadata (native gridding output itself is unchanged)
gx.mtx_reco = mtx_reco;
gx.interp_factor = interp_factor;
gx.delay = delay;
gx.lb = lb;
gx.mask_fac = mask_fac;
gx.coco = coco;
gx.plt = plt;
gx.grdwrp = grdwrp;
gx.dim = dim;
gx.ntimesteps = ntimesteps;
gx.ncoils = ncoils;
gx.bw = bw;
gx.fov = fov;
gx.fov_m = fov;

%% print execution time
ttv = toc(timerVal);
fprintf('recon_grid: runtime = %g [m] %02g [s]\n',...
    floor(ttv/60),round(mod(ttv,60)));

if nargout<1, clear bb; end

end      % recon_grid3d_gx_XIPline.m
